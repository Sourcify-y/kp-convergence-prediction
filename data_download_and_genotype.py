#!/usr/bin/env python3
"""
data_download_and_genotype.py

Phase 1, Step 2: pull every K. pneumoniae assembly from NCBI, get the raw
BioSample collection date (never backfilled from submission date), download
the genomes, run Kleborate, and merge everything into one metadata CSV
(kp_metadata_FULL.csv).

Needs: biopython, pandas, tqdm, the NCBI datasets CLI, and Kleborate
(plus blast+ and mash on PATH, see the Kleborate README).

Set these to avoid NCBI throttling:
    export NCBI_EMAIL="you@example.com"
    export NCBI_API_KEY="..."     # optional but much faster

Usage (each stage is checkpointed to --outdir, so you can stop and resume;
delete a stage's checkpoint file to redo it):

    python3 data_download_and_genotype.py list        --outdir ./kp_pull
    python3 data_download_and_genotype.py biosample   --outdir ./kp_pull
    python3 data_download_and_genotype.py download    --outdir ./kp_pull --limit 500
    python3 data_download_and_genotype.py kleborate   --outdir ./kp_pull
    python3 data_download_and_genotype.py merge       --outdir ./kp_pull

    python3 data_download_and_genotype.py all --outdir ./kp_pull --limit 500

--limit 0 (default) pulls everything, which is 20k+ assemblies, so budget
~5-10 GB disk per 1,000 genomes. Try --limit 500 first as a sanity check.
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

TAXON = "Klebsiella pneumoniae"
ENTREZ_BATCH = 200          # biosample esummary batch size
ENTREZ_SLEEP = 0.34         # ~3 req/sec without an API key
DOWNLOAD_BATCH = 100        # genomes per `datasets download` call

ASSEMBLY_LIST_FILE = "01_assembly_list.jsonl"
BIOSAMPLE_FILE = "02_biosample_metadata.jsonl"
DOWNLOAD_MANIFEST = "03_download_manifest.csv"
GENOME_DIR = "genomes"
KLEBORATE_OUT = "04_kleborate_full/klebsiella_pneumo_complex__mlst_output.txt"
FINAL_CSV = "kp_metadata_FULL.csv"


def run(cmd, **kwargs):
    print(f"[cmd] {' '.join(cmd)}")
    return subprocess.run(cmd, check=True, **kwargs)


def check_tool(name, install_hint):
    from shutil import which
    if which(name) is None:
        sys.exit(f"ERROR: required tool '{name}' not found on PATH.\n"
                  f"  Install it: {install_hint}")


def load_jsonl(path):
    if not Path(path).exists():
        return []
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


def append_jsonl(path, records):
    with open(path, "a") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")


# stage 1: list every K. pneumoniae assembly via the datasets CLI
def stage_list(outdir: Path, limit: int):
    check_tool("datasets", "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/")
    out_path = outdir / ASSEMBLY_LIST_FILE
    if out_path.exists():
        print(f"[list] {out_path} already exists, skipping (delete it to redo).")
        return

    print(f"[list] Querying NCBI datasets for all '{TAXON}' assemblies...")
    # one JSON record per assembly, includes the linked biosample accession
    cmd = ["datasets", "summary", "genome", "taxon", TAXON,
           "--as-json-lines"]
    if limit:
        cmd += ["--limit", str(limit)]

    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    print(f"[list] Got {len(lines)} assembly records.")

    records = []
    for line in lines:
        rec = json.loads(line)
        accession = rec.get("accession")
        biosample = (rec.get("assembly_info", {})
                        .get("biosample", {})
                        .get("accession"))
        organism = rec.get("organism", {}).get("organism_name", "")
        records.append({
            "assembly_accession": accession,
            "biosample_accession": biosample,
            "organism_name": organism,
            "submission_date": rec.get("assembly_info", {}).get("biosample", {}).get("submission_date"),
            "assembly_level": rec.get("assembly_info", {}).get("assembly_level"),
        })

    append_jsonl(out_path, records)
    print(f"[list] Wrote {len(records)} records -> {out_path}")


# stage 2: raw collection_date, country, strain etc. from BioSample
def stage_biosample(outdir: Path):
    from Bio import Entrez
    Entrez.email = os.environ.get("NCBI_EMAIL", "")
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key
    if not Entrez.email:
        print("[biosample] WARNING: NCBI_EMAIL not set. NCBI asks you to "
              "set this; requests will still work but may be throttled "
              "or blocked without it.")

    assemblies = load_jsonl(outdir / ASSEMBLY_LIST_FILE)
    if not assemblies:
        sys.exit("[biosample] No assembly list found — run 'list' stage first.")

    done_path = outdir / BIOSAMPLE_FILE
    done_ids = {r["biosample_accession"] for r in load_jsonl(done_path)}
    todo = [a for a in assemblies
            if a["biosample_accession"] and a["biosample_accession"] not in done_ids]
    print(f"[biosample] {len(todo)} biosamples left to fetch "
          f"({len(done_ids)} already done).")

    i = 0
    while i < len(todo):
        batch = todo[i:i + ENTREZ_BATCH]
        ids = [b["biosample_accession"] for b in batch]

        data = None
        for attempt in range(1, 6):
            try:
                # esummary wants numeric UIDs, not SAMN accessions, so resolve
                # them with one OR-joined esearch per batch
                term = " OR ".join(f"{_strip_prefix(x)}[Accession]" for x in ids)
                search_handle = Entrez.esearch(db="biosample", term=term,
                                               retmax=len(ids))
                search_result = Entrez.read(search_handle)
                search_handle.close()
                uids = search_result.get("IdList", [])

                if not uids:
                    data = {"DocumentSummarySet": {"DocumentSummary": []}}
                    break

                handle = Entrez.esummary(db="biosample", id=",".join(uids),
                                         retmode="xml")
                data = Entrez.read(handle, validate=False)
                handle.close()
                break
            except Exception as e:
                wait = 5 * attempt
                print(f"[biosample] batch {i} failed on attempt {attempt}/5 "
                      f"({e}), retrying in {wait}s...")
                time.sleep(wait)

        if data is None:
            sys.exit(f"[biosample] batch {i} failed 5 times in a row — "
                      f"stopping so nothing is silently skipped. Fix the "
                      f"underlying error (see messages above) and rerun "
                      f"this stage; already-fetched biosamples are cached "
                      f"in {done_path} and won't be re-fetched.")

        out_records = []
        found_accs = set()
        for doc in data.get("DocumentSummarySet", {}).get("DocumentSummary", []):
            sample_xml = doc.get("SampleData", "")
            parsed = _parse_biosample_xml(sample_xml)
            acc = doc.get("Accession", "")
            parsed["biosample_accession"] = acc
            found_accs.add(acc)
            out_records.append(parsed)

        # accessions esearch couldn't resolve (withdrawn/private/typos) still
        # get a blank record so they show up in the audit
        for acc in ids:
            if acc not in found_accs:
                out_records.append({
                    "biosample_accession": acc,
                    "collection_date_raw": None,
                    "country": None,
                    "strain": None,
                    "isolation_source": None,
                    "host": None,
                    "submission_date": None,
                })

        append_jsonl(done_path, out_records)
        print(f"[biosample] {i + len(batch)}/{len(todo)} done")
        time.sleep(ENTREZ_SLEEP)
        i += ENTREZ_BATCH

    print(f"[biosample] Complete -> {done_path}")


def _strip_prefix(acc):
    return acc


def _parse_biosample_xml(xml_blob: str) -> dict:
    # pull collection_date, geo_loc_name, strain, isolation_source, host and
    # submission_date out of a BioSample XML blob
    import xml.etree.ElementTree as ET
    out = {
        "collection_date_raw": None,
        "country": None,
        "strain": None,
        "isolation_source": None,
        "host": None,
        "submission_date": None,
    }
    if not xml_blob:
        return out
    try:
        root = ET.fromstring(xml_blob)
    except ET.ParseError:
        return out
    out["submission_date"] = root.get("submission_date")
    for attr in root.iter("Attribute"):
        name = (attr.get("attribute_name") or attr.get("harmonized_name") or "").lower()
        val = (attr.text or "").strip()
        if not val:
            continue
        if "collection_date" in name or name == "collection date":
            out["collection_date_raw"] = val
        elif "geo_loc_name" in name or "geographic location" in name:
            out["country"] = val.split(":")[0].strip()
        elif name == "strain":
            out["strain"] = val
        elif "isolation_source" in name:
            out["isolation_source"] = val
        elif name == "host":
            out["host"] = val
    return out


# stage 3: download FASTAs for everything that has a raw collection date
def stage_download(outdir: Path, limit: int):
    check_tool("datasets", "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/")
    genome_dir = outdir / GENOME_DIR
    genome_dir.mkdir(parents=True, exist_ok=True)

    assemblies = load_jsonl(outdir / ASSEMBLY_LIST_FILE)
    biosamples = {r["biosample_accession"]: r for r in load_jsonl(outdir / BIOSAMPLE_FILE)}

    manifest_path = outdir / DOWNLOAD_MANIFEST
    already = set()
    if manifest_path.exists():
        with open(manifest_path) as f:
            already = {row["assembly_accession"] for row in csv.DictReader(f)}

    # skip genomes with no raw collection date, they're unusable anyway
    targets = [a for a in assemblies
               if a["assembly_accession"] not in already
               and biosamples.get(a["biosample_accession"], {}).get("collection_date_raw")]
    if limit:
        targets = targets[:limit]
    print(f"[download] {len(targets)} genomes to download "
          f"({len(already)} already on disk).")

    write_header = not manifest_path.exists()
    with open(manifest_path, "a", newline="") as mf:
        writer = csv.DictWriter(mf, fieldnames=["assembly_accession", "path"])
        if write_header:
            writer.writeheader()

        for i in range(0, len(targets), DOWNLOAD_BATCH):
            batch = targets[i:i + DOWNLOAD_BATCH]
            accs = [b["assembly_accession"] for b in batch]
            zip_path = genome_dir / f"batch_{i}.zip"
            cmd = ["datasets", "download", "genome", "accession",
                   *accs, "--include", "genome", "--filename", str(zip_path)]
            try:
                run(cmd)
            except subprocess.CalledProcessError as e:
                print(f"[download] batch {i} failed: {e} — skipping batch")
                continue

            extract_dir = genome_dir / f"batch_{i}_extracted"
            run(["unzip", "-oq", str(zip_path), "-d", str(extract_dir)])
            zip_path.unlink()

            data_dir = extract_dir / "ncbi_dataset" / "data"
            for acc in accs:
                match = list(data_dir.glob(f"{acc}*/*.fna"))
                if match:
                    writer.writerow({"assembly_accession": acc, "path": str(match[0])})
            mf.flush()
            print(f"[download] {i + len(batch)}/{len(targets)} done")

    print(f"[download] Complete -> {manifest_path}")


# stage 4: Kleborate (ST/CC/AMR/virulence/capsule/plasmid calls)
def stage_kleborate(outdir: Path):
    check_tool("kleborate", "pip install kleborate  (also needs blastn, mash on PATH)")
    manifest_path = outdir / DOWNLOAD_MANIFEST
    if not manifest_path.exists():
        sys.exit("[kleborate] No download manifest — run 'download' stage first.")
    with open(manifest_path) as f:
        rows = list(csv.DictReader(f))
    fasta_paths = [r["path"] for r in rows if Path(r["path"]).exists()]
    print(f"[kleborate] Running Kleborate on {len(fasta_paths)} genomes...")
    out_path = outdir / KLEBORATE_OUT
    modules = (
        "klebsiella_pneumo_complex__mlst,klebsiella_pneumo_complex__amr,"
        "klebsiella_pneumo_complex__resistance_score,"
        "klebsiella_pneumo_complex__resistance_class_count,"
        "klebsiella_pneumo_complex__resistance_gene_count,"
        "klebsiella_pneumo_complex__virulence_score,"
        "klebsiella_pneumo_complex__kaptive,"
        "klebsiella__ybst,klebsiella__abst,klebsiella__smst,"
        "klebsiella__cbst,klebsiella__rmst,klebsiella__rmpa2"
    )
    BATCH_SIZE = 500
    batch_outputs = []
    for i in range(0, len(fasta_paths), BATCH_SIZE):
        batch = fasta_paths[i:i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        batch_out = outdir / f"{KLEBORATE_OUT}.batch{batch_num}.tmp"
        print(f"[kleborate] Batch {batch_num} ({len(batch)} genomes)...")
        cmd = ["kleborate", "-a", *batch, "-m", modules, "-o", str(batch_out)]
        run(cmd)
        batch_outputs.append(batch_out)
    # stitch the batch outputs together, keeping only the first header
    with open(out_path, "w") as outfile:
        for idx, bpath in enumerate(batch_outputs):
            with open(bpath) as infile:
                lines = infile.readlines()
            outfile.writelines(lines if idx == 0 else lines[1:])
            Path(bpath).unlink()
    print(f"[kleborate] Complete -> {out_path}")


# stage 5: merge into the final CSV
FINAL_FIELDS = [
    "genome_id", "assembly_accession", "biosample_accession",
    "collection_date_raw", "collection_date_quality", "collection_year",
    "submission_date", "date_suspicious_match",
    "country", "source", "strain",
    "ST", "gapA", "infB", "mdh", "pgi", "phoE", "rpoB", "tonB",
    "AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired",
    "Gly_acquired", "MLS_acquired", "Phe_acquired", "Rif_acquired",
    "Sul_acquired", "Tet_acquired", "Tgc_acquired", "Tmt_acquired",
    "Bla_acquired", "Bla_inhR_acquired", "Bla_ESBL_acquired",
    "Bla_ESBL_inhR_acquired", "Bla_Carb_acquired", "Bla_chr",
    "SHV_mutations", "Omp_mutations", "Col_mutations", "Flq_mutations",
    "truncated_resistance_hits", "spurious_resistance_hits",
    "resistance_score", "num_resistance_classes", "num_resistance_genes",
    "virulence_score", "spurious_virulence_hits",
    "K_locus", "K_type", "K_locus_confidence", "K_locus_problems",
    "K_locus_identity", "K_Missing_expected_genes",
    "O_locus", "O_type", "O_locus_confidence", "O_locus_problems",
    "O_locus_identity", "O_Missing_expected_genes",
    "YbST", "Yersiniabactin", "ybtS", "ybtX", "ybtQ", "ybtP", "ybtA",
    "irp2", "irp1", "ybtU", "ybtT", "ybtE", "fyuA", "spurious_ybt_hits",
    "AbST", "Aerobactin", "iucA", "iucB", "iucC", "iucD", "iutA",
    "spurious_abst_hits",
    "SmST", "Salmochelin", "iroB", "iroC", "iroD", "iroN",
    "spurious_smst_hits",
    "CbST", "Colibactin", "clbA", "clbB", "clbC", "clbD", "clbE",
    "clbF", "clbG", "clbH", "clbI", "clbL", "clbM", "clbN", "clbO",
    "clbP", "clbQ", "spurious_clb_hits",
    "RmST", "RmpADC", "rmpA", "rmpD", "rmpC", "spurious_rmst_hits",
    "rmpA2",
]


def _classify_date_quality(raw: str):
    if not raw:
        return "missing", None
    raw = raw.strip()
    m = re.fullmatch(r"\d{4}-\d{2}-\d{2}", raw)
    if m:
        return "exact", int(raw[:4])
    m = re.fullmatch(r"\d{4}-\d{2}", raw)
    if m:
        return "month-year", int(raw[:4])
    m = re.fullmatch(r"\d{4}", raw)
    if m:
        return "year-only", int(raw)
    if re.search(r"\d{4}", raw):
        return "inferred", int(re.search(r"\d{4}", raw).group())
    return "missing", None


def stage_merge(outdir: Path):
    assemblies = {a["assembly_accession"]: a for a in load_jsonl(outdir / ASSEMBLY_LIST_FILE)}
    biosamples = {r["biosample_accession"]: r for r in load_jsonl(outdir / BIOSAMPLE_FILE)}

    kleb_path = outdir / KLEBORATE_OUT
    kleb_rows = {}
    kb_prefix_len = None
    if kleb_path.exists():
        with open(kleb_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # strip Kleborate's module prefixes down to bare field names,
                # e.g. "klebsiella_pneumo_complex__kaptive__K_locus" -> "K_locus"
                bare_row = {}
                for k, v in row.items():
                    bare = k.split("__")[-1] if "__" in k else k
                    bare_row[bare] = v
                genome_key = bare_row.get("strain")
                kleb_rows[genome_key] = bare_row
    else:
        print("[merge] WARNING: no Kleborate output found — genotype "
              "columns will be blank. Run 'kleborate' stage first "
              "if you want them filled in.")

    out_path = outdir / FINAL_CSV
    n_written = 0
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FINAL_FIELDS)
        writer.writeheader()

        for acc, a in assemblies.items():
            bs = biosamples.get(a["biosample_accession"], {})
            raw_date = bs.get("collection_date_raw")
            quality, year = _classify_date_quality(raw_date)
            submission_date = bs.get("submission_date")
            suspicious = bool(raw_date and submission_date and
                              raw_date[:10] == str(submission_date)[:10])

            kleb = None
            for key, row in kleb_rows.items():
                if key and acc in key:
                    kleb = row
                    break

            record = {
                "genome_id": acc,
                "assembly_accession": acc,
                "biosample_accession": a.get("biosample_accession"),
                "collection_date_raw": raw_date,
                "collection_date_quality": quality,
                "collection_year": year,
                "submission_date": submission_date,
                "date_suspicious_match": suspicious,
                "country": bs.get("country"),
                "source": bs.get("isolation_source"),
                "strain": bs.get("strain"),
            }
            # everything else in FINAL_FIELDS comes straight from the Kleborate row
            already_set = set(record.keys())
            for field in FINAL_FIELDS:
                if field in already_set:
                    continue
                record[field] = kleb.get(field) if kleb else None

            writer.writerow(record)
            n_written += 1

    print(f"[merge] Wrote {n_written} rows -> {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stage", choices=["list", "biosample", "download",
                                          "kleborate", "merge", "all"])
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--limit", type=int, default=0,
                        help="Cap the number of genomes (0 = no limit). "
                             "Strongly recommended for a first test run.")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    stages = {
        "list": lambda: stage_list(args.outdir, args.limit),
        "biosample": lambda: stage_biosample(args.outdir),
        "download": lambda: stage_download(args.outdir, args.limit),
        "kleborate": lambda: stage_kleborate(args.outdir),
        "merge": lambda: stage_merge(args.outdir),
    }

    if args.stage == "all":
        for name, fn in stages.items():
            print(f"\n===== STAGE: {name} =====")
            fn()
    else:
        stages[args.stage]()


if __name__ == "__main__":
    main()