#!/usr/bin/env python3

import os
import csv
import sys
from pathlib import Path
from collections import Counter, defaultdict

# merged set (Mac 4,718 + Windows 3,600 = 8,318 genomes)
MOB_OUT = Path(os.environ.get("PHASE10_DIR", "phase10")) / "mob_out_all"
OUT = Path(os.environ.get("PHASE10_DIR", "phase10")) / "mob_consolidated"
OUT.mkdir(parents=True, exist_ok=True)

EXPECTED = 8318

CONTIG_OUT = OUT / "contig_to_cluster.tsv"
MOBTYPER_OUT = OUT / "mobtyper_results_all.tsv"
REPLICON_OUT = OUT / "replicon_matrix.csv"
SUMMARY_OUT = OUT / "plasmid_summary.csv"
QC_OUT = OUT / "code1_qc.txt"


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


def split_multi(value):
    if not value:
        return []
    value = value.strip()
    if value in {"", "-", "None", "none", "nan"}:
        return []
    return [x.strip() for x in value.split(";") if x.strip()]


# 1. find the per-genome result directories

if not MOB_OUT.exists():
    sys.exit(f"ERROR: MOB output directory does not exist: {MOB_OUT}")

genome_dirs = sorted(
    p for p in MOB_OUT.iterdir()
    if p.is_dir() and not p.name.startswith("_")
)

genome_ids = [p.name for p in genome_dirs]

if len(genome_ids) != EXPECTED:
    sys.exit(
        f"ERROR: found {len(genome_ids)} genome directories; "
        f"expected exactly {EXPECTED}. "
        f"Do not run consolidation yet."
    )

# 2. check the raw MOB files are there

missing_contig = []
missing_mobtyper = []
no_plasmid_genomes = []

for gid in genome_ids:
    gdir = MOB_OUT / gid

    if not (gdir / "contig_report.txt").is_file():
        missing_contig.append(gid)

    if not (gdir / "mobtyper_results.txt").is_file():
        # MOB-suite doesn't write mobtyper_results.txt when a genome has no
        # plasmid contigs (spot-checked 2026-09), so a missing file is only an
        # error if contig_report.txt shows plasmid contigs.
        contig_path = gdir / "contig_report.txt"
        has_plasmid_contig = False
        if contig_path.is_file():
            rows = read_tsv(contig_path)
            # molecule_type == "plasmid" is the real signal; rep_type(s) alone
            # isn't enough (GCF_001031965.1 has rep hits on a "chromosome" row)
            has_plasmid_contig = any(
                r.get("molecule_type", "").strip().lower() == "plasmid"
                for r in rows
            )
        if has_plasmid_contig:
            missing_mobtyper.append(gid)
        else:
            no_plasmid_genomes.append(gid)

if missing_contig or missing_mobtyper:
    print(f"ERROR: missing contig_report.txt: {len(missing_contig)}")
    print(f"ERROR: missing mobtyper_results.txt (genome HAS plasmid contigs, real error): {len(missing_mobtyper)}")

    if missing_contig:
        print("First missing contig reports:")
        print("\n".join(missing_contig[:20]))

    if missing_mobtyper:
        print("First genomes with plasmid contigs but no mobtyper file (real error):")
        print("\n".join(missing_mobtyper[:20]))

    sys.exit("QC FAILED — no consolidated files were generated.")

print(f"No-plasmid genomes (no mobtyper_results.txt, expected/valid): {len(no_plasmid_genomes)}")


# 3. consolidate contig_report.txt

contig_rows = []
contig_counts = Counter()

for gid in genome_ids:
    path = MOB_OUT / gid / "contig_report.txt"
    rows = read_tsv(path)

    if not rows:
        contig_counts[gid] = 0
        continue

    contig_counts[gid] = len(rows)

    for r in rows:
        new_row = dict(r)
        new_row["genome_id"] = gid
        contig_rows.append(new_row)


# keep MOB-suite's own column order, with genome_id first
contig_fields = ["genome_id"]

if contig_rows:
    for k in contig_rows[0].keys():
        if k != "genome_id":
            contig_fields.append(k)

write_tsv(CONTIG_OUT, contig_rows, contig_fields)


# 4. consolidate mobtyper_results.txt

mobtyper_rows = []

for gid in genome_ids:
    path = MOB_OUT / gid / "mobtyper_results.txt"
    if not path.is_file():
        continue  # no-plasmid genome, nothing to add
    rows = read_tsv(path)

    for r in rows:
        new_row = dict(r)
        new_row["genome_id"] = gid
        mobtyper_rows.append(new_row)

mobtyper_fields = ["genome_id"]

if mobtyper_rows:
    for k in mobtyper_rows[0].keys():
        if k != "genome_id":
            mobtyper_fields.append(k)

write_tsv(MOBTYPER_OUT, mobtyper_rows, mobtyper_fields)


# 5. replicon presence/absence matrix

replicons_by_genome = defaultdict(set)
all_replicons = set()

for r in contig_rows:
    gid = r["genome_id"]

    for rep in split_multi(r.get("rep_type(s)", "")):
        replicons_by_genome[gid].add(rep)
        all_replicons.add(rep)

replicons = sorted(all_replicons)

with open(REPLICON_OUT, "w", newline="", encoding="utf-8") as f:
    fieldnames = ["genome_id"] + replicons
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()

    for gid in genome_ids:
        row = {"genome_id": gid}
        present = replicons_by_genome.get(gid, set())

        for rep in replicons:
            row[rep] = 1 if rep in present else 0

        w.writerow(row)


# 6. per-genome summary

# group rows by genome once, rather than rescanning every row for every genome
# (that was O(genomes x rows) and took 9+ min at 8,318 genomes)
rows_by_genome = defaultdict(list)
for r in contig_rows:
    rows_by_genome[r["genome_id"]].append(r)

summary_rows = []

for gid in genome_ids:
    rows = rows_by_genome.get(gid, [])

    primary_clusters = {
        cid.strip()
        for r in rows
        for cid in [r.get("primary_cluster_id", "")]
        if cid.strip() and cid.strip() not in {"-", "None"}
    }

    secondary_clusters = {
        cid.strip()
        for r in rows
        for cid in [r.get("secondary_cluster_id", "")]
        if cid.strip() and cid.strip() not in {"-", "None"}
    }

    circular = sum(
        1 for r in rows
        if r.get("circularity_status", "").lower() == "circular"
    )

    total_bp = 0
    for r in rows:
        try:
            total_bp += int(float(r.get("size", 0)))
        except (ValueError, TypeError):
            pass

    replicon_hits = sum(
        1 for r in rows
        if split_multi(r.get("rep_type(s)", ""))
    )

    relaxase_hits = sum(
        1 for r in rows
        if split_multi(r.get("relaxase_type(s)", ""))
    )

    mpf_hits = sum(
        1 for r in rows
        if split_multi(r.get("mpf_type", ""))
    )

    orit_hits = sum(
        1 for r in rows
        if split_multi(r.get("orit_type(s)", ""))
    )

    mobility_counts = Counter(
        r.get("predicted_mobility", "").strip()
        for r in rows
        if r.get("predicted_mobility", "").strip()
    )

    summary_rows.append({
        "genome_id": gid,
        "mob_contig_count": len(rows),
        "primary_cluster_count": len(primary_clusters),
        "secondary_cluster_count": len(secondary_clusters),
        "circular_contig_count": circular,
        "total_mob_contig_bp": total_bp,
        "replicon_hit_contig_count": replicon_hits,
        "relaxase_hit_contig_count": relaxase_hits,
        "mpf_hit_contig_count": mpf_hits,
        "orit_hit_contig_count": orit_hits,
        "predicted_mobility_counts": ";".join(
            f"{k}:{v}" for k, v in sorted(mobility_counts.items())
        )
    })

summary_fields = [
    "genome_id",
    "mob_contig_count",
    "primary_cluster_count",
    "secondary_cluster_count",
    "circular_contig_count",
    "total_mob_contig_bp",
    "replicon_hit_contig_count",
    "relaxase_hit_contig_count",
    "mpf_hit_contig_count",
    "orit_hit_contig_count",
    "predicted_mobility_counts",
]

with open(SUMMARY_OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=summary_fields)
    w.writeheader()
    w.writerows(summary_rows)


# 7. QC report

total_contigs = len(contig_rows)
genomes_with_contigs = sum(
    1 for gid in genome_ids if contig_counts[gid] > 0
)

qc_lines = [
    "Phase 10 — Code 1 MOB consolidation QC",
    "=======================================",
    f"MOB output directory: {MOB_OUT}",
    f"Expected genomes: {EXPECTED}",
    f"Genome directories found: {len(genome_ids)}",
    f"Missing contig_report.txt: {len(missing_contig)}",
    f"Missing mobtyper_results.txt (real errors, genome has plasmid contigs): {len(missing_mobtyper)}",
    f"No-plasmid genomes (no mobtyper_results.txt, expected/valid): {len(no_plasmid_genomes)}",
    f"Total consolidated contig rows: {total_contigs}",
    f"Genomes with >=1 contig report row: {genomes_with_contigs}",
    f"Unique replicon types: {len(replicons)}",
    f"Unique primary clusters: {len({r.get('primary_cluster_id','') for r in contig_rows if r.get('primary_cluster_id','') not in {'','-'}})}",
    f"Unique secondary clusters: {len({r.get('secondary_cluster_id','') for r in contig_rows if r.get('secondary_cluster_id','') not in {'','-'}})}",
    "",
    "OUTPUT FILES",
    "------------",
    str(CONTIG_OUT),
    str(MOBTYPER_OUT),
    str(REPLICON_OUT),
    str(SUMMARY_OUT),
    str(QC_OUT),
    "",
    "STATUS",
    "------",
    "GREEN" if (
        len(genome_ids) == EXPECTED
        and not missing_contig
        and not missing_mobtyper
    ) else "RED",
]

OUT.mkdir(parents=True, exist_ok=True)

with open(QC_OUT, "w", encoding="utf-8") as f:
    f.write("\n".join(qc_lines) + "\n")


print("\n".join(qc_lines))
print("\nCode 1 completed successfully.")