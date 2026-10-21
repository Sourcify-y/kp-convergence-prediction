#!/usr/bin/env python3
"""
robustness_assembly_qc.py -- per-genome assembly QC computed directly from the .fna files
(the merged Kleborate table has no contig-stats or species columns).

Writes: genome_id, assembly_accession, contig_count, total_size, largest_contig, N50, ambiguous_bases
(ambiguous_bases = number of non-ACGT letters, i.e. Ns and IUPAC codes).

Time it first:  python3 robustness_assembly_qc.py --manifest <03_download_manifest.csv> --limit 200
Full run:       python3 robustness_assembly_qc.py --manifest <03_download_manifest.csv> --matrix phase3_out/character_matrix.csv
"""
import argparse
import gzip
import re
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

ACC_RE = re.compile(r"GC[AF]_\d{9}\.\d+")


def qc_one(path):
    try:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rb") as f:
            data = f.read()
        if not data.startswith(b">"):
            return path, None
        lengths, n_ambig = [], 0
        for rec in data[1:].split(b"\n>"):
            nl = rec.find(b"\n")
            seq = rec[nl + 1:] if nl >= 0 else b""
            seq = seq.replace(b"\r", b"").replace(b"\n", b"")
            lengths.append(len(seq))
            n_ambig += len(seq.translate(None, b"ACGTacgt"))
        total = sum(lengths)
        half, running, n50 = total / 2, 0, 0
        for contig_len in sorted(lengths, reverse=True):
            running += contig_len
            if running >= half:
                n50 = contig_len
                break
        return path, dict(contig_count=len(lengths), total_size=total, largest_contig=max(lengths), N50=n50, ambiguous_bases=n_ambig)
    except Exception as e:
        return path, None


def read_manifest(path):
    # first .fna/.fa/.fasta(.gz) field on each line
    paths = []
    for line in Path(path).read_text().splitlines():
        line = line.strip().replace("\r", "")
        if not line:
            continue
        for field in line.split(","):
            field = field.strip().strip('"')
            if re.search(r"\.(fna|fa|fasta)(\.gz)?$", field):
                paths.append(field)
                break
    return paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--matrix", default=None, help="restrict to genome_ids in this character_matrix.csv")
    parser.add_argument("--out", default="assembly_qc.csv")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int, default=None, help="only the first N genomes (timing test)")
    args = parser.parse_args()

    paths = read_manifest(args.manifest)
    if not paths:
        sys.exit(f"no .fna/.fa/.fasta paths found in {args.manifest} -- check the file (first lines: {Path(args.manifest).read_text()[:200]!r})")
    acc_of = {p: (ACC_RE.search(p).group(0) if ACC_RE.search(p) else None) for p in paths}
    if args.matrix:
        ids = set(pd.read_csv(args.matrix, usecols=["genome_id"], dtype=str)["genome_id"])
        paths = [p for p in paths if acc_of[p] in ids]
        print(f"restricted to {len(paths)} genomes present in {args.matrix} (matrix has {len(ids)})")
        if len(paths) < 0.95 * len(ids):
            sys.exit("fewer than 95% of matrix genomes were found in the manifest -- wrong manifest or genome_id != accession")
    if args.limit:
        paths = paths[: args.limit]
    missing = [p for p in paths if not Path(p).exists()]
    print(f"{len(paths)} genomes to process; {len(missing)} paths do not exist on disk")
    if len(missing) > 0.05 * len(paths):
        sys.exit(f"too many missing files, e.g. {missing[:3]}")
    paths = [p for p in paths if p not in set(missing)]
    t0, rows, n_fail = time.time(), [], 0
    with Pool(args.workers) as pool:
        for i, (p, res) in enumerate(pool.imap_unordered(qc_one, paths, chunksize=4), 1):
            if res is None:
                n_fail += 1
                continue
            res.update(genome_id=acc_of[p], assembly_accession=acc_of[p])
            rows.append(res)
            if i % 500 == 0:
                elapsed = time.time() - t0
                print(f"  {i}/{len(paths)}  {elapsed:.0f}s elapsed, ETA {elapsed / i * (len(paths) - i):.0f}s", flush=True)
    df = pd.DataFrame(rows)[["genome_id", "assembly_accession", "contig_count", "total_size", "largest_contig", "N50", "ambiguous_bases"]]
    df.to_csv(args.out, index=False)
    print(f"wrote {args.out}: {len(df)} genomes ({n_fail} unreadable) in {time.time() - t0:.0f}s")
    print(df[["contig_count", "total_size", "N50", "ambiguous_bases"]].describe().round(0).to_string())


if __name__ == "__main__":
    main()