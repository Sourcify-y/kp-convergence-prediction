#!/usr/bin/env python3
"""
data_build_character_matrix.py

Phase 3: build the genome x character presence/absence matrix that
HyperTraPS-CT runs on. One row per isolate (non-independence is handled
upstream by the Mash dedup and downstream by the Phase 9 lineage holdout).

Reads kp_metadata_FULL.csv and, optionally, the --manifest/--genome-list pair
to restrict to the Mash-dedup representative genomes. Writes:
    character_matrix.csv          metadata cols + 20 0/1 character cols
    character_matrix_summary.txt  prevalence of each character

Column order matters for HyperTraPS-CT, don't reorder between runs.
14 resistance + 6 virulence = 20 characters.

Resistance set locked 2026-09 (spec.md Section 4, applied to the 9,108-genome
dedup set): 11 non-beta-lactam classes over 1% prevalence, plus Bla_acquired
and Bla_ESBL_acquired (both over 1%) and Bla_Carb_acquired (kept regardless,
it's the definition of the primary resistance trait). Gly_acquired (0%),
Bla_inhR_acquired and Bla_ESBL_inhR_acquired (<0.5%) are excluded.

Usage:
    python3 data_build_character_matrix.py --csv kp_metadata_FULL.csv \\
        --outdir ./phase3_out \\
        [--manifest 03_download_manifest.csv --genome-list representative_genomes.txt]
"""

import argparse
import csv
import re
from pathlib import Path

import pandas as pd

RES_COLS = [
    "AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired",
    "MLS_acquired", "Phe_acquired", "Rif_acquired", "Sul_acquired",
    "Tet_acquired", "Tgc_acquired", "Tmt_acquired",
    "Bla_acquired", "Bla_ESBL_acquired", "Bla_Carb_acquired",
]

VIR_COLS = [
    "Yersiniabactin", "Aerobactin", "Salmochelin",
    "Colibactin", "rmpA", "rmpA2",
]

CHAR_COLS = RES_COLS + VIR_COLS

META_COLS = ["genome_id", "collection_year", "base_ST", "country"]


def load_rep_accessions(manifest_path: Path, genome_list_path: Path) -> set:
    # map genome file path -> assembly accession using the download manifest
    path2acc = {}
    with open(manifest_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            p = (row.get("path") or "").strip().rstrip("\r")
            acc = (row.get("assembly_accession") or "").strip()
            if p and acc:
                path2acc[p] = acc

    kept = set()
    n_unmatched = 0
    with open(genome_list_path) as f:
        for line in f:
            p = line.strip().rstrip("\r")
            if not p:
                continue
            acc = path2acc.get(p)
            if acc:
                kept.add(acc)
            else:
                n_unmatched += 1

    if n_unmatched:
        print(f"[dedup-filter] WARNING: {n_unmatched} paths in {genome_list_path} "
              f"had no matching row in {manifest_path} — check for a path-format "
              f"mismatch before trusting the filtered matrix.")
    print(f"[dedup-filter] Resolved {len(kept)} representative "
          f"accessions from {genome_list_path}")
    return kept


def _present(val) -> int:
    # Kleborate: '-' or blank = absent, anything else (a gene/allele call) = present
    s = str(val).strip()
    return 0 if s in ("", "-", "nan", "0") else 1


def build_matrix(df: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in CHAR_COLS if c not in df.columns]
    if missing:
        raise SystemExit(
            f"[error] These locked character columns are not present in the "
            f"input CSV: {missing}\nCheck RES_COLS at the top of "
            f"this script against the actual column names in kp_metadata_FULL.csv."
        )

    mat = df[META_COLS].copy()
    for col in CHAR_COLS:
        mat[col] = df[col].apply(_present)
    return mat


def write_summary(matrix: pd.DataFrame, outdir: Path):
    lines = []
    lines.append(f"Genome x character matrix — {len(matrix)} genomes, "
                  f"{len(CHAR_COLS)} characters "
                  f"({len(RES_COLS)} resistance + "
                  f"{len(VIR_COLS)} virulence)\n")
    lines.append("Character prevalence (fraction of genomes positive):")
    for col in CHAR_COLS:
        frac = matrix[col].mean()
        n = int(matrix[col].sum())
        tag = "resistance" if col in RES_COLS else "virulence"
        lines.append(f"  [{tag:10s}] {col:25s} {n:6d} / {len(matrix)}  ({frac:.1%})")
    lines.append("")
    lines.append("Sanity check: compare these prevalence numbers against Step "
                  "5's original character-set prevalence figures before "
                  "trusting this matrix — a big divergence likely means the "
                  "RES_COLS list doesn't actually match the "
                  "locked spec.md definition.")
    text = "\n".join(lines)
    (outdir / "character_matrix_summary.txt").write_text(text)
    print(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True, type=Path,
                        help="Path to kp_metadata_FULL.csv")
    parser.add_argument("--outdir", required=True, type=Path,
                        help="Directory to write the matrix + summary into")
    parser.add_argument("--manifest", type=Path, default=None,
                        help="Path to 03_download_manifest.csv. Required if "
                             "--genome-list is given.")
    parser.add_argument("--genome-list", type=Path, default=None,
                        help="Path to the Mash-dedup representative genome list "
                             "(e.g. mash_test/representative_genomes.txt). When "
                             "given (with --manifest), restricts the matrix to "
                             "just these genomes.")
    parser.add_argument("--usable-only", action="store_true", default=True,
                        help="Restrict to genomes with both a parsed collection "
                             "year and an ST call (default: on — HyperTraPS "
                             "needs a usable year for the temporal model).")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if bool(args.manifest) != bool(args.genome_list):
        parser.error("--manifest and --genome-list must be given together.")

    df = pd.read_csv(args.csv, dtype=str)
    df["collection_year_num"] = pd.to_numeric(df["collection_year"], errors="coerce")
    df["has_st"] = df["ST"].notna() & (df["ST"].astype(str).str.strip() != "")
    df["usable"] = df["collection_year_num"].notna() & df["has_st"]
    # drop the single-locus-variant suffix, e.g. "258-1LV" -> "258"
    df["base_ST"] = df["ST"].apply(
        lambda st: re.sub(r"-\d+LV$", "", str(st)) if pd.notna(st) else st
    )
    print(f"[load] {len(df)} total rows, {df['usable'].sum()} usable "
          f"(dated + genotyped) genomes")

    if args.genome_list:
        kept = load_rep_accessions(args.manifest, args.genome_list)
        before = len(df)
        df = df[df["genome_id"].isin(kept)].copy()
        print(f"[dedup-filter] {before} -> {len(df)} rows "
              f"({df['usable'].sum()} usable after filtering)")

    if args.usable_only:
        before = len(df)
        df = df[df["usable"]].copy()
        print(f"[usable-filter] {before} -> {len(df)} rows "
              f"(dropped genomes missing a year or ST call)")

    # HyperTraPS wants a clean int year
    df["collection_year"] = df["collection_year_num"].astype(int)

    matrix = build_matrix(df)

    out_path = args.outdir / "character_matrix.csv"
    matrix.to_csv(out_path, index=False)
    print(f"\n[write] {len(matrix)} genomes x {len(CHAR_COLS)} "
          f"characters -> {out_path}")

    write_summary(matrix, outdir=args.outdir)


if __name__ == "__main__":
    main()