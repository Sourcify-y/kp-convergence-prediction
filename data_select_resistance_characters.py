#!/usr/bin/env python3
"""
data_select_resistance_characters.py

Applies the spec.md Section 4 (Model B) decision rule to the real data: for
each candidate Kleborate AMR column, compute prevalence among usable genomes
(optionally only the Mash-dedup representatives), so the resistance
character list comes from numbers instead of guesses.

Candidates are the AMR fields from data_download_and_genotype.py's FINAL_FIELDS:
  - 12 non-beta-lactam "_acquired" columns
  - 5 beta-lactam related columns (Bla_Carb_acquired is mandatory regardless
    of prevalence, per spec.md Section 3)
  - 4 mutation-based columns (chromosomal, probably not part of the
    "acquired" set, shown so that can be checked)

Usage:
    python3 data_select_resistance_characters.py --csv kp_metadata_FULL.csv \\
        [--manifest ... --genome-list ...] [--threshold 0.01]
"""

import argparse
import csv
from pathlib import Path

import pandas as pd

NON_BLA_ACQUIRED = [
    "AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired",
    "Gly_acquired", "MLS_acquired", "Phe_acquired", "Rif_acquired",
    "Sul_acquired", "Tet_acquired", "Tgc_acquired", "Tmt_acquired",
]

BLA_RELATED = [
    "Bla_acquired", "Bla_inhR_acquired", "Bla_ESBL_acquired",
    "Bla_ESBL_inhR_acquired", "Bla_Carb_acquired",
]

MUTATION_BASED = [
    "SHV_mutations", "Omp_mutations", "Col_mutations", "Flq_mutations",
]

MANDATORY = {"Bla_Carb_acquired"}  # spec.md Section 3, keep regardless of prevalence


def load_rep_accessions(manifest_path: Path, genome_list_path: Path) -> set:
    path2acc = {}
    with open(manifest_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            p = (row.get("path") or "").strip().rstrip("\r")
            acc = (row.get("assembly_accession") or "").strip()
            if p and acc:
                path2acc[p] = acc
    kept = set()
    with open(genome_list_path) as f:
        for line in f:
            p = line.strip().rstrip("\r")
            if p and p in path2acc:
                kept.add(path2acc[p])
    return kept


def _present(val) -> int:
    s = str(val).strip()
    return 0 if s in ("", "-", "nan", "0") else 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--genome-list", type=Path, default=None)
    parser.add_argument("--threshold", type=float, default=0.01,
                        help="Prevalence threshold from spec.md Section 4 "
                             "(candidate default 1%%)")
    args = parser.parse_args()

    df = pd.read_csv(args.csv, dtype=str)
    df["collection_year_num"] = pd.to_numeric(df["collection_year"], errors="coerce")
    df["has_st"] = df["ST"].notna() & (df["ST"].astype(str).str.strip() != "")
    df["usable"] = df["collection_year_num"].notna() & df["has_st"]

    print(f"[load] {len(df)} total rows, {df['usable'].sum()} usable")

    if args.manifest and args.genome_list:
        kept = load_rep_accessions(args.manifest, args.genome_list)
        before = len(df)
        df = df[df["genome_id"].isin(kept)].copy()
        print(f"[dedup-filter] {before} -> {len(df)} rows")

    df = df[df["usable"]].copy()
    n = len(df)
    print(f"[usable] {n} genomes\n")

    def report(section_name, cols):
        print(f"--- {section_name} ---")
        for col in cols:
            if col not in df.columns:
                print(f"  {col:25s}  MISSING FROM CSV")
                continue
            present = df[col].apply(_present).sum()
            frac = present / n if n else 0
            flag = ""
            if col in MANDATORY:
                flag = "  [MANDATORY — Section 3 target trait, keep regardless]"
            elif frac < args.threshold:
                flag = f"  [BELOW {args.threshold:.1%} THRESHOLD]"
            print(f"  {col:25s}  {present:6d} / {n}  ({frac:.2%}){flag}")
        print()

    report("Non-beta-lactam acquired classes (12 candidates)", NON_BLA_ACQUIRED)
    report("Beta-lactam / carbapenem related (5 candidates)", BLA_RELATED)
    report("Mutation-based mechanisms (4 candidates, likely NOT part of "
           "the 'acquired' character set — shown for completeness)", MUTATION_BASED)

    print("Next: pick the final resistance character list as whichever "
          "columns clear the threshold, keeping Bla_Carb_acquired regardless, "
          "then edit RES_COLS in data_build_character_matrix.py to match.")


if __name__ == "__main__":
    main()