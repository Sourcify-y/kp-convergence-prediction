#!/usr/bin/env python3
"""
Phase 9: geographic-holdout TRAIN input.
Same logic as robustness_lineage_holdout_train_input.py, but excludes the lineages in the
geo holdout assignment instead of the lineage holdout assignment.
"""

import pandas as pd
import sys

CHAR_MATRIX_PATH = "phase3_out/character_matrix.csv"
HOLDOUT_PATH = "phase9_out/geo_holdout_assignment.csv"
OUT_DATA_PATH = "hypertraps-ct/Data/kp-geo-holdout-train-data-20char.txt"
OUT_LABELS_PATH = "hypertraps-ct/Data/kp-geo-holdout-train-labels-20char.txt"

CHARACTER_COLUMNS = [
    "AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired",
    "MLS_acquired", "Phe_acquired", "Rif_acquired", "Sul_acquired",
    "Tet_acquired", "Tgc_acquired", "Tmt_acquired", "Bla_acquired",
    "Bla_ESBL_acquired", "Bla_Carb_acquired", "Yersiniabactin",
    "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2",
]


def main():
    df = pd.read_csv(CHAR_MATRIX_PATH)
    holdout = pd.read_csv(HOLDOUT_PATH)

    required_cols = {"genome_id", "collection_year", "base_ST"} | set(CHARACTER_COLUMNS)
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        sys.exit(f"ERROR: character_matrix.csv missing expected columns: {missing_cols}")

    n_before = len(df)
    df = df[df["collection_year"] <= 2017].copy()
    n_after_year = len(df)

    held_out_lineages = set(holdout["base_ST"])
    df = df[~df["base_ST"].isin(held_out_lineages)].copy()
    n_after_holdout = len(df)

    print(f"Genomes total:                    {n_before}")
    print(f"After TRAIN-window filter (<=2017): {n_after_year}")
    print(f"After geo-holdout exclusion:        {n_after_holdout} "
          f"({n_after_year - n_after_holdout} genomes removed, "
          f"{len(held_out_lineages)} lineages excluded)")

    if df["genome_id"].duplicated().any():
        sys.exit("ERROR: duplicate genome_id after filtering — investigate before writing model input.")

    df[CHARACTER_COLUMNS].to_csv(OUT_DATA_PATH, sep=" ", header=False, index=False)
    df[["genome_id"]].to_csv(OUT_LABELS_PATH, sep=" ", header=False, index=False)

    print(f"\nWrote {OUT_DATA_PATH} ({n_after_holdout} rows)")
    print(f"Wrote {OUT_LABELS_PATH}")


if __name__ == "__main__":
    main()