#!/usr/bin/env python3
"""
Phase 9: lineage-holdout TRAIN input.

Same filtering as eval_build_train_input.py (collection_year <= 2017), plus one
addition: every genome in a lineage listed in phase9_out/lineage_holdout_assignment.csv
is dropped from training entirely, not just from the scored task tables. That's what
makes it a real lineage holdout (the model never sees these lineages) instead of
re-scoring the same model on a different subset.

Expects genome_id, collection_year, base_ST and the 20 character columns in
character_matrix.csv.
"""

import pandas as pd
import sys

CHAR_MATRIX_PATH = "phase3_out/character_matrix.csv"
HOLDOUT_PATH = "phase9_out/lineage_holdout_assignment.csv"
OUT_DATA_PATH = "hypertraps-ct/Data/kp-holdout-train-data-20char.txt"
OUT_LABELS_PATH = "hypertraps-ct/Data/kp-holdout-train-labels-20char.txt"

# Column order of character_matrix.csv (its header order), not the Phase 7/8 order. The holdout obs
# files are built in this order, so scoring must use it too: routes files encode characters by 0-indexed position.
# Checked against the real character_matrix.csv header (2026-09, via `head -1`).
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

    # TRAIN window (same as Phase 8)
    df = df[df["collection_year"] <= 2017].copy()
    n_after_year = len(df)

    # new for Phase 9: drop any genome whose lineage is held out, whatever its year
    held_out_lineages = set(holdout["base_ST"])
    df = df[~df["base_ST"].isin(held_out_lineages)].copy()
    n_after_holdout = len(df)

    print(f"Genomes total:                    {n_before}")
    print(f"After TRAIN-window filter (<=2017): {n_after_year}")
    print(f"After lineage-holdout exclusion:    {n_after_holdout} "
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