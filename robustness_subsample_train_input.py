import os
import pandas as pd
import numpy as np
import sys
from pathlib import Path

BASE_DIR = Path(os.environ.get("KP_BASE_DIR", "."))  # run from the repo root, or set KP_BASE_DIR
HYPERTRAPS_DATA = BASE_DIR / "hypertraps-ct" / "Data"

# Identical to eval_build_train_input.py — do not diverge from this list.
CHAR_ORDER = [
    "Bla_ESBL_acquired","Bla_Carb_acquired","Flq_acquired","Sul_acquired","AGly_acquired",
    "Col_acquired","Tet_acquired","MLS_acquired","Rif_acquired","Tmt_acquired",
    "Fcyn_acquired","Phe_acquired","Yersiniabactin","Aerobactin","Salmochelin",
    "Colibactin","rmpA","rmpA2","Bla_acquired","Tgc_acquired"
]

df = pd.read_csv(BASE_DIR / "phase3_out" / "character_matrix.csv")

needed = CHAR_ORDER + ["collection_year"]
missing_cols = [c for c in needed if c not in df.columns]
if missing_cols:
    print(f"FATAL: missing columns: {missing_cols}")
    print(f"Actual columns: {list(df.columns)}")
    sys.exit(1)

train = df[df["collection_year"] <= 2017].copy()
print(f"Full TRAIN population: {len(train)} genomes (collection_year<=2017)")

HYPERTRAPS_DATA.mkdir(parents=True, exist_ok=True)

REPLICATE_SEEDS = [101, 102, 103]  # fixed, unrelated to the model's own --seed
SUBSAMPLE_FRAC = 0.80

for i, seed in enumerate(REPLICATE_SEEDS, start=1):
    sub = train.sample(frac=SUBSAMPLE_FRAC, random_state=seed)
    print(f"Replicate {i} (pandas seed {seed}): {len(sub)} genomes "
          f"({len(sub)/len(train)*100:.1f}% of TRAIN)")
    out_path = HYPERTRAPS_DATA / f"kp-train-data-20char-sub{i}.txt"
    sub[CHAR_ORDER].to_csv(out_path, sep=" ", header=False, index=False)
    print(f"  Wrote {out_path}")

print("\nReusing existing kp-train-labels-20char.txt (character set unchanged) — no new labels file needed.")