import os
import pandas as pd
import sys
from pathlib import Path

BASE_DIR = Path(os.environ.get("KP_BASE_DIR", "."))  # run from the repo root, or set KP_BASE_DIR
HYPERTRAPS_DATA = BASE_DIR / "hypertraps-ct" / "Data"

RES_CHARS = [
    "Bla_ESBL_acquired","Bla_Carb_acquired","Flq_acquired","Sul_acquired","AGly_acquired",
    "Col_acquired","Tet_acquired","MLS_acquired","Rif_acquired","Tmt_acquired",
    "Fcyn_acquired","Phe_acquired","Bla_acquired","Tgc_acquired"
]
VIR_CHARS = ["Yersiniabactin","Aerobactin","Salmochelin","Colibactin","rmpA","rmpA2"]

# 14 resistance characters + one collapsed "any virulence locus" column = 15
COARSE_CHAR_ORDER = RES_CHARS + ["Hypervirulence_any"]

df = pd.read_csv(BASE_DIR / "phase3_out" / "character_matrix.csv")

needed = RES_CHARS + VIR_CHARS + ["collection_year"]
missing_cols = [c for c in needed if c not in df.columns]
if missing_cols:
    print(f"FATAL: missing columns: {missing_cols}")
    print(f"Actual columns: {list(df.columns)}")
    sys.exit(1)

train = df[df["collection_year"] <= 2017].copy()
print(f"TRAIN: {len(train)} genomes (collection_year<=2017)")

train["Hypervirulence_any"] = (train[VIR_CHARS].sum(axis=1) > 0).astype(int)
prevalence = train["Hypervirulence_any"].mean() * 100
print(f"Hypervirulence_any prevalence in TRAIN: {prevalence:.1f}%")

HYPERTRAPS_DATA.mkdir(parents=True, exist_ok=True)
train[COARSE_CHAR_ORDER].to_csv(
    HYPERTRAPS_DATA / "kp-train-data-15char-coarse.txt", sep=" ", header=False, index=False
)
with open(HYPERTRAPS_DATA / "kp-train-labels-15char-coarse.txt", "w") as f:
    f.write("\n".join(COARSE_CHAR_ORDER) + "\n")

print(f"Wrote kp-train-data-15char-coarse.txt / kp-train-labels-15char-coarse.txt to {HYPERTRAPS_DATA}")
print(f"Column order (0-indexed, matches hypertraps.ce output): {list(enumerate(COARSE_CHAR_ORDER))}")