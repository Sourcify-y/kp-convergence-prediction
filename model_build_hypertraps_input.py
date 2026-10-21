import pandas as pd

df = pd.read_csv('phase3_out/character_matrix.csv')

print(f"Full dataset: {len(df)} genomes")

# same column order as the pilot so runs are comparable
base_18 = [
    'Bla_ESBL_acquired', 'Bla_Carb_acquired', 'Flq_acquired', 'Sul_acquired',
    'AGly_acquired', 'Col_acquired', 'Tet_acquired', 'MLS_acquired',
    'Rif_acquired', 'Tmt_acquired', 'Fcyn_acquired', 'Phe_acquired',
    'Yersiniabactin', 'Aerobactin', 'Salmochelin', 'Colibactin',
    'rmpA', 'rmpA2'  # these two replace the old 'RmpADC' label (unverified)
]
extra_2 = ['Bla_acquired', 'Tgc_acquired']
CHAR_ORDER = base_18 + extra_2

missing_cols = [c for c in CHAR_ORDER if c not in df.columns]
if missing_cols:
    raise ValueError(f"Missing expected columns: {missing_cols}")

before = len(df)
df = df.dropna(subset=['collection_year'])
print(f"Dropped {before - len(df)} rows missing collection_year")

df[CHAR_ORDER].astype(int).to_csv(
    'hypertraps_data/kp-full-data-20char.txt',
    sep=' ', header=False, index=False
)

with open('hypertraps_data/kp-full-labels-20char.txt', 'w') as f:
    for c in CHAR_ORDER:
        f.write(c + '\n')

df['collection_year'].astype(int).to_csv(
    'hypertraps_data/kp-full-datatime-20char.txt',
    header=False, index=False
)

print(f"Wrote {len(df)} rows x {len(CHAR_ORDER)} characters")