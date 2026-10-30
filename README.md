# Klebsiella pneumoniae resistance and virulence convergence

Code for prospective genomic prediction of resistance and virulence convergence in 8,318 *Klebsiella pneumoniae* genomes: Mash deduplication, HyperTraPS evolutionary accumulation modelling, temporal validation, robustness analyses, and sequence-type-adjusted plasmid-context analysis.

## Requirements

Python 3.10+, R.

```bash
pip install -r requirements.txt
```

R packages: `pROC`, `data.table`.

| Tool | Version |
|---|---|
| NCBI datasets CLI | |
| Mash | 2.3 |
| HyperTraPS-CT | |
| abricate | 1.4.0 |
| CARD database | 2026-04-03 |
| PlasmidFinder database | |
| MOB-suite | 3.1.9 |

## Pipeline

**1. Download.** Assemblies and metadata from NCBI.
Output: `kp_pull_20k/kp_metadata_FULL.csv`, `kp_pull_20k/03_download_manifest.csv`.

**2. Mash deduplication.** 15,764 genomes reduced to 9,108 representatives (single linkage, distance 0.0003).

```bash
tail -n +2 kp_pull_20k/03_download_manifest.csv | cut -d',' -f2 | tr -d '\r' \
  | sed "s|^|$(pwd)/|" > mash_test/paths_all.txt
mash sketch -o mash_test/sketch_full -l mash_test/paths_all.txt -p 8
mash triangle mash_test/sketch_full.msh > mash_test/distances_full.tsv
python data_mash_dedup.py \
  --triangle mash_test/distances_full.tsv \
  --threshold 0.0003 \
  --out mash_test/representative_genomes.txt
```

**3. Deduplicated metadata.**

```python
import pandas as pd
df = pd.read_csv('kp_pull_20k/kp_metadata_FULL.csv')
manifest = pd.read_csv('kp_pull_20k/03_download_manifest.csv')
manifest['path'] = manifest['path'].str.strip()
with open('mash_test/representative_genomes.txt') as f:
    rep_paths = set(line.strip() for line in f)
rep_acc = set(manifest[manifest['path'].isin(rep_paths)]['assembly_accession'])
df[df['genome_id'].isin(rep_acc)].to_csv('kp_pull_20k/kp_metadata_FULL_dedup.csv', index=False)
```

**4. Character matrix.** 8,318 genomes, 20 characters (14 resistance classes, 6 virulence loci).

```bash
python data_build_character_matrix.py \
  --csv kp_pull_20k/kp_metadata_FULL_dedup.csv \
  --manifest kp_pull_20k/03_download_manifest.csv \
  --genome-list mash_test/representative_genomes.txt \
  --outdir phase3_out
```

**5. Tree and acquisition events.** Neighbour-joining tree from Mash distances, Fitch parsimony ancestral state reconstruction. `make_tree_figure.py`.

**6. Temporal split and baselines.** TRAIN ≤2017, VALIDATE 2018 to 2020, TEST 2021 to 2022. `data_temporal_split_analysis.py`, `phase6_baselines.py`.

**7. HyperTraPS.** Flags: `--length 4 --kernel 5 --model 1`, 3 seeds.
- Full data: `kp-full-data-20char.txt`, labels `kp_full_run_model1_s{0,1,2}`.
- Train only: `kp-train-data-20char.txt` from `eval_build_train_input.py`.
- `build_20char_full.py`, `compare_seeds.py`, `recompute_bubbles*.py`.

**8. Robustness.** Lineage holdout, geographic holdout, subsampling (`phase9_subsample_out/`), coarse 15-character model (`phase9_coarse_out/`), 3 seeds each.

**9. Genomic context.** abricate (CARD, PlasmidFinder), MOB-suite, sequence-type-adjusted association tests.

**10. Figures and tables.** `feasibility_plots.py`, `compute_resistance_prevalence.py`.

## Browser-based demonstration

This repository includes a small executable notebook demonstrating data validation and chronological evaluation concepts related to the research project.

Open the notebook in Google Colab: https://colab.research.google.com/github/Sourcify-y/kp-convergence-prediction/blob/main/demo/kp_convergence_temporal_demo.ipynb

**Run the demo**

1. Open the notebook using the Colab link above.

2. Run the first cell.

3. Run the remaining cells in order, or select Runtime → Run all.

**Scope and limitations**

The demonstration uses a small, entirely synthetic dataset. Its illustrative labels and calculated metrics are not research findings and do not represent measured outcomes from Klebsiella pneumoniae. The notebook demonstrates data checks, a chronological train/test split, and baseline model evaluation. It does not reproduce the full 8,318-genome analysis or run HyperTraPS-CT and the complete bioinformatics pipeline.

See demo/README_DEMO.md for full instructions and limitations.

## Data availability

Genome assemblies are available from NCBI. Accessions are listed in the Supplementary Materials.

## License

MIT.
