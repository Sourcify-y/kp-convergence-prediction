# kp-convergence-prediction

Code and derived data for a study of resistance and virulence convergence in 8,318 deduplicated
*Klebsiella pneumoniae* genomes. A HyperTraPS-CT model trained on genomes collected up to 2017 was
evaluated on 2018 to 2020 and 2021 to 2022. In the 2021 to 2022 window, AUROC was 0.686 (95% CI 0.559 to 0.813;
within-starting-state permutation p = 0.042).

## Contents

| Path | Description |
|---|---|
| `accessions_8318.txt` | NCBI assembly accessions of the 8,318 genomes |
| `data/phase4/kp_dedup_nj_tree.nwk` | Neighbour-joining tree from Mash distances, midpoint-rooted, 8,318 tips |
| `data/phase4/candidate_acquisition_events_flagged_v2.csv` | 7,393 acquisition events from Fitch parsimony, with exclusion flags |
| `scripts/` | Analysis scripts, one folder per phase |
| `figures/` | Figures |
| `requirements.txt` | Python packages |
| `software_versions.txt` | R packages and external tools with versions |
| `SHA256SUMS` | Checksums of repository files |

## Event table columns

`character`, `ancestral_node`, `descendant_node`, `descendant_is_tip`, `descendant_genome_id`,
`ancestral_state_confident`, `edge_length`, `collection_year`, `subtree_size`, `is_singleton`,
`is_reacquisition`, `subtree_reverts`.

Node numbers follow the ape and phangorn convention: tips are 1 to 8318 in file order, and internal nodes follow in preorder from the root.

## Tree figure

```
pip install -r requirements.txt
python3 scripts/phase4/make_tree_figure.py
```

Output: `figures/Phase4_NJ_tree_acquisition_events.png`. Markers show acquisition events at the descendant node of each event: carbapenem resistance (668), virulence-associated loci (881), other resistance classes (5,844). The script checks the tree against the event table before drawing.

## Data

Genome assemblies are available from NCBI under the accessions in `accessions_8318.txt`. They were retrieved in September 2026.

## License

MIT
