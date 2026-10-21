#!/usr/bin/env python3
"""
Phase 10 Track B fusion: build one genome-keyed master feature table from
  - character_matrix.csv                   Kleborate AMR + virulence (Phase 3), the anchor row set
  - card_summary_matrix.tsv                Abricate/CARD resistance genes
  - plasmidfinder_summary_matrix.tsv       Abricate/PlasmidFinder replicons
  - mob_consolidated/plasmid_summary.csv   MOB-suite genome-level plasmid summary

replicon_matrix.csv (MOB-suite contig-level replicon combinations, 532 cols) is
validated but not merged: its columns are combination strings rather than
independent features, so it stays a separate file, joinable on genome_id.

Edit the paths below to match your machine before running.
"""

import pandas as pd
import sys
import os

# ---- EDIT THESE PATHS ----
CHARACTER_MATRIX     = "phase3_out/character_matrix.csv"
CARD_MATRIX          = "phase10/abricate_out/card_summary_matrix.tsv"
PLASMIDFINDER_MATRIX = "phase10/abricate_out/plasmidfinder_summary_matrix.tsv"
MOB_PLASMID_SUMMARY  = "phase10/mob_consolidated/plasmid_summary.csv"
MOB_REPLICON_MATRIX  = "phase10/mob_consolidated/replicon_matrix.csv"  # validated only, not merged

OUTPUT_MASTER = "phase10/phase10_master_table_v2.csv"
OUTPUT_QC_LOG = "phase10/phase10_fusion_qc_v2.txt"
# ---------------------------


def clean_ids(paths):
    # '/some/dir/GCF_xxx.tab' -> 'GCF_xxx'
    return (
        paths
        .apply(os.path.basename)
        .str.replace(r"\.(tab|tsv)$", "", regex=True)
    )


def main():
    qc_lines = []

    def log(msg):
        print(msg)
        qc_lines.append(msg)

    log("=" * 60)
    log("Phase 10 Track B fusion — QC log")
    log("=" * 60)

    # anchor: character_matrix.csv
    if not os.path.exists(CHARACTER_MATRIX):
        sys.exit(f"ERROR: character_matrix.csv not found at {CHARACTER_MATRIX} — edit the path at top of script.")
    char_df = pd.read_csv(CHARACTER_MATRIX, dtype={"genome_id": str})
    if "genome_id" not in char_df.columns:
        sys.exit(f"ERROR: 'genome_id' column not found in character_matrix.csv. Columns present: {list(char_df.columns)}")
    n_anchor = char_df["genome_id"].nunique()
    log(f"Anchor (character_matrix.csv): {len(char_df)} rows, {n_anchor} unique genome_id")
    if len(char_df) != n_anchor:
        log(f"  WARNING: {len(char_df) - n_anchor} duplicate genome_id rows in character_matrix.csv")

    master = char_df.copy()

    # CARD
    card_df = pd.read_csv(CARD_MATRIX, sep="\t")
    card_df["genome_id"] = clean_ids(card_df["#FILE"])
    card_df = card_df.drop(columns=["#FILE"])
    n_dupe_card = card_df["genome_id"].duplicated().sum()
    if n_dupe_card:
        log(f"  WARNING: {n_dupe_card} duplicate genome_id in CARD matrix after normalization")
    card_df = card_df.add_prefix("CARD_").rename(columns={"CARD_genome_id": "genome_id"})

    n_before = len(master)
    master = master.merge(card_df, on="genome_id", how="left", indicator=True)
    n_unmatched = (master["_merge"] == "left_only").sum()
    log(f"CARD merge: {n_before} -> {len(master)} rows; {n_unmatched} anchor genomes had NO CARD match (should be 0)")
    master = master.drop(columns=["_merge"])

    # PlasmidFinder
    pf_df = pd.read_csv(PLASMIDFINDER_MATRIX, sep="\t")
    pf_df["genome_id"] = clean_ids(pf_df["#FILE"])
    pf_df = pf_df.drop(columns=["#FILE"])
    n_dupe_pf = pf_df["genome_id"].duplicated().sum()
    if n_dupe_pf:
        log(f"  WARNING: {n_dupe_pf} duplicate genome_id in PlasmidFinder matrix after normalization")
    pf_df = pf_df.add_prefix("PF_").rename(columns={"PF_genome_id": "genome_id"})

    n_before = len(master)
    master = master.merge(pf_df, on="genome_id", how="left", indicator=True)
    n_unmatched = (master["_merge"] == "left_only").sum()
    log(f"PlasmidFinder merge: {n_before} -> {len(master)} rows; {n_unmatched} anchor genomes had NO PlasmidFinder match (should be 0)")
    master = master.drop(columns=["_merge"])

    # MOB-suite plasmid_summary
    mob_df = pd.read_csv(MOB_PLASMID_SUMMARY, dtype={"genome_id": str})
    if "genome_id" not in mob_df.columns:
        sys.exit(f"ERROR: 'genome_id' column not found in plasmid_summary.csv. Columns present: {list(mob_df.columns)}")
    n_dupe_mob = mob_df["genome_id"].duplicated().sum()
    if n_dupe_mob:
        log(f"  WARNING: {n_dupe_mob} duplicate genome_id in MOB-suite plasmid_summary after normalization")
    mob_df = mob_df.add_prefix("MOB_").rename(columns={"MOB_genome_id": "genome_id"})

    n_before = len(master)
    master = master.merge(mob_df, on="genome_id", how="left", indicator=True)
    n_unmatched = (master["_merge"] == "left_only").sum()
    log(f"MOB-suite plasmid_summary merge: {n_before} -> {len(master)} rows; {n_unmatched} anchor genomes had NO MOB-suite match (should be 0)")
    master = master.drop(columns=["_merge"])

    # replicon_matrix.csv: validate only, don't merge
    if os.path.exists(MOB_REPLICON_MATRIX):
        rep_df = pd.read_csv(MOB_REPLICON_MATRIX, dtype={"genome_id": str}, nrows=0)  # header only, it's wide
        rep_ids = pd.read_csv(MOB_REPLICON_MATRIX, usecols=["genome_id"], dtype=str)["genome_id"]
        anchor_ids = set(master["genome_id"])
        rep_ids_set = set(rep_ids)
        missing_from_rep = anchor_ids - rep_ids_set
        extra_in_rep = rep_ids_set - anchor_ids
        log(f"replicon_matrix.csv: {len(rep_ids)} rows, {len(rep_df.columns)} columns — NOT merged into master (see script docstring)")
        log(f"  genomes in anchor but missing from replicon_matrix: {len(missing_from_rep)} (should be 0)")
        log(f"  genomes in replicon_matrix but not in anchor: {len(extra_in_rep)} (should be 0)")
    else:
        log(f"  NOTE: replicon_matrix.csv not found at {MOB_REPLICON_MATRIX}, skipped validation")

    # final QC
    log("-" * 60)
    log(f"FINAL MASTER TABLE: {len(master)} rows x {len(master.columns)} columns")
    log(f"  genome_id unique: {master['genome_id'].nunique()} (should equal row count)")
    na_counts = master.isna().sum()
    cols_with_na = na_counts[na_counts > 0]
    if len(cols_with_na) > 0:
        log(f"  WARNING: {len(cols_with_na)} columns have missing values after merge (top 10 shown):")
        for col, n in cols_with_na.sort_values(ascending=False).head(10).items():
            log(f"    {col}: {n} missing")
    else:
        log("  No missing values anywhere in the merged table.")

    master.to_csv(OUTPUT_MASTER, index=False)
    log(f"\nWrote master table to: {OUTPUT_MASTER}")

    with open(OUTPUT_QC_LOG, "w") as f:
        f.write("\n".join(qc_lines))
    log(f"Wrote QC log to: {OUTPUT_QC_LOG}")

    # summary report: hit counts, most prevalent genes/replicons, distributions
    summary_lines = []

    def slog(msg):
        print(msg)
        summary_lines.append(msg)

    slog("=" * 60)
    slog("Phase 10 Fusion — Summary Report")
    slog(f"Total genomes: {len(master)}")
    slog("=" * 60)

    def hit_summary(df, prefix, num_found_col, label, top_n=20):
        gene_cols = [c for c in df.columns if c.startswith(prefix) and c != num_found_col]
        if num_found_col in df.columns:
            nf = df[num_found_col]
            slog(f"\n--- {label} ---")
            slog(f"Genomes with 0 hits: {(nf == 0).sum()} ({(nf == 0).mean()*100:.1f}%)")
            slog(f"Genomes with >=1 hit: {(nf > 0).sum()} ({(nf > 0).mean()*100:.1f}%)")
            slog(f"Mean hits/genome: {nf.mean():.2f} | Median: {nf.median():.1f} | Max: {nf.max()}")
        # present = not null and not '.'
        presence_counts = {}
        for c in gene_cols:
            present = df[c].notna() & (df[c].astype(str) != ".")
            n = present.sum()
            if n > 0:
                presence_counts[c.replace(prefix, "", 1)] = n
        top = sorted(presence_counts.items(), key=lambda x: -x[1])[:top_n]
        slog(f"\nTop {min(top_n, len(top))} most prevalent {label} genes/replicons:")
        slog(f"{'name':<45} {'genomes':>8} {'pct':>7}")
        for name, n in top:
            slog(f"{name:<45} {n:>8} {n/len(df)*100:>6.1f}%")
        return presence_counts

    card_hits = hit_summary(master, "CARD_", "CARD_NUM_FOUND", "CARD (resistance genes)")
    pf_hits = hit_summary(master, "PF_", "PF_NUM_FOUND", "PlasmidFinder (replicons)")

    # MOB-suite numeric summary
    mob_numeric_cols = [
        "MOB_mob_contig_count", "MOB_primary_cluster_count", "MOB_secondary_cluster_count",
        "MOB_circular_contig_count", "MOB_total_mob_contig_bp", "MOB_replicon_hit_contig_count",
        "MOB_relaxase_hit_contig_count", "MOB_mpf_hit_contig_count", "MOB_orit_hit_contig_count",
    ]
    slog("\n--- MOB-suite (plasmid summary) ---")
    for c in mob_numeric_cols:
        if c in master.columns:
            col = master[c]
            slog(f"{c.replace('MOB_', ''):<28} mean={col.mean():>10.2f}  median={col.median():>8.1f}  max={col.max():>10}  zero={ (col==0).sum():>6} ({(col==0).mean()*100:.1f}%)")

    with open("phase10/phase10_fusion_summary.txt", "w") as f:
        f.write("\n".join(summary_lines))
    slog(f"\nWrote summary report to: phase10/phase10_fusion_summary.txt")


if __name__ == "__main__":
    main()