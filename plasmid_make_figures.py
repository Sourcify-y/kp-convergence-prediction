#!/usr/bin/env python3
"""
plasmid_make_figures.py -- figures and tables for the Phase 10 consolidation.

Reads phase10_out/localization.csv (CARD + PlasmidFinder hits joined to
MOB-suite contig assignments) and, optionally, mob_consolidated/
contig_to_cluster.tsv. Writes PNG + SVG figures and CSV tables.

Notes:
  * required columns are checked up front and the script stops if any are missing
  * each figure prints the numbers it plotted so you can check them against the QC file
  * the counting unit is genome x resistance class, not raw hit rows (a genome
    with five alleles of one class counts once)
  * "plasmid" always means a MOB-suite-reconstructed plasmid unit
    (primary_cluster_id), not an experimentally closed plasmid

Usage:
  python plasmid_make_figures.py \\
      --localization phase10_out/localization.csv \\
      --mob-table mob_consolidated/contig_to_cluster.tsv \\
      --outdir phase10_figs
"""
import argparse
import os
import sys
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LOCKED = ["AGly_acquired", "Bla_acquired", "Bla_ESBL_acquired", "Bla_Carb_acquired",
          "Col_acquired", "Fcyn_acquired", "Flq_acquired", "MLS_acquired",
          "Phe_acquired", "Rif_acquired", "Sul_acquired", "Tet_acquired",
          "Tgc_acquired", "Tmt_acquired"]

REQUIRED_COLS = ["genome_id", "database", "gene", "resistance_class", "contig_id",
                 "molecule_type", "primary_cluster_id"]

# rationale text for the supplementary label table
LABEL_NOTES = {
    "intrinsic_chromosomal_efflux_regulatory": "Intrinsic K. pneumoniae efflux/porin/regulator/LPS genes; near-universal, not acquired resistance",
    "oqxAB_efflux_not_in_locked_set": "OqxAB efflux pump; near-universal, split out of Flq",
    "fosA_chromosomal_not_in_locked_set": "FosA5/6/7; >99% chromosomal, split out of Fcyn",
    "disinfectant_qac_not_in_locked_set": "Quaternary ammonium resistance (class 1 integron marker); not a locked drug class",
    "bleomycin_not_in_locked_set": "BRP(MBL)/BLMT bleomycin resistance; NDM-linked marker, not a beta-lactamase",
    "streptothricin_not_in_locked_set": "SAT streptothricin acetyltransferases; distinct drug, not in locked set",
    "SHV_family_see_Kleborate_call": "SHV: narrow vs ESBL depends on mutations; use Kleborate call",
    "TEM_family_see_Kleborate_call": "TEM: narrow vs ESBL depends on mutations; use Kleborate call",
    "AmpC_acquired_not_in_locked_set": "Acquired AmpC (DHA/CMY/FOX/MOX/ACC/LAP); no slot in locked Bla_* split",
    "intrinsic_chromosomal_not_acquired": "LEN/OKP/KLUC intrinsic chromosomal beta-lactamases",
    "GES_family_needs_allele_check": "GES: ESBL vs carbapenemase depends on allele; not classified",
    "Bla_family_needs_allele_check": "SFO-1/FRI-1/RSA1-1: CARD drug class conflicts with expectation; not classified",
    "chromosomal_or_other_bla_not_in_locked_set": "SED-2/EC-18 chromosomal or other beta-lactamases",
}

# Okabe-Ito colour-blind-safe palette
C_PLASMID = "#D55E00"
C_BOTH = "#E69F00"
C_CHROM = "#0072B2"
C_GREY = "#999999"

plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 16, "axes.labelsize": 14,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 100, "savefig.bbox": "tight",
})
FIGSIZE = (12, 6.75)  # 16:9, fits a slide and scales to a paper column


def short(c):
    # readable class label ('Bla_acquired' alone is ambiguous)
    return {"Bla_acquired": "Bla (narrow)", "Bla_ESBL_acquired": "Bla ESBL",
            "Bla_Carb_acquired": "Bla carbapenemase"}.get(c, c.replace("_acquired", ""))


def die(msg):
    print(f"\nERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def save(fig, outdir, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), dpi=200)
    plt.close(fig)
    print(f"  wrote {name}.png / .svg")


def load_localization(path):
    if not os.path.exists(path):
        die(f"{path} not found")
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        die(f"missing columns {missing}. Columns present: {list(df.columns)}")
    df["database"] = df["database"].str.lower().str.strip()
    dbs = set(df["database"])
    if not dbs <= {"card", "plasmidfinder"}:
        die(f"unexpected database labels {dbs} (expected card / plasmidfinder)")
    df["is_plasmid"] = (df["molecule_type"].str.lower() == "plasmid")
    df["cluster"] = df["primary_cluster_id"].str.strip()
    df.loc[df["cluster"].isin(["", "-", "NA", "nan"]), "cluster"] = ""
    return df


def genome_class_table(card):
    # one row per (genome, class): plasmid only / chromosome only / both
    g = (card.groupby(["genome_id", "resistance_class"])["is_plasmid"]
         .agg(["any", "all"]).reset_index())
    g["location"] = np.select(
        [g["all"], ~g["any"]], ["plasmid_only", "chromosome_only"], default="both")
    return g


def table_summary(df, n_genomes, outdir):
    card = df[df.database == "card"]
    pf = df[df.database == "plasmidfinder"]
    rows = [
        ("genomes in cohort", n_genomes),
        ("genomes with >=1 hit", df.genome_id.nunique()),
        ("total hit rows", len(df)),
        ("CARD hit rows", len(card)),
        ("PlasmidFinder hit rows", len(pf)),
        ("hit rows on chromosome", int((~df.is_plasmid).sum())),
        ("hit rows on MOB plasmid contigs", int(df.is_plasmid.sum())),
        ("distinct MOB plasmid units with >=1 hit (genome x cluster)",
         df[df.is_plasmid & (df.cluster != "")][["genome_id", "cluster"]].drop_duplicates().shape[0]),
    ]
    t = pd.DataFrame(rows, columns=["metric", "value"])
    t.to_csv(os.path.join(outdir, "table0_cohort_summary.csv"), index=False)
    print(t.to_string(index=False))


def table_and_fig_classes(df, n_genomes, outdir):
    card = df[df.database == "card"]
    gc = genome_class_table(card)

    recs = []
    for cls in LOCKED:
        sub = gc[gc.resistance_class == cls]
        n = len(sub)
        counts = sub.location.value_counts()
        recs.append({
            "resistance_class": cls,
            "genomes_with_class": n,
            "pct_of_cohort": round(100 * n / n_genomes, 2),
            "plasmid_only": int(counts.get("plasmid_only", 0)),
            "both": int(counts.get("both", 0)),
            "chromosome_only": int(counts.get("chromosome_only", 0)),
        })
    t1 = pd.DataFrame(recs)
    t1["pct_any_plasmid"] = np.where(
        t1.genomes_with_class > 0,
        (100 * (t1.plasmid_only + t1.both) / t1.genomes_with_class.replace(0, np.nan)).round(1), np.nan)
    t1.to_csv(os.path.join(outdir, "table1_locked_classes.csv"), index=False)
    print("\nTable 1 (locked classes, genome x class):")
    print(t1.to_string(index=False))

    # figure 1: two panels sharing the y order
    plot = t1[t1.genomes_with_class > 0].sort_values("pct_any_plasmid")
    if plot.empty:
        print("  no locked-class hits found; skipping figure 1")
        return gc
    labels = [short(c) for c in plot.resistance_class]
    y = np.arange(len(plot))

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=FIGSIZE, sharey=True,
                                            gridspec_kw={"width_ratios": [1, 1.25]})
    ax_left.barh(y, plot.pct_of_cohort, color=C_GREY)
    ax_left.set_yticks(y)
    ax_left.set_yticklabels(labels)
    ax_left.set_xlabel("% of all genomes")
    ax_left.set_ylabel("Resistance class")
    for yi, v in zip(y, plot.pct_of_cohort):
        ax_left.text(v + 0.8, yi, f"{v:.0f}%", va="center", fontsize=11)

    n_carry = plot.genomes_with_class.values
    pct_plasmid = 100 * plot.plasmid_only.values / n_carry
    pct_both = 100 * plot.both.values / n_carry
    pct_chrom = 100 * plot.chromosome_only.values / n_carry
    ax_right.barh(y, pct_plasmid, color=C_PLASMID, label="plasmid only")
    ax_right.barh(y, pct_both, left=pct_plasmid, color=C_BOTH, label="plasmid + chromosome")
    ax_right.barh(y, pct_chrom, left=pct_plasmid + pct_both, color=C_CHROM, label="chromosome only")
    ax_right.set_xlim(0, 100)
    ax_right.set_xlabel("% of genomes carrying class")
    ax_right.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, frameon=False, fontsize=11)
    for yi, n in zip(y, n_carry):
        ax_right.text(101, yi, f"n={n}", va="center", fontsize=10, color="#444")
    fig.tight_layout()
    save(fig, outdir, "fig1_class_prevalence_and_location")
    return gc


def table_label_decisions(df, outdir):
    card = df[df.database == "card"]
    rows = []
    for lab, note in LABEL_NOTES.items():
        sub = card[card.resistance_class == lab]
        rows.append({"label": lab, "hit_rows": len(sub),
                     "genomes": sub.genome_id.nunique(),
                     "pct_hits_on_plasmid": round(100 * sub.is_plasmid.mean(), 1) if len(sub) else np.nan,
                     "rationale": note})
    t = pd.DataFrame(rows).sort_values("hit_rows", ascending=False)
    t.to_csv(os.path.join(outdir, "tableS1_non_locked_labels.csv"), index=False)
    print("\nTable S1 (non-locked labels):")
    print(t.drop(columns="rationale").to_string(index=False))
    other = set(card.resistance_class) - set(LOCKED) - set(LABEL_NOTES) - {""}
    if other:
        print(f"  NOTE: labels present but not described in LABEL_NOTES: {sorted(other)}")


def fig_replicons(df, n_genomes, outdir, top=20):
    pf = df[df.database == "plasmidfinder"]
    if pf.empty:
        print("  no PlasmidFinder rows; skipping replicon figure")
        return []
    c = pf.groupby("gene").genome_id.nunique().sort_values(ascending=False)
    t = (c / n_genomes * 100).round(2).rename("pct_of_cohort").to_frame()
    t["genomes"] = c
    t.to_csv(os.path.join(outdir, "table2_replicons.csv"))
    top_n = t.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.barh(np.arange(len(top_n)), top_n.pct_of_cohort, color=C_PLASMID)
    ax.set_yticks(np.arange(len(top_n)))
    ax.set_yticklabels([g[:32] for g in top_n.index], fontsize=11)
    ax.set_xlabel("% of genomes")
    ax.set_ylabel("PlasmidFinder replicon")
    for yi, v in enumerate(top_n.pct_of_cohort):
        ax.text(v + 0.4, yi, f"{v:.1f}%", va="center", fontsize=10)
    fig.tight_layout()
    save(fig, outdir, "fig2_top_replicons")
    print("  top replicons:", list(t.head(5).index))
    return list(t.head(10).index)


def plasmid_unit_table(df):
    # hits that sit on a MOB plasmid unit that has a cluster id
    return df[df.is_plasmid & (df.cluster != "")]


def fig_class_replicon(df, gc, top_reps, outdir):
    if not top_reps:
        return
    pl = plasmid_unit_table(df)
    card = pl[(pl.database == "card") & (pl.resistance_class.isin(LOCKED))]
    rep = pl[(pl.database == "plasmidfinder") & (pl.gene.isin(top_reps))]
    if card.empty or rep.empty:
        print("  no plasmid-localized class/replicon pairs; skipping fig 3")
        return
    merged = card[["genome_id", "cluster", "resistance_class"]].drop_duplicates().merge(
        rep[["genome_id", "cluster", "gene"]].drop_duplicates(), on=["genome_id", "cluster"])
    pair = merged.groupby(["resistance_class", "gene"]).genome_id.nunique().unstack(fill_value=0)
    carrying = gc.groupby("resistance_class").genome_id.nunique()
    classes = [c for c in LOCKED if c in pair.index]
    reps = [r for r in top_reps if r in pair.columns]
    mat = pair.reindex(index=classes, columns=reps, fill_value=0)
    pct = mat.div(carrying.reindex(classes).values, axis=0) * 100
    pct.to_csv(os.path.join(outdir, "table3_class_replicon_colocation_pct.csv"))
    mat.to_csv(os.path.join(outdir, "table3_class_replicon_colocation_counts.csv"))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    im = ax.imshow(pct.values, cmap="Oranges", aspect="auto", vmin=0)
    ax.set_xticks(range(len(reps)))
    ax.set_xticklabels([r[:30] for r in reps], rotation=60, ha="right", fontsize=11)
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels([short(c) for c in classes])
    for i in range(pct.shape[0]):
        for j in range(pct.shape[1]):
            v = pct.values[i, j]
            if v >= 1:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=10,
                        color="white" if v > pct.values.max() * 0.6 else "black")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("% of class-carrying genomes\non same plasmid unit")
    ax.set_xlabel("PlasmidFinder replicon")
    ax.set_ylabel("Resistance class")
    ax.spines[:].set_visible(False)
    fig.tight_layout()
    save(fig, outdir, "fig3_class_replicon_colocation")


def fig_class_class(df, gc, outdir, min_both=30):
    pl = plasmid_unit_table(df)
    card = pl[(pl.database == "card") & (pl.resistance_class.isin(LOCKED))]
    present = [c for c in LOCKED if (gc.resistance_class == c).any()]
    genome_sets = gc.groupby("resistance_class").genome_id.apply(set)
    unit_sets = (card.groupby(["genome_id", "cluster"]).resistance_class.apply(set).reset_index())
    same_unit = {}
    for _, r in unit_sets.iterrows():
        for a, b in combinations(sorted(r.resistance_class), 2):
            same_unit.setdefault((a, b), set()).add(r.genome_id)
    n = len(present)
    mat = np.full((n, n), np.nan)
    rows = []
    for i, a in enumerate(present):
        for j, b in enumerate(present):
            if i <= j:
                continue
            both = genome_sets[a] & genome_sets[b]
            if len(both) < min_both:
                continue
            k = len(same_unit.get(tuple(sorted((a, b))), set()))
            mat[i, j] = 100 * k / len(both)
            rows.append({"class_a": a, "class_b": b, "genomes_with_both": len(both),
                         "genomes_same_plasmid_unit": k, "pct_same_unit": round(mat[i, j], 1)})
    if not rows:
        print("  no class pairs with enough genomes; skipping fig 4")
        return
    pd.DataFrame(rows).to_csv(os.path.join(outdir, "table4_class_pairs_same_plasmid.csv"), index=False)
    # lower triangle only: drop the always-empty first row and last column
    sub = mat[1:, :-1]
    labels = [short(c) for c in present]
    row_labels, col_labels = labels[1:], labels[:-1]
    fig, ax = plt.subplots(figsize=(9, 7.5))
    im = ax.imshow(sub, cmap="Oranges", vmin=0, vmax=100)
    ax.set_xticks(range(len(col_labels))); ax.set_xticklabels(col_labels, rotation=60, ha="right")
    ax.set_yticks(range(len(row_labels))); ax.set_yticklabels(row_labels)
    for i in range(sub.shape[0]):
        for j in range(sub.shape[1]):
            if not np.isnan(sub[i, j]):
                ax.text(j, i, f"{sub[i, j]:.0f}", ha="center", va="center", fontsize=10,
                        color="white" if sub[i, j] > 60 else "black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046)
    cbar.set_label("% of genomes with both classes\non same plasmid unit")
    ax.set_xlabel("Resistance class")
    ax.set_ylabel("Resistance class")
    ax.spines[:].set_visible(False)
    fig.tight_layout()
    save(fig, outdir, "fig4_class_pairs_same_plasmid")


def fig_plasmid_burden(mob_path, all_genomes, outdir):
    if not mob_path or not os.path.exists(mob_path):
        print("  no MOB contig table given/found; skipping plasmid-burden figure")
        return
    mob = pd.read_csv(mob_path, sep="\t", dtype=str, keep_default_na=False)
    genome_col = next((c for c in ("genome_id", "sample_id") if c in mob.columns), None)
    type_col = "molecule_type" if "molecule_type" in mob.columns else None
    cluster_col = "primary_cluster_id" if "primary_cluster_id" in mob.columns else None
    if not (genome_col and type_col and cluster_col):
        print(f"  MOB table columns not recognised {list(mob.columns)}; skipping plasmid-burden figure")
        return
    plasmid_rows = mob[(mob[type_col].str.lower() == "plasmid") & (~mob[cluster_col].isin(["", "-", "NA"]))]
    per_genome = plasmid_rows.groupby(genome_col)[cluster_col].nunique()
    per_genome = per_genome.reindex(sorted(all_genomes), fill_value=0)
    if per_genome.sum() == 0:
        print("  genome IDs in MOB table did not match localization IDs (all zero); skipping. "
              f"Example MOB id: {mob[genome_col].iloc[0]!r}")
        return
    cap = 12
    capped = per_genome.clip(upper=cap)
    counts = capped.value_counts().sort_index()
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.bar(counts.index, counts.values, color=C_PLASMID)
    ax.set_xlabel("Plasmid units per genome")
    ax.set_xticks(range(0, cap + 1, 2))
    ax.set_xticklabels([str(t) if t < cap else f"{cap}+" for t in range(0, cap + 1, 2)])
    ax.set_ylabel("Genomes")
    fig.tight_layout()
    save(fig, outdir, "fig5_plasmid_burden")
    per_genome.rename("n_plasmid_units").to_csv(os.path.join(outdir, "table5_plasmid_units_per_genome.csv"))
    print(f"  plasmid units/genome: median {per_genome.median():.0f}, mean {per_genome.mean():.2f}, "
          f"genomes with 0: {(per_genome == 0).sum()}")


CAPTIONS = """Figure captions and caveats (kept out of the figures themselves)

Fig 1. Left: percentage of all genomes carrying each acquired resistance class (CARD/Abricate; genome x class counted once).
Right: for genomes carrying the class, where the gene is found: on a MOB-suite-reconstructed plasmid unit only, on both plasmid and chromosome
(different genes of the same class in different places), or chromosome only. n = genomes carrying the class.
Bla (narrow) and Bla ESBL exclude SHV/TEM, which are classified from Kleborate. Classes with no CARD hits (Tgc) are omitted.
'Plasmid' means a MOB-suite-reconstructed plasmid unit from draft assemblies, not an experimentally closed plasmid.

Fig 2. Most frequent PlasmidFinder replicons, as percentage of genomes (n = 8,318).

Fig 3. Percentage of class-carrying genomes in which a gene of that class and the given replicon are on the same reconstructed plasmid unit.
Empty cells are below 1%. Small classes (Col, Fcyn) give noisy percentages.

Fig 4. Among genomes carrying both classes, percentage in which the two classes are on the same reconstructed plasmid unit (pairs with >= 30 genomes).
No chance baseline is applied; interpret against the plasmid-unit burden (median 4 per genome).

Fig 5. Reconstructed plasmid units per genome; last bin pools 12 or more. Median and mean are in table5_plasmid_units_per_genome.csv.
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--localization", default="phase10_out/localization.csv")
    parser.add_argument("--mob-table", default="mob_consolidated/contig_to_cluster.tsv")
    parser.add_argument("--outdir", default="phase10_figs")
    parser.add_argument("--n-genomes", type=int, default=8318)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = load_localization(args.localization)
    ng = df.genome_id.nunique()
    print(f"Loaded {len(df)} hit rows, {ng} genomes")
    if ng != args.n_genomes:
        print(f"  WARNING: {ng} genomes with hits vs expected cohort {args.n_genomes}")

    print("\n== Table 0: cohort summary ==")
    table_summary(df, args.n_genomes, args.outdir)
    print("\n== Figure 1 / Table 1 ==")
    gc = table_and_fig_classes(df, args.n_genomes, args.outdir)
    table_label_decisions(df, args.outdir)
    print("\n== Figure 2 / Table 2: replicons ==")
    top_reps = fig_replicons(df, args.n_genomes, args.outdir)
    print("\n== Figure 3: class x replicon co-location ==")
    fig_class_replicon(df, gc, top_reps, args.outdir)
    print("\n== Figure 4: class pairs on the same plasmid unit ==")
    fig_class_class(df, gc, args.outdir)
    print("\n== Figure 5: plasmid burden ==")
    fig_plasmid_burden(args.mob_table, set(df.genome_id), args.outdir)
    with open(os.path.join(args.outdir, "figure_captions.txt"), "w") as f:
        f.write(CAPTIONS)
    print(f"\nDone. Outputs in {args.outdir}/ (captions in figure_captions.txt)")


if __name__ == "__main__":
    main()