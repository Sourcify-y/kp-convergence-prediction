#!/usr/bin/env python3
"""
data_cohort_overview_plots.py

Phase 1, Step 3: seven feasibility plots to check whether there's real
temporal depth in the data or if most of it sits in a narrow recent window.

Reads kp_metadata_FULL.csv and writes:
    (a) usable genomes by year
    (b) date completeness by year
    (c) usable genomes by ST
    (d) usable genomes by country
    (e) genomes with both AMR + virulence calls
    (f) longitudinal coverage of major lineages
    (g) dual-risk genomes over time (proxy only)

Usage:
    python3 data_cohort_overview_plots.py --csv kp_metadata_FULL.csv --outdir ./feasibility_out
"""

import argparse
import csv
import re
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_rep_accessions(manifest_path: Path, genome_list_path: Path) -> set:
    # map the Mash-dedup genome list (fasta paths) back to assembly accessions
    # via the download manifest. Both files get \r stripped since stray CRs
    # silently broke Mash before.
    path2acc = {}
    with open(manifest_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            p = (row.get("path") or "").strip().rstrip("\r")
            acc = (row.get("assembly_accession") or "").strip()
            if p and acc:
                path2acc[p] = acc

    kept = set()
    n_unmatched = 0
    with open(genome_list_path) as f:
        for line in f:
            p = line.strip().rstrip("\r")
            if not p:
                continue
            acc = path2acc.get(p)
            if acc:
                kept.add(acc)
            else:
                n_unmatched += 1

    if n_unmatched:
        print(f"[dedup-filter] WARNING: {n_unmatched} paths in {genome_list_path} "
              f"had no matching row in {manifest_path} — check for a path-format "
              f"mismatch (e.g. relative vs. absolute paths) before trusting the "
              f"filtered counts below.")
    print(f"[dedup-filter] Resolved {len(kept)} representative "
          f"accessions from {genome_list_path}")
    return kept


def load_data(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path, dtype=str)
    df["collection_year_num"] = pd.to_numeric(df["collection_year"], errors="coerce")
    # usable = has a parsed year AND an ST call from Kleborate
    df["has_st"] = df["ST"].notna() & (df["ST"].astype(str).str.strip() != "")
    df["usable"] = df["collection_year_num"].notna() & df["has_st"]
    # collapse locus-variant tags (e.g. "ST11-1LV" -> "ST11") so near-identical
    # variants don't count as separate lineages (208/1373 raw STs were LV-tagged)
    df["base_ST"] = df["ST"].apply(
        lambda st: re.sub(r"-\d+LV$", "", str(st)) if pd.notna(st) else st
    )
    return df


def save_fig(fig, outdir: Path, name: str):
    path = outdir / f"{name}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {path}")


# (a) usable genomes by year
def plot_a_usable_by_year(df: pd.DataFrame, outdir: Path):
    usable = df[df["usable"]]
    counts = usable["collection_year_num"].value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(counts.index.astype(int), counts.values, color="#3b6ea5")
    ax.set_xlabel("Collection year")
    ax.set_ylabel("Usable genomes (dated + genotyped)")
    ax.set_title("(a) Usable genomes by year")
    save_fig(fig, outdir, "a_usable_genomes_by_year")

    n_years = counts.shape[0]
    total_usable = int(counts.sum())
    print(f"  -> {total_usable} usable genomes spanning {n_years} distinct years "
          f"({int(counts.index.min())}-{int(counts.index.max())})" if n_years else
          "  -> NO usable genomes with both a year and an ST call")


# (b) date completeness by year
def plot_b_date_completeness(df: pd.DataFrame, outdir: Path):
    quality_order = ["exact", "month-year", "year-only", "inferred", "missing"]
    sub = df[df["collection_year_num"].notna()].copy()
    pivot = (
        sub.groupby(["collection_year_num", "collection_date_quality"])
        .size()
        .unstack(fill_value=0)
    )
    # only the tiers that are actually present, in stacking order
    present_order = [q for q in quality_order if q in pivot.columns]
    pivot = pivot[present_order]

    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = None
    colors = {"exact": "#2b8a3e", "month-year": "#82c91e",
              "year-only": "#f08c00", "inferred": "#e03131", "missing": "#868e96"}
    for q in present_order:
        vals = pivot[q].values
        ax.bar(pivot.index.astype(int), vals, bottom=bottom,
               label=q, color=colors.get(q, None))
        bottom = vals if bottom is None else bottom + vals
    ax.set_xlabel("Collection year")
    ax.set_ylabel("Genome count")
    ax.set_title("(b) Date completeness by year (stacked by quality tier)")
    ax.legend(title="Date quality")
    save_fig(fig, outdir, "b_date_completeness_by_year")


# (c) usable genomes by ST
def plot_c_usable_by_st(df: pd.DataFrame, outdir: Path, top_n: int = 20):
    usable = df[df["usable"]]
    counts = usable["base_ST"].value_counts().head(top_n)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(counts.index[::-1], counts.values[::-1], color="#5f3dc4")
    ax.set_xlabel("Usable genome count")
    ax.set_ylabel("Sequence type (base ST, LV-collapsed)")
    ax.set_title(f"(c) Usable genomes by base ST (top {top_n})")
    save_fig(fig, outdir, "c_usable_genomes_by_st")

    n_st_total = usable["base_ST"].nunique()
    n_raw_st_total = usable["ST"].nunique()
    print(f"  -> {n_st_total} distinct base STs (LV-collapsed) represented "
          f"among usable genomes ({n_raw_st_total} raw ST values before collapsing)")


# (d) usable genomes by country
def plot_d_usable_by_country(df: pd.DataFrame, outdir: Path, top_n: int = 20):
    usable = df[df["usable"] & df["country"].notna()]
    counts = usable["country"].value_counts().head(top_n)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(counts.index[::-1], counts.values[::-1], color="#0c8599")
    ax.set_xlabel("Usable genome count")
    ax.set_ylabel("Country")
    ax.set_title(f"(d) Usable genomes by country (top {top_n})")
    save_fig(fig, outdir, "d_usable_genomes_by_country")

    n_missing_country = int((df["usable"] & df["country"].isna()).sum())
    print(f"  -> {n_missing_country} usable genomes have no country recorded")


# (e) genomes with both AMR + virulence calls
def plot_e_amr_and_virulence(df: pd.DataFrame, outdir: Path):
    df = df.copy()
    df["has_amr"] = pd.to_numeric(df.get("num_resistance_genes"), errors="coerce").fillna(0) > 0
    # virulence = any of the main loci called (not "-")
    vir_cols = ["Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2"]
    vir_cols = [c for c in vir_cols if c in df.columns]
    df["has_virulence"] = df[vir_cols].apply(
        lambda row: any(str(v).strip() not in ("", "-", "nan") for v in row), axis=1
    ) if vir_cols else False

    categories = {
        "AMR only": ((df["has_amr"]) & (~df["has_virulence"])).sum(),
        "Virulence only": ((~df["has_amr"]) & (df["has_virulence"])).sum(),
        "Both AMR + virulence": ((df["has_amr"]) & (df["has_virulence"])).sum(),
        "Neither": ((~df["has_amr"]) & (~df["has_virulence"]) & df["has_st"]).sum(),
    }

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(categories.keys(), categories.values(),
           color=["#e8590c", "#37b24d", "#c92a2a", "#868e96"])
    ax.set_ylabel("Genome count")
    ax.set_title("(e) Genomes with AMR / virulence calls")
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    save_fig(fig, outdir, "e_amr_and_virulence_overlap")

    print(f"  -> {categories['Both AMR + virulence']} genomes carry BOTH AMR and "
          f"virulence-associated determinants — these are your convergence candidates")


# (f) longitudinal coverage of major lineages
def plot_f_lineage_longitudinal(df: pd.DataFrame, outdir: Path, top_n: int = 8):
    usable = df[df["usable"]]
    top_sts = usable["base_ST"].value_counts().head(top_n).index.tolist()

    fig, ax = plt.subplots(figsize=(10, 6))
    for st in top_sts:
        sub = usable[usable["base_ST"] == st]
        years = sorted(sub["collection_year_num"].dropna().astype(int).unique())
        ax.scatter(years, [st] * len(years), s=60)
    ax.set_xlabel("Collection year")
    ax.set_ylabel("Base ST (LV-collapsed)")
    ax.set_title(f"(f) Longitudinal coverage of top {top_n} lineages")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    save_fig(fig, outdir, "f_lineage_longitudinal_coverage")

    for st in top_sts:
        sub = usable[usable["base_ST"] == st]
        years = sub["collection_year_num"].dropna()
        span = int(years.max() - years.min()) if len(years) > 1 else 0
        print(f"  -> {st}: n={len(sub)}, years span={span} "
              f"({int(years.min()) if len(years) else '-'}-{int(years.max()) if len(years) else '-'})")


# (g) dual-risk genomes over time
def plot_g_candidate_acquisition_events(df: pd.DataFrame, outdir: Path, top_n: int = 6):
    # NOTE: this is NOT a real acquisition-event count. A real one needs
    # ancestral-state reconstruction (Phase 4) to separate new acquisition from
    # ancestral presence, loss/reappearance, HGT or genotyping artifacts. This
    # just plots genomes carrying both AMR and virulence by year for the top
    # lineages, to see if dual-risk genomes span enough of the timeline to make
    # Phase 4 worth doing. Don't use these counts as Phase 4/8 results.
    df = df.copy()
    df["has_amr"] = pd.to_numeric(df.get("num_resistance_genes"), errors="coerce").fillna(0) > 0
    vir_cols = ["Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2"]
    vir_cols = [c for c in vir_cols if c in df.columns]
    df["has_virulence"] = df[vir_cols].apply(
        lambda row: any(str(v).strip() not in ("", "-", "nan") for v in row), axis=1
    ) if vir_cols else False
    df["dual_risk"] = df["has_amr"] & df["has_virulence"]

    usable = df[df["usable"]]
    top_sts = usable["base_ST"].value_counts().head(top_n).index.tolist()

    fig, ax = plt.subplots(figsize=(10, 6))
    for st in top_sts:
        sub = usable[(usable["base_ST"] == st) & usable["dual_risk"]]
        years = sorted(sub["collection_year_num"].dropna().astype(int).unique())
        if years:
            ax.scatter(years, [st] * len(years), s=80, marker="D", color="#c92a2a")
    ax.set_xlabel("Collection year")
    ax.set_ylabel("ST")
    ax.set_title(f"(g) Candidate dual-risk genomes over time [PROXY — see docstring]\n"
                 f"NOT a validated acquisition-event count (Phase 4 required for that)")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    save_fig(fig, outdir, "g_candidate_acquisition_events_proxy")

    n_dual_total = int(usable["dual_risk"].sum())
    print(f"  -> {n_dual_total} usable genomes show dual AMR+virulence signal "
          f"(crude proxy, NOT validated acquisition events — see Phase 4)")


# summary numbers for the Go/No-Go table
def print_go_no_go_inputs(df: pd.DataFrame):
    usable = df[df["usable"]]
    print("\n" + "=" * 60)
    print("GO/NO-GO TABLE INPUTS (fill into Phase 1 Step 10)")
    print("=" * 60)
    print(f"Genomes with reliable collection year : {len(usable)}")
    print(f"Years represented                     : {usable['collection_year_num'].nunique()}")
    print(f"Lineages (base ST, LV-collapsed)       : {usable['base_ST'].nunique()}")
    print(f"Lineages (raw ST, pre-collapse)        : {usable['ST'].nunique()}")
    df2 = df.copy()
    df2["has_amr"] = pd.to_numeric(df2.get("num_resistance_genes"), errors="coerce").fillna(0) > 0
    vir_cols = ["Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2"]
    vir_cols = [c for c in vir_cols if c in df2.columns]
    df2["has_virulence"] = df2[vir_cols].apply(
        lambda row: any(str(v).strip() not in ("", "-", "nan") for v in row), axis=1
    ) if vir_cols else False
    n_both = int((df2["has_amr"] & df2["has_virulence"]).sum())
    print(f"Genomes with AMR + virulence calls     : {n_both}")
    print("=" * 60)
    print("NOTE: 'Candidate acquisition events' and 'lineages with >=N years "
          "of observations' require manual/further analysis beyond this script.")
    print("=" * 60 + "\n")


def dump_yearly_counts(df: pd.DataFrame, outdir: Path):
    # one row per year: usable and dual-risk counts, used to pick T1/T2 in Step 4
    df = df.copy()
    df["has_amr"] = pd.to_numeric(df.get("num_resistance_genes"), errors="coerce").fillna(0) > 0
    vir_cols = ["Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2"]
    vir_cols = [c for c in vir_cols if c in df.columns]
    df["has_virulence"] = df[vir_cols].apply(
        lambda row: any(str(v).strip() not in ("", "-", "nan") for v in row), axis=1
    ) if vir_cols else False
    df["dual_risk"] = df["has_amr"] & df["has_virulence"]

    usable = df[df["usable"]].copy()
    usable["year"] = usable["collection_year_num"].astype(int)

    by_year = usable.groupby("year").agg(
        usable_genomes=("year", "size"),
        dual_risk_genomes=("dual_risk", "sum"),
    ).reset_index().sort_values("year")
    by_year["cumulative_usable"] = by_year["usable_genomes"].cumsum()
    by_year["cumulative_dual_risk"] = by_year["dual_risk_genomes"].cumsum()

    out_path = outdir / "yearly_counts_for_t1_t2.csv"
    by_year.to_csv(out_path, index=False)
    print(f"[dump] wrote per-year usable/dual-risk counts -> {out_path}")
    return by_year


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True, type=Path,
                        help="Path to kp_metadata_FULL.csv")
    parser.add_argument("--outdir", required=True, type=Path,
                        help="Directory to write plot PNGs into")
    parser.add_argument("--manifest", type=Path, default=None,
                        help="Path to 03_download_manifest.csv (path <-> "
                             "assembly_accession mapping). Required if "
                             "--genome-list is given.")
    parser.add_argument("--genome-list", type=Path, default=None,
                        help="Path to a Mash-dedup representative genome list "
                             "(e.g. mash_test/representative_genomes.txt), one "
                             "fasta path per line. When given (with --manifest), "
                             "restricts every plot/count to just these genomes "
                             "instead of the full metadata CSV.")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if bool(args.manifest) != bool(args.genome_list):
        parser.error("--manifest and --genome-list must be given together.")

    df = load_data(args.csv)
    print(f"[load] {len(df)} total rows, {df['usable'].sum()} usable "
          f"(dated + genotyped) genomes\n")

    if args.genome_list:
        kept = load_rep_accessions(args.manifest, args.genome_list)
        before = len(df)
        df = df[df["genome_id"].isin(kept)].copy()
        print(f"[dedup-filter] Restricted metadata to representative genomes: "
              f"{before} -> {len(df)} rows "
              f"({df['usable'].sum()} usable after filtering)\n")

    print("(a) usable genomes by year")
    plot_a_usable_by_year(df, args.outdir)
    print("\n(b) date completeness by year")
    plot_b_date_completeness(df, args.outdir)
    print("\n(c) usable genomes by ST")
    plot_c_usable_by_st(df, args.outdir)
    print("\n(d) usable genomes by country")
    plot_d_usable_by_country(df, args.outdir)
    print("\n(e) genomes with AMR + virulence calls")
    plot_e_amr_and_virulence(df, args.outdir)
    print("\n(f) longitudinal coverage of major lineages")
    plot_f_lineage_longitudinal(df, args.outdir)
    print("\n(g) candidate dual-risk genomes over time [PROXY]")
    plot_g_candidate_acquisition_events(df, args.outdir)

    dump_yearly_counts(df, args.outdir)
    print_go_no_go_inputs(df)
    print(f"[done] All plots written to {args.outdir}")


if __name__ == "__main__":
    main()