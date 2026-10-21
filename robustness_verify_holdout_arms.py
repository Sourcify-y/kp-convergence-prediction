#!/usr/bin/env python3
"""
robustness_verify_holdout_arms.py -- independent re-scoring of the Phase 9 arms (lineage holdout,
geographic holdout, pooled sensitivity check) with the standardised LOOSE method, straight from
the raw -routes.txt files.

Why it exists: The project log (not included in this repository) has two incompatible descriptions of the
lineage-holdout result (clean re-run 0.677/0.729 with 0% NA vs "0.788/0.889 with 32% NA"). The
loose method can't produce NA for full-permutation routes except when truncation lands on the last
slot, so the NA-bearing CSVs are almost certainly stale outputs of the older strict-method run.
This script recomputes everything from scratch, prints the absolute path + mtime of every file it
touches, and (optionally) diffs its predictions against the CSVs the R script wrote.

Usage (from the project root):
  python3 robustness_verify_holdout_arms.py --arm holdout --routes-dir phase9_out
  python3 robustness_verify_holdout_arms.py --arm geo     --routes-dir phase9_out
  python3 robustness_verify_holdout_arms.py --arm sens
"""
import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd

import shared_lib as L

REPORTED = {   # numbers from the project log, for eyeballing only
    "holdout": "clean loose re-run: VALIDATE 0.677 (n=40) / TEST 0.729 (n=19); older STRICT run: 0.788 / 0.889 (13/40 and 6/19 NA)",
    "geo":     "clean loose re-run: VALIDATE 0.639 (n=24) / TEST 0.400 (n=12)",
    "sens":    "pooled 0.669 / 0.685 vs s0-only 0.664 / 0.686",
}


def mtime(p):
    p = Path(p)
    return dt.datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M") if p.exists() else "MISSING"


def report_window(name, df):
    y = df["became_dual_risk"].values
    s = df["predicted_risk"].values
    scored = ~np.isnan(s)
    print(f"\n  [{name}] table n={len(df)}  scored n={int(scored.sum())}  NA={int((~scored).sum())}  "
          f"positives(scored)={int(y[scored].sum())}  base rate={y[scored].mean():.3f}")
    for state, sub in df[scored].groupby("starting_state"):
        print(f"      {state:16s} n={len(sub):3d} pos={int(sub['became_dual_risk'].sum()):3d}")
    auc_val, lo, hi = L.delong_ci(y, s)
    boot_lo, boot_hi = L.bootstrap_ci(y, s, n_boot=5000)
    obs, p = L.perm_pvalue(y, s, df["starting_state"].values, n_perm=20000)
    print(f"      AUROC {auc_val:.3f}   DeLong 95% CI [{lo:.3f}, {hi:.3f}]   bootstrap CI [{boot_lo:.3f}, {boot_hi:.3f}]")
    print(f"      permutation p (labels shuffled within starting_state) = {p:.4f}   AUPRC = {L.avg_precision(y, s):.3f}")
    if int(y[scored].sum()) < 5:
        print("      WARNING: <5 positives among scored lineages -- treat AUROC as descriptive only.")
    return dict(window=name, n_table=len(df), n_scored=int(scored.sum()), n_na=int((~scored).sum()),
                n_pos=int(y[scored].sum()), auroc=auc_val, delong_lo=lo, delong_hi=hi, boot_lo=boot_lo, boot_hi=boot_hi, perm_p=p)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", required=True, choices=["holdout", "geo", "sens"])
    parser.add_argument("--matrix", default="phase3_out/character_matrix.csv")
    parser.add_argument("--tables-dir", default="phase8_out")
    parser.add_argument("--routes-dir", default=None,
                        help="dir with the arm's -routes.txt files (default: phase9_out for holdout/geo, phase8_out for sens)")
    parser.add_argument("--assign", default=None, help="holdout assignment csv (default per arm, in phase9_out/)")
    parser.add_argument("--r-output-dir", default=None, help="where the R script wrote *_scored.csv (default phase9_out / phase8_out)")
    parser.add_argument("--out", default="phase9_verify_out")
    parser.add_argument("--train-end", type=int, default=L.TRAIN_END)
    parser.add_argument("--seeds", default="0,1,2")
    args = parser.parse_args()

    out = L.ensure_dir(args.out)
    seeds = [int(x) for x in args.seeds.split(",")]
    arm = args.arm
    if arm == "holdout":
        routes_dir, pattern, order = args.routes_dir or "phase9_out", "kp_holdout_run_model1_s{s}-routes.txt", L.CHAR_ORDER_HDR
        assign = args.assign or "phase9_out/lineage_holdout_assignment.csv"
        r_csv_pattern = "phase9_{w}_holdout_scored.csv"; r_out_dir = args.r_output_dir or "phase9_out"
    elif arm == "geo":
        routes_dir, pattern, order = args.routes_dir or "phase9_out", "kp_geo_holdout_run_model1_s{s}-routes.txt", L.CHAR_ORDER_HDR
        assign = args.assign or "phase9_out/geo_holdout_assignment.csv"
        r_csv_pattern = "phase9_{w}_geo_scored.csv"; r_out_dir = args.r_output_dir or "phase9_out"
    else:
        routes_dir, pattern, order = args.routes_dir or "phase8_out", "kp_train_only_model1_s{s}-routes.txt", L.CHAR_ORDER_P8
        assign = None
        r_csv_pattern = "phase8_sensitivity_check_pooled_{w}.csv"; r_out_dir = args.r_output_dir or "phase8_out"

    print(f"=== robustness_verify_holdout_arms: {arm} ===")
    print(f"reported in project log: {REPORTED[arm]}")
    cm = L.load_matrix(args.matrix)
    print(f"matrix: {Path(args.matrix).resolve()}  [{mtime(args.matrix)}]  {len(cm)} genomes")
    tables = {w: L.load_task_table(Path(args.tables_dir) / f"phase8_{w}_task_table.csv") for w in ("validate", "test")}
    for w, tab in tables.items():
        print(f"task table {w}: {len(tab)} lineages [{mtime(Path(args.tables_dir) / f'phase8_{w}_task_table.csv')}]")

    if assign:
        asg = L.read_csv_str(assign)
        print(f"assignment: {Path(assign).resolve()} [{mtime(assign)}]  columns={list(asg.columns)}  unique base_ST={asg['base_ST'].nunique()}")
        held = set(asg["base_ST"])
        tables = {w: tab[tab["base_ST"].isin(held)].reset_index(drop=True) for w, tab in tables.items()}
        print("  -> if this file lists ALL lineages with a held-out flag column, the R filter (%in%) would have scored everything; "
              "unique count above should equal 228 (lineage) / 92 (geo).")

    paths = [Path(routes_dir) / pattern.format(s=s) for s in seeds]
    for p in paths:
        print(f"routes: {p.resolve()}  [{mtime(p)}]")
    states = L.lineage_states(cm, args.train_end, order)
    variants = {"pooled": L.load_pooled_routes(paths)}
    if arm == "sens":
        variants["s0-only"] = L.load_routes(paths[0])
    print(f"pooled trajectories: {variants['pooled'].shape[0]}")

    all_rows = []
    for vname, routes in variants.items():
        scorer = L.RouteScorer(routes)
        print(f"\n--- scoring variant: {vname} ({routes.shape[0]} routes) ---")
        for w, tab in tables.items():
            indexed = tab.set_index("base_ST")
            sc, nv, no_state = L.score_task(indexed, states, scorer, order)
            res = tab.copy()
            res["predicted_risk"] = sc
            res["n_valid_routes"] = nv
            if no_state:
                print(f"  WARNING: {no_state} lineages in {w} table have no TRAIN genomes in the matrix")
            res.to_csv(out / f"verify_{arm}_{vname}_{w}_scored.csv", index=False)
            row = report_window(f"{w.upper()} {vname}", res)
            row.update(arm=arm, variant=vname)
            all_rows.append(row)

            # diff against the R-written CSV (stale-file detector)
            if vname == "pooled":
                r_path = Path(r_out_dir) / r_csv_pattern.format(w=w)
                if r_path.exists():
                    r_df = L.read_csv_str(r_path)
                    merged = res.merge(r_df[["base_ST", "predicted_risk"]], on="base_ST", how="left", suffixes=("", "_R"))
                    both = merged["predicted_risk"].notna() & merged["predicted_risk_R"].notna()
                    max_diff = float((merged.loc[both, "predicted_risk"] - merged.loc[both, "predicted_risk_R"]).abs().max()) if both.any() else np.nan
                    print(f"      R file {r_path.resolve()} [{mtime(r_path)}]: rows={len(r_df)}, NA in R file={int(r_df['predicted_risk'].isna().sum())}, "
                          f"max |python - R| on rows scored by both = {max_diff:.2e}")
                    if r_df["predicted_risk"].isna().sum() > 0 and res["predicted_risk"].isna().sum() == 0:
                        print("      >>> R file has NA where the loose-method recomputation has none: that R file is most likely a stale output "
                              "of the older strict-method run. Do not use it; re-run the R script and re-check timestamps and NA counts.")
                    elif not np.isnan(max_diff) and max_diff > 1e-9:
                        print("      >>> R file disagrees with the recomputation: investigate before trusting either.")
                else:
                    print(f"      (no R output at {r_path} to compare)")
    pd.DataFrame(all_rows).to_csv(out / f"verify_{arm}_summary.csv", index=False)
    print(f"\nwrote {out.resolve()}/verify_{arm}_*.csv")


if __name__ == "__main__":
    main()