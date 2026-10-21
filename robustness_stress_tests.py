#!/usr/bin/env python3
"""
robustness_stress_tests.py -- Phase 12 adversarial pass on the frozen Phase 8 model.

Steps (all pre-specified in PHASE_11_12_DECISIONS.md before any output was seen):
  integrity   GATE. Rebuild the Phase 8 task tables from character_matrix.csv, rescore with the
              s0 routes, and require an exact match with phase8_*_task_table.csv and
              phase8_*_scored.csv. Nothing else runs if this fails.
  confounds   Is the score just lineage size / number of traits present / starting state?
  labels      Does the result survive a same-isolate outcome and a narrow hv/CRKP outcome?
  sampling    Size-matched outcomes: restrict to lineages with >=k outcome genomes; rarefy.
  relatives   Drop outcome genomes that have a close relative in TRAIN (Mash), recompute labels.
  qc          Re-derive TRAIN state and outcome from Kleborate-QC-passing genomes only.

Run from the project root:
  python3 robustness_stress_tests.py                       # everything that has its inputs
  python3 robustness_stress_tests.py --steps integrity,confounds
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import shared_lib as L

ORDER = L.CHAR_ORDER_P8          # Phase 7 and Phase 8 runs both use this column order


class Tee:
    # log to the console and to a file
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, *parts):
        msg = " ".join(str(x) for x in parts)
        print(msg)
        self.f.write(msg + "\n")
        self.f.flush()


# ---------------------------------------------------------------------------
# integrity gate
# ---------------------------------------------------------------------------
def integrity(cm, args, log):
    log("\n" + "=" * 78 + "\nSTEP integrity (gate)\n" + "=" * 78)
    routes = L.load_routes(args.routes_p8_s0)
    scorer = L.RouteScorer(routes)
    states = L.lineage_states(cm, L.TRAIN_END, ORDER)
    cov = L.lineage_covariates(cm, L.TRAIN_END, ORDER)
    log(f"Phase 8 s0 routes: {Path(args.routes_p8_s0).resolve()}  ({routes.shape[0]} routes)")
    tabs, windows = {}, {}
    for w, cands in [("validate", [L.VALIDATE_WINDOW]), ("test", L.TEST_WINDOW_CANDIDATES)]:
        orig = L.load_task_table(Path(args.tables_dir) / f"phase8_{w}_task_table.csv").set_index("base_ST")
        chosen = None
        for win in cands:
            rebuilt = L.build_task(cm, L.TRAIN_END, *win)
            same_set = set(rebuilt.index) == set(orig.index)
            common = rebuilt.index.intersection(orig.index)
            st_ok = (rebuilt.loc[common, "starting_state"] == orig.loc[common, "starting_state"]).all()
            lab_mis = int((rebuilt.loc[common, "label_union"].values != orig.loc[common, "became_dual_risk"].values).sum())
            log(f"  {w} window {win}: rebuilt {len(rebuilt)} vs table {len(orig)} lineages; same set={same_set}; "
                f"starting_state match={st_ok}; label mismatches={lab_mis}")
            if same_set and st_ok and lab_mis == 0:
                chosen = (win, rebuilt)
                break
        if chosen is None:
            L.die(f"integrity FAILED for {w}: task table could not be reproduced from character_matrix.csv under any "
                  f"candidate window {cands}. The label/window definition in the table differs from what this script "
                  "assumes -- fix that before interpreting ANY Phase 12 number.")
        windows[w] = chosen[0]
        tab = chosen[1].join(cov)
        sc, nv, no_state = L.score_task(tab, states, scorer, ORDER)
        tab["score"], tab["n_valid"] = sc, nv
        log(f"  {w}: window {chosen[0]} reproduces the table exactly ({len(tab)} lineages); NaN scores: {int(np.isnan(sc).sum())}")
        scored_path = Path(args.tables_dir) / f"phase8_{w}_scored.csv"
        if scored_path.exists():
            ref = L.read_csv_str(scored_path).set_index("base_ST")
            joined = tab.join(ref[["predicted_risk"]], how="inner")
            max_diff = float((joined["score"] - joined["predicted_risk"]).abs().max())
            log(f"  {w}: max |python score - phase8_{w}_scored.csv| = {max_diff:.2e} over {len(joined)} lineages")
            if not max_diff < 1e-9:
                L.die(f"integrity FAILED: rescoring {w} does not reproduce {scored_path}. Either the routes file, the character order, "
                      "or the state definition differs from what produced that CSV.")
        else:
            log(f"  (no {scored_path} to compare against -- scorer reproduction NOT verified against your R output)")
        log(f"  {w.upper()} AUROC (s0-only, primary): {L.fmt_auc(tab['label_union'].values, tab['score'].values)}")
        tab = tab.rename(columns={"label_union": "label"})
        tabs[w] = tab
    return dict(tabs=tabs, windows=windows, states=states, scorer=scorer, cov=cov)


# ---------------------------------------------------------------------------
# confounds
# ---------------------------------------------------------------------------
def _z(x):
    x = np.asarray(x, float)
    sd = x.std()
    return (x - x.mean()) / (sd if sd > 0 else 1.0)


def _fit_ll(X, y):
    from sklearn.linear_model import LogisticRegression
    model = LogisticRegression(C=1e6, max_iter=20000)
    model.fit(X, y)
    p = np.clip(model.predict_proba(X)[:, 1], 1e-12, 1 - 1e-12)
    return float((y * np.log(p) + (1 - y) * np.log(1 - p)).sum()), model


def confounds(ctx, args, out, log):
    log("\n" + "=" * 78 + "\nSTEP confounds\n" + "=" * 78)
    rows = []
    for w, tab in ctx["tabs"].items():
        tab = tab.dropna(subset=["score"]).copy()
        y = tab["label"].values
        score = tab["score"].values
        state = tab["starting_state"].values
        log(f"\n--- {w.upper()}  n={len(tab)}  pos={int(y.sum())}  base rate={y.mean():.3f} ---")
        # A. fair baselines (all knowable at prediction time) + one label-mechanism diagnostic
        base = {
            "log n_train genomes": np.log1p(tab["n_train"].values),
            "n traits present (of 20)": tab["n_present"].values.astype(float),
            "n countries in TRAIN": tab["n_countries"].values.astype(float),
            "n distinct TRAIN years": tab["n_train_years"].values.astype(float),
        }
        diag = {"[diagnostic, NOT knowable at prediction time] log n outcome genomes": np.log1p(tab["n_outcome"].values)}
        log("A. Model vs single-feature baselines (direction pre-specified: larger = riskier; no fitting)")
        log(f"   {'score':70s} AUROC [95% CI]           paired-DeLong diff(model-baseline)  p")
        auc0, lo0, hi0 = L.delong_ci(y, score)
        log(f"   {'HyperTraPS score (s0-only, primary)':70s} {auc0:.3f} [{lo0:.3f}, {hi0:.3f}]")
        rows.append(dict(window=w, section="A_baseline", item="model", auroc=auc0, lo=lo0, hi=hi0, n=len(y), pos=int(y.sum())))
        for name, vec in {**base, **diag}.items():
            auc_vec, lo_v, hi_v = L.delong_ci(y, vec)
            _, _, diff, pval = L.delong_paired(y, score, vec)
            log(f"   {name:70s} {auc_vec:.3f} [{lo_v:.3f}, {hi_v:.3f}]   diff={diff:+.3f}  p={pval:.3f}")
            rows.append(dict(window=w, section="A_baseline", item=name, auroc=auc_vec, lo=lo_v, hi=hi_v, n=len(y),
                             pos=int(y.sum()), diff_vs_model=diff, p_paired=pval))
        # B. correlations of the score with sampling covariates
        log("B. Spearman correlation of model score with covariates (overall / resistance_only / virulence_only)")
        for name, vec in base.items():
            r_all = stats.spearmanr(score, vec)[0]
            parts = []
            for grp in ("resistance_only", "virulence_only"):
                mask = state == grp
                parts.append(f"{stats.spearmanr(score[mask], vec[mask])[0]:+.2f}" if mask.sum() > 4 and np.std(score[mask]) > 0 and np.std(vec[mask]) > 0 else " n/a")
            log(f"   {name:40s} rho={r_all:+.2f}   [{parts[0]} / {parts[1]}]")
            rows.append(dict(window=w, section="B_spearman", item=name, auroc=r_all, n=len(y)))
        # C. stratified AUROC
        log("C. Stratified AUROC of the model score")
        for grp in ("resistance_only", "virulence_only"):
            mask = state == grp
            log(f"   starting_state={grp:16s} {L.fmt_auc(y[mask], score[mask])}   mean score {score[mask].mean():.3f}, conversion rate {y[mask].mean():.3f}")
            rows.append(dict(window=w, section="C_strata", item=grp, auroc=L.auroc(y[mask], score[mask]), n=int(mask.sum()), pos=int(y[mask].sum())))
        terc = pd.qcut(tab["n_train"].rank(method="first"), 3, labels=["small", "medium", "large"])
        for grp in ["small", "medium", "large"]:
            mask = (terc == grp).values
            log(f"   TRAIN-size tercile={grp:7s} (n_train {tab['n_train'][mask].min()}-{tab['n_train'][mask].max()})  {L.fmt_auc(y[mask], score[mask])}")
            rows.append(dict(window=w, section="C_strata", item=f"size_{grp}", auroc=L.auroc(y[mask], score[mask]), n=int(mask.sum()), pos=int(y[mask].sum())))
        # D. adjusted likelihood-ratio test
        within_pct = pd.Series(score).groupby(state).rank(pct=True).values
        X_base = np.column_stack([(state == "resistance_only").astype(float), _z(np.log1p(tab["n_train"])),
                                  _z(np.log1p(tab["n_outcome"])), _z(tab["n_present"])])
        ll_base, _ = _fit_ll(X_base, y)
        ll_full, model_full = _fit_ll(np.column_stack([X_base, within_pct]), y)
        chi2 = max(0.0, 2 * (ll_full - ll_base))
        p_lrt = float(stats.chi2.sf(chi2, 1))
        log("D. Does the model add information beyond starting_state + log n_train + log n_outcome + n_present?")
        log(f"   logistic LRT (1 df): chi2={chi2:.2f}  p={p_lrt:.3f}   coef on within-state score percentile={model_full.coef_[0][-1]:+.2f}")
        rows.append(dict(window=w, section="D_LRT", item="score beyond covariates", auroc=np.nan, p_paired=p_lrt, n=len(y), pos=int(y.sum())))
        # E. permutation + bootstrap
        obs, perm_p = L.perm_pvalue(y, score, state, n_perm=20000)
        boot_lo, boot_hi = L.bootstrap_ci(y, score, n_boot=5000)
        log(f"E. permutation p (labels shuffled within starting_state) = {perm_p:.4f};  bootstrap 95% CI [{boot_lo:.3f}, {boot_hi:.3f}]")
        rows.append(dict(window=w, section="E_perm", item="model", auroc=obs, lo=boot_lo, hi=boot_hi, p_paired=perm_p, n=len(y), pos=int(y.sum())))
    pd.DataFrame(rows).to_csv(out / "t1_confounds.csv", index=False)
    log(f"\nwrote {out / 't1_confounds.csv'}")


def _within(tab, lab, name, log):
    """Pooled AUROC mixes two starting states whose score scales aren't comparable;
    show each state on its own, plus a permutation p."""
    y, score, state = tab[lab].values, tab["score"].values, tab["starting_state"].values
    parts = []
    for grp in ("resistance_only", "virulence_only"):
        mask = state == grp
        parts.append(f"{grp}: {L.auroc(y[mask], score[mask]):.3f} (n={int(mask.sum())}, pos={int(y[mask].sum())}, mean score {score[mask].mean():.3f}, conv {y[mask].mean():.3f})")
    _, p = L.perm_pvalue(y, score, state, n_perm=10000)
    log(f"  {'':24s}  within-state -> " + " | ".join(parts) + f" | within-state permutation p={p:.3f}")


# ---------------------------------------------------------------------------
# label-definition sensitivity
# ---------------------------------------------------------------------------
def labels(cm, ctx, args, out, log):
    log("\n" + "=" * 78 + "\nSTEP labels\n" + "=" * 78)
    rows = []
    res_idx, vir_idx = L.class_indices(ORDER)
    hv_idx = [ORDER.index(c) for c in L.HV_MARKERS]
    carb_idx = [ORDER.index(L.CARB)]
    for w, tab in ctx["tabs"].items():
        win = ctx["windows"][w]
        scored = tab.dropna(subset=["score"])
        log(f"\n--- {w.upper()} window {win} ---")
        n_union = int(scored["label"].sum())
        n_iso = int(scored["label_isolate"].sum())
        log(f"Lineage-level ('union') positives: {n_union}; of these, lineages where a SINGLE isolate carries both classes: "
            f"{int(((scored['label'] == 1) & (scored['label_isolate'] == 1)).sum())}  "
            f"(=> {n_union - int(((scored['label'] == 1) & (scored['label_isolate'] == 1)).sum())} 'converted' only because different isolates carry the two classes)")
        log(f"  L0 union label       : {L.fmt_auc(scored['label'].values, scored['score'].values)}")
        log(f"  L1 same-isolate label: {L.fmt_auc(scored['label_isolate'].values, scored['score'].values)}   (same lineages, same scores)")
        _within(scored, "label_isolate", "L1", log)
        rows.append(dict(window=w, variant="L0_union_broad", n=len(scored), pos=n_union, auroc=L.auroc(scored['label'].values, scored['score'].values)))
        rows.append(dict(window=w, variant="L1_isolate_broad", n=len(scored), pos=n_iso, auroc=L.auroc(scored['label_isolate'].values, scored['score'].values)))
        # narrow hv/CRKP task: carbapenemase (Bla_Carb) x plasmid hv markers, frozen model, narrow target sets
        narrow = L.build_task(cm, L.TRAIN_END, *win, res_cols=[L.CARB], vir_cols=L.HV_MARKERS)
        targets = {"resistance_only": hv_idx, "virulence_only": carb_idx}   # 'resistance_only' here = carbapenemase-only
        sc, nv, _ = L.score_task(narrow, ctx["states"], ctx["scorer"], ORDER, targets=targets)
        narrow = narrow.assign(score=sc).join(ctx["cov"]).dropna(subset=["score"])
        log(f"  Narrow hv/CRKP task (carbapenemase [{L.CARB}] vs hv markers {L.HV_MARKERS}): {len(narrow)} lineages "
            f"({int((narrow['starting_state'] == 'resistance_only').sum())} carbapenemase-only, "
            f"{int((narrow['starting_state'] == 'virulence_only').sum())} hv-only at TRAIN)")
        for lab, nm in [("label_union", "N0 narrow union"), ("label_isolate", "N1 narrow same-isolate")]:
            yv = narrow[lab].values
            log(f"  {nm:24s}: {L.fmt_auc(yv, narrow['score'].values)}   base rate {yv.mean() if len(yv) else float('nan'):.3f}")
            size_auc = L.auroc(yv, np.log1p(narrow['n_train'].values))
            log(f"  {'':24s}  size-only baseline (log n_train) AUROC = {size_auc:.3f}")
            _within(narrow, lab, nm, log)
            rows.append(dict(window=w, variant=f"{nm}", n=len(narrow), pos=int(yv.sum()), auroc=L.auroc(yv, narrow['score'].values),
                             size_baseline=size_auc))
        if len(narrow) and narrow["label_isolate"].sum() < 5:
            log("  NOTE: <5 narrow same-isolate positives -- the narrow-definition numbers are descriptive only.")
    pd.DataFrame(rows).to_csv(out / "t2_label_variants.csv", index=False)
    log(f"\nwrote {out / 't2_label_variants.csv'}")


# ---------------------------------------------------------------------------
# sampling-effort matching
# ---------------------------------------------------------------------------
def sampling(cm, ctx, args, out, log, n_draw=500, seed=42):
    log("\n" + "=" * 78 + "\nSTEP sampling (size-matched outcomes)\n" + "=" * 78)
    rng = np.random.default_rng(seed)
    res_f, vir_f = L.genome_flags(cm, L.RESISTANCE, L.VIRULENCE)
    genomes = pd.DataFrame({"base_ST": cm["base_ST"].values, "res": res_f.values, "vir": vir_f.values, "yr": cm["collection_year"].values})
    rows = []
    for w, tab in ctx["tabs"].items():
        y0, y1 = ctx["windows"][w]
        tab = tab.dropna(subset=["score"])
        win = genomes[(genomes["yr"] >= y0) & (genomes["yr"] <= y1)]
        per_lineage = {st: (d["res"].values, d["vir"].values) for st, d in win.groupby("base_ST")}
        log(f"\n--- {w.upper()} ---  (outcome genomes per lineage: median {tab['n_outcome'].median():.0f}, "
            f"IQR {tab['n_outcome'].quantile(.25):.0f}-{tab['n_outcome'].quantile(.75):.0f}; "
            f"Spearman(n_outcome, label) = {stats.spearmanr(tab['n_outcome'], tab['label'])[0]:+.2f})")
        for k in (3, 5, 10):
            mask = (tab["n_outcome"] >= k).values
            if mask.sum() < 10 or tab["label"][mask].nunique() < 2:
                log(f"  restrict n_outcome>={k}: n={int(mask.sum())} -- too small / one class, skipped")
                continue
            log(f"  restrict n_outcome>={k}: {L.fmt_auc(tab['label'].values[mask], tab['score'].values[mask])}   base rate {tab['label'][mask].mean():.3f}")
            rows.append(dict(window=w, analysis=f"restrict_ge{k}", n=int(mask.sum()), pos=int(tab["label"][mask].sum()),
                             auroc=L.auroc(tab["label"].values[mask], tab["score"].values[mask])))
            # rarefy: every lineage contributes exactly k outcome genomes
            idx = np.nonzero(mask)[0]
            aucs, base_rates = [], []
            for _ in range(n_draw):
                lab = np.zeros(len(idx), int)
                for jj, i in enumerate(idx):
                    res_vec, vir_vec = per_lineage[tab.index[i]]
                    pick = rng.choice(len(res_vec), size=k, replace=False)
                    lab[jj] = int(res_vec[pick].any() and vir_vec[pick].any())
                if lab.sum() in (0, len(lab)):
                    continue
                aucs.append(L.auroc(lab, tab["score"].values[idx])); base_rates.append(lab.mean())
            if aucs:
                lo, hi = np.quantile(aucs, [.025, .975])
                log(f"     rarefied to exactly {k} genomes/lineage ({len(aucs)} draws): mean AUROC {np.mean(aucs):.3f} "
                    f"(draw-to-draw 2.5-97.5%: {lo:.3f}-{hi:.3f}); mean base rate {np.mean(base_rates):.3f}")
                rows.append(dict(window=w, analysis=f"rarefy_{k}", n=len(idx), pos=np.nan, auroc=float(np.mean(aucs)),
                                 lo=float(lo), hi=float(hi)))
    pd.DataFrame(rows).to_csv(out / "t5_sampling.csv", index=False)
    log(f"\nwrote {out / 't5_sampling.csv'}")


# ---------------------------------------------------------------------------
# close-relative removal (Mash)
# ---------------------------------------------------------------------------
def mash_min_dist_to_train(path, cm, train_end, log):
    """Stream a `mash triangle` lower-triangular file; return the min Mash distance from every
    matrix genome to any OTHER matrix genome collected <= train_end."""
    cm = cm.dropna(subset=["collection_year"])
    gid_year = dict(zip(cm["genome_id"], cm["collection_year"]))
    pos_of = {}
    names = []
    t0 = time.time()
    with open(path) as f:
        first = f.readline().strip()
        if not first.isdigit():
            L.die(f"{path}: first line should be the genome count from `mash triangle`, got {first[:60]!r}")
        n_total = int(first)
        # single pass over the file, arrays sized up front
        member = np.zeros(n_total, bool)
        train = np.zeros(n_total, bool)
        mind = np.full(n_total, np.inf, dtype=np.float64)
        gids = [None] * n_total
        for i, line in enumerate(f):
            line = line.rstrip("\r\n")
            name, _, rest = line.partition("\t")
            acc = L.accession_from_path(name)
            gids[i] = acc
            if acc in gid_year:
                member[i] = True
                train[i] = gid_year[acc] <= train_end
            if i == 0 or not rest:
                continue
            try:
                d = np.fromstring(rest, dtype=np.float64, sep="\t")
            except Exception:
                d = np.array(rest.split("\t"), dtype=np.float64)
            if d.size != i:
                L.die(f"{path}: row {i} has {d.size} distances, expected {i} (not a lower-triangular `mash triangle` file?)")
            if member[i]:
                train_j = train[:i]
                if train_j.any():
                    mind[i] = min(mind[i], d[train_j].min())
            if train[i]:
                cand = np.where(member[:i], d, np.inf)
                np.minimum(mind[:i], cand, out=mind[:i])
            if i % 2000 == 0:
                print(f"   parsed {i}/{n_total} rows ({time.time() - t0:.0f}s)", flush=True)
    found = set(g for g in gids if g)
    miss = [g for g in cm["genome_id"] if g not in found]
    log(f"  Mash file: {n_total} genomes; matrix genomes found in Mash names: {len(cm) - len(miss)}/{len(cm)}")
    if miss:
        L.die(f"{len(miss)} matrix genome_ids were not found in the Mash names (e.g. {miss[:3]}). "
              "genome_id must be the assembly accession that appears in the fasta path used for `mash sketch`.")
    return pd.Series(mind[member], index=[gids[i] for i in np.nonzero(member)[0]])


def relatives(cm, ctx, args, out, log):
    log("\n" + "=" * 78 + "\nSTEP relatives (Mash close-relative removal)\n" + "=" * 78)
    if "genome_id" not in cm.columns:
        L.die("character_matrix.csv has no genome_id column; cannot map to Mash names.")
    cache = out / "min_dist_to_train.csv"
    if cache.exists() and not args.recompute_mash:
        min_dist = pd.read_csv(cache, index_col=0).iloc[:, 0]
        log(f"  using cached {cache} ({len(min_dist)} genomes); pass --recompute-mash to redo")
    else:
        if not Path(args.mash).exists():
            L.die(f"Mash distance file not found: {args.mash} (use --mash)")
        min_dist = mash_min_dist_to_train(args.mash, cm, L.TRAIN_END, log)
        min_dist.rename("min_dist_to_train").to_csv(cache)
    res_f, vir_f = L.genome_flags(cm, L.RESISTANCE, L.VIRULENCE)
    genomes = pd.DataFrame({"gid": cm["genome_id"].values, "base_ST": cm["base_ST"].values, "res": res_f.values,
                            "vir": vir_f.values, "yr": cm["collection_year"].values})
    genomes["md"] = genomes["gid"].map(min_dist)
    rows = []
    for w, tab in ctx["tabs"].items():
        y0, y1 = ctx["windows"][w]
        tab = tab.dropna(subset=["score"])
        win = genomes[(genomes["yr"] >= y0) & (genomes["yr"] <= y1) & genomes["base_ST"].isin(tab.index)]
        log(f"\n--- {w.upper()} ---  outcome genomes in task lineages: {len(win)}; "
            f"median min-distance to a TRAIN genome = {win['md'].median():.4f}; "
            f"share within 0.001 / 0.003 of a TRAIN genome: "
            f"{(win['md'] <= .001).mean():.2f} / {(win['md'] <= .003).mean():.2f}")
        quantiles = win["md"].quantile([.1, .25, .5, .75, .9]).values
        log("   outcome-genome min-distance quantiles (10/25/50/75/90%): " + " / ".join(f"{x:.4f}" for x in quantiles))
        log(f"   unfiltered: {L.fmt_auc(tab['label'].values, tab['score'].values)}")
        for thr in args.mash_thresholds:
            keep = win[win["md"] > thr]
            agg = keep.groupby("base_ST").agg(wr=("res", "max"), wv=("vir", "max"), n=("res", "size"))
            lab = (agg["wr"].astype(bool) & agg["wv"].astype(bool)).astype(int)
            joined = tab.join(lab.rename("label_novel"), how="inner")
            log(f"   drop outcome genomes with a TRAIN relative within Mash {thr}: {len(keep)} genomes / {len(joined)} lineages remain -> "
                f"{L.fmt_auc(joined['label_novel'].values, joined['score'].values) if joined['label_novel'].nunique() == 2 else 'one class only'}   "
                f"base rate {joined['label_novel'].mean() if len(joined) else float('nan'):.3f}")
            rows.append(dict(window=w, mash_threshold=thr, n_lineages=len(joined), pos=int(joined["label_novel"].sum()) if len(joined) else 0,
                             auroc=L.auroc(joined["label_novel"].values, joined["score"].values) if len(joined) else np.nan))
    pd.DataFrame(rows).to_csv(out / "t3_relatives.csv", index=False)
    log(f"\nwrote {out / 't3_relatives.csv'}")


# ---------------------------------------------------------------------------
# QC filter
# ---------------------------------------------------------------------------
def _parse_amb(x):
    if pd.isna(x):
        return np.nan
    text = str(x).strip().lower()
    if text.startswith("no"):
        return 0.0
    if text.startswith("yes"):
        import re
        m = re.search(r"\d+", text)
        return float(m.group()) if m else 1.0
    try:
        return float(text)
    except ValueError:
        return np.nan


def qc(cm, ctx, args, out, log):
    log("\n" + "=" * 78 + "\nSTEP qc (assembly QC filter)\n" + "=" * 78)
    if "genome_id" not in cm.columns:
        L.die("character_matrix.csv has no genome_id column; cannot join QC metadata.")
    meta = L.read_csv_str(args.metadata, sep="\t") if str(args.metadata).endswith((".txt", ".tsv")) else L.read_csv_str(args.metadata)
    key = next((k for k in ("assembly_accession", "genome_id") if k in meta.columns), None)
    if key is None and "strain" in meta.columns:          # raw Kleborate output: strain = fasta basename
        meta["assembly_accession"] = meta["strain"].astype(str).map(L.accession_from_path)
        key = "assembly_accession"
    have_raw = all(c in meta.columns for c in ("N50", "ambiguous_bases", "total_size"))
    have_kleb = all(c in meta.columns for c in ("species_match", "QC_warnings"))
    if key is None or not (have_raw or have_kleb):
        L.die(f"{args.metadata}: need a genome key (assembly_accession / genome_id / strain) plus either "
              "N50+ambiguous_bases+total_size (robustness_assembly_qc.py output) or species_match+QC_warnings (Kleborate v3). "
              f"Actual columns: {list(meta.columns)}")
    meta = meta.drop_duplicates(key).set_index(key)
    variants = {}
    joined = None
    if have_raw:
        n50 = pd.to_numeric(cm["genome_id"].map(meta["N50"]), errors="coerce")
        tot = pd.to_numeric(cm["genome_id"].map(meta["total_size"]), errors="coerce")
        amb = cm["genome_id"].map(meta["ambiguous_bases"]).map(_parse_amb)
        joined = n50.notna()
        ok50 = (n50 >= 10000).fillna(False)
        clean = ok50 & (amb == 0)
        variants["permissive (N50 >= 10 kb)"] = ok50
        variants["strict (N50 >= 10 kb AND zero ambiguous bases; Kleborate's own warning rules)"] = clean
        variants["strict + total size 4.8-6.6 Mb (size window is an arbitrary plausibility range)"] = clean & tot.between(4.8e6, 6.6e6)
        # Added after looking at the QC distribution only (not at any AUROC): many scaffold-level assemblies
        # carry N gaps, so "zero ambiguous bases" removes most genomes. Moderate = N fraction <= 0.5% of the assembly.
        variants["moderate (N50 >= 10 kb, ambiguous <= 0.5% of assembly, size 4.8-6.6 Mb)"] = (
            ok50 & ((amb / tot) <= 0.005).fillna(False) & tot.between(4.8e6, 6.6e6))
    if have_kleb:
        species_match = cm["genome_id"].map(meta["species_match"])
        qc_warn = cm["genome_id"].map(meta["QC_warnings"])
        joined = species_match.notna() if joined is None else joined
        no_warn = qc_warn.fillna("-").astype(str).str.strip().isin(["-", "", "nan"])
        strong = species_match.astype(str).str.lower().eq("strong")
        variants["Kleborate strict (species_match=strong AND no QC_warnings)"] = strong & no_warn
    log(f"  QC metadata joined for {int(joined.sum())}/{len(cm)} matrix genomes ({Path(args.metadata).resolve()})")
    if joined.sum() < 0.9 * len(cm):
        L.die("fewer than 90% of matrix genomes matched the QC metadata -- wrong key column or wrong file.")
    rows = []
    for vname, keepmask in variants.items():
        cmq = cm[keepmask.values]
        log(f"\n--- QC variant: {vname}: keeps {len(cmq)}/{len(cm)} genomes ({len(cmq) / len(cm):.1%}) ---")
        states_q = L.lineage_states(cmq, L.TRAIN_END, ORDER)
        for wname, tab in ctx["tabs"].items():
            win = ctx["windows"][wname]
            tab_q = L.build_task(cmq, L.TRAIN_END, *win)
            sc, nv, _ = L.score_task(tab_q, states_q, ctx["scorer"], ORDER)
            tab_q = tab_q.assign(score=sc).dropna(subset=["score"])
            overlap = len(set(tab_q.index) & set(tab.index))
            log(f"  {wname.upper():8s} original {len(tab)} lineages -> QC-derived {len(tab_q)} (overlap {overlap}); "
                f"{L.fmt_auc(tab_q['label_union'].values, tab_q['score'].values)}   base rate {tab_q['label_union'].mean():.3f}")
            rows.append(dict(window=wname, variant=vname, genomes_kept=len(cmq), n=len(tab_q), pos=int(tab_q["label_union"].sum()),
                             overlap=overlap, auroc=L.auroc(tab_q["label_union"].values, tab_q["score"].values)))
    pd.DataFrame(rows).to_csv(out / "t4_qc.csv", index=False)
    log(f"\nwrote {out / 't4_qc.csv'}")
    log("NOTE: the frozen model was trained on ALL TRAIN genomes (including any that fail QC); only the state/label\n"
        "definitions are QC-filtered here. If results move materially, the fully rigorous version re-trains on QC-passing genomes.")


# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", default="integrity,confounds,labels,sampling,relatives,qc")
    parser.add_argument("--matrix", default="phase3_out/character_matrix.csv")
    parser.add_argument("--tables-dir", default="phase8_out")
    parser.add_argument("--routes-p8-s0", default="phase8_out/kp_train_only_model1_s0-routes.txt")
    parser.add_argument("--mash", default="mash_test/distances_full.tsv")
    parser.add_argument("--mash-thresholds", type=float, nargs="+", default=[0.001, 0.002, 0.003, 0.004])
    parser.add_argument("--recompute-mash", action="store_true")
    parser.add_argument("--metadata", default="kp_pull_20k/kp_metadata_FULL_dedup.csv")
    parser.add_argument("--out", default="phase12_out")
    args = parser.parse_args()

    out = L.ensure_dir(args.out)
    log = Tee(out / "phase12_log.txt")
    steps = [s.strip() for s in args.steps.split(",")]
    explicit = args.steps != parser.get_default("steps")
    cm = L.load_matrix(args.matrix, need_genome_id=False)
    log(f"matrix {Path(args.matrix).resolve()}: {len(cm)} genomes, {cm['base_ST'].nunique()} lineages")
    ctx = integrity(cm, args, log)          # always runs, it's the gate
    if "confounds" in steps:
        confounds(ctx, args, out, log)
    if "labels" in steps:
        labels(cm, ctx, args, out, log)
    if "sampling" in steps:
        sampling(cm, ctx, args, out, log)
    for name, fn, inp in [("relatives", relatives, args.mash), ("qc", qc, args.metadata)]:
        if name in steps:
            if not Path(inp).exists() and not explicit:
                log(f"\n(skipping {name}: {inp} not found; pass the path explicitly to run it)")
                continue
            if "genome_id" not in cm.columns:
                log(f"\n(skipping {name}: matrix has no genome_id column)")
                continue
            fn(cm, ctx, args, out, log)
    log("\nDONE. Everything above was pre-specified; report ALL of it, including anything unflattering.")


if __name__ == "__main__":
    main()