#!/usr/bin/env python3
"""
risk_current_lineage_shortlist.py -- Phase 11: currently at-risk, non-convergent lineages.

Rules fixed before any output was seen (see PHASE_11_12_DECISIONS.md):
  model         Phase 7 full-data --model 1 runs, 3 seeds pooled (2,370 routes); s0-only is a sensitivity column
  state         union of traits over ALL genomes in the matrix (through --state-through, default = last year present)
  population    lineages that are single-risk (resistance-only or virulence-only) AND observed in/after --recent-from (2021)
  score         loose Phase-8 truncation score: P(next trait acquired is in the complementary class)
  High tier     top tercile of score WITHIN starting_state (ties included; needs >=6 lineages and non-constant scores)
  size-adjusted top tercile within (starting_state x TRAIN-size half). Priority shortlist = High under BOTH rules
  evidence      the same tier rules are applied to the Phase 8 VALIDATE/TEST tables and the observed conversion rate
                in High vs not-High is reported (the ONLY performance evidence attached to the list)
  language      "genomic profile associated with elevated predicted probability of subsequent acquisition", never
                "will become". Scores are discrimination scores, NOT calibrated probabilities

Post-hoc amendment (made after seeing Phase 12 output, and labelled as such): only resistance_only
lineages are shortlisted.

Run from the project root:
  python3 risk_current_lineage_shortlist.py --p7-dir .
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import shared_lib as L
from robustness_stress_tests import integrity, Tee

ORDER = L.CHAR_ORDER_P8


def tier(scores: pd.Series, groups: pd.Series, min_group=6, q=2 / 3) -> pd.Series:
    out = pd.Series("n/a", index=scores.index, dtype=object)
    for grp in groups.dropna().unique():
        idx = groups.index[groups == grp]
        vals = scores.loc[idx].dropna()
        if len(vals) < min_group or vals.max() == vals.min():
            continue
        thr = vals.quantile(q)
        high = (vals >= thr) & (vals > vals.min())
        out.loc[vals.index] = np.where(high, "High", "Not high")
    return out


def size_half(n, states):
    # 'small'/'large' median split of n within each starting_state
    half = pd.Series("large", index=n.index, dtype=object)
    for state in states.unique():
        in_state = states == state
        half[in_state & (n < n[in_state].median())] = "small"
    return half


def next_trait_strings(scorer, present, target_idx, k=3):
    nxt = scorer.next_features(present)
    if nxt.size == 0:
        return "", ""
    freq = np.bincount(nxt, minlength=len(ORDER)) / nxt.size
    top_any = np.argsort(-freq)[:k]
    top_target = [i for i in np.argsort(-freq) if i in set(target_idx)][:k]
    fmt = lambda idxs: "; ".join(f"{ORDER[i]}:{freq[i]:.2f}" for i in idxs if freq[i] > 0)
    return fmt(top_any), fmt(top_target)


def analogue(ctx, log, out):
    log("\n" + "=" * 78 + "\nRETROSPECTIVE EVIDENCE: same tier rules applied to Phase 8 (frozen TRAIN<=2017 model, s0)\n" + "=" * 78)
    rows, lin_rows = [], []
    for win, tab in ctx["tabs"].items():
        tab = tab.dropna(subset=["score"]).copy()
        tab["tier_raw"] = tier(tab["score"], tab["starting_state"])
        grp = tab["starting_state"] + "|" + size_half(tab["n_train"], tab["starting_state"])
        tab["tier_size"] = tier(tab["score"], grp)
        tab["priority"] = np.where((tab["tier_raw"] == "High") & (tab["tier_size"] == "High"), "High", "Not high")
        for subset, mask in [("all single-risk", np.ones(len(tab), bool)),
                             ("resistance_only only", (tab["starting_state"] == "resistance_only").values)]:
            sub = tab[mask]
            log(f"\n--- {win.upper()} / {subset} (n={len(sub)}, base rate {sub['label'].mean():.3f}) ---")
            for col, name in [("tier_raw", "raw tercile"), ("tier_size", "size-adjusted tercile"), ("priority", "High under BOTH")]:
                high = sub[sub[col] == "High"]
                rest = sub[sub[col] != "High"]
                a, b = int(high["label"].sum()), len(high) - int(high["label"].sum())
                c, d = int(rest["label"].sum()), len(rest) - int(rest["label"].sum())
                rate_high = a / max(len(high), 1)
                rate_rest = c / max(len(rest), 1)
                pval = stats.fisher_exact([[a, b], [c, d]], alternative="greater")[1] if len(high) and len(rest) else np.nan
                lift = rate_high / sub["label"].mean() if sub["label"].mean() > 0 else np.nan
                log(f"  {name:24s}: High n={len(high):3d} converted {a:3d} ({rate_high:.2f})  | not-High n={len(rest):3d} converted {c:3d} ({rate_rest:.2f})"
                    f"  | lift vs base {lift:.2f}  Fisher one-sided p={pval:.3f}")
                rows.append(dict(window=win, subset=subset, rule=name, n_high=len(high), conv_high=a, rate_high=rate_high, n_not=len(rest),
                                 conv_not=c, rate_not=rate_rest, lift=lift, fisher_p=pval))
        keep = tab.reset_index()[["base_ST", "starting_state", "n_train", "score", "tier_raw", "tier_size", "priority", "label", "label_isolate"]]
        keep.insert(0, "window", win)
        lin_rows.append(keep)
    pd.DataFrame(rows).to_csv(out / "phase11_analogue_summary.csv", index=False)
    pd.concat(lin_rows).to_csv(out / "phase11_analogue_lineages.csv", index=False)
    log("\n'High and subsequently converged' (protocol category 4) = rows with tier High and label=1 in phase11_analogue_lineages.csv")


def directional(routes_by_seed, log, out):
    log("\n" + "=" * 78 + "\nDIRECTIONAL ROUTES (Phase 7 two-route comparison; post-hoc on -routes.txt, no new compute)\n" + "=" * 78)
    carb = ORDER.index(L.CARB)
    hv = [ORDER.index(c) for c in L.HV_MARKERS]
    rows = []
    for name, routes in list(routes_by_seed.items()) + [("pooled", np.vstack(list(routes_by_seed.values())))]:
        scorer = L.RouteScorer(routes)
        first_hv = scorer.pos[:, hv].min(axis=1)
        p_hv_first = float((first_hv < scorer.pos[:, carb]).mean())
        rows.append(dict(seed=name, n_routes=len(routes), P_hvmarker_before_carbapenemase=p_hv_first,
                         P_carbapenemase_before_hvmarker=1 - p_hv_first))
        log(f"  {name:7s} n={len(routes):5d}  P(first hv marker precedes Bla_Carb) = {p_hv_first:.3f}   P(Bla_Carb precedes all hv markers) = {1 - p_hv_first:.3f}")
    pd.DataFrame(rows).to_csv(out / "phase11_directional_routes.csv", index=False)
    pooled = np.vstack(list(routes_by_seed.values()))
    scorer = L.RouteScorer(pooled)
    # pairwise P(virulence locus before resistance class)
    pair = pd.DataFrame({L.RESISTANCE[j]: [float((scorer.pos[:, ORDER.index(v)] < scorer.pos[:, ORDER.index(L.RESISTANCE[j])]).mean())
                                           for v in L.VIRULENCE] for j in range(len(L.RESISTANCE))}, index=L.VIRULENCE)
    pair.to_csv(out / "phase11_pairwise_virulence_before_resistance.csv")
    mean_pos = pd.DataFrame({"trait": ORDER, "mean_route_position(0=earliest)": scorer.pos.mean(axis=0),
                             "sd": scorer.pos.std(axis=0)}).sort_values("mean_route_position(0=earliest)")
    mean_pos.to_csv(out / "phase11_mean_positions.csv", index=False)
    log("  mean route position (pooled, 0 = earliest):  " + ", ".join(f"{r.trait}:{r._2:.1f}" for r in mean_pos.rename(columns={"mean_route_position(0=earliest)": "_2"}).itertuples()))
    log("  CAVEAT: HyperTraPS orderings assume irreversible accumulation and reflect cross-sectional structure; the direction\n"
        "  statistic describes the model-inferred ordering, not observed within-lineage history or causality.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default="phase3_out/character_matrix.csv")
    parser.add_argument("--tables-dir", default="phase8_out")
    parser.add_argument("--routes-p8-s0", default="phase8_out/kp_train_only_model1_s0-routes.txt")
    parser.add_argument("--p7-dir", default=".")
    parser.add_argument("--p7-pattern", default="kp_full_run_model1_s{s}-routes.txt")
    parser.add_argument("--seeds", default="0,1,2")
    parser.add_argument("--state-through", type=int, default=None)
    parser.add_argument("--recent-from", type=int, default=2021)
    parser.add_argument("--out", default="phase11_out")
    args = parser.parse_args()

    out = L.ensure_dir(args.out)
    log = Tee(out / "phase11_log.txt")
    cm = L.load_matrix(args.matrix)
    last_year = int(cm["collection_year"].max())
    through = args.state_through or last_year
    log(f"matrix {Path(args.matrix).resolve()}: {len(cm)} genomes; state computed through {through}; 'recent' = last seen >= {args.recent_from}")

    # 0. integrity gate + retrospective evidence (needs the Phase 8 pieces)
    ctx = integrity(cm, args, log)
    analogue(ctx, log, out)

    # 1. Phase 7 routes
    seeds = [int(s) for s in args.seeds.split(",")]
    paths = [Path(args.p7_dir) / args.p7_pattern.format(s=s) for s in seeds]
    routes_by_seed = {f"s{s}": L.load_routes(p) for s, p in zip(seeds, paths)}
    for p in paths:
        log(f"Phase 7 routes: {p.resolve()}")
    pooled = np.vstack(list(routes_by_seed.values()))
    scorer_pool, scorer_s0 = L.RouteScorer(pooled), L.RouteScorer(routes_by_seed[f"s{seeds[0]}"])

    # 2. current lineage table
    sub = cm[cm["collection_year"] <= through]
    states = L.lineage_states(cm, through, ORDER)
    res_any = states[L.RESISTANCE].sum(axis=1) > 0
    vir_any = states[L.VIRULENCE].sum(axis=1) > 0
    res_flag, vir_flag = L.genome_flags(sub, L.RESISTANCE, L.VIRULENCE)
    same_isolate = pd.Series((res_flag & vir_flag).values, index=sub.index).groupby(sub["base_ST"]).any()
    lin = pd.DataFrame(index=states.index)
    lin["category"] = np.select([res_any & vir_any, res_any ^ vir_any], ["already_dual_risk", "single_risk"], "neither")
    lin["starting_state"] = np.where(res_any & ~vir_any, "resistance_only", np.where(vir_any & ~res_any, "virulence_only", ""))
    lin["convergence_evidence"] = np.where(same_isolate.reindex(lin.index).fillna(False), "same-isolate",
                                           np.where(lin["category"] == "already_dual_risk", "co-circulating (different isolates)", "none"))
    by_st = sub.groupby("base_ST")
    lin["n_genomes"] = by_st.size()
    lin["last_seen_year"] = by_st["collection_year"].max()
    lin["n_recent_genomes"] = sub[sub["collection_year"] >= args.recent_from].groupby("base_ST").size().reindex(lin.index).fillna(0).astype(int)
    lin["n_countries"] = by_st["country"].nunique()
    top_country = by_st["country"].agg(lambda x: x.value_counts().index[0] if x.notna().any() else "Unknown")
    top_share = by_st["country"].agg(lambda x: x.value_counts(normalize=True).iloc[0] if x.notna().any() else np.nan)
    lin["top_country"], lin["top_country_share"] = top_country, top_share
    lin["n_traits_present"] = states.sum(axis=1)
    lin["traits_present"] = states.apply(lambda r: ";".join(c for c in ORDER if r[c] == 1), axis=1)
    lin["recent"] = lin["last_seen_year"] >= args.recent_from
    log(f"\nlineages: {len(lin)} total | " + " | ".join(f"{k}: {v}" for k, v in lin["category"].value_counts().items())
        + f" | recent (last seen >= {args.recent_from}): {int(lin['recent'].sum())}")
    log(f"single-risk AND recent (the scored population): {int(((lin['category'] == 'single_risk') & lin['recent']).sum())}")

    # 3. score the population
    pop = lin[(lin["category"] == "single_risk") & lin["recent"]].copy()
    res_idx, vir_idx = L.class_indices(ORDER)
    scores_pooled, scores_s0, n_valid, next_any, next_target = [], [], [], [], []
    for st, start_state in zip(pop.index, pop["starting_state"]):
        present = np.nonzero(states.loc[st].values == 1)[0]
        target = vir_idx if start_state == "resistance_only" else res_idx
        s, nv = scorer_pool.score(present, target)
        s0, _ = scorer_s0.score(present, target)
        scores_pooled.append(s); scores_s0.append(s0); n_valid.append(nv)
        ta, tt = next_trait_strings(scorer_pool, present, target)
        next_any.append(ta); next_target.append(tt)
    pop["score_pooled"], pop["score_s0only"], pop["n_valid_routes"] = scores_pooled, scores_s0, n_valid
    pop["top3_next_trait_any"], pop["top3_next_trait_target_class"] = next_any, next_target
    ok = pop["score_pooled"].notna() & pop["score_s0only"].notna()
    if ok.sum() > 3:
        log(f"Spearman(pooled score, s0-only score) = {stats.spearmanr(pop.loc[ok, 'score_pooled'], pop.loc[ok, 'score_s0only'])[0]:.3f}")
    log(f"unscorable lineages (no valid route): {int(pop['score_pooled'].isna().sum())}")

    # 4. tiers
    pop["tier_raw"] = tier(pop["score_pooled"], pop["starting_state"])
    pop["size_half"] = size_half(pop["n_genomes"], pop["starting_state"])
    pop["tier_size_adjusted"] = tier(pop["score_pooled"], pop["starting_state"] + "|" + pop["size_half"])
    pop["high_under_both"] = (pop["tier_raw"] == "High") & (pop["tier_size_adjusted"] == "High")
    # Post-hoc amendment (after Phase 12): within virulence_only lineages the frozen model's AUROC was
    # 0.345 (VALIDATE, n=19) and 0.440 (TEST, n=13), i.e. no demonstrated skill, and for hv-only lineages
    # predicting carbapenemase acquisition it was 0.254. Their higher scores are a between-state scale
    # artefact. They are still tiered and listed, but NOT shortlisted.
    pop["shortlist_eligible"] = pop["starting_state"] == "resistance_only"
    pop["priority_shortlist"] = pop["high_under_both"] & pop["shortlist_eligible"]
    pop["thin_evidence"] = pop["n_genomes"] < 5
    pop["single_country"] = pop["n_countries"] <= 1
    log("\nTier counts (scored population):")
    log(pd.crosstab(pop["tier_raw"], pop["tier_size_adjusted"], margins=True).to_string())
    log(f"High under both rules: {int(pop['high_under_both'].sum())} lineages, of which "
        f"{int((pop['high_under_both'] & ~pop['shortlist_eligible']).sum())} are virulence_only (model has no demonstrated within-state skill there; NOT shortlisted)")
    log(f"PRIORITY SHORTLIST (resistance_only AND High under both): {int(pop['priority_shortlist'].sum())} lineages "
        f"({int((pop['priority_shortlist'] & pop['thin_evidence']).sum())} with <5 genomes; "
        f"{int((pop['priority_shortlist'] & ~pop['thin_evidence']).sum())} with >=5 genomes)")

    # 5. outputs
    lin = lin.join(pop[["score_pooled", "score_s0only", "n_valid_routes", "top3_next_trait_any", "top3_next_trait_target_class",
                        "tier_raw", "size_half", "tier_size_adjusted", "high_under_both", "shortlist_eligible", "priority_shortlist", "thin_evidence", "single_country"]])
    lin.reset_index(names="base_ST").to_csv(out / "phase11_current_lineages.csv", index=False)
    shortlist = pop[pop["priority_shortlist"]].sort_values("score_pooled", ascending=False).reset_index(names="base_ST")
    shortlist.to_csv(out / "phase11_priority_shortlist.csv", index=False)
    shortlist[~shortlist["thin_evidence"]].to_csv(out / "phase11_shortlist_well_sampled.csv", index=False)
    by_profile = (shortlist.groupby("traits_present")
                  .agg(n_lineages=("base_ST", "size"), score=("score_pooled", "first"), n_genomes_total=("n_genomes", "sum"),
                       lineages=("base_ST", lambda x: ";".join(x)))
                  .sort_values(["score", "n_lineages"], ascending=False).reset_index())
    by_profile.to_csv(out / "phase11_shortlist_by_profile.csv", index=False)
    log(f"\nShortlist collapses to {len(by_profile)} distinct trait profiles (score depends only on the profile; lineages sharing one are tied)")
    log("\nTop of the priority shortlist:")
    cols = ["base_ST", "starting_state", "score_pooled", "n_genomes", "n_recent_genomes", "top_country", "top3_next_trait_target_class"]
    log(shortlist[cols].head(15).to_string(index=False) if len(shortlist) else "  (empty)")

    directional(routes_by_seed, log, out)
    log(f"\nwrote {out.resolve()}/phase11_*.csv")
    log("\nREMINDER: report as 'genomic profile associated with elevated predicted probability of subsequent acquisition of the "
        "complementary trait'. Attach the retrospective High-vs-not-High table above, and the Phase 12 confound results.")


if __name__ == "__main__":
    main()