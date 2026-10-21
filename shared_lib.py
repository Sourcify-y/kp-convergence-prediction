#!/usr/bin/env python3
"""
shared_lib.py -- shared helpers for the Phase 11 / Phase 12 scripts.

Design rules (from this project's bug history):
  * Every input file is validated on load and fails loudly with the actual
    column list, never a silent wrong-column result (earlier bugs: Mash CRLF,
    .isdigit(), ASR tie-break, stale files).
  * The route scorer is a vectorised port of the loose method in
    eval_score_lineages.R / robustness_score_holdouts_pooled.R. It was checked
    against a literal line-by-line port (reference_score_lineage; the tests used
    during development are not included in this repository), and the `integrity`
    step in robustness_stress_tests.py re-checks it against the real
    phase8_*_scored.csv before anything else runs.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# --------------------------------------------------------------------------
# character sets / orders (copied from the R scripts, not re-derived)
# --------------------------------------------------------------------------
# order used by the Phase 7 full run, the Phase 8 train-only run and the 20-char build
CHAR_ORDER_P8 = [
    "Bla_ESBL_acquired", "Bla_Carb_acquired", "Flq_acquired", "Sul_acquired", "AGly_acquired",
    "Col_acquired", "Tet_acquired", "MLS_acquired", "Rif_acquired", "Tmt_acquired",
    "Fcyn_acquired", "Phe_acquired", "Yersiniabactin", "Aerobactin", "Salmochelin",
    "Colibactin", "rmpA", "rmpA2", "Bla_acquired", "Tgc_acquired",
]
# character_matrix.csv header order; used for the Phase 9 lineage/geo holdout runs
CHAR_ORDER_HDR = [
    "AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired",
    "MLS_acquired", "Phe_acquired", "Rif_acquired", "Sul_acquired",
    "Tet_acquired", "Tgc_acquired", "Tmt_acquired", "Bla_acquired",
    "Bla_ESBL_acquired", "Bla_Carb_acquired", "Yersiniabactin",
    "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2",
]
RESISTANCE = [
    "Bla_ESBL_acquired", "Bla_Carb_acquired", "Flq_acquired", "Sul_acquired",
    "AGly_acquired", "Col_acquired", "Tet_acquired", "MLS_acquired",
    "Rif_acquired", "Tmt_acquired", "Fcyn_acquired", "Phe_acquired",
    "Bla_acquired", "Tgc_acquired",
]
VIRULENCE = ["Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin", "rmpA", "rmpA2"]
# Plasmid-borne hypervirulence markers available in Kleborate output
# (Russo et al. 2018, J Clin Microbiol: iucA, iroB, prmpA, prmpA2; peg-344 isn't in Kleborate).
# Yersiniabactin (ICEKp) and Colibactin (chromosomal island) are deliberately excluded.
HV_MARKERS = ["Aerobactin", "Salmochelin", "rmpA", "rmpA2"]
CARB = "Bla_Carb_acquired"

assert set(CHAR_ORDER_P8) == set(CHAR_ORDER_HDR) == set(RESISTANCE + VIRULENCE)
assert len(RESISTANCE) == 14 and len(VIRULENCE) == 6

TRAIN_END = 2017
VALIDATE_WINDOW = (2018, 2020)
TEST_WINDOW_CANDIDATES = [(2021, 2022), (2018, 2022)]   # the integrity step picks the one that reproduces


# --------------------------------------------------------------------------
# loading with validation
# --------------------------------------------------------------------------
def die(msg: str):
    print("\nERROR: " + msg, file=sys.stderr)
    raise SystemExit(1)


def read_csv_str(path, **kw) -> pd.DataFrame:
    """read_csv that never lets base_ST / genome_id become floats ('11.0')."""
    p = Path(path)
    if not p.exists():
        die(f"file not found: {p.resolve()}")
    dtype = {"base_ST": str, "genome_id": str, "assembly_accession": str}
    df = pd.read_csv(p, dtype=dtype, low_memory=False, **kw)
    if "base_ST" in df.columns:
        df["base_ST"] = df["base_ST"].astype(str).str.strip()
    return df


def load_matrix(path, need_genome_id=False) -> pd.DataFrame:
    df = read_csv_str(path)
    need = RESISTANCE + VIRULENCE + ["base_ST", "collection_year", "country"]
    if need_genome_id:
        need = need + ["genome_id"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        die(f"{path}: missing columns {missing}\nActual columns: {list(df.columns)}")
    df["collection_year"] = pd.to_numeric(df["collection_year"], errors="coerce")
    for c in RESISTANCE + VIRULENCE:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_task_table(path) -> pd.DataFrame:
    df = read_csv_str(path)
    need = ["base_ST", "starting_state", "became_dual_risk"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        die(f"{path}: missing columns {missing}\nActual columns: {list(df.columns)}")
    bad = set(df["starting_state"].unique()) - {"resistance_only", "virulence_only"}
    if bad:
        die(f"{path}: unexpected starting_state values {bad}")
    return df


def load_routes(path, n_features=20) -> np.ndarray:
    p = Path(path)
    if not p.exists():
        die(f"routes file not found: {p.resolve()}")
    routes = np.loadtxt(p, dtype=np.int64, ndmin=2)
    if routes.shape[1] != n_features:
        die(f"{p}: expected {n_features} columns per route, got {routes.shape[1]}")
    if routes.min() == 1 and routes.max() == n_features:
        routes = routes - 1                         # same 1-indexed rescue as eval_score_lineages.R
    if routes.min() != 0 or routes.max() != n_features - 1:
        die(f"{p}: route values span [{routes.min()}, {routes.max()}], expected 0..{n_features-1}")
    bad = int((np.sort(routes, axis=1) != np.arange(n_features)[None, :]).any(axis=1).sum())
    if bad:
        die(f"{p}: {bad} rows are not full permutations of 0..{n_features-1}; "
            "the vectorised scorer assumes full orderings (as in every run so far).")
    return routes


def load_pooled_routes(paths) -> np.ndarray:
    return np.vstack([load_routes(p) for p in paths])


# --------------------------------------------------------------------------
# route scoring (loose Phase 8 method)
# --------------------------------------------------------------------------
class RouteScorer:
    """Vectorised port of score_lineage() (eval_score_lineages.R, loose method):
    for each route, truncate at the LATEST position of any already-present
    character, take the next character, and count a hit if it is in the target
    set. Routes where truncation lands on the last slot are skipped (not valid)."""

    def __init__(self, routes: np.ndarray):
        self.routes = routes
        n, length = routes.shape
        self.n, self.L = n, length
        # pos[r, c] = position of character c in route r
        pos = np.full((n, length), -1, dtype=np.int16)
        pos[np.arange(n)[:, None], routes] = np.arange(length, dtype=np.int16)[None, :]
        self.pos = pos

    def next_features(self, present_idx) -> np.ndarray:
        present_idx = np.asarray(present_idx, dtype=int)
        if present_idx.size == 0:
            return np.array([], dtype=int)
        sub = self.pos[:, present_idx]
        ok = (sub >= 0).all(axis=1)
        trunc = sub.max(axis=1).astype(int)
        ok &= trunc < self.L - 1
        idx = np.nonzero(ok)[0]
        return self.routes[idx, trunc[idx] + 1]

    def score(self, present_idx, target_idx):
        nxt = self.next_features(present_idx)
        if nxt.size == 0:
            return np.nan, 0
        return float(np.isin(nxt, target_idx).mean()), int(nxt.size)


def reference_score_lineage(route_lines, present_idx, target_idx):
    """LITERAL port of the R code (only used by tests / integrity checks).
    route_lines: whitespace-separated strings, exactly as readLines() returns them."""
    hits = 0
    valid = 0
    for line in route_lines:
        steps = [int(x) for x in line.strip().split()]
        try:
            pos = [steps.index(p) + 1 for p in present_idx]      # 1-based, like match()
        except ValueError:
            continue                                            # any(is.na(pos))
        if not pos:
            raise ValueError("empty present set (R would return -Inf here)")
        trunc_at = max(pos)
        if trunc_at >= len(steps):
            continue
        next_feat = steps[trunc_at]                              # steps[trunc_at + 1] in R (1-based)
        valid += 1
        if next_feat in target_idx:
            hits += 1
    return np.nan if valid == 0 else hits / valid


def class_indices(char_order):
    res = [char_order.index(c) for c in RESISTANCE]
    vir = [char_order.index(c) for c in VIRULENCE]
    return res, vir


# --------------------------------------------------------------------------
# task-table reconstruction (mirrors build_task in baseline_lineage_models.py)
# --------------------------------------------------------------------------
def genome_flags(cm, res_cols, vir_cols):
    res = (cm[res_cols] == 1).any(axis=1)
    vir = (cm[vir_cols] == 1).any(axis=1)
    return res, vir


def build_task(cm, train_end, y0, y1, res_cols=RESISTANCE, vir_cols=VIRULENCE) -> pd.DataFrame:
    """Single-risk lineages at TRAIN (<= train_end) that are observed again in [y0, y1].
    label_union   = lineage level: ANY genome in the window has resistance AND ANY genome has
                    virulence (what phase 6 / phase 8 used -> became_dual_risk)
    label_isolate = a SINGLE genome in the window carries both classes."""
    res, vir = genome_flags(cm, res_cols, vir_cols)
    g = pd.DataFrame({"base_ST": cm["base_ST"].values, "res": res.values, "vir": vir.values,
                      "yr": cm["collection_year"].values})
    g["both"] = g["res"] & g["vir"]
    train_agg = g[g["yr"] <= train_end].groupby("base_ST").agg(
        tr_res=("res", "max"), tr_vir=("vir", "max"), n_train=("res", "size"))
    window_agg = g[(g["yr"] >= y0) & (g["yr"] <= y1)].groupby("base_ST").agg(
        w_res=("res", "max"), w_vir=("vir", "max"), w_both=("both", "max"), n_outcome=("res", "size"))
    tab = train_agg.join(window_agg, how="inner")
    tab = tab[tab["tr_res"].astype(bool) ^ tab["tr_vir"].astype(bool)].copy()
    tab["starting_state"] = np.where(tab["tr_res"].astype(bool), "resistance_only", "virulence_only")
    tab["label_union"] = (tab["w_res"].astype(bool) & tab["w_vir"].astype(bool)).astype(int)
    tab["label_isolate"] = tab["w_both"].astype(bool).astype(int)
    return tab[["starting_state", "label_union", "label_isolate", "n_train", "n_outcome"]]


def lineage_states(cm, train_end, char_order) -> pd.DataFrame:
    """0/1 matrix, columns in char_order: 1 if any genome <= train_end in the lineage carries the trait."""
    sub = cm[cm["collection_year"] <= train_end]
    return (sub[char_order] == 1).groupby(sub["base_ST"]).any().astype(int)[char_order]


def lineage_covariates(cm, train_end, char_order) -> pd.DataFrame:
    sub = cm[cm["collection_year"] <= train_end]
    states = lineage_states(cm, train_end, char_order)
    out = pd.DataFrame(index=states.index)
    out["n_present"] = states.sum(axis=1)
    out["n_countries"] = sub.groupby("base_ST")["country"].nunique()
    out["n_train_years"] = sub.groupby("base_ST")["collection_year"].nunique()
    return out


def score_task(task: pd.DataFrame, states: pd.DataFrame, scorer: RouteScorer, char_order,
               targets=None):
    """Score every row of `task` (index = base_ST, column starting_state).
    targets: optional {starting_state: [target idx]} to override the default, which is
    the complementary CLASS (resistance_only -> virulence, virulence_only -> resistance)."""
    res_idx, vir_idx = class_indices(char_order)
    if targets is None:
        targets = {"resistance_only": vir_idx, "virulence_only": res_idx}
    scores, n_valid, no_state = [], [], 0
    for st, start_state in zip(task.index, task["starting_state"]):
        if st not in states.index:
            scores.append(np.nan); n_valid.append(0); no_state += 1
            continue
        present = np.nonzero(states.loc[st].values == 1)[0]
        s, nv = scorer.score(present, targets[start_state])
        scores.append(s); n_valid.append(nv)
    return np.array(scores, float), np.array(n_valid, int), no_state


# --------------------------------------------------------------------------
# statistics: AUROC, DeLong (single + paired), bootstrap, permutation
# --------------------------------------------------------------------------
def _clean(y, s):
    y = np.asarray(y).astype(int)
    s = np.asarray(s, float)
    ok = ~np.isnan(s)
    return y[ok], s[ok]


def auroc(y, s):
    y, s = _clean(y, s)
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    ranks = stats.rankdata(s)
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _placements(y, s):
    pos, neg = s[y == 1], s[y == 0]
    m, n = len(pos), len(neg)
    t_all = stats.rankdata(np.concatenate([pos, neg]))
    t_pos, t_neg = stats.rankdata(pos), stats.rankdata(neg)
    v10 = (t_all[:m] - t_pos) / n
    v01 = 1.0 - (t_all[m:] - t_neg) / m
    return v10, v01


def delong_ci(y, s, alpha=0.05):
    """DeLong CI, truncated to [0, 1] like pROC::ci.auc. Returns (auc, lo, hi)."""
    y, s = _clean(y, s)
    m, n = int((y == 1).sum()), int((y == 0).sum())
    if m < 2 or n < 2:
        return auroc(y, s), np.nan, np.nan
    v10, v01 = _placements(y, s)
    a = float(v10.mean())
    var = v10.var(ddof=1) / m + v01.var(ddof=1) / n
    if var <= 0:
        return a, a, a
    z = stats.norm.ppf(1 - alpha / 2)
    se = np.sqrt(var)
    return a, max(0.0, a - z * se), min(1.0, a + z * se)


def delong_paired(y, s1, s2):
    """Paired DeLong test for two correlated AUCs on the same labels.
    Returns (auc1, auc2, diff, p_two_sided)."""
    y = np.asarray(y).astype(int)
    s1, s2 = np.asarray(s1, float), np.asarray(s2, float)
    ok = ~(np.isnan(s1) | np.isnan(s2))
    y, s1, s2 = y[ok], s1[ok], s2[ok]
    m, n = int((y == 1).sum()), int((y == 0).sum())
    if m < 2 or n < 2:
        return np.nan, np.nan, np.nan, np.nan
    a10, a01 = _placements(y, s1)
    b10, b01 = _placements(y, s2)
    A, B = a10.mean(), b10.mean()
    S10 = np.cov(np.vstack([a10, b10]))
    S01 = np.cov(np.vstack([a01, b01]))
    S = S10 / m + S01 / n
    var = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    if var <= 1e-15:
        return float(A), float(B), float(A - B), np.nan
    z = (A - B) / np.sqrt(var)
    return float(A), float(B), float(A - B), float(2 * (1 - stats.norm.cdf(abs(z))))


def bootstrap_ci(y, s, n_boot=5000, seed=42, alpha=0.05):
    y, s = _clean(y, s)
    rng = np.random.default_rng(seed)
    n = len(y)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        v = auroc(y[idx], s[idx])
        if not np.isnan(v):
            vals.append(v)
    if not vals:
        return np.nan, np.nan
    return tuple(np.quantile(vals, [alpha / 2, 1 - alpha / 2]))


def perm_pvalue(y, s, strata, n_perm=20000, seed=42):
    """One-sided permutation p for AUROC >= observed, with labels permuted WITHIN strata
    (e.g. starting_state), so it asks whether the score discriminates beyond what the
    strata alone explain. Exact in the limit, no normal approximation."""
    y = np.asarray(y).astype(int)
    s = np.asarray(s, float)
    strata = np.asarray(strata)
    ok = ~np.isnan(s)
    y, s, strata = y[ok], s[ok], strata[ok]
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return np.nan, np.nan
    ranks = stats.rankdata(s)
    obs = (ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    rng = np.random.default_rng(seed)
    total = np.zeros(n_perm)
    for lab in np.unique(strata):
        idx = np.nonzero(strata == lab)[0]
        k = int(y[idx].sum())
        if k == 0:
            continue
        if k == len(idx):
            total += ranks[idx].sum()
            continue
        rand = rng.random((n_perm, len(idx)))
        pick = np.argpartition(rand, k - 1, axis=1)[:, :k]
        total += ranks[idx][pick].sum(axis=1)
    null = (total - n1 * (n1 + 1) / 2) / (n1 * n0)
    p = (1 + int((null >= obs - 1e-12).sum())) / (1 + n_perm)
    return float(obs), float(p)


def avg_precision(y, s):
    from sklearn.metrics import average_precision_score
    y, s = _clean(y, s)
    if len(np.unique(y)) < 2:
        return np.nan
    return float(average_precision_score(y, s))


def fmt_auc(y, s):
    """'0.664 [0.512, 0.816] (n=143, pos=24)'"""
    y, s = _clean(y, s)
    a, lo, hi = delong_ci(y, s)
    return f"{a:.3f} [{lo:.3f}, {hi:.3f}] (n={len(y)}, pos={int(y.sum())})"


# --------------------------------------------------------------------------
# misc
# --------------------------------------------------------------------------
ACC_RE = re.compile(r"GC[AF]_\d{9}\.\d+")


def accession_from_path(p: str):
    m = ACC_RE.search(p)
    return m.group(0) if m else None


def ensure_dir(d) -> Path:
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    return d