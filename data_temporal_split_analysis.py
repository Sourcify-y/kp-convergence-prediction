#!/usr/bin/env python3
"""
data_temporal_split_analysis.py

Phase 1, Step 4: work out the TRAIN / VALIDATE / TEST cutoffs (T1, T2) from
the actual per-year genome and dual-risk counts, instead of picking years
arbitrarily.

"Dual-risk" = at least one resistance class AND at least one virulence locus.
For each year and each candidate (T1, T2) split it prints genome and
dual-risk counts per window, also broken out for ST258 and ST11. It doesn't
pick T1/T2 for you, it just gives you the numbers to decide with.

Usage:
  python3 data_temporal_split_analysis.py --csv kp_pull_20k/kp_metadata_FULL_dedup.csv
  python3 data_temporal_split_analysis.py --csv ... --t1-candidates 2012 2013 2014 --t2-candidates 2015 2016 2017
"""

import argparse
import csv
import re


RES_COLS = [
    "Bla_ESBL_acquired", "Bla_Carb_acquired", "Flq_acquired", "Sul_acquired",
    "AGly_acquired", "Col_acquired", "Tet_acquired", "MLS_acquired",
    "Rif_acquired", "Tmt_acquired", "Fcyn_acquired", "Phe_acquired",
]
VIR_COLS = [
    "Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin",
    "RmpADC", "rmpA2",
]


def is_present(value):
    if value is None:
        return False
    v = value.strip()
    if v == "" or v == "-" or v.lower() in ("na", "n/a", "unknown"):
        return False
    return True


def has_any(row, cols):
    return any(is_present(row.get(c)) for c in cols)


def base_st(st):
    # collapse LV tags ("ST11-1LV" -> "ST11") so variants count toward their
    # real lineage and don't fall into OTHER
    if st is None:
        return st
    return re.sub(r"-\d+LV$", "", st.strip())


def load_usable(rows):
    # keep rows with a numeric collection_year, and tag each with year,
    # base ST and AMR / virulence / dual-risk flags
    usable = []
    for r in rows:
        yr = r.get("collection_year", "")
        st = r.get("ST", "")
        if not yr or not yr.strip():
            continue
        try:
            year_float = float(yr)
        except ValueError:
            continue
        r["_year"] = int(year_float)
        r["_base_st"] = base_st(st)
        r["_has_amr"] = has_any(r, RES_COLS)
        r["_has_vir"] = has_any(r, VIR_COLS)
        r["_dual_risk"] = r["_has_amr"] and r["_has_vir"]
        usable.append(r)
    return usable


def print_year_by_lineage_table(usable):
    print("=" * 100)
    print("Genomes and dual-risk genomes per year, by lineage")
    print("-" * 100)
    years = sorted(set(r["_year"] for r in usable))
    lineages = ["ST258", "ST11", "OTHER"]

    header = f"  {'year':<6}"
    for lin in lineages:
        header += f"{lin+'_n':<10}{lin+'_dual':<11}"
    header += f"{'TOTAL_n':<10}{'TOTAL_dual':<11}"
    print(header)

    for yr in years:
        yr_rows = [r for r in usable if r["_year"] == yr]
        row_str = f"  {yr:<6}"
        for lin in lineages:
            if lin == "OTHER":
                sub = [r for r in yr_rows if r["_base_st"] not in ("ST258", "ST11")]
            else:
                sub = [r for r in yr_rows if r["_base_st"] == lin]
            n = len(sub)
            dual = sum(1 for r in sub if r["_dual_risk"])
            row_str += f"{n:<10}{dual:<11}"
        total_n = len(yr_rows)
        total_dual = sum(1 for r in yr_rows if r["_dual_risk"])
        row_str += f"{total_n:<10}{total_dual:<11}"
        print(row_str)
    print()


def print_split_summary(usable, t1, t2):
    # TRAIN: year <= t1, VALIDATE: t1 < year <= t2, TEST: year > t2
    train = [r for r in usable if r["_year"] <= t1]
    validate = [r for r in usable if t1 < r["_year"] <= t2]
    test = [r for r in usable if r["_year"] > t2]

    def stats(subset, label):
        n = len(subset)
        dual = sum(1 for r in subset if r["_dual_risk"])
        st258 = sum(1 for r in subset if r["_base_st"] == "ST258")
        st258_dual = sum(1 for r in subset if r["_base_st"] == "ST258" and r["_dual_risk"])
        st11 = sum(1 for r in subset if r["_base_st"] == "ST11")
        st11_dual = sum(1 for r in subset if r["_base_st"] == "ST11" and r["_dual_risk"])
        print(f"  {label:<10} n={n:<6} dual_risk={dual:<6} "
              f"ST258: n={st258} dual={st258_dual}   ST11: n={st11} dual={st11_dual}")

    print(f"--- T1={t1}, T2={t2}  (TRAIN <= {t1} | VALIDATE {t1+1}-{t2} | TEST > {t2}) ---")
    stats(train, "TRAIN")
    stats(validate, "VALIDATE")
    stats(test, "TEST")

    warnings = []
    if len(validate) < 30:
        warnings.append(f"VALIDATE window has only {len(validate)} genomes — likely too thin")
    if len(test) < 30:
        warnings.append(f"TEST window has only {len(test)} genomes — likely too thin")
    test_dual = sum(1 for r in test if r["_dual_risk"])
    if test_dual < 10:
        warnings.append(f"TEST window has only {test_dual} dual-risk genomes — "
                        "may be too few observable transitions to validate against")
    if warnings:
        print("  [!] " + " | ".join(warnings))
    print()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--t1-candidates", nargs="+", type=int,
                        default=[2012, 2013, 2014, 2015])
    parser.add_argument("--t2-candidates", nargs="+", type=int,
                        default=[2015, 2016, 2017])
    args = parser.parse_args()

    with open(args.csv) as f:
        rows = list(csv.DictReader(f))
    usable = load_usable(rows)
    print(f"Loaded {len(rows)} genomes, {len(usable)} usable (have a numeric collection_year)\n")

    print_year_by_lineage_table(usable)

    print("=" * 100)
    print("CANDIDATE TRAIN / VALIDATE / TEST SPLITS")
    print("(pick T1 < T2 pairs from --t1-candidates x --t2-candidates)")
    print("=" * 100)
    for t1 in sorted(args.t1_candidates):
        for t2 in sorted(args.t2_candidates):
            if t2 <= t1:
                continue
            print_split_summary(usable, t1, t2)


if __name__ == "__main__":
    main()