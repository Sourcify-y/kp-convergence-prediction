#!/usr/bin/env python3
"""
model_synthetic_recovery_test.py (18-feature version)

Phase 1, Step 9: write two synthetic observation files for hypertraps.ce to
check that it (a) recovers a known acquisition order when there is one and
(b) doesn't invent an order when there isn't (null control).

Positive set: N independent "walks" through an 18-feature binary state space,
features acquired roughly in order 0..17 with noise. Features 16 and 17 are
coupled (17 almost always co-occurs with 16) to mimic a near-redundant locus
pair, so you can see whether coupled late features show extra posterior spread.

Null set: every feature present/absent independently with its own fixed
marginal, no ordering. The same 16/17 coupling is injected, so correlation-driven
spread can be told apart from ordering-driven spread.

Usage:
    python3 model_synthetic_recovery_test.py --outdir ./spike_data --n 500
"""

import argparse
import random
from pathlib import Path

N_FEATURES = 18

# correlated pair, stand-in for the RmpADC/rmpA pair from the earlier 18-character set
CORR_PARENT = 16
CORR_CHILD = 17
CORR_CHILD_GIVEN_PARENT = 0.9        # P(child=1 | parent=1), near-total redundancy
CORR_CHILD_GIVEN_NOT_PARENT = 0.05   # P(child=1 | parent=0), rare independent occurrence


def simulate_positive(n: int, seed: int = 42):
    # ordered accumulation: features come in order 0..17, each acquired with
    # some probability when its turn comes, plus noise so it isn't a rigid
    # staircase. The child feature depends on the parent instead of on order.
    rng = random.Random(seed)
    obs = []
    for _ in range(n):
        state = [0] * N_FEATURES
        for i in range(N_FEATURES):
            if i == CORR_CHILD:
                p = CORR_CHILD_GIVEN_PARENT if state[CORR_PARENT] == 1 else CORR_CHILD_GIVEN_NOT_PARENT
                if rng.random() < p:
                    state[i] = 1
                continue
            if rng.random() < 0.75:  # 75% chance of acquiring this feature on its turn
                state[i] = 1
            else:
                if rng.random() < 0.3:
                    continue  # skip it but keep walking
                else:
                    break  # walk stops here
            if rng.random() < 0.15:  # random early stopping
                break
        obs.append(state)
    return obs


def simulate_null(n: int, seed: int = 99):
    # independent features with fixed marginals, no ordering; only the 16/17
    # pair is coupled, same strength as the positive set
    rng = random.Random(seed)
    marginals = [rng.uniform(0.2, 0.8) for _ in range(N_FEATURES)]
    obs = []
    for _ in range(n):
        state = [1 if rng.random() < marginals[i] else 0 for i in range(N_FEATURES)]
        p = CORR_CHILD_GIVEN_PARENT if state[CORR_PARENT] == 1 else CORR_CHILD_GIVEN_NOT_PARENT
        state[CORR_CHILD] = 1 if rng.random() < p else 0
        obs.append(state)
    return obs


def write_obs(obs, outpath: Path):
    with open(outpath, "w") as f:
        for row in obs:
            f.write(" ".join(str(x) for x in row) + "\n")
    print(f"[write] {outpath} -- {len(obs)} observations x {N_FEATURES} features")


def summarize(obs, label):
    # mean presence per feature: positive data should show a gradient, null
    # shouldn't, and 16/17 should track each other in both
    n = len(obs)
    means = [sum(row[i] for row in obs) / n for i in range(N_FEATURES)]
    print(f"[{label}] per-feature presence rate (feature 0..{N_FEATURES-1}): "
          + ", ".join(f"{m:.2f}" for m in means))
    print(f"[{label}] correlated pair check -- feature {CORR_PARENT}: {means[CORR_PARENT]:.2f}, "
          f"feature {CORR_CHILD}: {means[CORR_CHILD]:.2f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--n", type=int, default=500,
                        help="Number of synthetic observations per dataset")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("Generating POSITIVE (known-order) synthetic dataset...")
    pos = simulate_positive(args.n)
    summarize(pos, "positive")
    write_obs(pos, args.outdir / "synthetic_positive_18.txt")

    print("\nGenerating NULL (no-order) synthetic dataset...")
    null = simulate_null(args.n)
    summarize(null, "null")
    write_obs(null, args.outdir / "synthetic_null_18.txt")

    print(f"\n[done] Both files written to {args.outdir}")
    print("\nNEXT STEPS:")
    print(f"  hypertraps.ce --obs {args.outdir}/synthetic_positive_18.txt --label positive_test_18")
    print(f"  hypertraps.ce --obs {args.outdir}/synthetic_null_18.txt --label null_test_18")
    print("\nINTERPRETATION GUIDE:")
    print("  - POSITIVE run: expect a roughly monotonic mean-presence gradient")
    print("    from feature 0 (highest) to feature 17 (lowest), and expect the")
    print("    LATE features (especially 15-17) to show much larger posterior")
    print("    spread than early features -- mirroring what was seen in the earlier real")
    print("    n=500, 18-character run (near-zero spread early, ~900-1000 late).")
    print("  - Feature 17 (the correlated child) may show even MORE spread than")
    print("    its neighbors, since it's both late AND dependent on feature 16 --")
    print("    if so, that isolates correlation as an amplifier of the")
    print("    late-acquisition sparsity effect, not a separate problem.")
    print("  - NULL run: expect flat/inconsistent spread across most features,")
    print("    EXCEPT possibly features 16/17, which are correlated by design")
    print("    even though there's no ordering. If 16/17 show elevated spread")
    print("    here too, that confirms it's the correlation itself (not order)")
    print("    driving some of the effect -- if they look normal here but huge")
    print("    in the positive run, that confirms it's specifically the")
    print("    late-order + sparse-data combination.")


if __name__ == "__main__":
    main()