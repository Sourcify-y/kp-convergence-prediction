#!/usr/bin/env python3
"""
Recompute the HyperTraPS-CT -bubbles.csv probabilities from -routes.txt after
dropping an early fraction of the posterior samples as burn-in.

-routes.txt has samples_per_row (default 10) simulated routes for each row of
-posterior.txt, in MCMC order, and the tool doesn't remove burn-in itself.
This drops the first --burnin fraction of posterior-sample BLOCKS (not raw
rows), then for each Time position and feature index computes the fraction of
the remaining routes with that feature at that position.

Usage:
    python3 model_burnin_trimmed_bubbles.py \\
        --routes kp_full_run_len4_s0-routes.txt \\
        --n-posterior 79 \\
        --samples-per-row 10 \\
        --burnin 0.2 \\
        --out kp_full_run_len4_s0-bubbles-trimmed.csv

--n-posterior is the row count of the matching -posterior.txt (check with wc -l).
"""
import argparse
import csv
from collections import defaultdict


def load_routes(path):
    # one list of feature indices (one per Time position) per route
    routes = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            routes.append([int(x) for x in line.split()])
    return routes


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--routes", required=True, help="-routes.txt file")
    parser.add_argument("--n-posterior", type=int, required=True,
                        help="row count of the matching -posterior.txt file")
    parser.add_argument("--samples-per-row", type=int, default=10,
                        help="routes simulated per posterior sample (tool default: 10)")
    parser.add_argument("--burnin", type=float, default=0.2,
                        help="fraction of posterior-sample BLOCKS to discard as burn-in (default 0.2)")
    parser.add_argument("--out", required=True, help="output CSV path")
    args = parser.parse_args()

    routes = load_routes(args.routes)
    expected = args.n_posterior * args.samples_per_row
    if len(routes) != expected:
        print(f"WARNING: routes.txt has {len(routes)} rows, expected "
              f"{args.n_posterior} x {args.samples_per_row} = {expected}. "
              f"Check --n-posterior/--samples-per-row before trusting the output.")

    n_chars = max(max(r) for r in routes) + 1
    n_time = len(routes[0])

    burnin_blocks = int(round(args.n_posterior * args.burnin))
    burnin_rows = burnin_blocks * args.samples_per_row
    kept = routes[burnin_rows:]

    print(f"Total posterior samples: {args.n_posterior}")
    print(f"Discarding first {burnin_blocks} blocks ({burnin_rows} of {len(routes)} routes) as burn-in")
    print(f"Retained {len(kept)} routes for recomputed bubbles")

    # counts[time][feature]
    counts = defaultdict(lambda: defaultdict(int))
    for route in kept:
        for t, feat in enumerate(route):
            counts[t][feat] += 1

    n_kept = len(kept)
    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Time", "OriginalIndex", "Probability"])
        for t in range(n_time):
            for idx in range(n_chars):
                prob = counts[t].get(idx, 0) / n_kept
                writer.writerow([t, idx, f"{prob:.15f}"])

    print(f"Wrote trimmed bubble probabilities to {args.out}")


if __name__ == "__main__":
    main()