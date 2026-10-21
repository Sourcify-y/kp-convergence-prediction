#!/usr/bin/env python3
"""
Compare the per-Time-slot acquisition-order probabilities across independent
HyperTraPS-CT seed runs, using their -bubbles.csv files.

Usage:
    python3 model_compare_seeds.py s0-bubbles.csv s1-bubbles.csv s2-bubbles.csv

For every pair of seeds and every Time slot (0-19) it takes the probability
vector over OriginalIndex (0-19), computes the Spearman correlation between the
two seeds, and flags where the top-1 (argmax) feature disagrees. It also prints
mean/min correlation, the number of top-1 disagreements, and each character's
peak Time position per seed.
"""
import sys
import csv
from collections import defaultdict
from itertools import combinations

try:
    from scipy.stats import spearmanr
except ImportError:
    print("scipy not found — install with: pip install scipy --break-system-packages")
    sys.exit(1)


def load_bubbles(path):
    # Time -> {OriginalIndex: Probability}
    data = defaultdict(dict)
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = int(row["Time"])
            idx = int(row["OriginalIndex"])
            prob = float(row["Probability"])
            data[t][idx] = prob
    return data


def compare_pair(name_a, data_a, name_b, data_b, n_chars=20):
    print(f"\n=== {name_a} vs {name_b} ===")
    corrs = []
    disagreements = 0
    for t in sorted(set(data_a) & set(data_b)):
        vec_a = [data_a[t].get(i, 0.0) for i in range(n_chars)]
        vec_b = [data_b[t].get(i, 0.0) for i in range(n_chars)]
        rho, _ = spearmanr(vec_a, vec_b)
        corrs.append(rho)
        top_a = max(range(n_chars), key=lambda i: vec_a[i])
        top_b = max(range(n_chars), key=lambda i: vec_b[i])
        flag = ""
        if top_a != top_b:
            disagreements += 1
            flag = f"  <-- TOP-1 DISAGREEMENT ({top_a} vs {top_b})"
        print(f"  Time={t:2d}  spearman_rho={rho:.3f}  top1_a={top_a:2d} ({vec_a[top_a]:.3f})  "
              f"top1_b={top_b:2d} ({vec_b[top_b]:.3f}){flag}")
    print(f"  --- mean rho = {sum(corrs)/len(corrs):.3f}, min rho = {min(corrs):.3f}, "
          f"top-1 disagreements = {disagreements}/{len(corrs)} Time slots ---")


def peak_times(data, n_chars=20):
    # for each feature index, the Time where its probability is highest
    peaks = {}
    for idx in range(n_chars):
        best_t, best_p = None, -1.0
        for t, char_probs in data.items():
            p = char_probs.get(idx, 0.0)
            if p > best_p:
                best_p, best_t = p, t
        peaks[idx] = (best_t, best_p)
    return peaks


def compare_peaks(names_and_data, char_names=None, n_chars=20):
    print("\n=== Per-character peak Time position (mode of probability across Time) ===")
    all_peaks = [(name, peak_times(data, n_chars)) for name, data in names_and_data]
    header = "Idx  " + ("Name".ljust(20) if char_names else "") + "  " + "  ".join(
        f"{name.split('/')[-1].split('-bubbles')[0]:>20}" for name, _ in all_peaks)
    print(header)
    for idx in range(n_chars):
        label = char_names[idx] if char_names else ""
        row = f"{idx:3d}  " + (label.ljust(20) if char_names else "") + "  "
        row += "  ".join(f"T={peaks[idx][0]:2d} (p={peaks[idx][1]:.3f})".rjust(20) for _, peaks in all_peaks)
        print(row)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    paths = sys.argv[1:]
    datasets = [(p, load_bubbles(p)) for p in paths]

    # index -> name, same order as CHAR_ORDER in the HyperTraPS input
    char_names = [
        "Bla_ESBL_acquired", "Bla_Carb_acquired", "Flq_acquired", "Sul_acquired",
        "AGly_acquired", "Col_acquired", "Tet_acquired", "MLS_acquired",
        "Rif_acquired", "Tmt_acquired", "Fcyn_acquired", "Phe_acquired",
        "Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin",
        "rmpA", "rmpA2", "Bla_acquired", "Tgc_acquired",
    ]

    for (name_a, data_a), (name_b, data_b) in combinations(datasets, 2):
        compare_pair(name_a, data_a, name_b, data_b)

    compare_peaks(datasets, char_names=char_names)


if __name__ == "__main__":
    main()