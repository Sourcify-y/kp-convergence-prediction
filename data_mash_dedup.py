#!/usr/bin/env python3
"""
data_mash_dedup.py

Reads a mash triangle output file (lower-triangular distance matrix), clusters
genomes below a distance threshold with single linkage (union-find), and
keeps one representative genome per cluster.

Usage:
  python3 data_mash_dedup.py --triangle mash_test/distances_500.tsv --threshold 0.001
"""

import argparse


def parse_triangle(path):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].strip())
    names = []
    dist_rows = []
    for line in lines[1:1 + n]:
        parts = line.split("\t")
        names.append(parts[0])
        dists = [float(x) for x in parts[1:]]
        dist_rows.append(dists)
    return names, dist_rows


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--triangle", required=True)
    parser.add_argument("--threshold", type=float, required=True,
                        help="Mash distance below which two genomes are considered near-identical")
    parser.add_argument("--out", default="cluster_representatives.txt")
    args = parser.parse_args()

    names, dist_rows = parse_triangle(args.triangle)
    n = len(names)
    uf = UnionFind(n)

    # dist_rows[i] = distances from genome i to genomes 0..i-1
    for i, row in enumerate(dist_rows):
        for j, d in enumerate(row):
            if d < args.threshold:
                uf.union(i, j)

    clusters = {}
    for i in range(n):
        root = uf.find(i)
        clusters.setdefault(root, []).append(i)

    n_clusters = len(clusters)
    n_singletons = sum(1 for members in clusters.values() if len(members) == 1)
    largest = max(len(members) for members in clusters.values())

    print(f"Genomes: {n}")
    print(f"Threshold: {args.threshold}")
    print(f"Clusters: {n_clusters}")
    print(f"Singleton clusters: {n_singletons}")
    print(f"Largest cluster size: {largest}")

    with open(args.out, "w") as f:
        for root, members in clusters.items():
            rep = names[members[0]]  # just the first member, no quality criterion
            f.write(f"{rep}\t{len(members)}\n")

    print(f"[wrote] {args.out} ({n_clusters} representative genomes)")


if __name__ == "__main__":
    main()