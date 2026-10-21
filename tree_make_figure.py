#!/usr/bin/env python3
"""Draws the Phase 4 NJ tree with acquisition events.

Input:  data/phase4/kp_dedup_nj_tree.nwk
        data/phase4/candidate_acquisition_events_flagged_v2.csv
Output: figures/Phase4_NJ_tree_acquisition_events.png

Markers sit at the descendant node of each acquisition event. Negative branch
lengths are drawn as zero, and branches beyond the 99.5th percentile of tip
depth are clipped in the drawing.
"""
import argparse, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

VIR = ["Yersiniabactin", "Aerobactin", "Salmochelin", "rmpA", "rmpA2", "Colibactin"]


def parse_newick(path):
    text = open(path).read().strip().rstrip(";")
    tokens = re.findall(r"\(|\)|,|:[^(),;]+|[^(),:;]+", text)
    parent, length, label, children = [], [], [], []

    def new_node(p):
        parent.append(p); length.append(0.0); label.append(None); children.append([])
        idx = len(parent) - 1
        if p >= 0: children[p].append(idx)
        return idx

    stack, last = [], None
    for tok in tokens:                  # iterative: the tree is very deep
        if tok == "(":
            stack.append(new_node(stack[-1] if stack else -1))
        elif tok == ")":
            last = stack.pop()
        elif tok == ",":
            pass
        elif tok[0] == ":":
            length[last] = float(tok[1:])
        else:
            last = new_node(stack[-1]); label[last] = tok
    # preorder traversal
    order, todo = [], [0]
    while todo:
        u = todo.pop(); order.append(u); todo.extend(reversed(children[u]))
    # node numbers: tips first (in preorder), then internal nodes
    num, n_tips = {}, 0
    for u in order:
        if label[u] is not None: n_tips += 1; num[u] = n_tips
    n_numbered = n_tips
    for u in order:
        if label[u] is None: n_numbered += 1; num[u] = n_numbered
    return parent, length, label, children, order, num, n_tips


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", default="data/phase4/kp_dedup_nj_tree.nwk")
    parser.add_argument("--events", default="data/phase4/candidate_acquisition_events_flagged_v2.csv")
    parser.add_argument("--out", default="figures/Phase4_NJ_tree_acquisition_events.png")
    args = parser.parse_args()

    parent, length, label, children, order, num, n_tips = parse_newick(args.tree)
    n_nodes = len(parent)
    node_of_num = {v: k for k, v in num.items()}
    ev = pd.read_csv(args.events)
    # integrity checks against the Phase 4 event table
    assert n_tips == 8318 and n_nodes == 16635, (n_tips, n_nodes)
    assert all(abs(length[node_of_num[d]] - e) < 1e-9 for d, e in zip(ev.descendant_node, ev.edge_length))
    assert all(num[parent[node_of_num[d]]] == p for d, p in zip(ev.descendant_node, ev.ancestral_node))
    tip_events = ev[ev.descendant_is_tip]
    assert all(label[node_of_num[d]] == g for d, g in zip(tip_events.descendant_node, tip_events.descendant_genome_id))

    # branch lengths (negatives drawn as zero) -> depth from the root
    blen = [max(0.0, x) for x in length]; blen[0] = 0.0
    depth = np.zeros(n_nodes)
    for u in order:
        if parent[u] >= 0: depth[u] = depth[parent[u]] + blen[u]
    # angles: tips evenly spaced around the circle, internal nodes midway between their children
    tips = [u for u in order if label[u] is not None]
    angle = np.zeros(n_nodes)
    for i, u in enumerate(tips): angle[u] = 2 * np.pi * i / len(tips)
    for u in reversed(order):
        if children[u]:
            kid_angles = [angle[c] for c in children[u]]
            angle[u] = (min(kid_angles) + max(kid_angles)) / 2
    clip_depth = np.quantile(depth[tips], 0.995)
    squash = lambda d: np.minimum(d / clip_depth, 1.04) * 0.96
    radius = squash(depth)
    segs, arcs = [], []
    for u in order:
        p = parent[u]
        if p >= 0:
            segs.append([(radius[p]*np.cos(angle[u]), radius[p]*np.sin(angle[u])),
                         (radius[u]*np.cos(angle[u]), radius[u]*np.sin(angle[u]))])
        if children[u]:
            kid_angles = [angle[c] for c in children[u]]
            theta = np.linspace(min(kid_angles), max(kid_angles), max(2, int((max(kid_angles)-min(kid_angles))*80)+2))
            arcs.append(np.c_[radius[u]*np.cos(theta), radius[u]*np.sin(theta)])
    kind = np.where(ev.character == "Bla_Carb_acquired", "carb", np.where(ev.character.isin(VIR), "vir", "res"))
    ev["kind"] = kind
    ev["node_idx"] = ev.descendant_node.map(node_of_num)
    # (colour, marker size, zorder, legend label)
    style = {"res": ("#9aa5b1", 3, 2, "Other resistance class events (n={})"),
             "carb": ("#1b6ca8", 9, 3, "Carbapenem resistance (Bla_Carb) events (n={})"),
             "vir": ("#d9480f", 9, 4, "Virulence-associated locus events (n={})")}
    fig, ax = plt.subplots(figsize=(11, 11), dpi=250)
    ax.add_collection(LineCollection(segs, colors="#4a4a4a", linewidths=0.2, alpha=.95))
    ax.add_collection(LineCollection(arcs, colors="#4a4a4a", linewidths=0.2, alpha=.95))
    for k, (color, size, z, text) in style.items():
        sub = ev[ev["kind"] == k]
        ax.scatter(radius[sub["node_idx"]]*np.cos(angle[sub["node_idx"]]),
                   radius[sub["node_idx"]]*np.sin(angle[sub["node_idx"]]),
                   s=size, c=color, alpha=.55 if k == "res" else .85, linewidths=0,
                   label=text.format(len(sub)), zorder=z)
    # dotted reference circles at fixed distances
    for dist in (0.005, 0.01, 0.015):
        theta = np.linspace(0, 2*np.pi, 400)
        ax.plot(squash(dist)*np.cos(theta), squash(dist)*np.sin(theta), color="#bbb", lw=.5, ls=":", zorder=0)
        ax.text(0, squash(dist), str(dist), fontsize=6, color="#888", ha="center", va="bottom")
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05)
    ax.legend(loc="upper left", frameon=False, fontsize=9, markerscale=1.8, bbox_to_anchor=(0, 1))
    plt.savefig(args.out, bbox_inches="tight", facecolor="white")
    print("wrote", args.out, ev["kind"].value_counts().to_dict())


if __name__ == "__main__":
    main()