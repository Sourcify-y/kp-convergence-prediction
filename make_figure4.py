import argparse, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.linear_model import LogisticRegression

BLUE, ORANGE, DARK = "#0072B2", "#D55E00", "#222222"
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 9, "legend.fontsize": 7.5,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": DARK, "axes.labelcolor": DARK, "text.color": DARK,
    "xtick.color": DARK, "ytick.color": DARK,
    "pdf.fonttype": 42, "savefig.facecolor": "white", "figure.facecolor": "white"})

parser = argparse.ArgumentParser()
parser.add_argument("--val", required=True)
parser.add_argument("--test", required=True)
parser.add_argument("--score", default="predicted_risk")
parser.add_argument("--label", default="became_dual_risk")
parser.add_argument("--out", default="Figure4_ROC_calibration_pub")
args = parser.parse_args()

panels = [("Validation (2018\u20132020)", args.val, BLUE), ("Test (2021\u20132022)", args.test, ORANGE)]
fig, axes = plt.subplots(2, 2, figsize=(6.6, 6.6))

for col, (title, path, color) in enumerate(panels):
    df = pd.read_csv(path)[[args.score, args.label]].dropna()
    y = (df[args.label].astype(float) > 0).astype(int).values
    s = df[args.score].astype(float).values
    auc = roc_auc_score(y, s)
    print("%s: n=%d, positives=%d, AUROC=%.3f" % (title, len(y), y.sum(), auc))

    # top row: ROC
    ax = axes[0, col]
    fpr, tpr, _ = roc_curve(y, s)
    ax.plot([0, 1], [0, 1], color="#9A9A9A", lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.plot(fpr, tpr, color=color, lw=1.6, zorder=3,
            label="HyperTraPS\nAUROC = %.3f\n(n = %d, pos = %d)" % (auc, len(y), y.sum()))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xticks([0, .2, .4, .6, .8, 1]); ax.set_yticks([0, .2, .4, .6, .8, 1])
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
    ax.set_title(title)
    ax.legend(loc="lower right", frameon=False, handlelength=1.6)

    # bottom row: calibration after in-sample Platt scaling, 4 quantile bins
    ax = axes[1, col]
    platt = LogisticRegression(C=1e6, max_iter=1000).fit(s.reshape(-1, 1), y)
    prob = platt.predict_proba(s.reshape(-1, 1))[:, 1]
    bins = pd.DataFrame({"p": prob, "y": y, "b": pd.qcut(prob, q=4, duplicates="drop")}).groupby("b", observed=True)
    mean_pred, obs_frac, bin_n = bins["p"].mean().values, bins["y"].mean().values, bins["y"].size().values
    lim = 0.6
    ax.plot([0, lim], [0, lim], color="#9A9A9A", lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.plot(mean_pred, obs_frac, "-", color=color, lw=1.4, zorder=2)
    ax.scatter(mean_pred, obs_frac, s=34, color=color, edgecolor="white", linewidth=0.6, zorder=3)
    for xi, yi, ni in zip(mean_pred, obs_frac, bin_n):
        ax.annotate("n=%d" % ni, (xi, yi), textcoords="offset points", xytext=(0, 7),
                    ha="center", fontsize=6.5, color="#555555")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
    ax.set_xticks(np.arange(0, lim + 0.01, 0.1)); ax.set_yticks(np.arange(0, lim + 0.01, 0.1))
    ax.set_xlabel("Predicted probability\n(in-sample Platt scaling)")
    ax.set_ylabel("Observed fraction becoming dual-risk")

for ax, tag in zip(axes.ravel(), "ABCD"):
    ax.text(-0.22, 1.05, tag, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")
fig.tight_layout(h_pad=1.6, w_pad=1.6)
fig.savefig(args.out + ".png", dpi=600, bbox_inches="tight")
fig.savefig(args.out + ".pdf", bbox_inches="tight")
print("saved %s.png and %s.pdf" % (args.out, args.out))