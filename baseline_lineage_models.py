# Phase 6 - simple baselines
# Among lineages (base_ST) that are single-risk (resistance-only or virulence-only)
# at TRAIN (<=2017), which ones become dual-risk by the end of VALIDATE (2018-2020)?
#
# run from the project root:  python3 baseline_lineage_models.py
# needs pandas + scikit-learn

import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import cross_val_predict, StratifiedKFold

MATRIX_PATH = "phase3_out/character_matrix.csv"

RES_COLS = [
    "Bla_acquired", "Bla_ESBL_acquired", "Bla_Carb_acquired",
    "AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired",
    "MLS_acquired", "Phe_acquired", "Rif_acquired", "Sul_acquired",
    "Tet_acquired", "Tgc_acquired", "Tmt_acquired",
]
VIR_COLS = [
    "Yersiniabactin", "Aerobactin", "Salmochelin", "Colibactin",
    "rmpA", "rmpA2",
]

TRAIN_CUTOFF = 2017
VALIDATE_END = 2020


def load_matrix(path):
    df = pd.read_csv(path)
    required = RES_COLS + VIR_COLS + ["base_ST", "collection_year", "country"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing expected columns: {missing}\n"
            f"Actual columns in file: {list(df.columns)}\n"
            "-> Check the real column names for genome_id/base_ST/collection_year/country "
            "in your character_matrix.csv and adjust this script's column list accordingly "
            "before trusting any output below."
        )
    return df


def risk_state(sub_df):
    # (has_resistance, has_virulence) for one lineage in one window:
    # True if any genome has any of the resistance / virulence characters
    has_res = (sub_df[RES_COLS].sum(axis=1) > 0).any()
    has_vir = (sub_df[VIR_COLS].sum(axis=1) > 0).any()
    return has_res, has_vir


def build_task(df):
    train = df[df["collection_year"] <= TRAIN_CUTOFF]
    val = df[(df["collection_year"] > TRAIN_CUTOFF) & (df["collection_year"] <= VALIDATE_END)]

    rows = []
    for st, st_train in train.groupby("base_ST"):
        has_res, has_vir = risk_state(st_train)
        if has_res and has_vir:
            continue  # already dual-risk at TRAIN
        if not has_res and not has_vir:
            continue  # neither trait, so not a candidate

        st_val = val[val["base_ST"] == st]
        if len(st_val) == 0:
            continue  # never seen again in VALIDATE, can't score it

        val_res, val_vir = risk_state(st_val)
        became_dual = int(val_res and val_vir)

        # most common country in TRAIN as a basic geography feature
        top_country = st_train["country"].mode().iloc[0] if st_train["country"].notna().any() else "Unknown"

        rows.append({
            "base_ST": st,
            "start_state": "resistance_only" if has_res else "virulence_only",
            "top_country": top_country,
            "n_train_genomes": len(st_train),
            "label": became_dual,
        })

    task_df = pd.DataFrame(rows)
    return task_df


def run_baselines(task_df):
    print(f"\nTask built: {len(task_df)} single-risk lineages scored")
    print(f"Positive rate (became dual-risk): {task_df['label'].mean():.3f}")
    print(task_df["start_state"].value_counts())
    print()

    if task_df["label"].nunique() < 2:
        print("WARNING: only one class present in labels -- AUROC/AUPRC undefined. "
              "Check TRAIN/VALIDATE window definitions before trusting anything downstream.")
        return

    y = task_df["label"].values
    cv = StratifiedKFold(n_splits=min(5, task_df["label"].sum()), shuffle=True, random_state=0)

    feature_sets = {}

    # 1. ST alone (one-hot per lineage)
    X1 = pd.get_dummies(task_df["base_ST"])
    feature_sets["ST alone"] = X1

    # 2. ST + geography
    X2 = pd.concat([X1, pd.get_dummies(task_df["top_country"], prefix="country")], axis=1)
    feature_sets["ST + geography"] = X2

    # 3. current state only, no ST identity
    X3 = pd.get_dummies(task_df["start_state"])
    feature_sets["Current state only"] = X3

    # 4. state + geography
    X4 = pd.concat([X3, pd.get_dummies(task_df["top_country"], prefix="country")], axis=1)
    feature_sets["State + geography"] = X4

    print(f"{'Model':<25} {'AUROC':>8} {'AUPRC':>8}")
    print("-" * 43)
    for name, X in feature_sets.items():
        X = X.astype(float)
        try:
            preds = cross_val_predict(
                LogisticRegression(max_iter=1000),
                X, y, cv=cv, method="predict_proba"
            )[:, 1]
            auroc = roc_auc_score(y, preds)
            auprc = average_precision_score(y, preds)
            print(f"{name:<25} {auroc:>8.3f} {auprc:>8.3f}")
        except Exception as e:
            print(f"{name:<25}  FAILED: {e}")

    print("\nCompare these numbers against Phase 7/8's HyperTraPS-based prediction "
          "once it's available -- HyperTraPS needs to clear the best of these to "
          "justify its added complexity, per Phase 6's own mandate.")


if __name__ == "__main__":
    df = load_matrix(MATRIX_PATH)
    task_df = build_task(df)
    run_baselines(task_df)
    task_df.to_csv("phase6_out_task_table.csv", index=False)
    print("\nWrote lineage-level task table to phase6_out_task_table.csv for inspection.")