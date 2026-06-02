#!/usr/bin/env python3
"""
Script 180: Semi-Supervised Learning for Fibrosis Prediction
=============================================================

Hypothesis: Leveraging unlabeled samples to learn better feature representations
can improve fibrosis (fib_ge3) prediction compared to supervised-only baselines.

Data structure (1,444 total samples):
  - 718 labeled in LOCO folds (6 cohort folds, used for evaluation)
  - 410 labeled but excluded from LOCO folds (GSE213621 subset; the earlier
    PRJNA512027 contribution was retired 2026-05-15)
  - 316 truly unlabeled (no fibrosis annotation)

Semi-supervised opportunity:
  For each LOCO test fold, training uses:
    Supervised-only: other 5 LOCO labeled folds (~600 samples)
    Semi-supervised:  other 5 LOCO labeled folds + 410 excluded-labeled + 316 unlabeled

Three semi-supervised methods compared:
  1. Self-Training: iterative pseudo-labeling with confidence threshold
  2. Label Spreading: sklearn graph-based propagation on kNN features
  3. Co-Training: two-view learning (cell-type vs expression views)

Baseline (from Script 161): M2_celltype supervised elastic net AUROC = 0.799

Validation:
  - LOCO-CV: test fold is NEVER used in any semi-supervised training step
  - Same 6 LOCO folds as Script 161 for direct comparison
  - Reports delta vs supervised baseline per fold

Outputs to results/novel_ml/semi_supervised/
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from copy import deepcopy

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, average_precision_score, brier_score_loss,
    precision_recall_curve
)
from sklearn.semi_supervised import LabelSpreading
from sklearn.neighbors import kneighbors_graph

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
try:
    from sklearn.exceptions import ConvergenceWarning
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
INTEGRATION = BASE / "RNA-seq/Human/Patient_Cohorts/analysis/integration"
RESULTS = INTEGRATION / "results/novel_ml/semi_supervised"

FEATURE_FILE = INTEGRATION / "results/prognosis_v2/feature_matrix_v2.csv"
LABEL_FILE = INTEGRATION / "results/prognosis_v2/labels_v2.csv"

RANDOM_STATE = 42
MAX_ITER = 5000
SELF_TRAIN_ROUNDS = 5
SELF_TRAIN_THRESHOLD = 0.9  # P > 0.9 or P < 0.1 for pseudo-labels
LABEL_SPREAD_K = 30
LABEL_SPREAD_ALPHA = 0.2  # clamping factor (higher = more weight to propagated labels)

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def compute_metrics(y_true, y_prob):
    """Compute AUROC, AUPRC, Brier score."""
    results = {}
    try:
        results["auroc"] = roc_auc_score(y_true, y_prob)
    except ValueError:
        results["auroc"] = np.nan
    try:
        results["auprc"] = average_precision_score(y_true, y_prob)
    except ValueError:
        results["auprc"] = np.nan
    try:
        results["brier"] = brier_score_loss(y_true, y_prob)
    except ValueError:
        results["brier"] = np.nan
    return results


def fit_supervised_baseline(X_train, y_train, X_test, feature_cols):
    """Fit supervised elastic net on labeled data only. Exact match to Script 161 M2."""
    scaler = StandardScaler().fit(X_train[feature_cols])
    X_tr = scaler.transform(X_train[feature_cols])
    X_te = scaler.transform(X_test[feature_cols])

    model = LogisticRegressionCV(
        penalty="elasticnet",
        l1_ratios=[0.5],
        solver="saga",
        max_iter=MAX_ITER,
        class_weight="balanced",
        cv=5,
        scoring="roc_auc",
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )
    model.fit(X_tr, y_train)
    y_prob = model.predict_proba(X_te)[:, 1]
    return y_prob, model


# ---------------------------------------------------------------------------
# Method 1: Self-Training (Iterative Pseudo-Labeling)
# ---------------------------------------------------------------------------
def self_training(X_labeled, y_labeled, X_unlabeled, X_test, feature_cols,
                  n_rounds=SELF_TRAIN_ROUNDS, threshold=SELF_TRAIN_THRESHOLD):
    """
    Self-training loop:
    1. Train elastic net on labeled data
    2. Predict on unlabeled data
    3. Add high-confidence pseudo-labels (P > threshold or P < 1-threshold)
    4. Retrain on expanded set
    5. Repeat n_rounds times

    Returns: test predictions, history of added samples per round
    """
    scaler = StandardScaler().fit(X_labeled[feature_cols])

    # Working copies
    X_lab = X_labeled[feature_cols].copy()
    y_lab = y_labeled.copy()
    X_unlab = X_unlabeled[feature_cols].copy()
    unlab_indices = X_unlabeled.index.tolist()

    history = []

    for round_i in range(n_rounds):
        # Scale all data using current labeled set statistics
        scaler = StandardScaler().fit(X_lab)
        X_tr_scaled = scaler.transform(X_lab)

        # Train model
        model = LogisticRegressionCV(
            penalty="elasticnet",
            l1_ratios=[0.5],
            solver="saga",
            max_iter=MAX_ITER,
            class_weight="balanced",
            cv=min(5, max(2, len(y_lab) // 10)),
            scoring="roc_auc",
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )

        # Check we have both classes
        if len(np.unique(y_lab)) < 2:
            break

        model.fit(X_tr_scaled, y_lab)

        # Predict on remaining unlabeled
        if len(X_unlab) == 0:
            break

        X_unlab_scaled = scaler.transform(X_unlab)
        probs = model.predict_proba(X_unlab_scaled)[:, 1]

        # Select high-confidence predictions
        confident_pos = probs > threshold
        confident_neg = probs < (1.0 - threshold)
        confident_mask = confident_pos | confident_neg

        n_added = confident_mask.sum()
        history.append({
            "round": round_i + 1,
            "n_labeled": len(y_lab),
            "n_unlabeled_remaining": len(X_unlab),
            "n_added": int(n_added),
            "n_confident_pos": int(confident_pos.sum()),
            "n_confident_neg": int(confident_neg.sum()),
            "threshold": threshold,
        })

        if n_added == 0:
            break

        # Create pseudo-labels
        pseudo_labels = (probs > 0.5).astype(int)

        # Add confident samples to labeled set
        confident_idx = X_unlab.index[confident_mask]
        X_lab = pd.concat([X_lab, X_unlab.loc[confident_idx]])
        y_lab = pd.concat([
            y_lab,
            pd.Series(pseudo_labels[confident_mask], index=confident_idx)
        ])

        # Remove from unlabeled
        X_unlab = X_unlab.loc[~confident_mask]

    # Final model: retrain on expanded labeled set
    scaler_final = StandardScaler().fit(X_lab)
    X_tr_final = scaler_final.transform(X_lab)
    X_te_final = scaler_final.transform(X_test[feature_cols])

    model_final = LogisticRegressionCV(
        penalty="elasticnet",
        l1_ratios=[0.5],
        solver="saga",
        max_iter=MAX_ITER,
        class_weight="balanced",
        cv=min(5, max(2, len(y_lab) // 10)),
        scoring="roc_auc",
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )

    if len(np.unique(y_lab)) < 2:
        # Fall back to majority class prediction
        y_prob = np.full(len(X_test), y_lab.mean())
    else:
        model_final.fit(X_tr_final, y_lab)
        y_prob = model_final.predict_proba(X_te_final)[:, 1]

    return y_prob, history


# ---------------------------------------------------------------------------
# Method 2: Label Spreading on kNN Graph
# ---------------------------------------------------------------------------
def label_spreading(X_labeled, y_labeled, X_unlabeled, X_test, feature_cols,
                    k=LABEL_SPREAD_K, alpha=LABEL_SPREAD_ALPHA):
    """
    Label propagation via sklearn LabelSpreading:
    1. Build kNN graph on ALL samples (labeled + unlabeled)
    2. Propagate labels from labeled to unlabeled
    3. Use propagated soft labels to train a final supervised model

    The key insight is that the graph structure captures manifold geometry
    from ALL 1,444 samples, potentially revealing clusters that align with
    fibrosis status.

    Returns: test predictions, propagation statistics
    """
    # Combine labeled + unlabeled for graph construction
    X_all = pd.concat([X_labeled[feature_cols], X_unlabeled[feature_cols]])
    n_labeled = len(X_labeled)
    n_unlabeled = len(X_unlabeled)

    # Scale features
    scaler = StandardScaler().fit(X_all)
    X_all_scaled = scaler.transform(X_all)

    # Prepare labels: -1 for unlabeled
    y_all = np.full(len(X_all), -1)
    y_all[:n_labeled] = y_labeled.values.astype(int)

    # Adjust k to not exceed number of samples
    effective_k = min(k, len(X_all) - 1)

    # LabelSpreading with RBF kernel
    label_spread = LabelSpreading(
        kernel="knn",
        n_neighbors=effective_k,
        alpha=alpha,
        max_iter=100,
        n_jobs=N_JOBS,
    )
    label_spread.fit(X_all_scaled, y_all)

    # Get soft labels (probabilities) for unlabeled samples
    proba = label_spread.label_distributions_
    unlabeled_proba = proba[n_labeled:, 1]  # P(class=1) for unlabeled

    # Statistics about propagation
    stats = {
        "n_labeled": n_labeled,
        "n_unlabeled": n_unlabeled,
        "k": effective_k,
        "alpha": alpha,
        "mean_propagated_prob": float(np.mean(unlabeled_proba)),
        "std_propagated_prob": float(np.std(unlabeled_proba)),
        "n_propagated_positive": int((unlabeled_proba > 0.5).sum()),
        "n_propagated_negative": int((unlabeled_proba <= 0.5).sum()),
    }

    # Now train a supervised model on labeled + pseudo-labeled data
    # Use hard pseudo-labels (thresholded at 0.5)
    pseudo_labels = (unlabeled_proba > 0.5).astype(int)

    X_expanded = pd.concat([X_labeled[feature_cols], X_unlabeled[feature_cols]])
    y_expanded = np.concatenate([y_labeled.values.astype(int), pseudo_labels])

    scaler_final = StandardScaler().fit(X_expanded)
    X_tr_final = scaler_final.transform(X_expanded)
    X_te_final = scaler_final.transform(X_test[feature_cols])

    model = LogisticRegressionCV(
        penalty="elasticnet",
        l1_ratios=[0.5],
        solver="saga",
        max_iter=MAX_ITER,
        class_weight="balanced",
        cv=5,
        scoring="roc_auc",
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )
    model.fit(X_tr_final, y_expanded)
    y_prob = model.predict_proba(X_te_final)[:, 1]

    return y_prob, stats


# ---------------------------------------------------------------------------
# Method 3: Co-Training (Two Views)
# ---------------------------------------------------------------------------
def co_training(X_labeled, y_labeled, X_unlabeled, X_test,
                view_a_cols, view_b_cols,
                n_rounds=SELF_TRAIN_ROUNDS, n_add_per_round=10):
    """
    Co-training with two independent feature views:
      View A: Cell-type features (17 dims) — from scRNA reference deconvolution
      View B: Divergence gene expression (100 dims) — bulk expression-derived

    The assumption: if View A is confident about a sample, that pseudo-label
    is informative for View B, and vice versa. This works when the two views
    are conditionally independent given the label.

    Algorithm:
    1. Train model_A on labeled data using View A
    2. Train model_B on labeled data using View B
    3. model_A predicts unlabeled -> top n confident added to model_B training
    4. model_B predicts unlabeled -> top n confident added to model_A training
    5. Repeat

    Returns: test predictions (average of both views), history
    """
    # Working labeled sets for each view (start identical)
    lab_idx_a = set(X_labeled.index.tolist())
    lab_idx_b = set(X_labeled.index.tolist())
    pseudo_labels = {}  # index -> pseudo-label

    # All labeled data
    y_dict = dict(zip(X_labeled.index, y_labeled.values.astype(int)))

    history = []

    for round_i in range(n_rounds):
        # Current training sets
        idx_a = sorted(lab_idx_a)
        idx_b = sorted(lab_idx_b)

        y_a = np.array([y_dict.get(i, pseudo_labels.get(i, 0)) for i in idx_a])
        y_b = np.array([y_dict.get(i, pseudo_labels.get(i, 0)) for i in idx_b])

        # Check both classes present
        if len(np.unique(y_a)) < 2 or len(np.unique(y_b)) < 2:
            break

        # Scale and train View A
        X_a_train = pd.concat([X_labeled, X_unlabeled]).loc[idx_a, view_a_cols]
        scaler_a = StandardScaler().fit(X_a_train)
        model_a = LogisticRegression(
            penalty="elasticnet", l1_ratio=0.5, solver="saga",
            max_iter=MAX_ITER, class_weight="balanced",
            C=1.0, random_state=RANDOM_STATE,
        )
        model_a.fit(scaler_a.transform(X_a_train), y_a)

        # Scale and train View B
        X_b_train = pd.concat([X_labeled, X_unlabeled]).loc[idx_b, view_b_cols]
        scaler_b = StandardScaler().fit(X_b_train)
        model_b = LogisticRegression(
            penalty="elasticnet", l1_ratio=0.5, solver="saga",
            max_iter=MAX_ITER, class_weight="balanced",
            C=1.0, random_state=RANDOM_STATE,
        )
        model_b.fit(scaler_b.transform(X_b_train), y_b)

        # Find remaining unlabeled for each view
        remaining_for_b = sorted(set(X_unlabeled.index) - lab_idx_b)
        remaining_for_a = sorted(set(X_unlabeled.index) - lab_idx_a)

        n_added_a = 0
        n_added_b = 0

        # View A predicts unlabeled -> add confident to View B's training set
        if len(remaining_for_b) > 0:
            X_rem_a = X_unlabeled.loc[remaining_for_b, view_a_cols]
            probs_a = model_a.predict_proba(scaler_a.transform(X_rem_a))[:, 1]

            # Select top confident from each class
            n_each = max(1, n_add_per_round // 2)
            # Most confident positives
            pos_order = np.argsort(-probs_a)
            neg_order = np.argsort(probs_a)

            added_this_round = []
            for idx_pos in pos_order[:n_each]:
                if probs_a[idx_pos] > 0.5:
                    real_idx = remaining_for_b[idx_pos]
                    pseudo_labels[real_idx] = 1
                    lab_idx_b.add(real_idx)
                    added_this_round.append(real_idx)
                    n_added_b += 1

            for idx_neg in neg_order[:n_each]:
                if probs_a[idx_neg] < 0.5:
                    real_idx = remaining_for_b[idx_neg]
                    pseudo_labels[real_idx] = 0
                    lab_idx_b.add(real_idx)
                    added_this_round.append(real_idx)
                    n_added_b += 1

        # View B predicts unlabeled -> add confident to View A's training set
        if len(remaining_for_a) > 0:
            X_rem_b = X_unlabeled.loc[remaining_for_a, view_b_cols]
            probs_b = model_b.predict_proba(scaler_b.transform(X_rem_b))[:, 1]

            n_each = max(1, n_add_per_round // 2)
            pos_order = np.argsort(-probs_b)
            neg_order = np.argsort(probs_b)

            for idx_pos in pos_order[:n_each]:
                if probs_b[idx_pos] > 0.5:
                    real_idx = remaining_for_a[idx_pos]
                    pseudo_labels[real_idx] = 1
                    lab_idx_a.add(real_idx)
                    n_added_a += 1

            for idx_neg in neg_order[:n_each]:
                if probs_b[idx_neg] < 0.5:
                    real_idx = remaining_for_a[idx_neg]
                    pseudo_labels[real_idx] = 0
                    lab_idx_a.add(real_idx)
                    n_added_a += 1

        history.append({
            "round": round_i + 1,
            "n_labeled_view_a": len(lab_idx_a),
            "n_labeled_view_b": len(lab_idx_b),
            "n_added_to_a": n_added_a,
            "n_added_to_b": n_added_b,
        })

        if n_added_a == 0 and n_added_b == 0:
            break

    # Final prediction: retrain both views on their expanded sets and average
    idx_a = sorted(lab_idx_a)
    idx_b = sorted(lab_idx_b)
    y_a = np.array([y_dict.get(i, pseudo_labels.get(i, 0)) for i in idx_a])
    y_b = np.array([y_dict.get(i, pseudo_labels.get(i, 0)) for i in idx_b])

    all_data = pd.concat([X_labeled, X_unlabeled])

    # View A final
    X_a_final = all_data.loc[idx_a, view_a_cols]
    scaler_a_f = StandardScaler().fit(X_a_final)
    model_a_f = LogisticRegressionCV(
        penalty="elasticnet", l1_ratios=[0.5], solver="saga",
        max_iter=MAX_ITER, class_weight="balanced",
        cv=min(5, max(2, len(y_a) // 10)),
        scoring="roc_auc", random_state=RANDOM_STATE, n_jobs=-1,
    )
    if len(np.unique(y_a)) >= 2:
        model_a_f.fit(scaler_a_f.transform(X_a_final), y_a)
        prob_a = model_a_f.predict_proba(
            scaler_a_f.transform(X_test[view_a_cols])
        )[:, 1]
    else:
        prob_a = np.full(len(X_test), np.mean(y_a))

    # View B final
    X_b_final = all_data.loc[idx_b, view_b_cols]
    scaler_b_f = StandardScaler().fit(X_b_final)
    model_b_f = LogisticRegressionCV(
        penalty="elasticnet", l1_ratios=[0.5], solver="saga",
        max_iter=MAX_ITER, class_weight="balanced",
        cv=min(5, max(2, len(y_b) // 10)),
        scoring="roc_auc", random_state=RANDOM_STATE, n_jobs=-1,
    )
    if len(np.unique(y_b)) >= 2:
        model_b_f.fit(scaler_b_f.transform(X_b_final), y_b)
        prob_b = model_b_f.predict_proba(
            scaler_b_f.transform(X_test[view_b_cols])
        )[:, 1]
    else:
        prob_b = np.full(len(X_test), np.mean(y_b))

    # Average both views
    y_prob = 0.5 * prob_a + 0.5 * prob_b

    return y_prob, history


# ---------------------------------------------------------------------------
# Method 1b: Self-Training with Excluded-Labeled + Truly Unlabeled
# ---------------------------------------------------------------------------
def self_training_with_excluded_labeled(
    X_fold_train, y_fold_train,
    X_excl_labeled, y_excl_labeled,
    X_truly_unlabeled,
    X_test, feature_cols,
    n_rounds=SELF_TRAIN_ROUNDS, threshold=SELF_TRAIN_THRESHOLD
):
    """
    Enhanced self-training that first incorporates excluded-but-labeled samples,
    then applies pseudo-labeling to truly unlabeled samples.

    This separates two effects:
    A) Adding more labeled data (excluded-labeled) — supervised augmentation
    B) Semi-supervised learning on truly unlabeled data
    """
    # Step 1: Combine LOCO-train + excluded-labeled (all have ground truth)
    X_all_labeled = pd.concat([X_fold_train[feature_cols], X_excl_labeled[feature_cols]])
    y_all_labeled = pd.concat([y_fold_train, y_excl_labeled])

    # Step 2: Self-train using truly unlabeled only
    y_prob, history = self_training(
        X_all_labeled.reset_index(drop=True),
        y_all_labeled.reset_index(drop=True),
        X_truly_unlabeled, X_test, feature_cols,
        n_rounds=n_rounds, threshold=threshold,
    )
    return y_prob, history


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    RESULTS.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("[180] Semi-Supervised Learning for Fibrosis Prediction")
    print("=" * 70)
    print(f"Output: {RESULTS}")
    print(f"Jobs:   {N_JOBS}")
    print()

    # Load data
    features = pd.read_csv(FEATURE_FILE)
    labels = pd.read_csv(LABEL_FILE)
    df = features.merge(labels, on="sample_id", how="inner")
    print(f"[180] Merged: {len(df)} samples, {df.shape[1]} columns")

    # Identify feature columns
    ct_cols = [c for c in features.columns if c.startswith("ct_")]
    div_cols = [c for c in features.columns if c.startswith("div_")]
    print(f"[180] Features: ct={len(ct_cols)}, div={len(div_cols)}")

    # -----------------------------------------------------------------------
    # Partition samples
    # -----------------------------------------------------------------------
    # LOCO fold samples (718 with labels)
    loco_mask = (df["loco_fold_fibrosis"] != "excluded") & df["loco_fold_fibrosis"].notna()
    loco_df = df[loco_mask].copy()
    fold_names = sorted(loco_df["loco_fold_fibrosis"].unique())

    # Excluded samples (726 total)
    excl_df = df[~loco_mask].copy()

    # Among excluded: those with labels vs truly unlabeled
    excl_labeled_mask = excl_df["fib_ge3"].notna()
    excl_labeled = excl_df[excl_labeled_mask].copy()
    truly_unlabeled = excl_df[~excl_labeled_mask].copy()

    print(f"\n[180] Sample partitions:")
    print(f"  LOCO fold labeled:   {len(loco_df)}")
    print(f"  Excluded w/ labels:  {len(excl_labeled)} (fib_ge3: {int(excl_labeled['fib_ge3'].sum())}/{len(excl_labeled)})")
    print(f"  Truly unlabeled:     {len(truly_unlabeled)}")
    print(f"  LOCO folds:          {fold_names}")
    print()

    # -----------------------------------------------------------------------
    # Methods to evaluate
    # -----------------------------------------------------------------------
    methods = [
        "supervised_only",             # Baseline: M2 on LOCO train only
        "supervised_augmented",        # +410 excluded-labeled (pure supervised)
        "self_training",               # Self-training on 726 unlabeled/excluded
        "self_training_split",         # Self-train with excluded-labeled + truly unlabeled
        "label_spreading",             # Graph-based propagation
        "co_training",                 # Two-view co-training
    ]

    all_results = []
    all_predictions = []
    all_histories = []

    target_col = "fib_ge3"

    for fold in fold_names:
        print(f"\n{'='*60}")
        print(f"[180] LOCO Fold: {fold}")
        print(f"{'='*60}")

        # Split LOCO
        test_mask = loco_df["loco_fold_fibrosis"] == fold
        train_loco = loco_df[~test_mask].copy()
        test_loco = loco_df[test_mask].copy()

        # Filter to valid fib_ge3
        train_loco = train_loco[train_loco[target_col].notna()].copy()
        test_loco = test_loco[test_loco[target_col].notna()].copy()

        if len(test_loco) == 0 or len(train_loco) == 0:
            print(f"  SKIP: no valid samples")
            continue

        y_train = train_loco[target_col].astype(int)
        y_test = test_loco[target_col].astype(int)

        if len(y_train.unique()) < 2 or len(y_test.unique()) < 2:
            print(f"  SKIP: single class in train or test")
            continue

        y_excl = excl_labeled[target_col].astype(int) if len(excl_labeled) > 0 else pd.Series(dtype=int)

        print(f"  Train (LOCO): {len(train_loco)} ({int(y_train.sum())} pos)")
        print(f"  Test:         {len(test_loco)} ({int(y_test.sum())} pos)")
        print(f"  Excl labeled: {len(excl_labeled)} ({int(y_excl.sum()) if len(y_excl) > 0 else 0} pos)")
        print(f"  Truly unlab:  {len(truly_unlabeled)}")

        # --- All unlabeled for self-training = excluded (labeled+unlabeled combined) ---
        # For self-training, treat ALL 726 excluded samples as unlabeled
        all_unlabeled_for_self_train = excl_df.copy()

        for method in methods:
            print(f"\n  --- {method} ---")
            t_method = time.time()

            try:
                if method == "supervised_only":
                    # Exact match to Script 161 M2_celltype
                    y_prob, _ = fit_supervised_baseline(
                        train_loco, y_train, test_loco, ct_cols
                    )

                elif method == "supervised_augmented":
                    # Add excluded-labeled to training (pure supervised, no semi-supervised)
                    X_aug = pd.concat([train_loco, excl_labeled])
                    y_aug = pd.concat([y_train, y_excl])
                    y_prob, _ = fit_supervised_baseline(
                        X_aug, y_aug, test_loco, ct_cols
                    )

                elif method == "self_training":
                    # Self-training: start with LOCO train, pseudo-label ALL 726 excluded
                    y_prob, hist = self_training(
                        train_loco, y_train, all_unlabeled_for_self_train,
                        test_loco, ct_cols,
                    )
                    for h in hist:
                        h.update({"method": method, "fold": fold})
                    all_histories.extend(hist)

                elif method == "self_training_split":
                    # Self-training: use excluded-labeled as extra labeled, self-train on truly unlabeled
                    y_prob, hist = self_training_with_excluded_labeled(
                        train_loco, y_train,
                        excl_labeled, y_excl,
                        truly_unlabeled,
                        test_loco, ct_cols,
                    )
                    for h in hist:
                        h.update({"method": method, "fold": fold})
                    all_histories.extend(hist)

                elif method == "label_spreading":
                    # Label propagation on kNN graph
                    # Labeled pool = LOCO train + excluded-labeled
                    X_lab_all = pd.concat([train_loco, excl_labeled])
                    y_lab_all = pd.concat([y_train, y_excl])
                    y_prob, stats = label_spreading(
                        X_lab_all, y_lab_all,
                        truly_unlabeled,
                        test_loco, ct_cols,
                    )
                    stats.update({"method": method, "fold": fold})
                    all_histories.append(stats)

                elif method == "co_training":
                    # Co-training with two views
                    y_prob, hist = co_training(
                        train_loco, y_train, all_unlabeled_for_self_train,
                        test_loco, view_a_cols=ct_cols, view_b_cols=div_cols,
                    )
                    for h in hist:
                        h.update({"method": method, "fold": fold})
                    all_histories.extend(hist)

                else:
                    raise ValueError(f"Unknown method: {method}")

            except Exception as e:
                print(f"    ERROR: {e}")
                y_prob = np.full(len(test_loco), 0.5)

            # Compute metrics
            metrics = compute_metrics(y_test.values, y_prob)
            elapsed_method = time.time() - t_method

            print(
                f"    AUROC={metrics['auroc']:.3f}, "
                f"AUPRC={metrics['auprc']:.3f}, "
                f"Brier={metrics['brier']:.3f} "
                f"({elapsed_method:.1f}s)"
            )

            all_results.append({
                "method": method,
                "fold": fold,
                "n_train_labeled": len(train_loco),
                "n_excl_labeled": len(excl_labeled),
                "n_truly_unlabeled": len(truly_unlabeled),
                "n_test": len(test_loco),
                "n_pos_test": int(y_test.sum()),
                "auroc": metrics["auroc"],
                "auprc": metrics["auprc"],
                "brier": metrics["brier"],
                "elapsed_s": elapsed_method,
            })

            for sid, prob, true_label in zip(
                test_loco["sample_id"].values, y_prob, y_test.values
            ):
                all_predictions.append({
                    "sample_id": sid,
                    "method": method,
                    "predicted_prob": prob,
                    "true_label": int(true_label),
                    "fold": fold,
                })

    # -----------------------------------------------------------------------
    # Save outputs
    # -----------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("[180] Saving outputs...")
    print(f"{'='*70}")

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(RESULTS / "semi_supervised_results.csv", index=False)
    print(f"  semi_supervised_results.csv: {len(results_df)} rows")

    pred_df = pd.DataFrame(all_predictions)
    pred_df.to_csv(RESULTS / "semi_supervised_predictions.csv", index=False)
    print(f"  semi_supervised_predictions.csv: {len(pred_df)} rows")

    if all_histories:
        hist_df = pd.DataFrame(all_histories)
        hist_df.to_csv(RESULTS / "semi_supervised_histories.csv", index=False)
        print(f"  semi_supervised_histories.csv: {len(hist_df)} rows")

    # -----------------------------------------------------------------------
    # Summary: mean ± std across folds, delta vs supervised baseline
    # -----------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("[180] SUMMARY")
    print(f"{'='*70}")
    print(f"{'Method':<28s} {'AUROC':>14s} {'AUPRC':>14s} {'Brier':>14s} {'Delta AUROC':>12s}")
    print("-" * 85)

    summary_rows = []
    supervised_mean = None

    for method in methods:
        sub = results_df[results_df["method"] == method]
        if len(sub) == 0:
            continue

        mean_auroc = sub["auroc"].mean()
        std_auroc = sub["auroc"].std()
        mean_auprc = sub["auprc"].mean()
        std_auprc = sub["auprc"].std()
        mean_brier = sub["brier"].mean()
        std_brier = sub["brier"].std()

        if method == "supervised_only":
            supervised_mean = mean_auroc
            delta = 0.0
        else:
            delta = mean_auroc - supervised_mean if supervised_mean is not None else np.nan

        summary_rows.append({
            "method": method,
            "n_folds": len(sub),
            "mean_auroc": mean_auroc,
            "std_auroc": std_auroc,
            "mean_auprc": mean_auprc,
            "std_auprc": std_auprc,
            "mean_brier": mean_brier,
            "std_brier": std_brier,
            "delta_auroc_vs_supervised": delta,
        })

        sign = "+" if delta > 0 else ""
        print(
            f"  {method:<26s} {mean_auroc:.3f}±{std_auroc:.3f}  "
            f"{mean_auprc:.3f}±{std_auprc:.3f}  "
            f"{mean_brier:.3f}±{std_brier:.3f}  "
            f"{sign}{delta:.3f}"
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(RESULTS / "semi_supervised_summary.csv", index=False)
    print(f"\n  semi_supervised_summary.csv: {len(summary_df)} rows")

    # -----------------------------------------------------------------------
    # Interpretation guidance
    # -----------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("[180] INTERPRETATION")
    print(f"{'='*70}")

    if supervised_mean is not None:
        best_row = summary_df.loc[summary_df["mean_auroc"].idxmax()]
        best_method = best_row["method"]
        best_auroc = best_row["mean_auroc"]
        best_delta = best_row["delta_auroc_vs_supervised"]

        print(f"  Supervised-only (M2 celltype): AUROC = {supervised_mean:.3f}")
        print(f"  Best method: {best_method} (AUROC = {best_auroc:.3f}, delta = {best_delta:+.3f})")
        print()

        if best_delta > 0.02:
            print("  CONCLUSION: Semi-supervised learning provides meaningful improvement.")
            print("  The 726 unlabeled samples help learn better feature geometry.")
        elif best_delta > 0.005:
            print("  CONCLUSION: Marginal improvement from semi-supervised learning.")
            print("  Unlabeled data provides slight benefit, but effect is small.")
        elif best_delta > -0.005:
            print("  CONCLUSION: Unlabeled data does not improve fibrosis prediction.")
            print("  Consistent with the transcriptome-wide signal ceiling — the 17")
            print("  cell-type features already capture the decision boundary well.")
        else:
            print("  CONCLUSION: Semi-supervised learning HURTS performance.")
            print("  Pseudo-labels may introduce noise that degrades the classifier.")
            print("  The labeled sample size (718) is sufficient for this feature space.")

        # Check supervised_augmented specifically
        aug_row = summary_df[summary_df["method"] == "supervised_augmented"]
        if len(aug_row) > 0:
            aug_delta = aug_row.iloc[0]["delta_auroc_vs_supervised"]
            print(f"\n  Supervised augmented (+ excluded labeled): delta = {aug_delta:+.3f}")
            if aug_delta > 0.01:
                print("  -> Adding more labeled data helps; improvement is from data quantity, not SSL.")
            elif aug_delta < -0.01:
                print("  -> Adding excluded cohorts HURTS: domain shift between LOCO and excluded datasets.")
            else:
                print("  -> No benefit from adding excluded labeled data: sample size is not the bottleneck.")

    elapsed = time.time() - t0
    print(f"\n[180] Done in {elapsed/60:.1f} min")


if __name__ == "__main__":
    main()
