#!/usr/bin/env python3
"""
Script 181: Ordinal CORAL Model in the v2 Framework
=====================================================
Predicts F0-F4 fibrosis as ordinal (not binary) using the same LOCO-CV
structure as Script 161.

Two approaches:
  Option A: Cumulative Binary Classifiers (sklearn)
    - Fit 4 binary classifiers: P(F>=1), P(F>=2), P(F>=3), P(F>=4)
    - Stack predictions -> class probabilities via differencing
    - Monotonicity enforced via isotonic clipping

  Option B: CORAL Neural Network (PyTorch, from-scratch implementation)
    - Small MLP: input -> 64 -> 32 -> CoralLayer(32, 4)
    - CORAL loss enforces ordinal consistency
    - Early stopping on inner-CV QWK

Feature configs (matching Script 161):
  M2_celltype    — bp_*/ct_* cell-type fractions (17 features)
  M3_expression  — div_* divergence genes (inner-CV-selected)
  M4_combined    — M2 + inner-CV-selected expression genes

Primary metrics:
  - QWK (quadratic weighted kappa) — primary ordinal metric
  - Adjacent accuracy (predicted within +/-1 stage)
  - Per-class accuracy
  - Macro-averaged AUROC (one-vs-rest)

Comparison targets:
  - v1 concept bottleneck QWK=0.624
  - v1 ordinal elastic net QWK=0.550
  - v2 binary fib_ge3 AUROC from Script 161

Outputs to results/novel_ml/ordinal_coral_v2/
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import cohen_kappa_score, accuracy_score
from joblib import Parallel, delayed

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
INTEGRATION = BASE / "RNA-seq/Human/Patient_Cohorts/analysis/integration"
V2_RESULTS = INTEGRATION / "results/prognosis_v2"
OUT_DIR = INTEGRATION / "results/novel_ml/ordinal_coral_v2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FEATURE_FILE = V2_RESULTS / "feature_matrix_v2.csv"
LABEL_FILE = V2_RESULTS / "labels_v2.csv"

N_INNER_FOLDS = 5
RANDOM_STATE = 42
MAX_ITER = 5000
N_CLASSES = 5          # F0, F1, F2, F3, F4
N_THRESHOLDS = 4       # K-1 cumulative thresholds

# Inner-CV feature selection grid (same as Script 161)
N_GENES_GRID = [10, 25, 50, 100]
L1_RATIO_GRID = [0.1, 0.5, 0.9]

# CORAL neural net hyperparams
CORAL_HIDDEN1 = 64
CORAL_HIDDEN2 = 32
CORAL_LR = 1e-3
CORAL_EPOCHS = 200
CORAL_PATIENCE = 20
CORAL_BATCH_SIZE = 64
CORAL_WEIGHT_DECAY = 1e-4

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
print(f"[181] Using {N_JOBS} parallel jobs")


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def compute_qwk(y_true, y_pred):
    """Quadratic weighted kappa."""
    return cohen_kappa_score(y_true, y_pred, weights="quadratic")


def compute_adjacent_accuracy(y_true, y_pred):
    """Fraction of predictions within +/-1 of true stage."""
    return np.mean(np.abs(y_true - y_pred) <= 1)


def compute_per_class_accuracy(y_true, y_pred, n_classes=N_CLASSES):
    """Per-class accuracy dict."""
    acc = {}
    for c in range(n_classes):
        mask = y_true == c
        if mask.sum() == 0:
            acc[f"F{c}"] = np.nan
        else:
            acc[f"F{c}"] = np.mean(y_pred[mask] == c)
    return acc


def compute_macro_auroc(y_true, y_prob_matrix, n_classes=N_CLASSES):
    """Macro-averaged one-vs-rest AUROC from probability matrix."""
    from sklearn.metrics import roc_auc_score
    aurocs = []
    for c in range(n_classes):
        binary = (y_true == c).astype(int)
        if binary.sum() == 0 or binary.sum() == len(binary):
            continue
        try:
            aurocs.append(roc_auc_score(binary, y_prob_matrix[:, c]))
        except ValueError:
            pass
    return np.mean(aurocs) if aurocs else np.nan


def compute_ordinal_metrics(y_true, y_pred, y_prob_matrix=None):
    """Compute all ordinal metrics."""
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)

    qwk = compute_qwk(y_true, y_pred)
    adj_acc = compute_adjacent_accuracy(y_true, y_pred)
    acc = accuracy_score(y_true, y_pred)
    per_class = compute_per_class_accuracy(y_true, y_pred)
    mae = np.mean(np.abs(y_true - y_pred))

    macro_auroc = np.nan
    if y_prob_matrix is not None:
        macro_auroc = compute_macro_auroc(y_true, y_prob_matrix)

    # Also compute binary fib_ge3 AUROC for direct comparison with Script 161
    fib_ge3_auroc = np.nan
    if y_prob_matrix is not None:
        from sklearn.metrics import roc_auc_score
        binary_true = (y_true >= 3).astype(int)
        if binary_true.sum() > 0 and binary_true.sum() < len(binary_true):
            # P(F>=3) = sum of P(F3) + P(F4)
            prob_ge3 = y_prob_matrix[:, 3:].sum(axis=1)
            try:
                fib_ge3_auroc = roc_auc_score(binary_true, prob_ge3)
            except ValueError:
                pass

    return {
        "qwk": qwk,
        "adjacent_accuracy": adj_acc,
        "accuracy": acc,
        "mae": mae,
        "macro_auroc": macro_auroc,
        "fib_ge3_auroc": fib_ge3_auroc,
        **{f"acc_{k}": v for k, v in per_class.items()},
    }


# ---------------------------------------------------------------------------
# Feature selection (from Script 161)
# ---------------------------------------------------------------------------
def compute_univariate_auroc(X, y, feature_cols):
    """Compute AUROC for each feature vs binary F>=3 (ordinal proxy).
    Returns sorted Series."""
    from sklearn.metrics import roc_auc_score
    y_binary = (y >= 3).astype(int)
    aurocs = {}
    for col in feature_cols:
        vals = X[col].values
        if np.std(vals) == 0:
            aurocs[col] = 0.5
            continue
        try:
            auc = roc_auc_score(y_binary, vals)
            aurocs[col] = max(auc, 1 - auc)
        except ValueError:
            aurocs[col] = 0.5
    return pd.Series(aurocs).sort_values(ascending=False)


# ---------------------------------------------------------------------------
# Option A: Cumulative Binary Classifiers
# ---------------------------------------------------------------------------
def cumulative_probs_to_class_probs(cum_probs):
    """
    Convert cumulative probabilities P(Y>=k) to class probabilities P(Y=k).

    cum_probs: array of shape (n_samples, K-1) where K-1=4
               cum_probs[:, j] = P(Y >= j+1)

    Returns: class_probs of shape (n_samples, K) where K=5
    """
    n = cum_probs.shape[0]
    K = cum_probs.shape[1] + 1

    # Enforce monotonicity: P(Y>=1) >= P(Y>=2) >= P(Y>=3) >= P(Y>=4)
    cum_sorted = np.copy(cum_probs)
    for i in range(n):
        for j in range(1, K - 1):
            if cum_sorted[i, j] > cum_sorted[i, j - 1]:
                cum_sorted[i, j] = cum_sorted[i, j - 1]

    # Clip to [0, 1]
    cum_sorted = np.clip(cum_sorted, 0.0, 1.0)

    # P(Y=0) = 1 - P(Y>=1)
    # P(Y=k) = P(Y>=k) - P(Y>=k+1) for k=1,...,K-2
    # P(Y=K-1) = P(Y>=K-1)
    class_probs = np.zeros((n, K))
    class_probs[:, 0] = 1.0 - cum_sorted[:, 0]
    for k in range(1, K - 1):
        class_probs[:, k] = cum_sorted[:, k - 1] - cum_sorted[:, k]
    class_probs[:, K - 1] = cum_sorted[:, K - 2]

    # Safety: clip negative probs and renormalize
    class_probs = np.maximum(class_probs, 0.0)
    row_sums = class_probs.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1e-10)
    class_probs = class_probs / row_sums

    return class_probs


def fit_cumulative_binary(X_train, y_train, X_test, feature_cols,
                          l1_ratio=0.5, C=1.0):
    """
    Fit 4 binary logistic regression classifiers for cumulative thresholds.
    Returns class probabilities (n_test, 5) and coefficients.
    """
    scaler = StandardScaler().fit(X_train[feature_cols])
    X_tr = scaler.transform(X_train[feature_cols])
    X_te = scaler.transform(X_test[feature_cols])

    cum_probs_train = np.zeros((X_tr.shape[0], N_THRESHOLDS))
    cum_probs_test = np.zeros((X_te.shape[0], N_THRESHOLDS))
    all_coefs = {}

    for k in range(N_THRESHOLDS):
        # Binary target: Y >= k+1
        y_binary = (y_train >= (k + 1)).astype(int)

        # Skip if only one class present
        if len(np.unique(y_binary)) < 2:
            # All below threshold or all above — assign constant
            p = y_binary.mean()
            cum_probs_train[:, k] = p
            cum_probs_test[:, k] = p
            continue

        model = LogisticRegression(
            penalty="elasticnet",
            l1_ratio=l1_ratio,
            solver="saga",
            max_iter=MAX_ITER,
            class_weight="balanced",
            C=C,
            random_state=RANDOM_STATE,
        )
        model.fit(X_tr, y_binary)
        cum_probs_train[:, k] = model.predict_proba(X_tr)[:, 1]
        cum_probs_test[:, k] = model.predict_proba(X_te)[:, 1]

        for feat_name, coef_val in zip(feature_cols, model.coef_[0]):
            key = f"threshold_ge{k+1}_{feat_name}"
            all_coefs[key] = coef_val

    class_probs = cumulative_probs_to_class_probs(cum_probs_test)
    return class_probs, all_coefs


def fit_cumulative_cv(X_train, y_train, X_test, feature_cols):
    """
    Fit cumulative binary classifiers with inner-CV for C selection.
    """
    scaler = StandardScaler().fit(X_train[feature_cols])
    X_tr = scaler.transform(X_train[feature_cols])
    X_te = scaler.transform(X_test[feature_cols])

    cum_probs_test = np.zeros((X_te.shape[0], N_THRESHOLDS))
    all_coefs = {}

    for k in range(N_THRESHOLDS):
        y_binary = (y_train >= (k + 1)).astype(int)

        if len(np.unique(y_binary)) < 2:
            p = y_binary.mean()
            cum_probs_test[:, k] = p
            continue

        model = LogisticRegressionCV(
            penalty="elasticnet",
            l1_ratios=[0.1, 0.5, 0.9],
            solver="saga",
            max_iter=MAX_ITER,
            class_weight="balanced",
            cv=min(N_INNER_FOLDS, min(y_binary.sum(), (1 - y_binary).sum())),
            scoring="roc_auc",
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )
        model.fit(X_tr, y_binary)
        cum_probs_test[:, k] = model.predict_proba(X_te)[:, 1]

        for feat_name, coef_val in zip(feature_cols, model.coef_[0]):
            key = f"threshold_ge{k+1}_{feat_name}"
            all_coefs[key] = coef_val

    class_probs = cumulative_probs_to_class_probs(cum_probs_test)
    return class_probs, all_coefs


def fit_cumulative_expression_inner_cv(X_train, y_train, X_test,
                                       expression_cols, extra_cols=None):
    """
    Nested inner CV for expression-based cumulative models.
    Feature selection tuned on inner folds using QWK (not AUROC).
    """
    best_inner_qwk = -999
    best_params = (25, 0.5)

    inner_cv = StratifiedKFold(
        n_splits=N_INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE
    )

    # Use fib_ge3 (binary) for stratified splitting
    y_binary_strat = (y_train >= 3).astype(int)

    for n_genes in N_GENES_GRID:
        for l1_ratio in L1_RATIO_GRID:
            inner_qwks = []

            for tr_idx, val_idx in inner_cv.split(X_train, y_binary_strat):
                X_inner_tr = X_train.iloc[tr_idx]
                X_inner_val = X_train.iloc[val_idx]
                y_inner_tr = y_train.iloc[tr_idx]
                y_inner_val = y_train.iloc[val_idx]

                # Feature selection on inner train only
                gene_aurocs = compute_univariate_auroc(
                    X_inner_tr, y_inner_tr, expression_cols
                )
                selected = gene_aurocs.head(
                    min(n_genes, len(gene_aurocs))
                ).index.tolist()

                if extra_cols is not None:
                    use_cols = list(extra_cols) + selected
                else:
                    use_cols = selected

                # Fit cumulative binary on inner train
                class_probs, _ = fit_cumulative_binary(
                    X_inner_tr, y_inner_tr, X_inner_val,
                    use_cols, l1_ratio=l1_ratio,
                )
                y_pred = class_probs.argmax(axis=1)

                try:
                    qwk = compute_qwk(y_inner_val.values, y_pred)
                except Exception:
                    qwk = 0.0
                inner_qwks.append(qwk)

            mean_qwk = np.mean(inner_qwks)
            if mean_qwk > best_inner_qwk:
                best_inner_qwk = mean_qwk
                best_params = (n_genes, l1_ratio)

    # Retrain on full outer training set with best params
    best_n, best_l1 = best_params
    gene_aurocs = compute_univariate_auroc(X_train, y_train, expression_cols)
    selected = gene_aurocs.head(min(best_n, len(gene_aurocs))).index.tolist()

    if extra_cols is not None:
        use_cols = list(extra_cols) + selected
    else:
        use_cols = selected

    class_probs, coefs = fit_cumulative_binary(
        X_train, y_train, X_test, use_cols, l1_ratio=best_l1,
    )

    params = {
        "n_genes": best_n,
        "l1_ratio": best_l1,
        "inner_qwk": best_inner_qwk,
    }
    return class_probs, coefs, selected, params


# ---------------------------------------------------------------------------
# Option B: CORAL Neural Network (from-scratch PyTorch implementation)
# ---------------------------------------------------------------------------
class CoralMLP(nn.Module):
    """
    MLP with CORAL ordinal output layer.

    CORAL (Consistent Rank Logits) uses shared weights with
    per-threshold biases to enforce ordinal consistency.
    Architecture: input -> hidden1 -> hidden2 -> coral_logits
    """

    def __init__(self, input_dim, hidden1=CORAL_HIDDEN1,
                 hidden2=CORAL_HIDDEN2, n_classes=N_CLASSES):
        super().__init__()
        self.n_thresholds = n_classes - 1

        self.shared = nn.Sequential(
            nn.Linear(input_dim, hidden1),
            nn.BatchNorm1d(hidden1),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden1, hidden2),
            nn.BatchNorm1d(hidden2),
            nn.ReLU(),
            nn.Dropout(0.2),
        )

        # CORAL: single linear layer + per-threshold bias
        self.coral_weight = nn.Linear(hidden2, 1, bias=False)
        self.coral_bias = nn.Parameter(torch.zeros(self.n_thresholds))

    def forward(self, x):
        h = self.shared(x)
        # Shared weight projection: (batch, 1)
        logit = self.coral_weight(h)
        # Broadcast + per-threshold bias: (batch, n_thresholds)
        logits = logit + self.coral_bias.unsqueeze(0)
        return logits

    def predict_proba(self, x):
        """Return class probabilities (batch, n_classes)."""
        with torch.no_grad():
            logits = self.forward(x)
            cum_probs = torch.sigmoid(logits)  # P(Y > k)
            # Convert to numpy for the differencing
            cum_np = cum_probs.cpu().numpy()
            return cumulative_probs_to_class_probs(cum_np)


def coral_loss_fn(logits, labels, n_classes=N_CLASSES):
    """
    CORAL loss: sum of binary cross-entropy losses for each threshold.

    logits: (batch, K-1) raw logits for P(Y > k)
    labels: (batch,) integer class labels 0..K-1
    """
    n_thresholds = n_classes - 1
    # Create binary targets for each threshold
    # For threshold k: target = 1 if label > k, else 0
    levels = torch.zeros(labels.size(0), n_thresholds,
                         device=labels.device, dtype=torch.float32)
    for k in range(n_thresholds):
        levels[:, k] = (labels > k).float()

    loss = nn.functional.binary_cross_entropy_with_logits(logits, levels)
    return loss


def train_coral_model(X_train_np, y_train_np, X_val_np, y_val_np,
                      input_dim, seed=RANDOM_STATE):
    """
    Train a CORAL MLP with early stopping on validation QWK.
    Returns trained model.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)

    device = torch.device("cpu")  # CPU partition

    model = CoralMLP(input_dim).to(device)
    optimizer = optim.Adam(model.parameters(), lr=CORAL_LR,
                           weight_decay=CORAL_WEIGHT_DECAY)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="max", factor=0.5, patience=10
    )

    # Compute class weights for balanced training
    class_counts = np.bincount(y_train_np, minlength=N_CLASSES).astype(float)
    class_counts = np.maximum(class_counts, 1.0)
    class_weights = 1.0 / class_counts
    class_weights = class_weights / class_weights.sum() * N_CLASSES
    sample_weights = class_weights[y_train_np]

    # Create data loaders
    X_tr_t = torch.FloatTensor(X_train_np).to(device)
    y_tr_t = torch.LongTensor(y_train_np).to(device)
    w_tr_t = torch.FloatTensor(sample_weights).to(device)

    X_val_t = torch.FloatTensor(X_val_np).to(device)
    y_val_t = torch.LongTensor(y_val_np).to(device)

    train_ds = TensorDataset(X_tr_t, y_tr_t, w_tr_t)
    train_loader = DataLoader(train_ds, batch_size=CORAL_BATCH_SIZE,
                              shuffle=True, drop_last=False)

    best_qwk = -999
    best_state = None
    patience_counter = 0

    for epoch in range(CORAL_EPOCHS):
        model.train()
        epoch_loss = 0.0

        for X_batch, y_batch, w_batch in train_loader:
            optimizer.zero_grad()
            logits = model(X_batch)

            # Weighted CORAL loss
            n_thresholds = N_CLASSES - 1
            levels = torch.zeros(y_batch.size(0), n_thresholds,
                                 device=device, dtype=torch.float32)
            for k in range(n_thresholds):
                levels[:, k] = (y_batch > k).float()

            # Per-sample BCE
            per_sample_loss = nn.functional.binary_cross_entropy_with_logits(
                logits, levels, reduction="none"
            ).mean(dim=1)
            loss = (per_sample_loss * w_batch).mean()

            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        # Validation
        model.eval()
        with torch.no_grad():
            val_probs = model.predict_proba(X_val_t)
            val_pred = val_probs.argmax(axis=1)
            val_qwk = compute_qwk(y_val_np, val_pred)

        scheduler.step(val_qwk)

        if val_qwk > best_qwk:
            best_qwk = val_qwk
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= CORAL_PATIENCE:
                break

    # Restore best model
    if best_state is not None:
        model.load_state_dict(best_state)

    return model, best_qwk


def fit_coral_nn(X_train, y_train, X_test, feature_cols, seed=RANDOM_STATE):
    """
    Fit CORAL neural network with inner-CV for early stopping.
    Splits training into train/val for early stopping.
    """
    scaler = StandardScaler().fit(X_train[feature_cols])
    X_tr_np = scaler.transform(X_train[feature_cols]).astype(np.float32)
    X_te_np = scaler.transform(X_test[feature_cols]).astype(np.float32)
    y_tr_np = y_train.values.astype(np.int64)

    # Use stratified split for early-stopping validation
    y_binary_strat = (y_tr_np >= 3).astype(int)
    inner_cv = StratifiedKFold(
        n_splits=N_INNER_FOLDS, shuffle=True, random_state=seed
    )

    # Train N_INNER_FOLDS models and ensemble their predictions
    all_test_probs = []

    for fold_idx, (tr_idx, val_idx) in enumerate(
        inner_cv.split(X_tr_np, y_binary_strat)
    ):
        model, val_qwk = train_coral_model(
            X_tr_np[tr_idx], y_tr_np[tr_idx],
            X_tr_np[val_idx], y_tr_np[val_idx],
            input_dim=X_tr_np.shape[1],
            seed=seed + fold_idx,
        )

        # Predict on test
        model.eval()
        with torch.no_grad():
            X_te_t = torch.FloatTensor(X_te_np)
            test_probs = model.predict_proba(X_te_t)
        all_test_probs.append(test_probs)

    # Average ensemble
    class_probs = np.mean(all_test_probs, axis=0)
    return class_probs


def fit_coral_expression_inner_cv(X_train, y_train, X_test,
                                  expression_cols, extra_cols=None):
    """
    Nested inner CV for feature selection + CORAL neural network.
    Feature selection on inner folds, CORAL ensemble on outer training set.
    """
    best_inner_qwk = -999
    best_n_genes = 25

    inner_cv = StratifiedKFold(
        n_splits=N_INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE
    )
    y_binary_strat = (y_train >= 3).astype(int)

    # Only tune n_genes for CORAL (l1_ratio not applicable)
    for n_genes in N_GENES_GRID:
        inner_qwks = []

        for tr_idx, val_idx in inner_cv.split(X_train, y_binary_strat):
            X_inner_tr = X_train.iloc[tr_idx]
            X_inner_val = X_train.iloc[val_idx]
            y_inner_tr = y_train.iloc[tr_idx]
            y_inner_val = y_train.iloc[val_idx]

            gene_aurocs = compute_univariate_auroc(
                X_inner_tr, y_inner_tr, expression_cols
            )
            selected = gene_aurocs.head(
                min(n_genes, len(gene_aurocs))
            ).index.tolist()

            if extra_cols is not None:
                use_cols = list(extra_cols) + selected
            else:
                use_cols = selected

            # Quick CORAL fit for feature selection (single train/val split
            # within this inner fold, not full ensemble)
            scaler = StandardScaler().fit(X_inner_tr[use_cols])
            X_in_tr = scaler.transform(X_inner_tr[use_cols]).astype(np.float32)
            X_in_val = scaler.transform(X_inner_val[use_cols]).astype(np.float32)
            y_in_tr = y_inner_tr.values.astype(np.int64)

            # Split inner train further for early stopping
            n_es = max(int(0.2 * len(y_in_tr)), 10)
            es_idx = np.random.RandomState(RANDOM_STATE).permutation(len(y_in_tr))
            es_val_idx = es_idx[:n_es]
            es_tr_idx = es_idx[n_es:]

            model, _ = train_coral_model(
                X_in_tr[es_tr_idx], y_in_tr[es_tr_idx],
                X_in_tr[es_val_idx], y_in_tr[es_val_idx],
                input_dim=X_in_tr.shape[1],
                seed=RANDOM_STATE,
            )

            model.eval()
            with torch.no_grad():
                val_probs = model.predict_proba(
                    torch.FloatTensor(X_in_val)
                )
            y_pred = val_probs.argmax(axis=1)
            try:
                qwk = compute_qwk(y_inner_val.values, y_pred)
            except Exception:
                qwk = 0.0
            inner_qwks.append(qwk)

        mean_qwk = np.mean(inner_qwks)
        if mean_qwk > best_inner_qwk:
            best_inner_qwk = mean_qwk
            best_n_genes = n_genes

    # Retrain on full outer training set with best n_genes
    gene_aurocs = compute_univariate_auroc(X_train, y_train, expression_cols)
    selected = gene_aurocs.head(
        min(best_n_genes, len(gene_aurocs))
    ).index.tolist()

    if extra_cols is not None:
        use_cols = list(extra_cols) + selected
    else:
        use_cols = selected

    class_probs = fit_coral_nn(
        X_train, y_train, X_test, use_cols, seed=RANDOM_STATE
    )

    params = {"n_genes": best_n_genes, "inner_qwk": best_inner_qwk}
    return class_probs, selected, params


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print(f"[181] Ordinal CORAL v2 — Loading data")

    features = pd.read_csv(FEATURE_FILE)
    labels = pd.read_csv(LABEL_FILE)
    df = features.merge(labels, on="sample_id", how="inner")
    print(f"[181] Merged: {len(df)} samples, {df.shape[1]} columns")

    # Identify feature columns
    div_cols = [c for c in features.columns if c.startswith("div_")]
    ct_cols = [c for c in features.columns
               if c.startswith("ct_") or c.startswith("bp_")]
    print(f"[181] Feature groups: div={len(div_cols)}, ct={len(ct_cols)}")

    # Filter to valid LOCO folds + valid fib_ordinal
    valid = df[
        (df["loco_fold_fibrosis"] != "excluded")
        & df["loco_fold_fibrosis"].notna()
        & df["fib_ordinal"].notna()
    ].copy()
    valid["fib_ordinal"] = valid["fib_ordinal"].astype(int)
    fold_names = sorted(valid["loco_fold_fibrosis"].unique())

    print(f"[181] Valid samples: {len(valid)}, folds: {fold_names}")
    print(f"[181] Class distribution:")
    for stage in range(N_CLASSES):
        n = (valid["fib_ordinal"] == stage).sum()
        print(f"  F{stage}: {n} ({100*n/len(valid):.1f}%)")

    # Feature configs
    config_specs = {
        "M2_celltype": {
            "type": "fixed", "cols": ct_cols,
        },
        "M3_expression": {
            "type": "expression", "expr_cols": div_cols, "extra_cols": None,
        },
        "M4_combined": {
            "type": "expression", "expr_cols": div_cols, "extra_cols": ct_cols,
        },
    }

    # Methods
    methods = ["cumulative_binary", "coral_nn"]

    # Storage
    all_results = []
    all_predictions = []
    all_importance = []

    # -----------------------------------------------------------------------
    # Main LOCO loop
    # -----------------------------------------------------------------------
    for method in methods:
        print(f"\n{'='*70}")
        print(f"[181] METHOD: {method}")
        print(f"{'='*70}")

        for config_name, config in config_specs.items():
            print(f"\n  --- Config: {config_name} ---")
            config_importance = defaultdict(list)

            for fold in fold_names:
                test_mask = valid["loco_fold_fibrosis"] == fold
                train_df = valid[~test_mask].copy()
                test_df = valid[test_mask].copy()

                if len(test_df) == 0 or len(train_df) == 0:
                    print(f"    Fold {fold}: SKIP (empty)")
                    continue

                y_train = train_df["fib_ordinal"].astype(int)
                y_test = test_df["fib_ordinal"].astype(int)

                # Check all classes represented in training
                n_train_classes = len(y_train.unique())
                if n_train_classes < 2:
                    print(f"    Fold {fold}: SKIP (single class in train)")
                    continue

                # -----------------------------------------------------------
                # Fit model
                # -----------------------------------------------------------
                coefs = {}
                selected_genes = []
                params = {}

                if method == "cumulative_binary":
                    if config["type"] == "fixed":
                        class_probs, coefs = fit_cumulative_cv(
                            train_df, y_train, test_df, config["cols"]
                        )
                        params = {"selected_n": len(config["cols"])}
                    else:
                        class_probs, coefs, selected_genes, params = \
                            fit_cumulative_expression_inner_cv(
                                train_df, y_train, test_df,
                                config["expr_cols"],
                                extra_cols=config.get("extra_cols"),
                            )

                elif method == "coral_nn":
                    if config["type"] == "fixed":
                        class_probs = fit_coral_nn(
                            train_df, y_train, test_df, config["cols"]
                        )
                        params = {"selected_n": len(config["cols"])}
                    else:
                        class_probs, selected_genes, params = \
                            fit_coral_expression_inner_cv(
                                train_df, y_train, test_df,
                                config["expr_cols"],
                                extra_cols=config.get("extra_cols"),
                            )

                y_pred = class_probs.argmax(axis=1)

                # Compute metrics
                metrics = compute_ordinal_metrics(
                    y_test.values, y_pred, class_probs
                )

                print(
                    f"    Fold {fold}: n_train={len(train_df)}, "
                    f"n_test={len(test_df)}, "
                    f"QWK={metrics['qwk']:.3f}, "
                    f"AdjAcc={metrics['adjacent_accuracy']:.3f}, "
                    f"MacroAUROC={metrics['macro_auroc']:.3f}, "
                    f"FibGe3_AUROC={metrics['fib_ge3_auroc']:.3f}"
                )

                # Record results
                result_row = {
                    "method": method,
                    "config": config_name,
                    "fold": fold,
                    "n_train": len(train_df),
                    "n_test": len(test_df),
                    **metrics,
                    **{f"param_{k}": v for k, v in params.items()},
                }
                all_results.append(result_row)

                # Record per-sample predictions
                for i, (sid, true_stage) in enumerate(
                    zip(test_df["sample_id"].values, y_test.values)
                ):
                    pred_row = {
                        "sample_id": sid,
                        "method": method,
                        "config": config_name,
                        "fold": fold,
                        "true_stage": int(true_stage),
                        "predicted_stage": int(y_pred[i]),
                    }
                    for c in range(N_CLASSES):
                        pred_row[f"prob_F{c}"] = float(class_probs[i, c])
                    all_predictions.append(pred_row)

                # Record feature importance (cumulative binary only)
                for feat, coef_val in coefs.items():
                    config_importance[feat].append(coef_val)

            # Summarize feature importance
            for feat, coef_list in config_importance.items():
                all_importance.append({
                    "feature": feat,
                    "method": method,
                    "config": config_name,
                    "n_folds_selected": len(coef_list),
                    "mean_coef": np.mean(coef_list),
                    "mean_abs_coef": np.mean(np.abs(coef_list)),
                    "std_coef": np.std(coef_list) if len(coef_list) > 1 else 0,
                })

    # -----------------------------------------------------------------------
    # Save outputs
    # -----------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("[181] Saving outputs...")
    print(f"{'='*70}")

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT_DIR / "ordinal_results.csv", index=False)
    print(f"  ordinal_results.csv: {len(results_df)} rows")

    pred_df = pd.DataFrame(all_predictions)
    pred_df.to_csv(OUT_DIR / "ordinal_predictions.csv", index=False)
    print(f"  ordinal_predictions.csv: {len(pred_df)} rows")

    if all_importance:
        imp_df = pd.DataFrame(all_importance)
        imp_df.to_csv(OUT_DIR / "ordinal_feature_importance.csv", index=False)
        print(f"  ordinal_feature_importance.csv: {len(imp_df)} rows")

    # -----------------------------------------------------------------------
    # Summary table
    # -----------------------------------------------------------------------
    summary_rows = []
    for method in methods:
        for config_name in config_specs:
            sub = results_df[
                (results_df["method"] == method)
                & (results_df["config"] == config_name)
            ]
            if len(sub) == 0:
                continue

            row = {
                "method": method,
                "config": config_name,
                "n_folds": len(sub),
            }
            for metric in ["qwk", "adjacent_accuracy", "accuracy", "mae",
                           "macro_auroc", "fib_ge3_auroc"]:
                if metric in sub.columns:
                    row[f"mean_{metric}"] = sub[metric].mean()
                    row[f"std_{metric}"] = sub[metric].std()

            # Per-class accuracy means
            for c in range(N_CLASSES):
                col = f"acc_F{c}"
                if col in sub.columns:
                    row[f"mean_{col}"] = sub[col].mean()

            summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / "ordinal_summary.csv", index=False)
    print(f"  ordinal_summary.csv: {len(summary_df)} rows")

    # -----------------------------------------------------------------------
    # Comparison with v1 benchmarks
    # -----------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("[181] SUMMARY — Comparison with v1 benchmarks")
    print(f"{'='*70}")
    print(f"  v1 concept bottleneck QWK: 0.624")
    print(f"  v1 ordinal elastic net QWK: 0.550")
    print(f"  v2 binary fib_ge3 AUROC (Script 161 M2): ~0.799")
    print(f"  v2 binary fib_ge3 AUROC (Script 161 M3): ~0.853")
    print()

    for _, row in summary_df.iterrows():
        qwk_str = f"QWK={row['mean_qwk']:.3f}+/-{row['std_qwk']:.3f}"
        adj_str = f"AdjAcc={row['mean_adjacent_accuracy']:.3f}"
        auroc_str = f"MacroAUROC={row['mean_macro_auroc']:.3f}"
        ge3_str = f"FibGe3_AUROC={row['mean_fib_ge3_auroc']:.3f}"

        print(
            f"  {row['method']:20s} | {row['config']:15s} | "
            f"{qwk_str} | {adj_str} | {auroc_str} | {ge3_str}"
        )

    # Confusion matrix from pooled predictions (best method)
    print(f"\n  --- Confusion matrix (pooled, best QWK method) ---")
    if len(summary_df) > 0:
        best_row = summary_df.loc[summary_df["mean_qwk"].idxmax()]
        best_method = best_row["method"]
        best_config = best_row["config"]
        print(f"  Best: {best_method} / {best_config} "
              f"(QWK={best_row['mean_qwk']:.3f})")

        best_preds = pred_df[
            (pred_df["method"] == best_method)
            & (pred_df["config"] == best_config)
        ]
        if len(best_preds) > 0:
            from sklearn.metrics import confusion_matrix
            cm = confusion_matrix(
                best_preds["true_stage"],
                best_preds["predicted_stage"],
                labels=list(range(N_CLASSES)),
            )
            print(f"  {'':>5s}  Pred_F0  Pred_F1  Pred_F2  Pred_F3  Pred_F4")
            for i in range(N_CLASSES):
                row_str = f"  F{i:>1d}  "
                for j in range(N_CLASSES):
                    row_str += f"  {cm[i, j]:>5d}  "
                print(row_str)

    elapsed = time.time() - t0
    print(f"\n[181] Done in {elapsed/60:.1f} min")


if __name__ == "__main__":
    main()
