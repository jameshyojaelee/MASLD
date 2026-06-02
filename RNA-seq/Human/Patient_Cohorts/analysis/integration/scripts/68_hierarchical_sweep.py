#!/usr/bin/env python3
"""
68_hierarchical_sweep.py
Multi-task neural network + SVM baselines + cascade error propagation analysis
for the hierarchical multi-resolution staging classifier (Phase 3, Plan 2).

Architecture:
  - Shared backbone: Linear(3000->500)->ReLU->Dropout->Linear(500->100)->ReLU->Dropout
  - 4 task heads: Disease(2), Severity(4), NAS ordinal(4), Fibrosis ordinal(5)
  - Masked cross-entropy loss (handles missing labels per task)
  - LOCO cross-validation per task

Also runs:
  - SVM-RBF baseline per tier
  - Cascade error propagation analysis (Tier 1 -> Tier 2 -> Tier 3)

Inputs:
  - results/staging_classifier/prepared_data.h5

Outputs:
  - hierarchical_neural_results.csv: Per-task neural net metrics
  - cascade_error_analysis.csv: Cascaded vs oracle accuracy per tier
  - svm_tier_results.csv: SVM baseline per tier
  - hierarchical_neural_training_curves.csv: Loss per epoch

Usage: python 68_hierarchical_sweep.py
SLURM: gpu, 1 GPU, 8 CPUs, 64GB RAM, 48h
"""

import os
import sys
import numpy as np
import pandas as pd
import h5py
from collections import defaultdict

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, accuracy_score, f1_score,
    confusion_matrix, classification_report
)

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# --- Paths ---
BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
OUTDIR = os.path.join(RDIR, "staging_classifier")
os.makedirs(OUTDIR, exist_ok=True)

# --- Device ---
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"=== 68: Hierarchical Sweep (Multi-Task Neural Net) ===")
print(f"Device: {device}")
print(f"Output: {OUTDIR}\n")

# ============================================================
# Load Data
# ============================================================
print("Loading prepared data from HDF5...")
h5_path = os.path.join(OUTDIR, "prepared_data.h5")

with h5py.File(h5_path, "r") as h5:
    # Expression matrix (genes x samples) -> transpose to (samples x genes)
    X_rank = h5["rank_expression"][:].T.astype(np.float32)
    gene_names = [g.decode() for g in h5["gene_names"][:]]
    sample_ids = [s.decode() for s in h5["sample_ids"][:]]

    # Metadata
    meta = {}
    for key in h5["metadata"].keys():
        vals = h5["metadata"][key][:]
        if vals.dtype.kind == "S":  # byte string
            meta[key] = [v.decode() for v in vals]
        else:
            meta[key] = vals.astype(np.float32)

print(f"Expression: {X_rank.shape[0]} samples x {X_rank.shape[1]} genes")

# Build metadata DataFrame
meta_df = pd.DataFrame({k: v for k, v in meta.items()})
meta_df["sample_id"] = sample_ids

# Extract labels
y_disease = meta["is_disease"].astype(np.int64)           # 0/1
y_severity = meta["severity4"].astype(np.int64)            # 0-3, -1=missing
y_nas = meta["nas_group4"].astype(np.int64)                # 0-3, -1=missing
y_fib = meta["fib_stage"].astype(np.int64)                 # 0-4, -1=missing

datasets = meta["dataset"] if isinstance(meta["dataset"], list) else \
    [d.decode() if isinstance(d, bytes) else d for d in meta["dataset"]]
datasets = np.array(datasets)

MISSING = -1
N_SAMPLES, N_GENES = X_rank.shape

print(f"\nLabel availability:")
print(f"  Disease: {np.sum(y_disease >= 0)} (0={np.sum(y_disease == 0)}, 1={np.sum(y_disease == 1)})")
print(f"  Severity: {np.sum(y_severity >= 0)} (distribution: {dict(zip(*np.unique(y_severity[y_severity >= 0], return_counts=True)))})")
print(f"  NAS: {np.sum(y_nas >= 0)}")
print(f"  Fibrosis: {np.sum(y_fib >= 0)}")

# --- Datasets for LOCO ---
# Disease LOCO: cohorts with both disease and control
DISEASE_LOCO = ["GSE130970", "GSE135251", "GSE162694", "GSE174478",
                "GSE193066", "GSE240729", "GSE213621", "GSE126848"]
# NAS LOCO
NAS_LOCO = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066"]
# Fibrosis LOCO
FIB_LOCO = ["GSE130970", "GSE135251", "GSE162694", "GSE174478",
            "GSE193066", "GSE240729"]

# ============================================================
# Custom Dataset
# ============================================================
class MultiTaskDataset(Dataset):
    """Dataset returning features and a dict of task labels.
    Missing labels are stored as -1 and masked in loss computation."""

    def __init__(self, X, y_disease, y_severity, y_nas, y_fib):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y_disease = torch.tensor(y_disease, dtype=torch.long)
        self.y_severity = torch.tensor(y_severity, dtype=torch.long)
        self.y_nas = torch.tensor(y_nas, dtype=torch.long)
        self.y_fib = torch.tensor(y_fib, dtype=torch.long)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return (
            self.X[idx],
            {
                "disease": self.y_disease[idx],
                "severity": self.y_severity[idx],
                "nas": self.y_nas[idx],
                "fib": self.y_fib[idx],
            }
        )


# ============================================================
# Multi-Task Model
# ============================================================
class MultiTaskNet(nn.Module):
    """Shared backbone with task-specific heads."""

    def __init__(self, n_features, hidden1=500, hidden2=100, dropout=0.3):
        super().__init__()
        self.backbone = nn.Sequential(
            nn.Linear(n_features, hidden1),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden1, hidden2),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        # Task heads
        self.head_disease = nn.Linear(hidden2, 2)
        self.head_severity = nn.Linear(hidden2, 4)
        self.head_nas = nn.Linear(hidden2, 4)
        self.head_fib = nn.Linear(hidden2, 5)

    def forward(self, x):
        z = self.backbone(x)
        return {
            "disease": self.head_disease(z),
            "severity": self.head_severity(z),
            "nas": self.head_nas(z),
            "fib": self.head_fib(z),
        }


def compute_class_weights(y, n_classes, device):
    """Inverse-frequency class weights, ignoring missing (-1) labels."""
    valid = y[y >= 0]
    if len(valid) == 0:
        return torch.ones(n_classes, device=device)
    counts = np.bincount(valid, minlength=n_classes).astype(np.float32)
    counts = np.maximum(counts, 1.0)
    weights = len(valid) / (n_classes * counts)
    return torch.tensor(weights, dtype=torch.float32, device=device)


def masked_cross_entropy(logits, targets, weight=None):
    """Cross-entropy loss ignoring targets == -1."""
    mask = targets >= 0
    if mask.sum() == 0:
        return torch.tensor(0.0, device=logits.device, requires_grad=True)
    return nn.functional.cross_entropy(logits[mask], targets[mask], weight=weight)


# ============================================================
# Training Loop
# ============================================================
def train_multitask(X_train, labels_train, X_val, labels_val,
                    n_epochs=200, lr=1e-3, batch_size=64, patience=20,
                    task_weights=None):
    """Train multi-task model and return predictions + training curves."""

    if task_weights is None:
        task_weights = {"disease": 1.0, "severity": 1.0, "nas": 1.0, "fib": 1.0}

    n_classes = {"disease": 2, "severity": 4, "nas": 4, "fib": 5}

    # Build datasets
    train_ds = MultiTaskDataset(X_train, labels_train["disease"],
                                 labels_train["severity"],
                                 labels_train["nas"], labels_train["fib"])
    val_ds = MultiTaskDataset(X_val, labels_val["disease"],
                               labels_val["severity"],
                               labels_val["nas"], labels_val["fib"])

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True,
                              drop_last=False, num_workers=0)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False,
                            num_workers=0)

    # Compute class weights per task (on training data)
    cw = {}
    for task in n_classes:
        cw[task] = compute_class_weights(labels_train[task], n_classes[task], device)

    # Model
    model = MultiTaskNet(X_train.shape[1]).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=10,
                                                       factor=0.5, min_lr=1e-6)

    best_val_loss = float("inf")
    best_state = None
    patience_counter = 0
    history = []

    for epoch in range(n_epochs):
        # --- Train ---
        model.train()
        train_losses = defaultdict(float)
        n_batches = 0

        for X_batch, y_batch in train_loader:
            X_batch = X_batch.to(device)
            y_batch = {k: v.to(device) for k, v in y_batch.items()}

            logits = model(X_batch)

            loss = torch.tensor(0.0, device=device)
            for task in n_classes:
                task_loss = masked_cross_entropy(logits[task], y_batch[task], cw[task])
                loss = loss + task_weights[task] * task_loss
                train_losses[task] += task_loss.item()

            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            n_batches += 1

        # --- Validate ---
        model.eval()
        val_losses = defaultdict(float)
        n_val = 0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch = X_batch.to(device)
                y_batch = {k: v.to(device) for k, v in y_batch.items()}
                logits = model(X_batch)
                for task in n_classes:
                    val_losses[task] += masked_cross_entropy(
                        logits[task], y_batch[task], cw[task]
                    ).item()
                n_val += 1

        total_val = sum(task_weights[t] * val_losses[t] / max(n_val, 1)
                        for t in n_classes)
        scheduler.step(total_val)

        history.append({
            "epoch": epoch + 1,
            "train_loss": sum(train_losses[t] / max(n_batches, 1) for t in n_classes),
            "val_loss": total_val,
            **{f"train_{t}": train_losses[t] / max(n_batches, 1) for t in n_classes},
            **{f"val_{t}": val_losses[t] / max(n_val, 1) for t in n_classes},
        })

        # Early stopping
        if total_val < best_val_loss:
            best_val_loss = total_val
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1

        if patience_counter >= patience:
            print(f"    Early stopping at epoch {epoch + 1}")
            break

        if (epoch + 1) % 50 == 0:
            print(f"    Epoch {epoch+1}: train_loss={history[-1]['train_loss']:.4f} "
                  f"val_loss={total_val:.4f}")

    # Load best model
    if best_state is not None:
        model.load_state_dict(best_state)
    model.to(device)

    # --- Predict on validation ---
    model.eval()
    all_logits = {t: [] for t in n_classes}
    with torch.no_grad():
        for X_batch, _ in val_loader:
            X_batch = X_batch.to(device)
            logits = model(X_batch)
            for t in n_classes:
                all_logits[t].append(logits[t].cpu())

    predictions = {}
    for t in n_classes:
        logit_cat = torch.cat(all_logits[t], dim=0)
        probs = torch.softmax(logit_cat, dim=1).numpy()
        predictions[t] = probs

    return predictions, pd.DataFrame(history)


# ============================================================
# Evaluation Helpers
# ============================================================
def eval_task(y_true, prob_mat, task_name, class_labels=None):
    """Compute classification metrics for a single task."""
    valid = y_true >= 0
    if valid.sum() < 5:
        return {"task": task_name, "n_valid": int(valid.sum()),
                "accuracy": np.nan, "macro_auroc": np.nan, "macro_f1": np.nan}

    y_v = y_true[valid]
    p_v = prob_mat[valid]
    y_pred = p_v.argmax(axis=1)

    if class_labels is not None:
        y_pred = np.array([class_labels[i] for i in y_pred])

    acc = accuracy_score(y_v, y_pred)
    f1_macro = f1_score(y_v, y_pred, average="macro", zero_division=0)

    # Macro AUROC (one-vs-rest)
    try:
        if len(np.unique(y_v)) > 1:
            auroc = roc_auc_score(y_v, p_v, multi_class="ovr", average="macro")
        else:
            auroc = np.nan
    except ValueError:
        auroc = np.nan

    # QWK for ordinal tasks
    qwk = np.nan
    if task_name in ["nas", "fib"]:
        qwk = _quadratic_weighted_kappa(y_v, y_pred)

    # Adjacent accuracy for ordinal tasks
    adj_acc = np.nan
    if task_name in ["nas", "fib"]:
        adj_acc = np.mean(np.abs(y_v - y_pred) <= 1)

    return {
        "task": task_name, "n_valid": int(valid.sum()),
        "accuracy": acc, "macro_auroc": auroc, "macro_f1": f1_macro,
        "qwk": qwk, "adj_acc": adj_acc
    }


def _quadratic_weighted_kappa(y_true, y_pred):
    """Compute QWK manually."""
    classes = np.union1d(y_true, y_pred)
    n_classes = len(classes)
    if n_classes < 2:
        return 0.0
    class_to_idx = {c: i for i, c in enumerate(classes)}
    n = len(y_true)

    cm = np.zeros((n_classes, n_classes), dtype=float)
    for t, p in zip(y_true, y_pred):
        cm[class_to_idx[t], class_to_idx[p]] += 1

    w = np.zeros((n_classes, n_classes))
    for i in range(n_classes):
        for j in range(n_classes):
            w[i, j] = (i - j) ** 2 / (n_classes - 1) ** 2

    hist_t = cm.sum(axis=1)
    hist_p = cm.sum(axis=0)
    e = np.outer(hist_t, hist_p) / n

    num = np.sum(w * cm)
    den = np.sum(w * e)
    return 1 - num / den if den > 0 else 0.0


# ============================================================
# 1. Multi-Task Neural Network with LOCO
# ============================================================
print("\n=== 1. Multi-Task Neural Network (LOCO) ===\n")

# Use disease LOCO datasets (largest set) as the primary fold set
# For each fold, evaluate all tasks on the held-out cohort
LOCO_ALL = sorted(set(DISEASE_LOCO) | set(NAS_LOCO) | set(FIB_LOCO))

neural_results = []
all_training_curves = []

for fold_ds in LOCO_ALL:
    print(f"--- Fold: {fold_ds} ---")

    test_mask = datasets == fold_ds
    train_mask = ~test_mask

    if test_mask.sum() == 0:
        print(f"  SKIP: no samples for {fold_ds}")
        continue

    X_train = X_rank[train_mask]
    X_test = X_rank[test_mask]

    labels_train = {
        "disease": y_disease[train_mask],
        "severity": y_severity[train_mask],
        "nas": y_nas[train_mask],
        "fib": y_fib[train_mask],
    }
    labels_test = {
        "disease": y_disease[test_mask],
        "severity": y_severity[test_mask],
        "nas": y_nas[test_mask],
        "fib": y_fib[test_mask],
    }

    # Standardize features (fit on train)
    scaler = StandardScaler()
    X_train_sc = scaler.fit_transform(X_train)
    X_test_sc = scaler.transform(X_test)

    # Train
    predictions, history = train_multitask(
        X_train_sc, labels_train, X_test_sc, labels_test,
        n_epochs=200, lr=1e-3, batch_size=64, patience=20
    )
    history["fold"] = fold_ds
    all_training_curves.append(history)

    # Evaluate per task
    for task, n_cls in [("disease", 2), ("severity", 4), ("nas", 4), ("fib", 5)]:
        m = eval_task(labels_test[task], predictions[task], task)
        m["fold"] = fold_ds
        m["model"] = "multitask_neural"
        neural_results.append(m)
        if not np.isnan(m["accuracy"]):
            print(f"  {task}: acc={m['accuracy']:.3f} auroc={m['macro_auroc']:.3f} "
                  f"f1={m['macro_f1']:.3f}")

    print()

neural_df = pd.DataFrame(neural_results)
neural_df.to_csv(os.path.join(OUTDIR, "hierarchical_neural_results.csv"), index=False)
print(f"Neural results saved: {len(neural_df)} rows")

# Training curves
if all_training_curves:
    curves_df = pd.concat(all_training_curves, ignore_index=True)
    curves_df.to_csv(os.path.join(OUTDIR, "hierarchical_neural_training_curves.csv"),
                     index=False)
    print(f"Training curves saved: {len(curves_df)} rows")

# Summary
print("\nNeural net mean metrics by task:")
nn_summary = neural_df.groupby("task").agg(
    mean_acc=("accuracy", "mean"),
    mean_auroc=("macro_auroc", "mean"),
    mean_f1=("macro_f1", "mean"),
    mean_qwk=("qwk", "mean"),
    n_folds=("accuracy", "count"),
).round(4)
print(nn_summary.to_string())

# ============================================================
# 2. SVM-RBF Baselines
# ============================================================
print("\n\n=== 2. SVM-RBF Baselines ===\n")

svm_results = []

TASK_CONFIG = {
    "disease": {"y": y_disease, "n_classes": 2, "loco_ds": DISEASE_LOCO},
    "severity": {"y": y_severity, "n_classes": 4, "loco_ds": DISEASE_LOCO},
    "nas": {"y": y_nas, "n_classes": 4, "loco_ds": NAS_LOCO},
    "fib": {"y": y_fib, "n_classes": 5, "loco_ds": FIB_LOCO},
}

for task_name, cfg in TASK_CONFIG.items():
    print(f"--- SVM: {task_name} ---")
    y_task = cfg["y"]
    valid_mask = y_task >= 0

    for fold_ds in cfg["loco_ds"]:
        test_mask = (datasets == fold_ds) & valid_mask
        train_mask = (datasets != fold_ds) & valid_mask

        if test_mask.sum() == 0 or train_mask.sum() < 10:
            continue

        X_tr = X_rank[train_mask]
        y_tr = y_task[train_mask]
        X_te = X_rank[test_mask]
        y_te = y_task[test_mask]

        # Need at least 2 classes
        if len(np.unique(y_tr)) < 2:
            continue

        # Standardize
        scaler = StandardScaler()
        X_tr_sc = scaler.fit_transform(X_tr)
        X_te_sc = scaler.transform(X_te)

        # SVM-RBF
        try:
            svm = SVC(kernel="rbf", probability=True, class_weight="balanced",
                       C=1.0, gamma="scale", random_state=42, max_iter=5000)
            svm.fit(X_tr_sc, y_tr)
            y_pred = svm.predict(X_te_sc)
            y_prob = svm.predict_proba(X_te_sc)

            acc = accuracy_score(y_te, y_pred)
            f1_macro = f1_score(y_te, y_pred, average="macro", zero_division=0)

            try:
                auroc = roc_auc_score(y_te, y_prob, multi_class="ovr", average="macro")
            except ValueError:
                auroc = np.nan

            qwk = _quadratic_weighted_kappa(y_te, y_pred) if task_name in ["nas", "fib"] else np.nan
            adj_acc = np.mean(np.abs(y_te - y_pred) <= 1) if task_name in ["nas", "fib"] else np.nan

            svm_results.append({
                "task": task_name, "fold": fold_ds, "model": "svm_rbf",
                "n_test": int(test_mask.sum()),
                "accuracy": acc, "macro_auroc": auroc, "macro_f1": f1_macro,
                "qwk": qwk, "adj_acc": adj_acc
            })
            print(f"  {fold_ds}: acc={acc:.3f} auroc={auroc:.3f}")

        except Exception as e:
            print(f"  {fold_ds}: ERROR: {e}")
            continue

    print()

svm_df = pd.DataFrame(svm_results)
svm_df.to_csv(os.path.join(OUTDIR, "svm_tier_results.csv"), index=False)
print(f"SVM results saved: {len(svm_df)} rows")

print("\nSVM mean metrics by task:")
if len(svm_df) > 0:
    svm_summary = svm_df.groupby("task").agg(
        mean_acc=("accuracy", "mean"),
        mean_auroc=("macro_auroc", "mean"),
        mean_f1=("macro_f1", "mean"),
        mean_qwk=("qwk", "mean"),
        n_folds=("accuracy", "count"),
    ).round(4)
    print(svm_summary.to_string())

# ============================================================
# 3. Cascade Error Propagation Analysis
# ============================================================
print("\n\n=== 3. Cascade Error Propagation Analysis ===\n")
print("Quantifying how Tier 1 errors degrade Tier 2 and Tier 3 performance.\n")

# We need per-sample predictions from all tiers.
# Retrain a single multi-task model on a 80/20 random split for cascade analysis,
# OR use the LOCO predictions collected above.
#
# Strategy: collect all LOCO predictions into per-sample arrays, then simulate
# the cascade on samples that have predictions from all tiers.

# Reconstruct per-sample neural predictions from LOCO folds
# We need to retrain per fold and collect, but we already did this above.
# Instead, do a simpler approach: single 80/20 stratified split for cascade.

print("Running cascade analysis on 80/20 stratified split...\n")
np.random.seed(42)

# Use samples that have disease label + at least one staging label
has_staging = (y_severity >= 0) | (y_nas >= 0) | (y_fib >= 0)
cascade_mask = (y_disease >= 0) & has_staging
cascade_idx = np.where(cascade_mask)[0]
np.random.shuffle(cascade_idx)

split_point = int(0.8 * len(cascade_idx))
train_idx = cascade_idx[:split_point]
test_idx = cascade_idx[split_point:]

print(f"Cascade split: train={len(train_idx)}, test={len(test_idx)}")

# Standardize
scaler_cascade = StandardScaler()
X_train_c = scaler_cascade.fit_transform(X_rank[train_idx])
X_test_c = scaler_cascade.transform(X_rank[test_idx])

labels_train_c = {
    "disease": y_disease[train_idx],
    "severity": y_severity[train_idx],
    "nas": y_nas[train_idx],
    "fib": y_fib[train_idx],
}
labels_test_c = {
    "disease": y_disease[test_idx],
    "severity": y_severity[test_idx],
    "nas": y_nas[test_idx],
    "fib": y_fib[test_idx],
}

# Train multi-task model
cascade_preds, _ = train_multitask(
    X_train_c, labels_train_c, X_test_c, labels_test_c,
    n_epochs=200, lr=1e-3, batch_size=64, patience=20
)

# Tier 1 predictions (Disease vs Control)
tier1_pred = cascade_preds["disease"].argmax(axis=1)  # 0=control, 1=disease
tier1_true = labels_test_c["disease"]

# Oracle: use ground-truth Tier 1
tier1_oracle = tier1_true.copy()

cascade_rows = []

# --- Tier 1 metrics ---
t1_valid = tier1_true >= 0
t1_acc = accuracy_score(tier1_true[t1_valid], tier1_pred[t1_valid])
print(f"Tier 1 (Disease Detection) accuracy: {t1_acc:.3f}")
cascade_rows.append({
    "tier": "Tier1_Disease", "scenario": "standalone",
    "accuracy": t1_acc, "n_samples": int(t1_valid.sum())
})

# --- Tier 2: Severity (cascaded vs oracle) ---
sev_valid = labels_test_c["severity"] >= 0
if sev_valid.sum() > 0:
    sev_true = labels_test_c["severity"]
    sev_pred = cascade_preds["severity"].argmax(axis=1)

    # Oracle: assume perfect Tier 1 -> evaluate severity on all disease samples
    oracle_disease = sev_true >= 0  # has severity label = was deemed "disease"
    oracle_acc = accuracy_score(sev_true[oracle_disease], sev_pred[oracle_disease])

    # Cascaded: only evaluate on samples where Tier 1 predicted "disease"
    # Samples misclassified as "control" by Tier 1 are lost
    cascaded_mask = (tier1_pred == 1) & sev_valid
    if cascaded_mask.sum() > 0:
        cascaded_acc = accuracy_score(sev_true[cascaded_mask], sev_pred[cascaded_mask])
        # Also count how many disease samples were lost (Tier 1 false negatives)
        missed_disease = sev_valid & (tier1_pred == 0)
        n_missed = missed_disease.sum()
    else:
        cascaded_acc = np.nan
        n_missed = 0

    print(f"\nTier 2 (Severity):")
    print(f"  Oracle accuracy: {oracle_acc:.3f} (n={oracle_disease.sum()})")
    print(f"  Cascaded accuracy: {cascaded_acc:.3f} (n={cascaded_mask.sum()}, "
          f"missed={n_missed})")

    cascade_rows.append({
        "tier": "Tier2_Severity", "scenario": "oracle",
        "accuracy": oracle_acc, "n_samples": int(oracle_disease.sum()),
        "n_missed_by_tier1": 0
    })
    cascade_rows.append({
        "tier": "Tier2_Severity", "scenario": "cascaded",
        "accuracy": cascaded_acc, "n_samples": int(cascaded_mask.sum()),
        "n_missed_by_tier1": int(n_missed)
    })

# --- Tier 3a: NAS (cascaded vs oracle) ---
nas_valid = labels_test_c["nas"] >= 0
if nas_valid.sum() > 0:
    nas_true = labels_test_c["nas"]
    nas_pred = cascade_preds["nas"].argmax(axis=1)

    oracle_mask = nas_valid
    oracle_acc = accuracy_score(nas_true[oracle_mask], nas_pred[oracle_mask])

    # Cascaded: must pass Tier 1 (predicted disease)
    cascaded_mask = (tier1_pred == 1) & nas_valid
    if cascaded_mask.sum() > 0:
        cascaded_acc = accuracy_score(nas_true[cascaded_mask], nas_pred[cascaded_mask])
        n_missed = (nas_valid & (tier1_pred == 0)).sum()
    else:
        cascaded_acc = np.nan
        n_missed = 0

    print(f"\nTier 3a (NAS Ordinal):")
    print(f"  Oracle accuracy: {oracle_acc:.3f} (n={oracle_mask.sum()})")
    print(f"  Cascaded accuracy: {cascaded_acc:.3f} (n={cascaded_mask.sum()}, "
          f"missed={n_missed})")

    cascade_rows.append({
        "tier": "Tier3a_NAS", "scenario": "oracle",
        "accuracy": oracle_acc, "n_samples": int(oracle_mask.sum()),
        "n_missed_by_tier1": 0
    })
    cascade_rows.append({
        "tier": "Tier3a_NAS", "scenario": "cascaded",
        "accuracy": cascaded_acc, "n_samples": int(cascaded_mask.sum()),
        "n_missed_by_tier1": int(n_missed)
    })

# --- Tier 3b: Fibrosis (cascaded vs oracle) ---
fib_valid = labels_test_c["fib"] >= 0
if fib_valid.sum() > 0:
    fib_true = labels_test_c["fib"]
    fib_pred = cascade_preds["fib"].argmax(axis=1)

    oracle_mask = fib_valid
    oracle_acc = accuracy_score(fib_true[oracle_mask], fib_pred[oracle_mask])

    cascaded_mask = (tier1_pred == 1) & fib_valid
    if cascaded_mask.sum() > 0:
        cascaded_acc = accuracy_score(fib_true[cascaded_mask], fib_pred[cascaded_mask])
        n_missed = (fib_valid & (tier1_pred == 0)).sum()
    else:
        cascaded_acc = np.nan
        n_missed = 0

    print(f"\nTier 3b (Fibrosis Ordinal):")
    print(f"  Oracle accuracy: {oracle_acc:.3f} (n={oracle_mask.sum()})")
    print(f"  Cascaded accuracy: {cascaded_acc:.3f} (n={cascaded_mask.sum()}, "
          f"missed={n_missed})")

    cascade_rows.append({
        "tier": "Tier3b_Fibrosis", "scenario": "oracle",
        "accuracy": oracle_acc, "n_samples": int(oracle_mask.sum()),
        "n_missed_by_tier1": 0
    })
    cascade_rows.append({
        "tier": "Tier3b_Fibrosis", "scenario": "cascaded",
        "accuracy": cascaded_acc, "n_samples": int(cascaded_mask.sum()),
        "n_missed_by_tier1": int(n_missed)
    })

# --- Two-stage cascade: Tier 1 -> Tier 2 -> Tier 3 ---
# Tier 2 predicts severity; only MASH samples (severity >= 2) proceed to Tier 3
if sev_valid.sum() > 0 and (nas_valid.sum() > 0 or fib_valid.sum() > 0):
    print("\n--- Two-stage cascade: Tier1 -> Tier2 -> Tier3 ---")

    # Samples predicted as MASH by Tier 2
    tier2_pred = cascade_preds["severity"].argmax(axis=1)
    mash_pred = tier2_pred >= 2

    for t3_name, t3_valid, t3_true, t3_pred_arr in [
        ("Tier3a_NAS_2stage", nas_valid, labels_test_c["nas"],
         cascade_preds["nas"].argmax(axis=1)),
        ("Tier3b_Fib_2stage", fib_valid, labels_test_c["fib"],
         cascade_preds["fib"].argmax(axis=1)),
    ]:
        if t3_valid.sum() == 0:
            continue

        # Oracle 2-stage: perfect Tier 1 + perfect Tier 2
        oracle_2s = t3_valid
        oracle_acc_2s = accuracy_score(t3_true[oracle_2s], t3_pred_arr[oracle_2s])

        # Cascaded 2-stage: Tier 1 says disease AND Tier 2 says MASH
        cascaded_2s = (tier1_pred == 1) & mash_pred & t3_valid
        if cascaded_2s.sum() > 0:
            cascaded_acc_2s = accuracy_score(t3_true[cascaded_2s], t3_pred_arr[cascaded_2s])
        else:
            cascaded_acc_2s = np.nan

        n_lost_t1 = (t3_valid & (tier1_pred == 0)).sum()
        n_lost_t2 = (t3_valid & (tier1_pred == 1) & (~mash_pred)).sum()

        print(f"  {t3_name}: oracle={oracle_acc_2s:.3f} cascaded={cascaded_acc_2s:.3f} "
              f"(lost T1={n_lost_t1}, lost T2={n_lost_t2})")

        cascade_rows.append({
            "tier": t3_name, "scenario": "oracle_2stage",
            "accuracy": oracle_acc_2s, "n_samples": int(oracle_2s.sum()),
            "n_missed_by_tier1": 0, "n_missed_by_tier2": 0
        })
        cascade_rows.append({
            "tier": t3_name, "scenario": "cascaded_2stage",
            "accuracy": cascaded_acc_2s, "n_samples": int(cascaded_2s.sum()),
            "n_missed_by_tier1": int(n_lost_t1),
            "n_missed_by_tier2": int(n_lost_t2)
        })

cascade_df = pd.DataFrame(cascade_rows)
cascade_df.to_csv(os.path.join(OUTDIR, "cascade_error_analysis.csv"), index=False)
print(f"\nCascade error analysis saved: {len(cascade_df)} rows")

# ============================================================
# Final Summary
# ============================================================
print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)

print("\n--- Multi-Task Neural Net (mean LOCO) ---")
for task in ["disease", "severity", "nas", "fib"]:
    task_data = neural_df[neural_df["task"] == task]
    if len(task_data) > 0:
        print(f"  {task}: acc={task_data['accuracy'].mean():.3f} "
              f"auroc={task_data['macro_auroc'].mean():.3f} "
              f"f1={task_data['macro_f1'].mean():.3f}")

print("\n--- SVM-RBF (mean LOCO) ---")
for task in ["disease", "severity", "nas", "fib"]:
    task_data = svm_df[svm_df["task"] == task]
    if len(task_data) > 0:
        print(f"  {task}: acc={task_data['accuracy'].mean():.3f} "
              f"auroc={task_data['macro_auroc'].mean():.3f} "
              f"f1={task_data['macro_f1'].mean():.3f}")

print("\n--- Cascade Error Propagation ---")
print(cascade_df.to_string(index=False))

print(f"\n=== 68_hierarchical_sweep.py completed ===")
