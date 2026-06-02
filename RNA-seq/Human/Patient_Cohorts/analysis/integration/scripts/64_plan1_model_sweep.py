#!/usr/bin/env python3
"""
64_plan1_model_sweep.py
Python model sweep on rank-transformed features with LOCO cross-validation.

Plan 1 (Ordinal Evidence Convergence) -- Phase 2, Script 2

Models:
  1. Random Forest (n_estimators=500, max_features=sqrt)
  2. XGBoost (n_estimators=300, max_depth=6, lr=0.1)
  3. LightGBM (n_estimators=300, num_leaves=31, lr=0.1)
  4. SVM-RBF (C=1.0, gamma=scale, probability=True)
  5. Simple MLP (hidden=[500, 100], dropout=0.3, epochs=100, lr=1e-3)

Targets:
  - Fibrosis 5-class (F0-F4): 6-fold LOCO
  - NAS 4-group: 5-fold LOCO
  - Binary F>=3: 6-fold LOCO
  - Binary NAS>=5: 5-fold LOCO

Outputs:
  - plan1_sweep_fibrosis.csv
  - plan1_sweep_nas.csv
  - plan1_sweep_binary.csv
  - plan1_confusion_matrices.json
  - plan1_feature_importances.csv

Usage: python 64_plan1_model_sweep.py
SLURM: cpu, 8 CPUs, 64GB RAM, 48h
"""

import os
import sys
import json
import time
import warnings
import numpy as np
import pandas as pd
import h5py
from collections import defaultdict

# Sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, mean_absolute_error,
    roc_auc_score, f1_score, confusion_matrix, cohen_kappa_score
)
from sklearn.utils.class_weight import compute_class_weight, compute_sample_weight

# XGBoost / LightGBM
import xgboost as xgb
import lightgbm as lgb

# PyTorch for MLP
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

warnings.filterwarnings("ignore", category=UserWarning)

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)
torch.manual_seed(RANDOM_STATE)

# =============================================================================
# PATHS
# =============================================================================
BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
OUTDIR = os.path.join(RDIR, "staging_classifier")
LOCO_GENE_DIR = os.path.join(OUTDIR, "loco_gene_lists")
os.makedirs(OUTDIR, exist_ok=True)

print("=== 64: Plan 1 Model Sweep ===")
print(f"Output dir: {OUTDIR}")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")

# =============================================================================
# LOAD DATA
# =============================================================================
h5_path = os.path.join(OUTDIR, "prepared_data.h5")
if not os.path.exists(h5_path):
    sys.exit(f"ERROR: {h5_path} not found. Run Script 62 first.")

print(f"\nLoading HDF5: {h5_path}")
with h5py.File(h5_path, "r") as h5:
    # Expression matrix (genes x samples) -> transpose to (samples x genes)
    rank_expr = h5["rank_expression"][:].T  # Now samples x genes
    gene_names = np.array([g.decode() if isinstance(g, bytes) else g
                           for g in h5["gene_names"][:]])
    sample_ids = np.array([s.decode() if isinstance(s, bytes) else s
                           for s in h5["sample_ids"][:]])

    # Metadata
    meta = {}
    for key in h5["metadata"].keys():
        vals = h5["metadata"][key][:]
        if vals.dtype.kind == 'S' or vals.dtype.kind == 'O':
            meta[key] = np.array([v.decode() if isinstance(v, bytes) else v
                                  for v in vals])
        else:
            meta[key] = vals.astype(np.float32)

print(f"  Expression: {rank_expr.shape[0]} samples x {rank_expr.shape[1]} genes")
print(f"  Gene names: {len(gene_names)}")
print(f"  Sample IDs: {len(sample_ids)}")

# Build gene_name -> column index mapping for fast subsetting
gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}

# Load global feature candidates as fallback
_global_fc_path = os.path.join(OUTDIR, "feature_candidates_3000.csv")
if os.path.exists(_global_fc_path):
    _global_fc = pd.read_csv(_global_fc_path)
    feature_genes_global = _global_fc["gene"].tolist()
else:
    feature_genes_global = list(gene_names)


def load_fold_gene_indices(fold_id):
    """Return column indices into rank_expr for the per-fold gene list.

    Falls back to the global feature_candidates_3000.csv if the fold-specific
    file does not exist.
    """
    fold_gene_file = os.path.join(
        LOCO_GENE_DIR, f"fold_{fold_id}_feature_candidates_3000.csv")
    if os.path.exists(fold_gene_file):
        fold_genes = pd.read_csv(fold_gene_file)["gene"].tolist()
        print(f"    Using per-fold gene list for fold: {fold_id}")
    else:
        fold_genes = feature_genes_global
        print(f"    Fold gene list not found for {fold_id} -- using global fallback")
    # Return indices of genes present in the expression matrix
    return np.array([gene_name_to_idx[g] for g in fold_genes
                     if g in gene_name_to_idx], dtype=np.int64)


# Build metadata DataFrame
meta_df = pd.DataFrame(meta)
meta_df["sample_id"] = sample_ids
print(f"  Metadata columns: {list(meta_df.columns)}")

# =============================================================================
# DEFINE TARGETS AND FOLDS
# =============================================================================
# Fibrosis ordinal: fib_stage in {0,1,2,3,4}, LOCO by loco_fold_fibrosis
# NAS ordinal: nas_group4 in {0,1,2,3}, LOCO by loco_fold_nas
# Binary F>=3: fib_ge3 in {0,1}, LOCO by loco_fold_fibrosis
# Binary NAS>=5: nas_ge5 in {0,1}, LOCO by loco_fold_nas

def get_valid_mask(meta_df, target_col, fold_col):
    """Return boolean mask of samples with valid target and non-excluded fold."""
    target_valid = meta_df[target_col].values >= 0
    fold_valid = meta_df[fold_col].values != "excluded"
    # fold_col might be numeric array stored as float — handle that
    if meta_df[fold_col].dtype in [np.float32, np.float64]:
        fold_valid = np.ones(len(meta_df), dtype=bool)
    return target_valid & fold_valid


def get_folds(meta_df, fold_col, valid_mask):
    """Return sorted unique fold names, excluding 'excluded'."""
    vals = meta_df[fold_col].values[valid_mask]
    if vals.dtype in [np.float32, np.float64]:
        return sorted(set(vals.astype(int)))
    folds = sorted(set(v for v in vals if v != "excluded" and v != "NA"))
    return folds


# =============================================================================
# QUADRATIC WEIGHTED KAPPA
# =============================================================================
def quadratic_weighted_kappa(y_true, y_pred):
    """Compute QWK using sklearn's cohen_kappa_score with quadratic weights."""
    return cohen_kappa_score(y_true, y_pred, weights="quadratic")


# =============================================================================
# TORCH MLP
# =============================================================================
class SimpleMLP(nn.Module):
    def __init__(self, input_dim, hidden_dims, n_classes, dropout=0.3):
        super().__init__()
        layers = []
        prev_dim = input_dim
        for h in hidden_dims:
            layers.append(nn.Linear(prev_dim, h))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
            prev_dim = h
        layers.append(nn.Linear(prev_dim, n_classes))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


def train_mlp(X_train, y_train, X_test, n_classes, class_weights,
              hidden_dims=(500, 100), dropout=0.3, epochs=100, lr=1e-3,
              patience=10, batch_size=64):
    """Train a simple MLP and return predicted class + probabilities on test set."""
    device = torch.device("cpu")  # CPU-only for this pipeline

    # Prepare tensors
    X_tr = torch.tensor(X_train, dtype=torch.float32).to(device)
    y_tr = torch.tensor(y_train, dtype=torch.long).to(device)
    X_te = torch.tensor(X_test, dtype=torch.float32).to(device)

    # Class weight tensor for CrossEntropyLoss
    cw_tensor = torch.tensor(class_weights, dtype=torch.float32).to(device)

    model = SimpleMLP(X_train.shape[1], hidden_dims, n_classes, dropout).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    criterion = nn.CrossEntropyLoss(weight=cw_tensor)

    dataset = TensorDataset(X_tr, y_tr)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True,
                        drop_last=False)

    best_loss = float("inf")
    best_state = None
    no_improve = 0

    model.train()
    for epoch in range(epochs):
        epoch_loss = 0.0
        n_batches = 0
        for xb, yb in loader:
            optimizer.zero_grad()
            logits = model(xb)
            loss = criterion(logits, yb)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
            n_batches += 1

        avg_loss = epoch_loss / max(n_batches, 1)
        if avg_loss < best_loss - 1e-4:
            best_loss = avg_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            no_improve = 0
        else:
            no_improve += 1

        if no_improve >= patience:
            break

    # Restore best
    if best_state is not None:
        model.load_state_dict(best_state)

    # Predict
    model.eval()
    with torch.no_grad():
        logits = model(X_te)
        probs = torch.softmax(logits, dim=1).cpu().numpy()
        preds = np.argmax(probs, axis=1)

    return preds, probs


# =============================================================================
# MODEL DEFINITIONS
# =============================================================================
def get_models(n_classes, is_binary=False):
    """Return dict of model_name -> (model_or_None, needs_scaling)."""
    models = {}

    # 1. Random Forest
    models["RandomForest"] = (
        RandomForestClassifier(
            n_estimators=500,
            max_features="sqrt",
            class_weight="balanced",
            random_state=RANDOM_STATE,
            n_jobs=-1
        ),
        False  # no scaling needed
    )

    # 2. XGBoost
    # XGBoost doesn't take class_weight directly — we pass sample_weight later
    models["XGBoost"] = (
        xgb.XGBClassifier(
            n_estimators=300,
            max_depth=6,
            learning_rate=0.1,
            random_state=RANDOM_STATE,
            n_jobs=-1,
            eval_metric="mlogloss" if not is_binary else "logloss",
            tree_method="hist",
            verbosity=0
        ),
        False
    )

    # 3. LightGBM
    models["LightGBM"] = (
        lgb.LGBMClassifier(
            n_estimators=300,
            num_leaves=31,
            learning_rate=0.1,
            class_weight="balanced",
            random_state=RANDOM_STATE,
            n_jobs=-1,
            verbose=-1
        ),
        False
    )

    # 4. SVM-RBF
    models["SVM_RBF"] = (
        SVC(
            C=1.0,
            gamma="scale",
            kernel="rbf",
            probability=True,
            class_weight="balanced",
            random_state=RANDOM_STATE
        ),
        True  # SVM needs scaling
    )

    # 5. MLP — handled separately via train_mlp; placeholder
    models["MLP"] = (None, True)  # None = custom training

    return models


# =============================================================================
# RUN SINGLE FOLD
# =============================================================================
def run_fold(model_name, model_obj, needs_scaling, X_train, y_train,
             X_test, y_test, n_classes, is_binary=False):
    """Train model on fold, return predictions dict."""
    # Scaling
    if needs_scaling:
        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)
    else:
        X_train_s = X_train
        X_test_s = X_test

    # Compute sample weights for XGBoost
    if model_name == "XGBoost":
        sw = compute_sample_weight("balanced", y_train)
    else:
        sw = None

    y_pred = None
    y_prob = None
    feature_importance = None

    if model_name == "MLP":
        # Compute class weights for MLP loss
        classes_present = np.unique(y_train)
        cw = compute_class_weight("balanced", classes=classes_present, y=y_train)
        # Map to full class range
        cw_full = np.ones(n_classes, dtype=np.float32)
        for i, c in enumerate(classes_present):
            cw_full[int(c)] = cw[i]

        y_pred, y_prob = train_mlp(
            X_train_s, y_train, X_test_s,
            n_classes=n_classes,
            class_weights=cw_full,
            hidden_dims=(500, 100),
            dropout=0.3,
            epochs=100,
            lr=1e-3,
            patience=10,
            batch_size=min(64, len(y_train))
        )
    else:
        if sw is not None:
            model_obj.fit(X_train_s, y_train, sample_weight=sw)
        else:
            model_obj.fit(X_train_s, y_train)

        y_pred = model_obj.predict(X_test_s)
        if hasattr(model_obj, "predict_proba"):
            y_prob = model_obj.predict_proba(X_test_s)

        # Feature importances (RF, XGB, LGBM)
        if hasattr(model_obj, "feature_importances_"):
            feature_importance = model_obj.feature_importances_

    # Metrics
    metrics = {}
    metrics["accuracy"] = accuracy_score(y_test, y_pred)
    metrics["balanced_accuracy"] = balanced_accuracy_score(y_test, y_pred)
    metrics["mae"] = mean_absolute_error(y_test, y_pred)

    if not is_binary:
        metrics["qwk"] = quadratic_weighted_kappa(y_test, y_pred)
        metrics["adjacent_accuracy"] = np.mean(np.abs(y_test - y_pred) <= 1)
    else:
        # Binary-specific
        if y_prob is not None and y_prob.ndim == 2 and y_prob.shape[1] >= 2:
            pos_prob = y_prob[:, 1]
        elif y_prob is not None and y_prob.ndim == 1:
            pos_prob = y_prob
        else:
            pos_prob = None

        if pos_prob is not None and len(np.unique(y_test)) == 2:
            metrics["auroc"] = roc_auc_score(y_test, pos_prob)
        else:
            metrics["auroc"] = np.nan

    # Per-class F1
    labels = sorted(np.unique(np.concatenate([y_test, y_pred])))
    f1_per_class = f1_score(y_test, y_pred, labels=labels, average=None,
                            zero_division=0)
    for i, lab in enumerate(labels):
        metrics[f"f1_class_{lab}"] = f1_per_class[i]
    metrics["f1_macro"] = f1_score(y_test, y_pred, average="macro",
                                    zero_division=0)

    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred, labels=list(range(n_classes)))

    return {
        "y_pred": y_pred,
        "y_prob": y_prob,
        "metrics": metrics,
        "confusion_matrix": cm.tolist(),
        "feature_importance": feature_importance,
        "n_train": len(y_train),
        "n_test": len(y_test)
    }


# =============================================================================
# RUN MODEL SWEEP FOR A TARGET
# =============================================================================
def run_sweep(X, y, folds, fold_labels, target_name, n_classes, is_binary=False):
    """Run all models across all LOCO folds for a given target."""
    results = []
    confusion_matrices = {}
    feature_importances = defaultdict(list)

    unique_folds = sorted(set(fold_labels))
    print(f"\n{'='*60}")
    print(f"TARGET: {target_name} | Classes: {n_classes} | "
          f"Folds: {len(unique_folds)} | Samples: {len(y)}")
    print(f"Class distribution: {dict(zip(*np.unique(y, return_counts=True)))}")
    print(f"{'='*60}")

    models_def = get_models(n_classes, is_binary)

    for fold_name in unique_folds:
        test_mask = fold_labels == fold_name
        train_mask = ~test_mask

        # Load per-fold gene list and subset features (removes leakage)
        fold_gene_idx = load_fold_gene_indices(fold_name)
        X_fold = X[:, fold_gene_idx]

        X_train = X_fold[train_mask]
        y_train = y[train_mask]
        X_test = X_fold[test_mask]
        y_test = y[test_mask]

        if len(X_test) < 5 or len(X_train) < 20:
            print(f"  Skipping fold {fold_name} (train={len(X_train)}, test={len(X_test)})")
            continue

        if len(np.unique(y_train)) < 2:
            print(f"  Skipping fold {fold_name} (< 2 classes in training)")
            continue

        print(f"\n  Fold: {fold_name} | Train: {len(X_train)} | Test: {len(X_test)} | Features: {X_fold.shape[1]}")

        for model_name, (model_obj, needs_scaling) in models_def.items():
            t0 = time.time()
            # Clone the model for each fold to avoid state leakage
            if model_obj is not None:
                from sklearn.base import clone
                model_clone = clone(model_obj)
            else:
                model_clone = None

            try:
                result = run_fold(
                    model_name, model_clone, needs_scaling,
                    X_train, y_train, X_test, y_test,
                    n_classes, is_binary
                )
            except Exception as e:
                print(f"    {model_name}: FAILED - {e}")
                continue

            elapsed = time.time() - t0
            m = result["metrics"]
            if is_binary:
                print(f"    {model_name:15s}: acc={m['accuracy']:.3f} "
                      f"bal_acc={m['balanced_accuracy']:.3f} "
                      f"auroc={m.get('auroc', np.nan):.3f} "
                      f"f1={m['f1_macro']:.3f} ({elapsed:.1f}s)")
            else:
                print(f"    {model_name:15s}: acc={m['accuracy']:.3f} "
                      f"bal_acc={m['balanced_accuracy']:.3f} "
                      f"qwk={m.get('qwk', np.nan):.3f} "
                      f"mae={m['mae']:.3f} "
                      f"adj={m.get('adjacent_accuracy', np.nan):.3f} ({elapsed:.1f}s)")

            # Collect results
            row = {
                "target": target_name,
                "model": model_name,
                "fold": fold_name,
                "n_train": result["n_train"],
                "n_test": result["n_test"],
            }
            row.update(result["metrics"])
            results.append(row)

            # Confusion matrix
            cm_key = f"{target_name}__{model_name}__{fold_name}"
            confusion_matrices[cm_key] = result["confusion_matrix"]

            # Feature importances: store (gene_names_for_fold, importance_vector)
            if result["feature_importance"] is not None:
                feature_importances[model_name].append(
                    (gene_names[fold_gene_idx], result["feature_importance"]))

    # Aggregate feature importances across folds (gene-level mean, fold-aware)
    fi_records = []
    for model_name, imp_pairs in feature_importances.items():
        if len(imp_pairs) == 0:
            continue
        # Accumulate per-gene importance across folds
        gene_imp_sum = defaultdict(float)
        gene_imp_count = defaultdict(int)
        for fold_gene_arr, imp_vec in imp_pairs:
            for g, v in zip(fold_gene_arr, imp_vec):
                gene_imp_sum[g] += float(v)
                gene_imp_count[g] += 1
        for g, total in gene_imp_sum.items():
            mean_imp = total / gene_imp_count[g]
            if mean_imp > 0:
                fi_records.append({
                    "target": target_name,
                    "model": model_name,
                    "gene": g,
                    "importance_mean": mean_imp,
                    "importance_std": np.nan,  # std not tracked per-gene across folds
                    "n_folds": gene_imp_count[g]
                })

    return results, confusion_matrices, fi_records


# =============================================================================
# PREPARE TARGETS
# =============================================================================
# Fibrosis ordinal
fib_target = meta_df["fib_stage"].values.astype(int)
fib_fold = meta_df["loco_fold_fibrosis"].values
fib_valid = (fib_target >= 0) & (fib_fold != "excluded") & (fib_fold != "NA")
fib_X = rank_expr[fib_valid]
fib_y = fib_target[fib_valid]
fib_folds = fib_fold[fib_valid]
print(f"\nFibrosis ordinal: {fib_X.shape[0]} samples, {len(np.unique(fib_folds))} folds")

# NAS ordinal (4 groups)
nas_target = meta_df["nas_group4"].values.astype(int)
nas_fold = meta_df["loco_fold_nas"].values
nas_valid = (nas_target >= 0) & (nas_fold != "excluded") & (nas_fold != "NA")
nas_X = rank_expr[nas_valid]
nas_y = nas_target[nas_valid]
nas_folds = nas_fold[nas_valid]
print(f"NAS 4-group: {nas_X.shape[0]} samples, {len(np.unique(nas_folds))} folds")

# Binary F>=3
fib_bin_target = meta_df["fib_ge3"].values.astype(int)
fib_bin_valid = (fib_bin_target >= 0) & (fib_fold != "excluded") & (fib_fold != "NA")
fib_bin_X = rank_expr[fib_bin_valid]
fib_bin_y = fib_bin_target[fib_bin_valid]
fib_bin_folds = fib_fold[fib_bin_valid]
print(f"Binary F>=3: {fib_bin_X.shape[0]} samples")

# Binary NAS>=5
nas_bin_target = meta_df["nas_ge5"].values.astype(int)
nas_bin_valid = (nas_bin_target >= 0) & (nas_fold != "excluded") & (nas_fold != "NA")
nas_bin_X = rank_expr[nas_bin_valid]
nas_bin_y = nas_bin_target[nas_bin_valid]
nas_bin_folds = nas_fold[nas_bin_valid]
print(f"Binary NAS>=5: {nas_bin_X.shape[0]} samples")

# =============================================================================
# RUN SWEEPS
# =============================================================================
all_results = []
all_cm = {}
all_fi = []

# 1. Fibrosis ordinal (F0-F4 = 5 classes)
res, cm, fi = run_sweep(fib_X, fib_y, None, fib_folds,
                         "fibrosis_5class", n_classes=5, is_binary=False)
all_results.extend(res)
all_cm.update(cm)
all_fi.extend(fi)

# 2. NAS 4-group
res, cm, fi = run_sweep(nas_X, nas_y, None, nas_folds,
                         "nas_4group", n_classes=4, is_binary=False)
all_results.extend(res)
all_cm.update(cm)
all_fi.extend(fi)

# 3. Binary F>=3
res, cm, fi = run_sweep(fib_bin_X, fib_bin_y, None, fib_bin_folds,
                         "fib_ge3_binary", n_classes=2, is_binary=True)
all_results.extend(res)
all_cm.update(cm)
all_fi.extend(fi)

# 4. Binary NAS>=5
res, cm, fi = run_sweep(nas_bin_X, nas_bin_y, None, nas_bin_folds,
                         "nas_ge5_binary", n_classes=2, is_binary=True)
all_results.extend(res)
all_cm.update(cm)
all_fi.extend(fi)

# =============================================================================
# SAVE OUTPUTS
# =============================================================================
print(f"\n{'='*60}")
print("SAVING OUTPUTS")
print(f"{'='*60}")

results_df = pd.DataFrame(all_results)

# Split by target type
fib_results = results_df[results_df["target"] == "fibrosis_5class"]
nas_results = results_df[results_df["target"] == "nas_4group"]
binary_results = results_df[results_df["target"].str.contains("binary")]

# 1. Fibrosis sweep
fib_out = os.path.join(OUTDIR, "plan1_sweep_fibrosis.csv")
fib_results.to_csv(fib_out, index=False)
print(f"Saved: {fib_out} ({len(fib_results)} rows)")

# 2. NAS sweep
nas_out = os.path.join(OUTDIR, "plan1_sweep_nas.csv")
nas_results.to_csv(nas_out, index=False)
print(f"Saved: {nas_out} ({len(nas_results)} rows)")

# 3. Binary sweep
bin_out = os.path.join(OUTDIR, "plan1_sweep_binary.csv")
binary_results.to_csv(bin_out, index=False)
print(f"Saved: {bin_out} ({len(binary_results)} rows)")

# 4. Confusion matrices
cm_out = os.path.join(OUTDIR, "plan1_confusion_matrices.json")
with open(cm_out, "w") as f:
    json.dump(all_cm, f, indent=2)
print(f"Saved: {cm_out} ({len(all_cm)} matrices)")

# 5. Feature importances
fi_df = pd.DataFrame(all_fi)
if len(fi_df) > 0:
    fi_df = fi_df.sort_values(["target", "model", "importance_mean"],
                               ascending=[True, True, False])
fi_out = os.path.join(OUTDIR, "plan1_feature_importances.csv")
fi_df.to_csv(fi_out, index=False)
print(f"Saved: {fi_out} ({len(fi_df)} rows)")

# =============================================================================
# AGGREGATE SUMMARY
# =============================================================================
print(f"\n{'='*60}")
print("AGGREGATE SUMMARY")
print(f"{'='*60}")

for target_name in results_df["target"].unique():
    print(f"\n--- {target_name} ---")
    sub = results_df[results_df["target"] == target_name]
    is_binary = "binary" in target_name

    for model_name in sub["model"].unique():
        msub = sub[sub["model"] == model_name]

        acc_mean = msub["accuracy"].mean()
        acc_std = msub["accuracy"].std()
        bal_mean = msub["balanced_accuracy"].mean()
        mae_mean = msub["mae"].mean()

        if is_binary:
            auroc_mean = msub["auroc"].mean() if "auroc" in msub else np.nan
            print(f"  {model_name:15s}: acc={acc_mean:.3f}+/-{acc_std:.3f} "
                  f"bal_acc={bal_mean:.3f} auroc={auroc_mean:.3f} "
                  f"({len(msub)} folds)")
        else:
            qwk_mean = msub["qwk"].mean() if "qwk" in msub else np.nan
            adj_mean = msub["adjacent_accuracy"].mean() if "adjacent_accuracy" in msub else np.nan
            print(f"  {model_name:15s}: acc={acc_mean:.3f}+/-{acc_std:.3f} "
                  f"bal_acc={bal_mean:.3f} qwk={qwk_mean:.3f} "
                  f"mae={mae_mean:.3f} adj={adj_mean:.3f} "
                  f"({len(msub)} folds)")

# Top features across all tree-based models
if len(fi_df) > 0:
    print(f"\n--- Top 20 Features (averaged across tree-based models & folds) ---")
    top_fi = (fi_df.groupby("gene")["importance_mean"]
              .mean().sort_values(ascending=False).head(20))
    for gene, imp in top_fi.items():
        print(f"  {gene:20s}: {imp:.6f}")

print(f"\n=== 64_plan1_model_sweep.py completed: {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
