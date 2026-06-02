#!/usr/bin/env python3
"""
70_classical_on_embeddings.py
Classical ML models on VAE embeddings vs raw features comparison.

Models (all with class_weight="balanced"):
  1. Logistic elastic net (l1_ratio=0.5, saga solver)
  2. Random Forest (n_estimators=500)
  3. XGBoost (n_estimators=300)
  4. SVM-RBF (probability=True)
  5. LightGBM (n_estimators=300)

Experiments:
  A) Train each model on 64-dim VAE embeddings
  B) Train each model on top 3000 rank-transformed genes (raw)
  Compare A vs B on same LOCO folds

Targets: Fibrosis 5-class, NAS 4-group, F>=3 binary, NAS>=5 binary

Input:
  - results/staging_classifier/embeddings_loco/embeddings_loco_{fold}.csv (per-fold LOCO embeddings from Script 69)
  - results/staging_classifier/embeddings_all_samples.csv (full-data embeddings, exploratory only)
  - results/staging_classifier/prepared_data.h5 (raw features + metadata)

Output:
  - embedding_model_results.csv (model x target x fold x metric)
  - embedding_vs_raw_comparison.csv (paired comparison per fold)
  - embedding_feature_summary.csv (most predictive embedding dimensions)

SLURM: cpu partition, 8 CPUs, 32GB RAM, 48h
Env:    micromamba activate rapids_singlecell
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import h5py
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, balanced_accuracy_score, cohen_kappa_score,
    f1_score, classification_report,
)
import xgboost as xgb
import lightgbm as lgb

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR = os.path.join(INT, "results/staging_classifier")
os.makedirs(OUTDIR, exist_ok=True)

EMB_PATH = os.path.join(OUTDIR, "embeddings_all_samples.csv")  # exploratory only
EMB_LOCO_DIR = os.path.join(OUTDIR, "embeddings_loco")  # per-fold LOCO embeddings
H5_PATH = os.path.join(OUTDIR, "prepared_data.h5")

print("=" * 60)
print("70: Classical ML on VAE Embeddings vs Raw Features")
print("=" * 60)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("\nLoading data...")

# VAE embeddings
emb_df = pd.read_csv(EMB_PATH, index_col="sample_id")
emb_cols = [c for c in emb_df.columns if c.startswith("z")]
print(f"  Embeddings: {emb_df.shape[0]} samples x {len(emb_cols)} dims")

# Raw expression + metadata from HDF5
with h5py.File(H5_PATH, "r") as h5:
    raw_expr = h5["rank_expression"][:].T.astype(np.float32)  # samples x genes
    gene_names = [g.decode() for g in h5["gene_names"][:]]
    sample_ids = [s.decode() for s in h5["sample_ids"][:]]

    meta = {}
    for key in h5["metadata"].keys():
        raw = h5["metadata"][key][:]
        if raw.dtype.kind == "S":
            meta[key] = np.array([v.decode() for v in raw])
        else:
            meta[key] = raw.astype(np.float32)

N_SAMPLES, N_GENES = raw_expr.shape
print(f"  Raw features: {N_SAMPLES} samples x {N_GENES} genes")

# Build metadata DataFrame aligned with sample order
meta_df = pd.DataFrame(meta, index=sample_ids)
meta_df.index.name = "sample_id"

# Ensure embedding and raw data alignment
assert list(emb_df.index) == sample_ids, "Embedding/raw sample order mismatch"


# ---------------------------------------------------------------------------
# Target definitions
# ---------------------------------------------------------------------------
TARGETS = {
    "fib_5class": {
        "col": "fib_stage",
        "n_classes": 5,
        "task": "multiclass",
        "description": "Fibrosis F0-F4",
    },
    "nas_4group": {
        "col": "nas_group4",
        "n_classes": 4,
        "task": "multiclass",
        "description": "NAS 4-group (Low/Mod/High/VeryHigh)",
    },
    "fib_ge3": {
        "col": "fib_ge3",
        "n_classes": 2,
        "task": "binary",
        "description": "F>=3 (advanced fibrosis)",
    },
    "nas_ge5": {
        "col": "nas_ge5",
        "n_classes": 2,
        "task": "binary",
        "description": "NAS>=5 (definite NASH)",
    },
}


# ---------------------------------------------------------------------------
# Model factory
# ---------------------------------------------------------------------------
def get_models(n_classes, random_state=42):
    """Return dict of model_name -> fitted model constructor."""
    models = {}

    # 1. Logistic elastic net
    models["LogisticElasticNet"] = LogisticRegression(
        penalty="elasticnet", l1_ratio=0.5, solver="saga",
        class_weight="balanced", max_iter=5000, random_state=random_state,
        C=1.0, n_jobs=-1,
    )

    # 2. Random Forest
    models["RandomForest"] = RandomForestClassifier(
        n_estimators=500, class_weight="balanced",
        random_state=random_state, n_jobs=-1,
    )

    # 3. XGBoost
    xgb_params = dict(
        n_estimators=300, learning_rate=0.05, max_depth=6,
        random_state=random_state, n_jobs=-1,
        eval_metric="mlogloss" if n_classes > 2 else "logloss",
        use_label_encoder=False,
    )
    # XGBoost handles class_weight via sample_weight in fit
    models["XGBoost"] = xgb.XGBClassifier(**xgb_params)

    # 4. SVM-RBF
    models["SVM_RBF"] = SVC(
        kernel="rbf", probability=True, class_weight="balanced",
        random_state=random_state,
    )

    # 5. LightGBM
    models["LightGBM"] = lgb.LGBMClassifier(
        n_estimators=300, learning_rate=0.05, max_depth=-1,
        class_weight="balanced", random_state=random_state,
        n_jobs=-1, verbose=-1,
    )

    return models


# ---------------------------------------------------------------------------
# Compute sample weights for XGBoost (which doesn't take class_weight)
# ---------------------------------------------------------------------------
def compute_sample_weights(y):
    """Compute balanced sample weights for XGBoost."""
    classes, counts = np.unique(y, return_counts=True)
    n = len(y)
    n_classes = len(classes)
    weights = np.zeros(n)
    for cls, cnt in zip(classes, counts):
        weights[y == cls] = n / (n_classes * cnt)
    return weights


# ---------------------------------------------------------------------------
# AUROC helper (handles multiclass via OVR)
# ---------------------------------------------------------------------------
def safe_auroc(y_true, y_prob, n_classes):
    """Compute AUROC, handling edge cases."""
    try:
        if n_classes == 2:
            if y_prob.ndim == 2:
                return roc_auc_score(y_true, y_prob[:, 1])
            return roc_auc_score(y_true, y_prob)
        else:
            return roc_auc_score(y_true, y_prob, multi_class="ovr", average="weighted")
    except (ValueError, IndexError):
        return np.nan


# ---------------------------------------------------------------------------
# Run a single experiment: model x target x fold
# ---------------------------------------------------------------------------
def evaluate_fold(model, X_train, y_train, X_test, y_test, n_classes, model_name):
    """Fit model on train, evaluate on test. Returns dict of metrics."""
    scaler = StandardScaler()
    X_train_sc = scaler.fit_transform(X_train)
    X_test_sc = scaler.transform(X_test)

    # Handle XGBoost sample weights
    fit_kwargs = {}
    if model_name == "XGBoost":
        fit_kwargs["sample_weight"] = compute_sample_weights(y_train)

    model.fit(X_train_sc, y_train, **fit_kwargs)
    y_pred = model.predict(X_test_sc)

    # Probability estimates
    try:
        y_prob = model.predict_proba(X_test_sc)
    except Exception:
        y_prob = None

    metrics = {
        "balanced_accuracy": balanced_accuracy_score(y_test, y_pred),
        "cohen_kappa": cohen_kappa_score(y_test, y_pred,
                                         weights="quadratic" if n_classes > 2 else None),
        "f1_macro": f1_score(y_test, y_pred, average="macro", zero_division=0),
    }

    if y_prob is not None:
        metrics["auroc"] = safe_auroc(y_test, y_prob, n_classes)
    else:
        metrics["auroc"] = np.nan

    return metrics


# ---------------------------------------------------------------------------
# LOCO cross-validation
# ---------------------------------------------------------------------------
def run_loco_experiments(X_raw, meta_df, target_info, target_name, sample_ids):
    """Run all models on embedding and raw features for one target.

    Embeddings are loaded per-fold from EMB_LOCO_DIR to avoid data leakage
    (the full-data VAE saw all samples, including the test fold).
    """
    col = target_info["col"]
    n_classes = target_info["n_classes"]
    task = target_info["task"]

    y = meta_df[col].values.astype(float)
    valid = y >= 0
    folds_fib = meta_df["loco_fold_fibrosis"].values
    folds_nas = meta_df["loco_fold_nas"].values

    # Choose appropriate folds based on target
    if "nas" in target_name:
        folds = folds_nas
    else:
        folds = folds_fib

    unique_folds = sorted(set(f for f in folds if f != "excluded" and f != "NA"))
    if len(unique_folds) == 0:
        print(f"    No LOCO folds for {target_name}, skipping")
        return []

    results = []

    for fold_name in unique_folds:
        test_mask = (folds == fold_name) & valid
        train_mask = (folds != fold_name) & (folds != "excluded") & (folds != "NA") & valid

        if test_mask.sum() < 5 or train_mask.sum() < 10:
            print(f"    Fold {fold_name}: too few samples "
                  f"(train={train_mask.sum()}, test={test_mask.sum()}), skipping")
            continue

        y_train = y[train_mask].astype(int)
        y_test = y[test_mask].astype(int)

        # Check class coverage
        train_classes = set(np.unique(y_train))
        test_classes = set(np.unique(y_test))
        if len(train_classes) < 2:
            print(f"    Fold {fold_name}: <2 classes in train, skipping")
            continue

        # Load fold-specific LOCO embeddings (VAE trained without this fold)
        loco_emb_path = os.path.join(EMB_LOCO_DIR, f"embeddings_loco_{fold_name}.csv")
        if not os.path.exists(loco_emb_path):
            print(f"    Fold {fold_name}: LOCO embedding file not found at "
                  f"{loco_emb_path}, skipping")
            continue
        loco_emb_df = pd.read_csv(loco_emb_path, index_col="sample_id")
        loco_emb_cols = [c for c in loco_emb_df.columns if c.startswith("z")]
        # Align LOCO embeddings with master sample order
        loco_emb_aligned = loco_emb_df.reindex(sample_ids)
        if loco_emb_aligned[loco_emb_cols].isna().all(axis=None):
            print(f"    Fold {fold_name}: LOCO embeddings have no matching samples, skipping")
            continue
        X_emb_fold = loco_emb_aligned[loco_emb_cols].values

        X_emb_train = X_emb_fold[train_mask]
        X_emb_test = X_emb_fold[test_mask]
        X_raw_train = X_raw[train_mask]
        X_raw_test = X_raw[test_mask]

        models_emb = get_models(n_classes)
        models_raw = get_models(n_classes)

        for model_name in models_emb:
            # A) Embeddings
            try:
                m_emb = evaluate_fold(
                    models_emb[model_name], X_emb_train, y_train,
                    X_emb_test, y_test, n_classes, model_name,
                )
                m_emb.update({
                    "model": model_name, "target": target_name, "fold": fold_name,
                    "input": "embedding",
                    "n_train": int(train_mask.sum()), "n_test": int(test_mask.sum()),
                    "n_train_classes": len(train_classes),
                    "n_test_classes": len(test_classes),
                })
                results.append(m_emb)
            except Exception as e:
                print(f"    {model_name} embedding fold={fold_name}: {e}")

            # B) Raw features
            try:
                m_raw = evaluate_fold(
                    models_raw[model_name], X_raw_train, y_train,
                    X_raw_test, y_test, n_classes, model_name,
                )
                m_raw.update({
                    "model": model_name, "target": target_name, "fold": fold_name,
                    "input": "raw",
                    "n_train": int(train_mask.sum()), "n_test": int(test_mask.sum()),
                    "n_train_classes": len(train_classes),
                    "n_test_classes": len(test_classes),
                })
                results.append(m_raw)
            except Exception as e:
                print(f"    {model_name} raw fold={fold_name}: {e}")

    return results


# ---------------------------------------------------------------------------
# Embedding dimension importance analysis
# ---------------------------------------------------------------------------
def embedding_importance(X_emb, meta_df, emb_cols):
    """Rank embedding dimensions by their predictive power for fibrosis."""
    y = meta_df["fib_stage"].values.astype(float)
    valid = y >= 0
    if valid.sum() < 20:
        return pd.DataFrame()

    X = X_emb[valid]
    yv = y[valid].astype(int)

    # Use per-dimension Spearman correlation with fibrosis stage as a simple metric
    from scipy.stats import spearmanr
    importance = []
    for i, col in enumerate(emb_cols):
        rho, pval = spearmanr(X[:, i], yv)
        importance.append({
            "dimension": col,
            "spearman_rho": rho,
            "abs_rho": abs(rho),
            "pvalue": pval,
        })

    # Also fit a quick Random Forest to get feature importance
    try:
        rf = RandomForestClassifier(n_estimators=200, class_weight="balanced",
                                     random_state=42, n_jobs=-1)
        scaler = StandardScaler()
        X_sc = scaler.fit_transform(X)
        rf.fit(X_sc, yv)
        rf_imp = rf.feature_importances_
        for i, row in enumerate(importance):
            row["rf_importance"] = rf_imp[i]
    except Exception:
        pass

    df = pd.DataFrame(importance).sort_values("abs_rho", ascending=False)
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()

    X_emb_full = emb_df[emb_cols].values  # full-data embeddings (exploratory only)
    X_raw = raw_expr  # samples x genes (already rank-transformed)

    all_results = []

    for target_name, target_info in TARGETS.items():
        print(f"\n--- Target: {target_name} ({target_info['description']}) ---")
        results = run_loco_experiments(
            X_raw, meta_df, target_info, target_name, sample_ids,
        )
        all_results.extend(results)
        print(f"    {len(results)} result rows")

    # -----------------------------------------------------------------------
    # Save model results
    # -----------------------------------------------------------------------
    results_df = pd.DataFrame(all_results)
    results_path = os.path.join(OUTDIR, "embedding_model_results.csv")
    results_df.to_csv(results_path, index=False)
    print(f"\nModel results saved: {results_path} ({len(results_df)} rows)")

    # -----------------------------------------------------------------------
    # Paired embedding vs raw comparison
    # -----------------------------------------------------------------------
    print("\n--- Embedding vs Raw Comparison ---")
    comparison_rows = []

    for target_name in TARGETS:
        for model_name in ["LogisticElasticNet", "RandomForest", "XGBoost",
                           "SVM_RBF", "LightGBM"]:
            emb_rows = results_df[
                (results_df["target"] == target_name) &
                (results_df["model"] == model_name) &
                (results_df["input"] == "embedding")
            ]
            raw_rows = results_df[
                (results_df["target"] == target_name) &
                (results_df["model"] == model_name) &
                (results_df["input"] == "raw")
            ]

            if len(emb_rows) == 0 or len(raw_rows) == 0:
                continue

            # Merge by fold for paired comparison
            merged = emb_rows.merge(raw_rows, on="fold", suffixes=("_emb", "_raw"))

            for metric in ["auroc", "balanced_accuracy", "cohen_kappa", "f1_macro"]:
                emb_vals = merged[f"{metric}_emb"].dropna().values
                raw_vals = merged[f"{metric}_raw"].dropna().values

                if len(emb_vals) == 0:
                    continue

                comparison_rows.append({
                    "target": target_name,
                    "model": model_name,
                    "metric": metric,
                    "emb_mean": np.mean(emb_vals),
                    "emb_std": np.std(emb_vals),
                    "raw_mean": np.mean(raw_vals),
                    "raw_std": np.std(raw_vals),
                    "delta": np.mean(emb_vals) - np.mean(raw_vals),
                    "n_folds": len(emb_vals),
                    "emb_wins": int(np.sum(emb_vals > raw_vals)),
                })

    comparison_df = pd.DataFrame(comparison_rows)
    comp_path = os.path.join(OUTDIR, "embedding_vs_raw_comparison.csv")
    comparison_df.to_csv(comp_path, index=False)
    print(f"Comparison saved: {comp_path} ({len(comparison_df)} rows)")

    # Print summary
    if len(comparison_df) > 0:
        print("\n  Summary (mean AUROC, embedding vs raw):")
        auroc_comp = comparison_df[comparison_df["metric"] == "auroc"]
        for _, row in auroc_comp.iterrows():
            delta_str = f"+{row['delta']:.3f}" if row["delta"] >= 0 else f"{row['delta']:.3f}"
            print(f"    {row['target']:15s} {row['model']:20s}: "
                  f"emb={row['emb_mean']:.3f} raw={row['raw_mean']:.3f} "
                  f"delta={delta_str} (emb wins {row['emb_wins']}/{row['n_folds']})")

    # -----------------------------------------------------------------------
    # Embedding dimension importance
    # -----------------------------------------------------------------------
    # NOTE: Dimension importance uses full-data embeddings (exploratory/descriptive only,
    # not used for model evaluation). This is acceptable since it's not part of LOCO-CV.
    print("\n--- Embedding Dimension Importance (exploratory, full-data embeddings) ---")
    imp_df = embedding_importance(X_emb_full, meta_df, emb_cols)
    if len(imp_df) > 0:
        imp_path = os.path.join(OUTDIR, "embedding_feature_summary.csv")
        imp_df.to_csv(imp_path, index=False)
        print(f"Feature summary saved: {imp_path}")

        print("\n  Top 10 dimensions by |Spearman rho| with fibrosis:")
        for _, row in imp_df.head(10).iterrows():
            rf_str = f"  RF_imp={row['rf_importance']:.4f}" if "rf_importance" in row else ""
            print(f"    {row['dimension']:5s}: rho={row['spearman_rho']:+.3f} "
                  f"p={row['pvalue']:.2e}{rf_str}")

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    dt = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"70_classical_on_embeddings.py completed in {dt / 60:.1f} min")
    print(f"  Results:    {results_path}")
    print(f"  Comparison: {comp_path}")
    if len(imp_df) > 0:
        print(f"  Features:   {imp_path}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
