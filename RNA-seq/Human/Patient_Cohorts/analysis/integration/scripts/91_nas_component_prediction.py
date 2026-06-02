#!/usr/bin/env python3
import os as _os
_os.environ["RAPIDS_NO_INITIALIZE"] = "1"
"""
91_nas_component_prediction.py
Predict NAS components (steatosis, inflammation, ballooning) separately,
then sum to get NAS. Compare component-sum vs direct ordinal prediction.

Strategy:
  For each component (steatosis 0-3, inflammation 0-2, ballooning 0-2):
    - Leave-one-out CV across 78 GSE130970 samples
    - Models: ordinal logistic, RandomForest, XGBoost
    - Feature sets: (a) NAS-VAE 64-dim, (b) top 500 genes by variance, (c) combined
  Then NAS_predicted = steatosis_pred + inflammation_pred + ballooning_pred

Input:
  - results/staging_classifier/prepared_data.h5 (expression + metadata)
  - results/staging_classifier/nas_embeddings_all_samples.csv (NAS-VAE, fallback to embeddings_all_samples.csv)
  - results/staging_classifier/nas_component_data.csv (78 samples with steatosis, inflammation, ballooning)

Output:
  - results/staging_classifier/component_prediction_results.csv
  - results/staging_classifier/component_vs_direct_comparison.csv

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
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    cohen_kappa_score, balanced_accuracy_score, mean_absolute_error,
    f1_score,
)
import xgboost as xgb

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

H5_PATH = os.path.join(OUTDIR, "prepared_data.h5")
NAS_EMB_PATH = os.path.join(OUTDIR, "nas_embeddings_all_samples.csv")
FALLBACK_EMB_PATH = os.path.join(OUTDIR, "embeddings_all_samples.csv")
COMPONENT_PATH = os.path.join(OUTDIR, "nas_component_data.csv")

print("=" * 60)
print("91: NAS Component Prediction (LOO on 78 GSE130970 samples)")
print("=" * 60)
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading data...")

# 1) NAS component labels (78 samples from GSE130970)
comp_df = pd.read_csv(COMPONENT_PATH)
comp_samples = set(comp_df["sample_id"].values)
print(f"  Component labels: {len(comp_df)} samples")
print(f"    Steatosis range: {comp_df['steatosis_grade'].min()}-{comp_df['steatosis_grade'].max()}")
print(f"    Inflammation range: {comp_df['lobular_inflammation_grade'].min()}-{comp_df['lobular_inflammation_grade'].max()}")
print(f"    Ballooning range: {comp_df['ballooning_grade'].min()}-{comp_df['ballooning_grade'].max()}")

# 2) VAE embeddings (prefer NAS-VAE, fallback to fibrosis-VAE)
if os.path.exists(NAS_EMB_PATH):
    emb_path_used = NAS_EMB_PATH
    print(f"  Using NAS-VAE embeddings: {NAS_EMB_PATH}")
else:
    emb_path_used = FALLBACK_EMB_PATH
    print(f"  NAS-VAE not found, using fibrosis-VAE: {FALLBACK_EMB_PATH}")

emb_df = pd.read_csv(emb_path_used, index_col="sample_id")
emb_cols = [c for c in emb_df.columns if c.startswith("z")]
print(f"  Embeddings: {emb_df.shape[0]} samples x {len(emb_cols)} dims")

# 3) Raw expression from HDF5 (for top-variance gene features)
with h5py.File(H5_PATH, "r") as h5:
    raw_expr = h5["rank_expression"][:].T.astype(np.float32)  # samples x genes
    gene_names = [g.decode() for g in h5["gene_names"][:]]
    sample_ids = [s.decode() for s in h5["sample_ids"][:]]

N_SAMPLES, N_GENES = raw_expr.shape
print(f"  Raw expression: {N_SAMPLES} samples x {N_GENES} genes")

# Build sample-indexed DataFrames for alignment
expr_df = pd.DataFrame(raw_expr, index=sample_ids, columns=gene_names)

# ---------------------------------------------------------------------------
# Subset to 78 component-labeled samples
# ---------------------------------------------------------------------------
# Ensure string consistency for sample ID matching
comp_sample_ids = comp_df["sample_id"].astype(str).str.strip().values
expr_df.index = expr_df.index.astype(str).str.strip()
emb_df.index = emb_df.index.astype(str).str.strip()

# Find matching samples (some may be QC-filtered out)
comp_in_expr = [s for s in comp_sample_ids if s in expr_df.index]
comp_in_emb = [s for s in comp_sample_ids if s in emb_df.index]
comp_sample_ids = np.array([s for s in comp_sample_ids if s in expr_df.index and s in emb_df.index])
print(f"  Component samples in expression: {len(comp_in_expr)}/{len(comp_df)}")
print(f"  Component samples in embeddings: {len(comp_in_emb)}/{len(comp_df)}")
print(f"  Component samples in both: {len(comp_sample_ids)}")
assert len(comp_sample_ids) >= 50, f"Too few matched component samples: {len(comp_sample_ids)}"

N_COMP = len(comp_sample_ids)
print(f"\n  Working with {N_COMP} component-labeled samples")

# Extract aligned feature matrices for the 78 samples
X_emb = emb_df.loc[comp_sample_ids, emb_cols].values.astype(np.float32)
X_expr_full = expr_df.loc[comp_sample_ids].values.astype(np.float32)

# Select top 500 genes by variance across these 78 samples
gene_var = X_expr_full.var(axis=0)
top500_idx = np.argsort(gene_var)[::-1][:500]
X_top500 = X_expr_full[:, top500_idx]
top500_genes = [gene_names[i] for i in top500_idx]
print(f"  Top 500 genes by variance selected (var range: {gene_var[top500_idx[0]]:.2f} - {gene_var[top500_idx[-1]]:.2f})")

# Combined features: embeddings + top 500 genes
X_combined = np.hstack([X_emb, X_top500])
print(f"  Feature sets: embedding={X_emb.shape[1]}, top500={X_top500.shape[1]}, combined={X_combined.shape[1]}")

# Component labels — filter comp_df to matched samples only (same order as X arrays)
comp_df_matched = comp_df.set_index(comp_df["sample_id"].astype(str).str.strip()).loc[comp_sample_ids].reset_index(drop=True)
assert len(comp_df_matched) == N_COMP, f"comp_df_matched has {len(comp_df_matched)} rows, expected {N_COMP}"

y_steatosis = comp_df_matched["steatosis_grade"].values.astype(int)
y_inflammation = comp_df_matched["lobular_inflammation_grade"].values.astype(int)
y_ballooning = comp_df_matched["ballooning_grade"].values.astype(int)
y_nas = comp_df_matched["nas_score_original"].values.astype(int)

COMPONENTS = {
    "steatosis": {"y": y_steatosis, "n_classes": 4, "range": "0-3"},
    "inflammation": {"y": y_inflammation, "n_classes": 3, "range": "0-2"},
    "ballooning": {"y": y_ballooning, "n_classes": 3, "range": "0-2"},
}

FEATURE_SETS = {
    "embedding": X_emb,
    "top500_genes": X_top500,
    "combined": X_combined,
}


# ---------------------------------------------------------------------------
# Model factory
# ---------------------------------------------------------------------------
def get_models(n_classes, random_state=42):
    """Return dict of model_name -> model instance."""
    models = {}

    # 1. Ordinal logistic (approximated via LogisticRegression with ordinal encoding)
    # sklearn doesn't have native ordinal regression, so we use multinomial
    # which is the standard ML approximation
    models["LogisticOrdinal"] = LogisticRegression(
        penalty="l2", solver="lbfgs", class_weight="balanced",
        max_iter=5000, random_state=random_state, C=1.0,
    )

    # 2. Random Forest
    models["RandomForest"] = RandomForestClassifier(
        n_estimators=500, class_weight="balanced",
        random_state=random_state, n_jobs=-1,
    )

    # 3. XGBoost
    models["XGBoost"] = xgb.XGBClassifier(
        n_estimators=300, learning_rate=0.05, max_depth=4,
        random_state=random_state, n_jobs=-1,
        eval_metric="mlogloss" if n_classes > 2 else "logloss",
    )

    return models


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
# LOO CV for a single component
# ---------------------------------------------------------------------------
def loo_cv_component(X, y, n_classes, component_name, feature_name, random_state=42):
    """Leave-one-out CV for a single component. Returns per-sample predictions."""
    n = len(y)
    models = get_models(n_classes, random_state)
    results = []

    for model_name, model_template in models.items():
        predictions = np.full(n, -1, dtype=int)

        for i in range(n):
            # LOO split
            test_idx = [i]
            train_idx = [j for j in range(n) if j != i]

            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # Check: at least 2 classes in training
            if len(np.unique(y_train)) < 2:
                predictions[i] = int(np.median(y_train))
                continue

            # Scale features
            scaler = StandardScaler()
            X_train_sc = scaler.fit_transform(X_train)
            X_test_sc = scaler.transform(X_test)

            # Fit
            # Need fresh model each fold
            fold_model = get_models(n_classes, random_state)[model_name]

            fit_kwargs = {}
            if model_name == "XGBoost":
                fit_kwargs["sample_weight"] = compute_sample_weights(y_train)

            try:
                fold_model.fit(X_train_sc, y_train, **fit_kwargs)
                pred = fold_model.predict(X_test_sc)
                predictions[i] = int(pred[0])
            except Exception as e:
                # Fallback: predict mode
                from scipy.stats import mode as scipy_mode
                mode_result = scipy_mode(y_train, keepdims=True)
                predictions[i] = int(mode_result.mode[0])

        # Compute metrics for this model
        valid = predictions >= 0
        y_valid = y[valid]
        p_valid = predictions[valid]

        qwk = cohen_kappa_score(y_valid, p_valid, weights="quadratic")
        ba = balanced_accuracy_score(y_valid, p_valid)
        mae = mean_absolute_error(y_valid, p_valid)
        f1 = f1_score(y_valid, p_valid, average="macro", zero_division=0)

        results.append({
            "component": component_name,
            "feature_set": feature_name,
            "model": model_name,
            "n_samples": int(valid.sum()),
            "qwk": qwk,
            "balanced_accuracy": ba,
            "mae": mae,
            "f1_macro": f1,
        })

        print(f"    {component_name}/{feature_name}/{model_name}: "
              f"QWK={qwk:.3f}, BA={ba:.3f}, MAE={mae:.3f}")

    return results, predictions


# ---------------------------------------------------------------------------
# LOO CV for direct NAS prediction
# ---------------------------------------------------------------------------
def loo_cv_direct_nas(X, y_nas, feature_name, random_state=42):
    """LOO CV predicting NAS directly as ordinal."""
    n = len(y_nas)
    n_classes = len(np.unique(y_nas))
    models = get_models(n_classes, random_state)
    results = []

    for model_name, model_template in models.items():
        predictions = np.full(n, -1, dtype=int)

        for i in range(n):
            train_idx = [j for j in range(n) if j != i]
            test_idx = [i]

            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y_nas[train_idx], y_nas[test_idx]

            if len(np.unique(y_train)) < 2:
                predictions[i] = int(np.median(y_train))
                continue

            scaler = StandardScaler()
            X_train_sc = scaler.fit_transform(X_train)
            X_test_sc = scaler.transform(X_test)

            fold_model = get_models(n_classes, random_state)[model_name]

            fit_kwargs = {}
            if model_name == "XGBoost":
                fit_kwargs["sample_weight"] = compute_sample_weights(y_train)

            try:
                fold_model.fit(X_train_sc, y_train, **fit_kwargs)
                pred = fold_model.predict(X_test_sc)
                predictions[i] = int(pred[0])
            except Exception:
                from scipy.stats import mode as scipy_mode
                mode_result = scipy_mode(y_train, keepdims=True)
                predictions[i] = int(mode_result.mode[0])

        valid = predictions >= 0
        y_valid = y_nas[valid]
        p_valid = predictions[valid]

        qwk = cohen_kappa_score(y_valid, p_valid, weights="quadratic")
        ba = balanced_accuracy_score(y_valid, p_valid)
        mae = mean_absolute_error(y_valid, p_valid)
        f1 = f1_score(y_valid, p_valid, average="macro", zero_division=0)

        results.append({
            "approach": "direct_nas",
            "feature_set": feature_name,
            "model": model_name,
            "n_samples": int(valid.sum()),
            "qwk": qwk,
            "balanced_accuracy": ba,
            "mae": mae,
            "f1_macro": f1,
        })

        print(f"    direct_nas/{feature_name}/{model_name}: "
              f"QWK={qwk:.3f}, BA={ba:.3f}, MAE={mae:.3f}")

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    np.random.seed(42)

    all_component_results = []
    # Store per-sample component predictions for the component-sum approach
    # Key: (feature_set, model) -> {component: predictions_array}
    component_preds = {}

    # -----------------------------------------------------------------------
    # STEP 1: Per-component LOO CV
    # -----------------------------------------------------------------------
    print("\n=== STEP 1: Per-Component LOO CV ===")

    for feat_name, X in FEATURE_SETS.items():
        print(f"\n  Feature set: {feat_name} ({X.shape[1]} features)")

        for comp_name, comp_info in COMPONENTS.items():
            res, preds = loo_cv_component(
                X, comp_info["y"], comp_info["n_classes"],
                comp_name, feat_name,
            )
            all_component_results.extend(res)

            # Store predictions per model for component-sum later
            for r in res:
                key = (feat_name, r["model"])
                if key not in component_preds:
                    component_preds[key] = {}
                # We need per-sample predictions — rerun a bit differently
                # Actually, we need to keep the last predictions array.
                # The loo_cv_component returns predictions for the LAST model.
                # Fix: store per-model predictions separately.

    # Re-run to capture per-sample predictions for component-sum computation
    # (refactored to store predictions per model)
    print("\n=== STEP 2: Component-Sum Predictions ===")
    component_preds_all = {}  # (feat, model, comp) -> predictions

    for feat_name, X in FEATURE_SETS.items():
        for comp_name, comp_info in COMPONENTS.items():
            n = len(comp_info["y"])
            n_classes = comp_info["n_classes"]
            models_dict = get_models(n_classes, random_state=42)

            for model_name in models_dict:
                predictions = np.full(n, -1, dtype=int)

                for i in range(n):
                    train_idx = [j for j in range(n) if j != i]
                    X_train, X_test = X[train_idx], X[[i]]
                    y_train = comp_info["y"][train_idx]

                    if len(np.unique(y_train)) < 2:
                        predictions[i] = int(np.median(y_train))
                        continue

                    scaler = StandardScaler()
                    X_train_sc = scaler.fit_transform(X_train)
                    X_test_sc = scaler.transform(X_test)

                    fold_model = get_models(n_classes, 42)[model_name]
                    fit_kwargs = {}
                    if model_name == "XGBoost":
                        fit_kwargs["sample_weight"] = compute_sample_weights(y_train)

                    try:
                        fold_model.fit(X_train_sc, y_train, **fit_kwargs)
                        pred = fold_model.predict(X_test_sc)
                        predictions[i] = int(pred[0])
                    except Exception:
                        from scipy.stats import mode as scipy_mode
                        mode_result = scipy_mode(y_train, keepdims=True)
                        predictions[i] = int(mode_result.mode[0])

                component_preds_all[(feat_name, model_name, comp_name)] = predictions

    # Compute component-sum NAS predictions
    component_sum_results = []

    for feat_name in FEATURE_SETS:
        for model_name in ["LogisticOrdinal", "RandomForest", "XGBoost"]:
            p_steat = component_preds_all.get((feat_name, model_name, "steatosis"))
            p_inflam = component_preds_all.get((feat_name, model_name, "inflammation"))
            p_balloon = component_preds_all.get((feat_name, model_name, "ballooning"))

            if p_steat is None or p_inflam is None or p_balloon is None:
                continue

            nas_pred_sum = p_steat + p_inflam + p_balloon
            valid = (p_steat >= 0) & (p_inflam >= 0) & (p_balloon >= 0)

            y_valid = y_nas[valid]
            p_valid = nas_pred_sum[valid]

            qwk = cohen_kappa_score(y_valid, p_valid, weights="quadratic")
            ba = balanced_accuracy_score(y_valid, p_valid)
            mae = mean_absolute_error(y_valid, p_valid)

            component_sum_results.append({
                "approach": "component_sum",
                "feature_set": feat_name,
                "model": model_name,
                "n_samples": int(valid.sum()),
                "qwk": qwk,
                "balanced_accuracy": ba,
                "mae": mae,
            })

            print(f"  component_sum/{feat_name}/{model_name}: "
                  f"QWK={qwk:.3f}, BA={ba:.3f}, MAE={mae:.3f}")

    # -----------------------------------------------------------------------
    # STEP 3: Direct NAS prediction (for comparison)
    # -----------------------------------------------------------------------
    print("\n=== STEP 3: Direct NAS LOO CV (for comparison) ===")

    direct_nas_results = []
    for feat_name, X in FEATURE_SETS.items():
        print(f"\n  Feature set: {feat_name}")
        res = loo_cv_direct_nas(X, y_nas, feat_name)
        direct_nas_results.extend(res)

    # -----------------------------------------------------------------------
    # STEP 4: Save results
    # -----------------------------------------------------------------------
    print("\n=== STEP 4: Saving Results ===")

    # Component-level results
    comp_results_df = pd.DataFrame(all_component_results)
    comp_out = os.path.join(OUTDIR, "component_prediction_results.csv")
    comp_results_df.to_csv(comp_out, index=False)
    print(f"  Component results: {comp_out} ({len(comp_results_df)} rows)")

    # Component-sum vs direct comparison
    comparison_rows = []

    for feat_name in FEATURE_SETS:
        for model_name in ["LogisticOrdinal", "RandomForest", "XGBoost"]:
            # Find component-sum result
            cs_match = [r for r in component_sum_results
                        if r["feature_set"] == feat_name and r["model"] == model_name]
            # Find direct result
            dr_match = [r for r in direct_nas_results
                        if r["feature_set"] == feat_name and r["model"] == model_name]

            if cs_match and dr_match:
                cs = cs_match[0]
                dr = dr_match[0]
                comparison_rows.append({
                    "feature_set": feat_name,
                    "model": model_name,
                    "component_sum_qwk": cs["qwk"],
                    "direct_nas_qwk": dr["qwk"],
                    "delta_qwk": cs["qwk"] - dr["qwk"],
                    "component_sum_mae": cs["mae"],
                    "direct_nas_mae": dr["mae"],
                    "delta_mae": cs["mae"] - dr["mae"],
                    "component_sum_ba": cs["balanced_accuracy"],
                    "direct_nas_ba": dr["balanced_accuracy"],
                    "delta_ba": cs["balanced_accuracy"] - dr["balanced_accuracy"],
                    "n_samples": cs["n_samples"],
                })

    comparison_df = pd.DataFrame(comparison_rows)
    comparison_out = os.path.join(OUTDIR, "component_vs_direct_comparison.csv")
    comparison_df.to_csv(comparison_out, index=False)
    print(f"  Comparison: {comparison_out} ({len(comparison_df)} rows)")

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    dt = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"91_nas_component_prediction.py completed in {dt / 60:.1f} min")
    print(f"  Component predictions: {comp_out}")
    print(f"  Component-sum vs direct: {comparison_out}")

    if len(comparison_df) > 0:
        best = comparison_df.loc[comparison_df["component_sum_qwk"].idxmax()]
        print(f"\n  Best component-sum: {best['feature_set']}/{best['model']}")
        print(f"    Component-sum QWK = {best['component_sum_qwk']:.3f}")
        print(f"    Direct NAS QWK    = {best['direct_nas_qwk']:.3f}")
        print(f"    Delta QWK         = {best['delta_qwk']:+.3f}")

    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
