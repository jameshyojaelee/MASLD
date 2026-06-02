#!/usr/bin/env python3
"""
45b_ml_target_prediction.py
============================================================================
ML-based MASLD target prioritization using multi-evidence features.

Compares 4 models: Elastic Net, Random Forest, LightGBM, XGBoost
Evaluates with stratified 5-fold CV (10 repeats).
Performs source-group ablation and SHAP analysis.

Input:  RNA-seq/results/multi_evidence/multi_evidence_atlas.csv (112 cols)
Output: RNA-seq/results/multi_evidence/ml_target_predictions.csv
        RNA-seq/results/multi_evidence/ml_model_comparison.csv
        RNA-seq/results/multi_evidence/ml_feature_importance.csv
        RNA-seq/results/multi_evidence/ml_ablation.csv
============================================================================
"""

import os
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score,
    precision_score, recall_score
)
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import xgboost as xgb
import lightgbm as lgb
import shap

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================================
# Configuration
# ============================================================================
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
ATLAS_FILE = f"{BASE}/RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
POS_CTRL_FILE = f"{BASE}/RNA-seq/results/validation/positive_control_validation.csv"
DRUG_FILE = f"{BASE}/RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"
OUT_DIR = f"{BASE}/RNA-seq/results/multi_evidence"

SEED = 42
N_FOLDS = 5
N_REPEATS = 10

# ============================================================================
# Feature definitions grouped by independent source
# ============================================================================
# Each source group maps to a list of (atlas_column, transform, new_name)
# transform: None = use as-is, "abs" = absolute value, "neglog10" = -log10,
#            "binary" = convert to 0/1, "max_cols" = take max of cols

SOURCE_FEATURES = {
    "S1_human_bulk": [
        ("dream_logFC", "abs", "abs_dream_logFC"),
        ("dream_padj", "neglog10", "neglog10_dream_padj"),
        ("dream_tstat", None, "dream_tstat"),
    ],
    "S2_mouse_bulk": [
        ("mouse_meta_logFC", "abs", "abs_mouse_meta_logFC"),
        ("mouse_meta_padj", "neglog10", "neglog10_mouse_meta_padj"),
        ("n_diets_sig", None, "n_diets_sig"),
    ],
    "S3_genetic_causal": [
        # MR (mr_beta, mr_pval) and sc-TWAS MR (sceqtl_twas_*) features removed
        # 2026-04-22 — MR ditched from paper; TWAS + COLOC + INTACT is the
        # causal framework.
        ("twas_z", "abs", "abs_twas_z"),
        ("twas_pval", "neglog10", "neglog10_twas_pval"),
        # Best COLOC across Ghodsian + Broadaway
        ("coloc_pp4", None, "coloc_pp4"),
        ("broadaway_coloc_pp4", None, "broadaway_coloc_pp4"),
        ("n_liver_enzyme_coloc", None, "n_liver_enzyme_coloc"),
        ("best_liver_enzyme_pp4", None, "best_liver_enzyme_pp4"),
        ("ieqtl_interaction_pval", "neglog10", "neglog10_ieqtl_pval"),
        ("pleiotropy_n_traits", None, "pleiotropy_n_traits"),
        ("n_coloc_sources", None, "n_coloc_sources"),
    ],
    "S4_essentiality": [
        ("essentiality_chronos", None, "essentiality_chronos"),
    ],
    "S5_epigenomic": [
        ("mouse_da_logFC", None, "mouse_da_logFC"),
        ("mouse_da_padj", "neglog10", "neglog10_mouse_da_padj"),
        ("hepatocyte_da_logFC", None, "hepatocyte_da_logFC"),
        ("hepatocyte_da_padj", "neglog10", "neglog10_hepatocyte_da_padj"),
        ("scenic_regulon_activity_diff", None, "scenic_regulon_activity_diff"),
        ("cross_species_promoter_conserved", "binary", "cross_species_promoter_conserved"),
    ],
    "S6_spatial": [
        ("spatial_morans_i", None, "spatial_morans_i"),
        ("spatial_is_svg", "binary", "spatial_is_svg"),
        ("spatial_hep_spatial_fc", None, "spatial_hep_spatial_fc"),
        ("spatial_niche_fc", None, "spatial_niche_fc"),
    ],
    "S7_singlecell": [
        ("sc_hepatocyte_logFC", None, "sc_hepatocyte_logFC"),
        ("sc_hepatocyte_padj", "neglog10", "neglog10_sc_hepatocyte_padj"),
        ("sc_n_celltypes_sig", None, "sc_n_celltypes_sig"),
        ("sc_max_abs_logFC", None, "sc_max_abs_logFC"),
        ("liana_n_diff_interactions", None, "liana_n_diff_interactions"),
    ],
    "Derived_concordance": [
        ("translatability_score", None, "translatability_score"),
        ("is_conserved", "binary", "is_conserved"),
    ],
    "Derived_sex": [
        ("dream_logFC_M", None, "dream_logFC_M"),
        ("dream_logFC_F", None, "dream_logFC_F"),
        ("sex_interaction_padj", "neglog10", "neglog10_sex_interaction_padj"),
    ],
    "Derived_pathway": [
        ("n_leading_edge_pathways", None, "n_leading_edge_pathways"),
    ],
}


def transform_feature(series, transform):
    """Apply transformation to a feature column."""
    if transform is None:
        return series.astype(float)
    elif transform == "abs":
        return series.astype(float).abs()
    elif transform == "neglog10":
        s = series.astype(float)
        # Clip to avoid -log10(0) = inf
        s = s.clip(lower=1e-300)
        return -np.log10(s)
    elif transform == "binary":
        return series.map({True: 1, False: 0, "TRUE": 1, "FALSE": 0}).fillna(0).astype(float)
    else:
        raise ValueError(f"Unknown transform: {transform}")


def build_feature_matrix(atlas):
    """Build the ML feature matrix from the atlas DataFrame."""
    feature_names = []
    source_groups = []
    X_cols = []

    for source, features in SOURCE_FEATURES.items():
        for col, tfm, name in features:
            if col not in atlas.columns:
                print(f"  WARNING: Column {col} not found, skipping")
                continue
            X_cols.append(transform_feature(atlas[col], tfm))
            feature_names.append(name)
            source_groups.append(source)

    X = pd.DataFrame(dict(zip(feature_names, X_cols)), index=atlas.index)

    # Replace inf with NaN
    X = X.replace([np.inf, -np.inf], np.nan)

    return X, feature_names, source_groups


def define_labels(atlas, pos_ctrl, drug_targets):
    """Define positive/negative labels."""
    # Expression-driven positive controls
    expr_genes = set(pos_ctrl[pos_ctrl["control_type"] == "Expression_driven"]["gene"])

    # Drug targets with Strong or Moderate support
    strong_mod = drug_targets[drug_targets["atlas_support"].isin(["Strong", "Moderate"])]
    drug_genes = set(strong_mod["target_gene"])

    positives = expr_genes | drug_genes
    positives_in_atlas = positives & set(atlas["human_symbol"])

    print(f"Positive controls: {len(positives_in_atlas)} "
          f"({len(expr_genes & set(atlas['human_symbol']))} expression-driven, "
          f"{len(drug_genes & set(atlas['human_symbol']))} drug targets)")

    y = atlas["human_symbol"].isin(positives_in_atlas).astype(int)
    print(f"Label distribution: {y.sum()} positives, {(~y.astype(bool)).sum()} negatives")

    return y, positives_in_atlas


def precision_at_k(y_true, y_score, k):
    """Precision at top-k predictions."""
    top_k_idx = np.argsort(y_score)[-k:]
    return y_true.iloc[top_k_idx].mean() if hasattr(y_true, 'iloc') else y_true[top_k_idx].mean()


# ============================================================================
# Main
# ============================================================================
def main():
    print("=" * 70)
    print("ML Target Prediction Pipeline")
    print("=" * 70)

    # Load data
    atlas = pd.read_csv(ATLAS_FILE)
    print(f"Atlas: {atlas.shape[0]} genes x {atlas.shape[1]} cols")

    pos_ctrl = pd.read_csv(POS_CTRL_FILE)
    drug_targets = pd.read_csv(DRUG_FILE)

    # Build features
    print("\nBuilding feature matrix...")
    X, feature_names, source_groups = build_feature_matrix(atlas)
    print(f"Feature matrix: {X.shape[0]} x {X.shape[1]}")
    print(f"Feature groups: {len(set(source_groups))}")
    for sg in dict.fromkeys(source_groups):
        n = source_groups.count(sg)
        print(f"  {sg}: {n} features")

    # Define labels
    print("\nDefining labels...")
    y, pos_genes = define_labels(atlas, pos_ctrl, drug_targets)

    # ========================================================================
    # Model definitions
    # ========================================================================
    pos_weight = (y == 0).sum() / max(y.sum(), 1)
    print(f"\nClass weight ratio: {pos_weight:.1f}:1")

    models = {
        "ElasticNet": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(
                penalty="elasticnet", solver="saga", l1_ratio=0.5,
                class_weight="balanced", max_iter=5000, random_state=SEED,
                C=1.0
            ))
        ]),
        "RandomForest": RandomForestClassifier(
            n_estimators=500, max_depth=10, min_samples_leaf=5,
            class_weight="balanced", random_state=SEED, n_jobs=-1
        ),
        "LightGBM": lgb.LGBMClassifier(
            n_estimators=500, max_depth=6, learning_rate=0.05,
            scale_pos_weight=pos_weight, random_state=SEED,
            verbose=-1, n_jobs=-1, num_leaves=31,
            min_child_samples=10
        ),
        "XGBoost": xgb.XGBClassifier(
            n_estimators=500, max_depth=6, learning_rate=0.05,
            scale_pos_weight=pos_weight, random_state=SEED,
            eval_metric="logloss", verbosity=0, n_jobs=-1
        ),
    }

    # ========================================================================
    # Cross-validation
    # ========================================================================
    print(f"\nRunning {N_FOLDS}-fold CV x {N_REPEATS} repeats...")

    cv = RepeatedStratifiedKFold(n_splits=N_FOLDS, n_repeats=N_REPEATS,
                                  random_state=SEED)

    # Fill NaN for models that need it (ElasticNet, RandomForest)
    X_filled = X.fillna(0)

    results = {name: {"auroc": [], "auprc": [], "f1": [], "prec50": [], "prec100": []}
               for name in models}

    for fold_i, (train_idx, test_idx) in enumerate(cv.split(X_filled, y)):
        X_train, X_test = X_filled.iloc[train_idx], X_filled.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        for name, model in models.items():
            model_clone = __import__("sklearn.base", fromlist=["clone"]).clone(model)
            model_clone.fit(X_train, y_train)
            y_prob = model_clone.predict_proba(X_test)[:, 1]

            results[name]["auroc"].append(roc_auc_score(y_test, y_prob))
            results[name]["auprc"].append(average_precision_score(y_test, y_prob))

            y_pred = (y_prob >= 0.5).astype(int)
            results[name]["f1"].append(f1_score(y_test, y_pred, zero_division=0))
            results[name]["prec50"].append(precision_at_k(y_test.values, y_prob, 50))
            results[name]["prec100"].append(precision_at_k(y_test.values, y_prob, 100))

        if (fold_i + 1) % N_FOLDS == 0:
            rep = (fold_i + 1) // N_FOLDS
            print(f"  Repeat {rep}/{N_REPEATS} done "
                  f"(LightGBM AUROC: {np.mean(results['LightGBM']['auroc'][-N_FOLDS:]):.3f})")

    # ========================================================================
    # Single-source baselines
    # ========================================================================
    print("\nRunning single-source baselines...")

    baselines = {}

    # Baseline 1: L1 only (dream padj ranking)
    dream_score = -np.log10(atlas["dream_padj"].clip(lower=1e-300).astype(float))
    dream_score = dream_score.fillna(0)
    baselines["L1_dream_only"] = {
        "auroc": roc_auc_score(y, dream_score),
        "auprc": average_precision_score(y, dream_score),
    }

    # Baseline 2: Genetic causal only (n_coloc_sources)
    genetic_score = atlas["n_coloc_sources"].fillna(0).astype(float)
    baselines["L4_genetic_only"] = {
        "auroc": roc_auc_score(y, genetic_score),
        "auprc": average_precision_score(y, genetic_score),
    }

    # Baseline 3: layers_active count
    layers_score = atlas["layers_active"].fillna(0).astype(float)
    baselines["layers_active_count"] = {
        "auroc": roc_auc_score(y, layers_score),
        "auprc": average_precision_score(y, layers_score),
    }

    # Baseline 4: sources_active count
    sources_score = atlas["sources_active"].fillna(0).astype(float)
    baselines["sources_active_count"] = {
        "auroc": roc_auc_score(y, sources_score),
        "auprc": average_precision_score(y, sources_score),
    }

    print("\nBaseline results:")
    for bname, bvals in baselines.items():
        print(f"  {bname}: AUROC={bvals['auroc']:.3f}, AUPRC={bvals['auprc']:.3f}")

    # ========================================================================
    # Model comparison summary
    # ========================================================================
    print("\n" + "=" * 70)
    print("Model Comparison (mean ± std across CV folds)")
    print("=" * 70)

    comparison_rows = []
    for name in models:
        row = {
            "model": name,
            "auroc_mean": np.mean(results[name]["auroc"]),
            "auroc_std": np.std(results[name]["auroc"]),
            "auprc_mean": np.mean(results[name]["auprc"]),
            "auprc_std": np.std(results[name]["auprc"]),
            "f1_mean": np.mean(results[name]["f1"]),
            "f1_std": np.std(results[name]["f1"]),
            "prec50_mean": np.mean(results[name]["prec50"]),
            "prec50_std": np.std(results[name]["prec50"]),
            "prec100_mean": np.mean(results[name]["prec100"]),
            "prec100_std": np.std(results[name]["prec100"]),
        }
        comparison_rows.append(row)
        print(f"  {name:15s}: AUROC={row['auroc_mean']:.3f}±{row['auroc_std']:.3f}  "
              f"AUPRC={row['auprc_mean']:.3f}±{row['auprc_std']:.3f}  "
              f"F1={row['f1_mean']:.3f}  P@50={row['prec50_mean']:.3f}  "
              f"P@100={row['prec100_mean']:.3f}")

    # Add baselines to comparison
    for bname, bvals in baselines.items():
        comparison_rows.append({
            "model": bname,
            "auroc_mean": bvals["auroc"],
            "auroc_std": 0.0,
            "auprc_mean": bvals["auprc"],
            "auprc_std": 0.0,
            "f1_mean": np.nan,
            "f1_std": np.nan,
            "prec50_mean": np.nan,
            "prec50_std": np.nan,
            "prec100_mean": np.nan,
            "prec100_std": np.nan,
        })

    comparison_df = pd.DataFrame(comparison_rows)
    comparison_df.to_csv(f"{OUT_DIR}/ml_model_comparison.csv", index=False)
    print(f"\nSaved: ml_model_comparison.csv")

    # ========================================================================
    # Select best model and train on full data
    # ========================================================================
    # WARNING: This model is retrained on ALL data INCLUDING positive controls.
    # Its predictions are useful for RANKING genes (e.g., prioritizing targets)
    # but must NOT be used to report prediction metrics (AUROC, AUPRC, etc.)
    # in the paper. Only the CV metrics computed above are honest/unbiased.
    # ========================================================================
    best_model_name = max(models.keys(),
                          key=lambda n: np.mean(results[n]["auroc"]))
    print(f"\nBest model: {best_model_name} "
          f"(AUROC={np.mean(results[best_model_name]['auroc']):.3f})")

    best_model = __import__("sklearn.base", fromlist=["clone"]).clone(
        models[best_model_name])
    best_model.fit(X_filled, y)

    # Predictions
    y_prob_full = best_model.predict_proba(X_filled)[:, 1]
    predictions = atlas[["human_symbol", "ensembl_id"]].copy()
    predictions["ml_score"] = y_prob_full
    predictions["ml_rank"] = predictions["ml_score"].rank(ascending=False).astype(int)
    predictions["is_positive_control"] = y.values
    predictions = predictions.sort_values("ml_rank")
    predictions.to_csv(f"{OUT_DIR}/ml_target_predictions.csv", index=False)
    print(f"Saved: ml_target_predictions.csv")

    # Top 20
    print("\nTop 20 predicted targets:")
    for _, row in predictions.head(20).iterrows():
        flag = " *PC*" if row["is_positive_control"] else ""
        print(f"  {row['ml_rank']:4d}. {row['human_symbol']:15s} "
              f"score={row['ml_score']:.4f}{flag}")

    # ========================================================================
    # Source-group ablation
    # ========================================================================
    print("\n" + "=" * 70)
    print("Source-Group Ablation")
    print("=" * 70)

    unique_sources = list(dict.fromkeys(source_groups))
    ablation_rows = []

    # Full model AUROC (from CV)
    full_auroc = np.mean(results[best_model_name]["auroc"])

    for source in unique_sources:
        # Columns to drop for this source
        drop_cols = [feature_names[i] for i in range(len(feature_names))
                     if source_groups[i] == source]

        # Build ablated feature matrix
        X_ablated = X_filled.drop(columns=drop_cols)

        # CV on ablated features
        ablated_aurocs = []
        for train_idx, test_idx in RepeatedStratifiedKFold(
                n_splits=N_FOLDS, n_repeats=3, random_state=SEED).split(X_ablated, y):
            model_clone = __import__("sklearn.base", fromlist=["clone"]).clone(
                models[best_model_name])
            model_clone.fit(X_ablated.iloc[train_idx], y.iloc[train_idx])
            y_prob = model_clone.predict_proba(X_ablated.iloc[test_idx])[:, 1]
            ablated_aurocs.append(roc_auc_score(y.iloc[test_idx], y_prob))

        ablated_mean = np.mean(ablated_aurocs)
        delta = full_auroc - ablated_mean

        ablation_rows.append({
            "source_group": source,
            "n_features": len(drop_cols),
            "auroc_without": ablated_mean,
            "auroc_with": full_auroc,
            "auroc_delta": delta,
        })
        print(f"  Without {source:25s} ({len(drop_cols):2d} feats): "
              f"AUROC={ablated_mean:.3f}  delta={delta:+.3f}")

    ablation_df = pd.DataFrame(ablation_rows)
    ablation_df.to_csv(f"{OUT_DIR}/ml_ablation.csv", index=False)
    print(f"Saved: ml_ablation.csv")

    # ========================================================================
    # SHAP analysis
    # ========================================================================
    # NOTE: SHAP values are computed on the full-data model (trained on all
    # samples including positive controls). They are valid for interpreting
    # which features drive model predictions, but NOT for validating
    # predictive performance. Use CV metrics for prediction validation.
    # ========================================================================
    print("\n" + "=" * 70)
    print("SHAP Feature Importance")
    print("=" * 70)

    # Use TreeExplainer for tree-based models
    if best_model_name in ("LightGBM", "XGBoost", "RandomForest"):
        explainer = shap.TreeExplainer(best_model)
        shap_values = explainer.shap_values(X_filled)
        # For binary classification, shap_values may be a list [neg, pos]
        # or a 3D array (n_samples, n_features, n_classes)
        if isinstance(shap_values, list):
            shap_values = shap_values[1]  # positive class
        elif shap_values.ndim == 3:
            shap_values = shap_values[:, :, 1]  # positive class
    else:
        # For logistic regression, use a sample-based explainer
        explainer = shap.LinearExplainer(best_model.named_steps["clf"],
                                          best_model.named_steps["scaler"].transform(X_filled))
        shap_values = explainer.shap_values(
            best_model.named_steps["scaler"].transform(X_filled))

    # Mean absolute SHAP per feature — ensure 1D
    shap_arr = np.abs(shap_values)
    if shap_arr.ndim > 1:
        mean_shap = shap_arr.mean(axis=0)
    else:
        mean_shap = shap_arr
    # Flatten in case still multi-dimensional
    mean_shap = np.asarray(mean_shap).flatten()
    assert len(mean_shap) == len(feature_names), \
        f"SHAP shape mismatch: {len(mean_shap)} vs {len(feature_names)}"

    importance_rows = []
    for i, fname in enumerate(feature_names):
        importance_rows.append({
            "feature": fname,
            "source_group": source_groups[i],
            "mean_abs_shap": mean_shap[i],
        })

    importance_df = pd.DataFrame(importance_rows).sort_values(
        "mean_abs_shap", ascending=False)
    importance_df.to_csv(f"{OUT_DIR}/ml_feature_importance.csv", index=False)
    print(f"Saved: ml_feature_importance.csv")

    print("\nTop 15 features by SHAP:")
    for _, row in importance_df.head(15).iterrows():
        print(f"  {row['feature']:35s} [{row['source_group']:25s}] "
              f"SHAP={row['mean_abs_shap']:.4f}")

    # Aggregate SHAP by source group
    print("\nSource-group aggregate SHAP:")
    source_shap = importance_df.groupby("source_group")["mean_abs_shap"].sum()
    source_shap = source_shap.sort_values(ascending=False)
    for sg, val in source_shap.items():
        print(f"  {sg:30s}: {val:.4f}")

    # ========================================================================
    # Final summary
    # ========================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Best model: {best_model_name}")
    print(f"  CV AUROC: {np.mean(results[best_model_name]['auroc']):.3f} "
          f"± {np.std(results[best_model_name]['auroc']):.3f}")
    print(f"  CV AUPRC: {np.mean(results[best_model_name]['auprc']):.3f} "
          f"± {np.std(results[best_model_name]['auprc']):.3f}")
    print(f"Beats L1-only baseline: "
          f"{'YES' if full_auroc > baselines['L1_dream_only']['auroc'] else 'NO'} "
          f"({full_auroc:.3f} vs {baselines['L1_dream_only']['auroc']:.3f})")
    print(f"Beats layers_active: "
          f"{'YES' if full_auroc > baselines['layers_active_count']['auroc'] else 'NO'} "
          f"({full_auroc:.3f} vs {baselines['layers_active_count']['auroc']:.3f})")
    # WARNING: These top-K counts are from the full-data model (trained on
    # all samples including positives) and are OPTIMISTICALLY BIASED.
    # They must NOT be reported as validation metrics in the paper.
    # Only CV AUROC/AUPRC (above) are honest performance estimates.
    print(f"\nPositive controls in top 50 (BIASED — full-data model): "
          f"{predictions.head(50)['is_positive_control'].sum()}")
    print(f"Positive controls in top 100 (BIASED — full-data model): "
          f"{predictions.head(100)['is_positive_control'].sum()}")
    print(f"Positive controls in top 200 (BIASED — full-data model): "
          f"{predictions.head(200)['is_positive_control'].sum()}")

    print(f"\nOutputs in: {OUT_DIR}/")
    print("Done.")


if __name__ == "__main__":
    main()
