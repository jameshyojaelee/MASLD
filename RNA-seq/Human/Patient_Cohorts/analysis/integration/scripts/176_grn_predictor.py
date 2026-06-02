#!/usr/bin/env python3
"""
Script 176: Cell-Type-Resolved GRN Predictor
=============================================
Uses TF activity x cell-type proportion interactions to identify which
regulatory programs drive fibrosis in each patient.

Approach (simple, robust at N=718):
  1. Load TF activity (296 dims) + cell-type proportions (13 cell types)
  2. Identify top 10 TFs by univariate AUROC with fib_ge3
  3. Create interaction terms: top10 TFs x 4 key cell types = 40 features
  4. Three feature sets compared via LOCO-CV elastic net:
       A. TF-only:      296 TFs -> PCA 50
       B. CellType-only: 13 cell-type proportions
       C. GRN (full):   13 proportions + 50 TF-PCs + 40 interactions
  5. Report which TF x cell-type interactions are most predictive

Key question: Does TF x cell-type interaction improve over either alone?

Inputs:
  - results/staging_classifier/tf_activity_features.csv (1,444 x 296)
  - results/staging_classifier/modeling_metadata.csv    (sample ordering)
  - results/progression/cibersortx_celltype_expression/bayesprism_proportions.csv
  - results/prognosis_v2/labels_v2.csv                 (targets + LOCO folds)

Outputs to results/novel_ml/grn_predictor/:
  - grn_loco_results.csv         -- AUROC per fold per feature set
  - grn_top_interactions.csv     -- top 20 TF x cell-type interactions by |coef|
  - grn_comparison.csv           -- summary: mean AUROC, std across feature sets
  - grn_feature_importance.csv   -- all feature importances from GRN model
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss

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
STAGING = INTEGRATION / "results/staging_classifier"
PROGRESSION = INTEGRATION / "results/progression/cibersortx_celltype_expression"
PROGNOSIS = INTEGRATION / "results/prognosis_v2"
OUTDIR = INTEGRATION / "results/novel_ml/grn_predictor"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Input files
TF_FILE = STAGING / "tf_activity_features.csv"
META_FILE = STAGING / "modeling_metadata.csv"
PROP_FILE = PROGRESSION / "bayesprism_proportions.csv"
LABEL_FILE = PROGNOSIS / "labels_v2.csv"

# Parameters
RANDOM_STATE = 42
MAX_ITER = 5000
N_TF_PCA = 50       # PCA components for TF features
N_TOP_TFS = 10       # Top TFs for interaction terms
KEY_CELL_TYPES = ["Hepatocyte", "Stellate", "Macrophage", "Endothelial"]
L1_RATIOS = [0.1, 0.5, 0.9]
N_INNER_FOLDS = 5
TARGETS = ["fib_ge3", "loco_s2"]


def log(msg):
    print(f"[176] {msg}", flush=True)


def compute_metrics(y_true, y_prob):
    """Compute AUROC, AUPRC, Brier."""
    metrics = {}
    try:
        metrics["auroc"] = roc_auc_score(y_true, y_prob)
    except ValueError:
        metrics["auroc"] = np.nan
    try:
        metrics["auprc"] = average_precision_score(y_true, y_prob)
    except ValueError:
        metrics["auprc"] = np.nan
    try:
        metrics["brier"] = brier_score_loss(y_true, y_prob)
    except ValueError:
        metrics["brier"] = np.nan
    return metrics


def fit_elastic_net(X_train, y_train, X_test, feature_names):
    """Fit LogisticRegressionCV (elastic net), return predictions + coefficients."""
    scaler = StandardScaler().fit(X_train)
    X_tr = scaler.transform(X_train)
    X_te = scaler.transform(X_test)

    model = LogisticRegressionCV(
        penalty="elasticnet",
        solver="saga",
        l1_ratios=L1_RATIOS,
        Cs=10,
        cv=min(N_INNER_FOLDS, min(y_train.sum(), (1 - y_train).sum())),
        scoring="roc_auc",
        max_iter=MAX_ITER,
        random_state=RANDOM_STATE,
        n_jobs=-1,
        class_weight="balanced",
    )
    model.fit(X_tr, y_train)
    y_prob = model.predict_proba(X_te)[:, 1]
    coefs = dict(zip(feature_names, model.coef_[0]))
    return y_prob, coefs, model


def main():
    t0 = time.time()
    log("Loading data...")

    # ---- 1. Load TF activity (last column is sample_id) ----
    tf_raw = pd.read_csv(TF_FILE)
    if "sample_id" in tf_raw.columns:
        tf_raw = tf_raw.set_index("sample_id")
    else:
        # Fallback: align with modeling_metadata row order
        meta = pd.read_csv(META_FILE)
        tf_raw.index = meta["sample_id"].values

    sample_ids = tf_raw.index.values
    tf_df = tf_raw.copy()
    tf_cols = list(tf_df.columns)
    log(f"  TF activity: {tf_df.shape[0]} samples x {tf_df.shape[1]} TFs")

    # ---- 2. Load cell-type proportions ----
    prop_df = pd.read_csv(PROP_FILE)
    prop_df = prop_df.set_index("sample_id")
    # Align to same sample order
    prop_df = prop_df.reindex(sample_ids).fillna(0.0)
    ct_cols = list(prop_df.columns)
    log(f"  Cell-type proportions: {prop_df.shape[1]} cell types")

    # ---- 3. Load labels + LOCO folds ----
    labels = pd.read_csv(LABEL_FILE)
    labels = labels.set_index("sample_id")
    labels = labels.reindex(sample_ids)

    # Filter to non-excluded samples with valid fib_ge3
    valid_mask = (labels["loco_fold_fibrosis"] != "excluded") & labels["fib_ge3"].notna()
    valid_ids = sample_ids[valid_mask.values]
    log(f"  Valid samples for LOCO: {len(valid_ids)}")

    labels_valid = labels.loc[valid_ids]
    tf_valid = tf_df.loc[valid_ids]
    prop_valid = prop_df.loc[valid_ids]

    # ---- 4. Identify top TFs by univariate AUROC (on full dataset for feature engineering) ----
    y_fib = labels_valid["fib_ge3"].values.astype(int)
    tf_aurocs = {}
    for col in tf_cols:
        vals = tf_valid[col].values
        if np.std(vals) == 0:
            tf_aurocs[col] = 0.5
            continue
        try:
            auc = roc_auc_score(y_fib, vals)
            tf_aurocs[col] = max(auc, 1 - auc)
        except ValueError:
            tf_aurocs[col] = 0.5

    tf_auroc_series = pd.Series(tf_aurocs).sort_values(ascending=False)
    top_tfs = tf_auroc_series.head(N_TOP_TFS).index.tolist()
    log(f"  Top {N_TOP_TFS} TFs (by univariate AUROC): {top_tfs}")
    for tf_name in top_tfs:
        log(f"    {tf_name}: AUROC={tf_auroc_series[tf_name]:.3f}")

    # ---- 5. Build feature matrices ----
    # A. TF-only: PCA to N_TF_PCA dims
    # B. CellType-only: 13 proportions
    # C. GRN: proportions + TF-PCA + interactions (top TFs x key cell types)

    # Build interaction features (on valid samples)
    interaction_cols = []
    interaction_data = {}
    for tf_name in top_tfs:
        for ct in KEY_CELL_TYPES:
            col_name = f"int_{tf_name}_x_{ct}"
            interaction_cols.append(col_name)
            interaction_data[col_name] = (
                tf_valid[tf_name].values * prop_valid[ct].values
            )
    interaction_df = pd.DataFrame(interaction_data, index=valid_ids)
    log(f"  Interaction features: {len(interaction_cols)}")

    # ---- 6. LOCO-CV ----
    folds = labels_valid["loco_fold_fibrosis"].values
    unique_folds = sorted(set(folds))
    log(f"  LOCO folds: {unique_folds}")

    all_results = []
    all_importances = []

    for target_name in TARGETS:
        y_col = labels_valid[target_name].values
        # Skip if target has missing values
        target_valid = ~np.isnan(y_col)
        if target_valid.sum() < 50:
            log(f"  Skipping target {target_name}: only {target_valid.sum()} valid samples")
            continue

        log(f"\n=== Target: {target_name} (n={target_valid.sum()}, "
            f"pos={int(y_col[target_valid].sum())}) ===")

        for fold_name in unique_folds:
            test_mask = (folds == fold_name) & target_valid
            train_mask = (folds != fold_name) & target_valid

            if test_mask.sum() < 5 or train_mask.sum() < 20:
                log(f"  Fold {fold_name}: too few samples (train={train_mask.sum()}, "
                    f"test={test_mask.sum()}), skipping")
                continue

            y_train = y_col[train_mask].astype(int)
            y_test = y_col[test_mask].astype(int)

            # Skip if test set has only one class
            if len(np.unique(y_test)) < 2:
                log(f"  Fold {fold_name}: single class in test set, skipping")
                continue

            # Skip if train set has fewer than 2 of minority class
            if min(y_train.sum(), (1 - y_train).sum()) < 2:
                log(f"  Fold {fold_name}: <2 minority in train, skipping")
                continue

            # --- Feature set A: TF-only (PCA) ---
            tf_train = tf_valid.values[train_mask]
            tf_test = tf_valid.values[test_mask]

            scaler_tf = StandardScaler().fit(tf_train)
            tf_train_s = scaler_tf.transform(tf_train)
            tf_test_s = scaler_tf.transform(tf_test)

            n_components = min(N_TF_PCA, tf_train_s.shape[0], tf_train_s.shape[1])
            pca_tf = PCA(n_components=n_components, random_state=RANDOM_STATE)
            tf_train_pca = pca_tf.fit_transform(tf_train_s)
            tf_test_pca = pca_tf.transform(tf_test_s)

            tf_pca_cols = [f"tf_pc{i}" for i in range(n_components)]

            try:
                y_prob_a, coefs_a, _ = fit_elastic_net(
                    tf_train_pca, y_train, tf_test_pca, tf_pca_cols
                )
                m_a = compute_metrics(y_test, y_prob_a)
                log(f"  Fold {fold_name} TF-only:     AUROC={m_a['auroc']:.3f}")
            except Exception as e:
                log(f"  Fold {fold_name} TF-only FAILED: {e}")
                m_a = {"auroc": np.nan, "auprc": np.nan, "brier": np.nan}

            # --- Feature set B: CellType-only ---
            ct_train = prop_valid.values[train_mask]
            ct_test = prop_valid.values[test_mask]
            try:
                y_prob_b, coefs_b, _ = fit_elastic_net(
                    ct_train, y_train, ct_test, ct_cols
                )
                m_b = compute_metrics(y_test, y_prob_b)
                log(f"  Fold {fold_name} CellType:    AUROC={m_b['auroc']:.3f}")
            except Exception as e:
                log(f"  Fold {fold_name} CellType FAILED: {e}")
                m_b = {"auroc": np.nan, "auprc": np.nan, "brier": np.nan}

            # --- Feature set C: GRN (full) ---
            # Combine: proportions + TF-PCA + interactions
            int_train = interaction_df.values[train_mask]
            int_test = interaction_df.values[test_mask]

            grn_train = np.hstack([ct_train, tf_train_pca, int_train])
            grn_test = np.hstack([ct_test, tf_test_pca, int_test])
            grn_cols = ct_cols + tf_pca_cols + interaction_cols

            try:
                y_prob_c, coefs_c, _ = fit_elastic_net(
                    grn_train, y_train, grn_test, grn_cols
                )
                m_c = compute_metrics(y_test, y_prob_c)
                log(f"  Fold {fold_name} GRN (full):  AUROC={m_c['auroc']:.3f}")

                # Store feature importances for GRN model
                for feat, coef in coefs_c.items():
                    all_importances.append({
                        "target": target_name,
                        "fold": fold_name,
                        "feature": feat,
                        "coefficient": coef,
                        "abs_coefficient": abs(coef),
                    })
            except Exception as e:
                log(f"  Fold {fold_name} GRN FAILED: {e}")
                m_c = {"auroc": np.nan, "auprc": np.nan, "brier": np.nan}

            # Store results
            for config, metrics in [
                ("A_TF_only", m_a),
                ("B_CellType_only", m_b),
                ("C_GRN_full", m_c),
            ]:
                all_results.append({
                    "target": target_name,
                    "fold": fold_name,
                    "config": config,
                    "auroc": metrics["auroc"],
                    "auprc": metrics["auprc"],
                    "brier": metrics["brier"],
                    "n_train": int(train_mask.sum()),
                    "n_test": int(test_mask.sum()),
                    "train_pos_rate": float(y_train.mean()),
                    "test_pos_rate": float(y_test.mean()),
                })

    # ---- 7. Save results ----
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUTDIR / "grn_loco_results.csv", index=False)
    log(f"\nSaved per-fold results: {len(results_df)} rows")

    # Summary comparison
    summary_rows = []
    for target_name in TARGETS:
        for config in ["A_TF_only", "B_CellType_only", "C_GRN_full"]:
            subset = results_df[
                (results_df["target"] == target_name) & (results_df["config"] == config)
            ]
            if len(subset) == 0:
                continue
            summary_rows.append({
                "target": target_name,
                "config": config,
                "n_folds": len(subset),
                "mean_auroc": subset["auroc"].mean(),
                "std_auroc": subset["auroc"].std(),
                "mean_auprc": subset["auprc"].mean(),
                "std_auprc": subset["auprc"].std(),
                "mean_brier": subset["brier"].mean(),
                "std_brier": subset["brier"].std(),
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUTDIR / "grn_comparison.csv", index=False)
    log("\n=== GRN Comparison Summary ===")
    for _, row in summary_df.iterrows():
        log(f"  {row['target']:12s} | {row['config']:20s} | "
            f"AUROC={row['mean_auroc']:.3f} +/- {row['std_auroc']:.3f} | "
            f"AUPRC={row['mean_auprc']:.3f}")

    # Feature importance (GRN model)
    if all_importances:
        imp_df = pd.DataFrame(all_importances)
        imp_df.to_csv(OUTDIR / "grn_feature_importance.csv", index=False)

        # Top interactions
        interaction_imp = imp_df[imp_df["feature"].str.startswith("int_")]
        if len(interaction_imp) > 0:
            mean_imp = (
                interaction_imp.groupby(["target", "feature"])["abs_coefficient"]
                .mean()
                .reset_index()
                .sort_values("abs_coefficient", ascending=False)
            )
            top_int = mean_imp.head(20)
            top_int.to_csv(OUTDIR / "grn_top_interactions.csv", index=False)
            log("\n=== Top TF x Cell-Type Interactions (fib_ge3) ===")
            fib_top = mean_imp[mean_imp["target"] == "fib_ge3"].head(10)
            for _, row in fib_top.iterrows():
                log(f"  {row['feature']:40s} | mean |coef| = {row['abs_coefficient']:.4f}")

            # Also report: do interactions add value?
            log("\n=== Key Question: Do interactions improve over components alone? ===")
            for target_name in TARGETS:
                s = summary_df[summary_df["target"] == target_name]
                if len(s) < 3:
                    continue
                tf_auroc = s[s["config"] == "A_TF_only"]["mean_auroc"].values[0]
                ct_auroc = s[s["config"] == "B_CellType_only"]["mean_auroc"].values[0]
                grn_auroc = s[s["config"] == "C_GRN_full"]["mean_auroc"].values[0]
                best_single = max(tf_auroc, ct_auroc)
                delta = grn_auroc - best_single
                verdict = "YES" if delta > 0.01 else ("MARGINAL" if delta > 0 else "NO")
                log(f"  {target_name}: TF={tf_auroc:.3f}, CT={ct_auroc:.3f}, "
                    f"GRN={grn_auroc:.3f} -> delta={delta:+.3f} ({verdict})")
    else:
        log("WARNING: No feature importances collected")

    elapsed = time.time() - t0
    log(f"\nDone in {elapsed:.0f}s ({elapsed / 60:.1f} min)")


if __name__ == "__main__":
    main()
