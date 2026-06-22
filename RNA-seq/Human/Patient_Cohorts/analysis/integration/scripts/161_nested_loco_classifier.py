#!/usr/bin/env python3
"""
Script 161: Nested LOCO Classifier
===================================
Properly nested cross-validation for MASLD progression prediction.

Outer loop: 6 LOCO cohort folds (leave-one-cohort-out)
Inner loop: 5-fold stratified CV for feature selection + hyperparameter tuning

Targets:
  1. fib_ge3    — binary fibrosis >= F3 (primary)
  2. loco_s2    — binary LOCO-NMF subtype (S2 progressor)

Feature configs:
  M1_clinical    — sex, age, fibrosis_stage (2-3 features)
  M2_celltype    — bp_*/ct_* cell-type fractions (17 features)
  M3_expression  — div_* divergence genes (FULL candidate pool, top-N chosen
                   fold-internally by inner-CV univariate-AUROC; A7 leakage fix
                   2026-06-13 — Script 150 no longer bakes a full-cohort
                   |cohens_d| top-100 screen into the static matrix)
  M4_combined    — M2 + inner-CV-selected expression genes

Baselines:
  - Permutation null (100 permutations × M2_celltype)
  - Random gene baseline (100 draws × 25 random genes)

Outputs to results/prognosis_v2/:
  - nested_loco_results.csv
  - nested_loco_predictions.csv
  - nested_loco_feature_importance.csv
  - nested_loco_baselines.csv
  - nested_loco_summary.csv
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from joblib import Parallel, delayed

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
RESULTS = INTEGRATION / "results/prognosis_v2"

FEATURE_FILE = RESULTS / "feature_matrix_v2.csv"
LABEL_FILE = RESULTS / "labels_v2.csv"

N_INNER_FOLDS = 5
RANDOM_STATE = 42
N_PERMUTATIONS = 100
N_RANDOM_GENE_DRAWS = 100
N_RANDOM_GENES = 25
MAX_ITER = 5000

# Hyperparameter grid for expression-based models
N_GENES_GRID = [10, 25, 50, 100]
L1_RATIO_GRID = [0.1, 0.5, 0.9]

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))
print(f"[161] Using {N_JOBS} parallel jobs")


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def compute_univariate_auroc(X, y, feature_cols):
    """Compute AUROC for each feature vs binary target. Returns sorted Series."""
    aurocs = {}
    for col in feature_cols:
        vals = X[col].values
        # Skip constant features
        if np.std(vals) == 0:
            aurocs[col] = 0.5
            continue
        try:
            auc = roc_auc_score(y, vals)
            # Flip if < 0.5 (direction doesn't matter for ranking)
            aurocs[col] = max(auc, 1 - auc)
        except ValueError:
            aurocs[col] = 0.5
    return pd.Series(aurocs).sort_values(ascending=False)


def compute_metrics(y_true, y_prob):
    """Compute AUROC, AUPRC, Brier score."""
    try:
        auroc = roc_auc_score(y_true, y_prob)
    except ValueError:
        auroc = np.nan
    try:
        auprc = average_precision_score(y_true, y_prob)
    except ValueError:
        auprc = np.nan
    try:
        brier = brier_score_loss(y_true, y_prob)
    except ValueError:
        brier = np.nan
    return auroc, auprc, brier


def fit_fixed_features(X_train, y_train, X_test, feature_cols):
    """Fit LogisticRegressionCV on fixed feature set (M1/M2). Returns predictions and model."""
    scaler = StandardScaler().fit(X_train[feature_cols])
    X_tr = scaler.transform(X_train[feature_cols])
    X_te = scaler.transform(X_test[feature_cols])

    model = LogisticRegressionCV(
        penalty="elasticnet",
        l1_ratios=[0.5],
        solver="saga",
        max_iter=MAX_ITER,
        class_weight="balanced",
        cv=N_INNER_FOLDS,
        scoring="roc_auc",
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )
    model.fit(X_tr, y_train)
    y_prob = model.predict_proba(X_te)[:, 1]

    # Extract coefficients
    coefs = dict(zip(feature_cols, model.coef_[0]))
    return y_prob, coefs, feature_cols, {"l1_ratio": model.l1_ratio_[0], "C": model.C_[0]}


def fit_expression_model_inner_cv(X_train, y_train, X_test, expression_cols,
                                  extra_cols=None):
    """
    Nested inner CV for expression-based models (M3/M4).
    Feature selection and hyperparameter tuning happen ONLY on training data.
    Returns predictions, coefficients, selected features, best params.
    """
    best_inner_auroc = -1
    best_params = (25, 0.5)  # defaults

    inner_cv = StratifiedKFold(
        n_splits=N_INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE
    )

    for n_genes in N_GENES_GRID:
        for l1_ratio in L1_RATIO_GRID:
            inner_aurocs = []

            for train_idx, val_idx in inner_cv.split(X_train, y_train):
                X_inner_train = X_train.iloc[train_idx]
                X_inner_val = X_train.iloc[val_idx]
                y_inner_train = y_train.iloc[train_idx]
                y_inner_val = y_train.iloc[val_idx]

                # Feature selection on inner train ONLY
                gene_aurocs = compute_univariate_auroc(
                    X_inner_train, y_inner_train, expression_cols
                )
                selected_genes = gene_aurocs.head(min(n_genes, len(gene_aurocs))).index.tolist()

                # Build feature set
                if extra_cols is not None:
                    use_cols = list(extra_cols) + selected_genes
                else:
                    use_cols = selected_genes

                # Scale on inner train
                scaler = StandardScaler().fit(X_inner_train[use_cols])
                X_tr = scaler.transform(X_inner_train[use_cols])
                X_val = scaler.transform(X_inner_val[use_cols])

                # Fit
                model = LogisticRegression(
                    penalty="elasticnet",
                    l1_ratio=l1_ratio,
                    solver="saga",
                    max_iter=MAX_ITER,
                    class_weight="balanced",
                    C=1.0,
                    random_state=RANDOM_STATE,
                )
                model.fit(X_tr, y_inner_train)

                try:
                    auc = roc_auc_score(
                        y_inner_val, model.predict_proba(X_val)[:, 1]
                    )
                except ValueError:
                    auc = 0.5
                inner_aurocs.append(auc)

            mean_inner = np.mean(inner_aurocs)
            if mean_inner > best_inner_auroc:
                best_inner_auroc = mean_inner
                best_params = (n_genes, l1_ratio)

    # Retrain on full outer training set with best params
    best_n, best_l1 = best_params
    gene_aurocs = compute_univariate_auroc(X_train, y_train, expression_cols)
    selected_genes = gene_aurocs.head(min(best_n, len(gene_aurocs))).index.tolist()

    if extra_cols is not None:
        use_cols = list(extra_cols) + selected_genes
    else:
        use_cols = selected_genes

    scaler = StandardScaler().fit(X_train[use_cols])
    X_tr = scaler.transform(X_train[use_cols])
    X_te = scaler.transform(X_test[use_cols])

    model = LogisticRegression(
        penalty="elasticnet",
        l1_ratio=best_l1,
        solver="saga",
        max_iter=MAX_ITER,
        class_weight="balanced",
        C=1.0,
        random_state=RANDOM_STATE,
    )
    model.fit(X_tr, y_train)
    y_prob = model.predict_proba(X_te)[:, 1]

    coefs = dict(zip(use_cols, model.coef_[0]))
    params = {
        "n_genes": best_n,
        "l1_ratio": best_l1,
        "inner_auroc": best_inner_auroc,
    }
    return y_prob, coefs, selected_genes, params


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print(f"[161] Loading data from {RESULTS}")

    features = pd.read_csv(FEATURE_FILE)
    labels = pd.read_csv(LABEL_FILE)

    # Merge on sample_id
    df = features.merge(labels, on="sample_id", how="inner")
    print(f"[161] Merged: {len(df)} samples, {df.shape[1]} columns")

    # Identify feature columns by prefix
    div_cols = [c for c in features.columns if c.startswith("div_")]
    ct_cols = [c for c in features.columns if c.startswith("ct_") or c.startswith("bp_")]
    clin_cols_feat = [c for c in features.columns if c.startswith("clin_")]

    print(f"[161] Feature groups: div={len(div_cols)}, ct={len(ct_cols)}, clin={len(clin_cols_feat)}")

    # Build clinical feature set from labels + features
    # Map sex to numeric: F=0, M=1
    df["_clin_sex"] = df["sex"].map({"F": 0, "M": 1})
    df["_clin_age"] = df["age"]
    # fibrosis_stage from feature matrix (clin_fibrosis_stage)
    df["_clin_fib_stage"] = df["clin_fibrosis_stage"] if "clin_fibrosis_stage" in df.columns else np.nan

    clinical_cols = ["_clin_sex", "_clin_age", "_clin_fib_stage"]

    # Filter to samples with valid LOCO folds
    valid_folds = df[
        (df["loco_fold_fibrosis"] != "excluded") & df["loco_fold_fibrosis"].notna()
    ].copy()
    fold_names = sorted(valid_folds["loco_fold_fibrosis"].unique())
    print(f"[161] Valid LOCO samples: {len(valid_folds)}, folds: {fold_names}")

    # Define targets
    targets = {
        "fib_ge3": "fib_ge3",
        "loco_s2": "loco_s2",
    }

    # Define feature configs
    config_specs = {
        "M1_clinical": {"type": "fixed", "cols": clinical_cols},
        "M2_celltype": {"type": "fixed", "cols": ct_cols},
        "M3_expression": {"type": "expression", "expr_cols": div_cols, "extra_cols": None},
        "M4_combined": {"type": "expression", "expr_cols": div_cols, "extra_cols": ct_cols},
    }

    # Storage
    all_results = []
    all_predictions = []
    all_importance = []

    # -----------------------------------------------------------------------
    # Main nested CV loop
    # -----------------------------------------------------------------------
    for target_name, target_col in targets.items():
        print(f"\n{'='*70}")
        print(f"[161] TARGET: {target_name}")
        print(f"{'='*70}")

        for config_name, config in config_specs.items():
            print(f"\n  --- Config: {config_name} ---")
            config_importance = {}  # gene -> list of coefs across folds

            for fold in fold_names:
                # Split outer fold
                test_mask = valid_folds["loco_fold_fibrosis"] == fold
                train_mask = ~test_mask

                train_df = valid_folds[train_mask].copy()
                test_df = valid_folds[test_mask].copy()

                # Filter to valid target
                train_df = train_df[train_df[target_col].notna()].copy()
                test_df = test_df[test_df[target_col].notna()].copy()

                if len(test_df) == 0 or len(train_df) == 0:
                    print(f"    Fold {fold}: SKIP (no valid samples)")
                    continue

                y_train = train_df[target_col].astype(int)
                y_test = test_df[target_col].astype(int)

                # Check that both classes are present
                if len(y_train.unique()) < 2 or len(y_test.unique()) < 2:
                    print(f"    Fold {fold}: SKIP (single class)")
                    continue

                # Handle clinical NaN: median imputation WITHIN training fold
                if config["type"] == "fixed" and config["cols"] == clinical_cols:
                    for col in clinical_cols:
                        median_val = train_df[col].median()
                        train_df[col] = train_df[col].fillna(median_val)
                        test_df[col] = test_df[col].fillna(median_val)
                    # Drop columns that are still all NaN or constant
                    usable_clin = [
                        c for c in clinical_cols
                        if train_df[c].notna().all() and train_df[c].std() > 0
                    ]
                    if len(usable_clin) == 0:
                        print(f"    Fold {fold}: SKIP (no usable clinical features)")
                        continue

                    # NOTE: For fib_ge3 target, exclude fibrosis_stage to avoid circularity
                    if target_name == "fib_ge3" and "_clin_fib_stage" in usable_clin:
                        usable_clin = [c for c in usable_clin if c != "_clin_fib_stage"]
                    if len(usable_clin) == 0:
                        print(f"    Fold {fold}: SKIP (no usable clinical features after circularity check)")
                        continue

                    y_prob, coefs, used_features, params = fit_fixed_features(
                        train_df, y_train, test_df, usable_clin
                    )
                    selected_n = len(usable_clin)
                    best_l1 = params.get("l1_ratio", 0.5)

                elif config["type"] == "fixed":
                    y_prob, coefs, used_features, params = fit_fixed_features(
                        train_df, y_train, test_df, config["cols"]
                    )
                    selected_n = len(config["cols"])
                    best_l1 = params.get("l1_ratio", 0.5)

                elif config["type"] == "expression":
                    y_prob, coefs, selected_genes, params = fit_expression_model_inner_cv(
                        train_df, y_train, test_df,
                        config["expr_cols"],
                        extra_cols=config["extra_cols"],
                    )
                    used_features = selected_genes
                    selected_n = params.get("n_genes", len(selected_genes))
                    best_l1 = params.get("l1_ratio", 0.5)

                # Compute metrics
                auroc, auprc, brier = compute_metrics(y_test.values, y_prob)

                print(
                    f"    Fold {fold}: n_train={len(train_df)}, n_test={len(test_df)}, "
                    f"AUROC={auroc:.3f}, AUPRC={auprc:.3f}, Brier={brier:.3f}"
                )

                # Record results
                all_results.append({
                    "target": target_name,
                    "config": config_name,
                    "fold": fold,
                    "n_train": len(train_df),
                    "n_test": len(test_df),
                    "n_pos_train": int(y_train.sum()),
                    "n_pos_test": int(y_test.sum()),
                    "auroc": auroc,
                    "auprc": auprc,
                    "brier": brier,
                    "selected_n_genes": selected_n,
                    "best_l1_ratio": best_l1,
                    "inner_auroc": params.get("inner_auroc", np.nan)
                        if config["type"] == "expression" else np.nan,
                })

                # Record per-sample predictions
                for sid, prob, true_label in zip(
                    test_df["sample_id"].values, y_prob, y_test.values
                ):
                    all_predictions.append({
                        "sample_id": sid,
                        "target": target_name,
                        "config": config_name,
                        "predicted_prob": prob,
                        "true_label": int(true_label),
                        "fold": fold,
                    })

                # Accumulate feature importance
                for feat, coef in coefs.items():
                    if feat not in config_importance:
                        config_importance[feat] = []
                    config_importance[feat].append(coef)

            # Summarize feature importance for this config × target
            for feat, coef_list in config_importance.items():
                all_importance.append({
                    "gene": feat,
                    "target": target_name,
                    "config": config_name,
                    "n_folds_selected": len(coef_list),
                    "mean_coef": np.mean(coef_list),
                    "mean_abs_coef": np.mean(np.abs(coef_list)),
                    "std_coef": np.std(coef_list) if len(coef_list) > 1 else 0,
                })

    # -----------------------------------------------------------------------
    # Baselines
    # -----------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("[161] Running baselines...")
    print(f"{'='*70}")

    all_baselines = []

    def run_permutation(perm_idx, target_name, target_col):
        """Single permutation: shuffle labels, run M2_celltype across all folds."""
        rng = np.random.RandomState(RANDOM_STATE + perm_idx)
        perm_aurocs = []

        for fold in fold_names:
            test_mask = valid_folds["loco_fold_fibrosis"] == fold
            train_mask = ~test_mask

            train_df = valid_folds[train_mask].copy()
            test_df = valid_folds[test_mask].copy()

            train_df = train_df[train_df[target_col].notna()].copy()
            test_df = test_df[test_df[target_col].notna()].copy()

            if len(test_df) == 0 or len(train_df) == 0:
                continue

            y_train = train_df[target_col].astype(int)
            y_test = test_df[target_col].astype(int)

            if len(y_train.unique()) < 2 or len(y_test.unique()) < 2:
                continue

            # Permute training labels
            y_train_perm = pd.Series(
                rng.permutation(y_train.values), index=y_train.index
            )

            scaler = StandardScaler().fit(train_df[ct_cols])
            X_tr = scaler.transform(train_df[ct_cols])
            X_te = scaler.transform(test_df[ct_cols])

            model = LogisticRegression(
                penalty="elasticnet",
                l1_ratio=0.5,
                solver="saga",
                max_iter=MAX_ITER,
                class_weight="balanced",
                C=1.0,
                random_state=RANDOM_STATE,
            )
            model.fit(X_tr, y_train_perm)
            y_prob = model.predict_proba(X_te)[:, 1]
            try:
                auc = roc_auc_score(y_test, y_prob)
            except ValueError:
                auc = 0.5
            perm_aurocs.append(auc)

        return {
            "baseline_type": "permutation",
            "target": target_name,
            "draw_idx": perm_idx,
            "mean_auroc": np.mean(perm_aurocs) if perm_aurocs else np.nan,
            "fold_aurocs": perm_aurocs,
        }

    def run_random_genes(draw_idx, target_name, target_col, all_gene_cols):
        """Single random gene draw: pick 25 random non-div genes, run across folds."""
        rng = np.random.RandomState(RANDOM_STATE + 10000 + draw_idx)
        random_genes = rng.choice(all_gene_cols, size=min(N_RANDOM_GENES, len(all_gene_cols)), replace=False).tolist()
        draw_aurocs = []

        for fold in fold_names:
            test_mask = valid_folds["loco_fold_fibrosis"] == fold
            train_mask = ~test_mask

            train_df = valid_folds[train_mask].copy()
            test_df = valid_folds[test_mask].copy()

            train_df = train_df[train_df[target_col].notna()].copy()
            test_df = test_df[test_df[target_col].notna()].copy()

            if len(test_df) == 0 or len(train_df) == 0:
                continue

            y_train = train_df[target_col].astype(int)
            y_test = test_df[target_col].astype(int)

            if len(y_train.unique()) < 2 or len(y_test.unique()) < 2:
                continue

            scaler = StandardScaler().fit(train_df[random_genes])
            X_tr = scaler.transform(train_df[random_genes])
            X_te = scaler.transform(test_df[random_genes])

            model = LogisticRegression(
                penalty="elasticnet",
                l1_ratio=0.5,
                solver="saga",
                max_iter=MAX_ITER,
                class_weight="balanced",
                C=1.0,
                random_state=RANDOM_STATE,
            )
            model.fit(X_tr, y_train)
            y_prob = model.predict_proba(X_te)[:, 1]
            try:
                auc = roc_auc_score(y_test, y_prob)
            except ValueError:
                auc = 0.5
            draw_aurocs.append(auc)

        return {
            "baseline_type": "random_genes",
            "target": target_name,
            "draw_idx": draw_idx,
            "mean_auroc": np.mean(draw_aurocs) if draw_aurocs else np.nan,
            "fold_aurocs": draw_aurocs,
        }

    # For random gene baseline, draw from raw expression features ONLY (div_*)
    # NEVER use learned embeddings (nasvae_*, bf_*, trans_*, hcc_*) — they are
    # trained on ALL samples including test folds, causing information leakage
    # (same issue as retracted Script 108 stacking)
    random_pool = div_cols
    print(f"[161] Random gene pool: {len(random_pool)} div_* expression features (no embeddings)")

    for target_name, target_col in targets.items():
        print(f"\n  Permutation null for {target_name}...")
        perm_results = Parallel(n_jobs=N_JOBS, verbose=1)(
            delayed(run_permutation)(i, target_name, target_col)
            for i in range(N_PERMUTATIONS)
        )
        for r in perm_results:
            all_baselines.append({
                "baseline_type": r["baseline_type"],
                "target": r["target"],
                "draw_idx": r["draw_idx"],
                "mean_auroc": r["mean_auroc"],
            })

        print(f"  Random gene baseline for {target_name}...")
        rand_results = Parallel(n_jobs=N_JOBS, verbose=1)(
            delayed(run_random_genes)(i, target_name, target_col, random_pool)
            for i in range(N_RANDOM_GENE_DRAWS)
        )
        for r in rand_results:
            all_baselines.append({
                "baseline_type": r["baseline_type"],
                "target": r["target"],
                "draw_idx": r["draw_idx"],
                "mean_auroc": r["mean_auroc"],
            })

    # -----------------------------------------------------------------------
    # Save outputs
    # -----------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("[161] Saving outputs...")
    print(f"{'='*70}")

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(RESULTS / "nested_loco_results.csv", index=False)
    print(f"  nested_loco_results.csv: {len(results_df)} rows")

    pred_df = pd.DataFrame(all_predictions)
    pred_df.to_csv(RESULTS / "nested_loco_predictions.csv", index=False)
    print(f"  nested_loco_predictions.csv: {len(pred_df)} rows")

    imp_df = pd.DataFrame(all_importance)
    imp_df.to_csv(RESULTS / "nested_loco_feature_importance.csv", index=False)
    print(f"  nested_loco_feature_importance.csv: {len(imp_df)} rows")

    base_df = pd.DataFrame(all_baselines)
    base_df.to_csv(RESULTS / "nested_loco_baselines.csv", index=False)
    print(f"  nested_loco_baselines.csv: {len(base_df)} rows")

    # -----------------------------------------------------------------------
    # Summary table with permutation p-values
    # -----------------------------------------------------------------------
    summary_rows = []
    for target_name in targets:
        for config_name in config_specs:
            sub = results_df[
                (results_df["target"] == target_name)
                & (results_df["config"] == config_name)
            ]
            if len(sub) == 0:
                continue

            mean_auroc = sub["auroc"].mean()
            std_auroc = sub["auroc"].std()
            mean_auprc = sub["auprc"].mean()
            std_auprc = sub["auprc"].std()
            mean_brier = sub["brier"].mean()
            std_brier = sub["brier"].std()

            # Permutation p-value
            perm_sub = base_df[
                (base_df["baseline_type"] == "permutation")
                & (base_df["target"] == target_name)
            ]
            if len(perm_sub) > 0:
                perm_aurocs = perm_sub["mean_auroc"].dropna().values
                p_perm = (np.sum(perm_aurocs >= mean_auroc) + 1) / (len(perm_aurocs) + 1)
            else:
                p_perm = np.nan

            # Random gene p-value (compare to M3/M4)
            rand_sub = base_df[
                (base_df["baseline_type"] == "random_genes")
                & (base_df["target"] == target_name)
            ]
            if len(rand_sub) > 0 and config_name in ("M3_expression", "M4_combined"):
                rand_aurocs = rand_sub["mean_auroc"].dropna().values
                p_rand = (np.sum(rand_aurocs >= mean_auroc) + 1) / (len(rand_aurocs) + 1)
            else:
                p_rand = np.nan

            summary_rows.append({
                "target": target_name,
                "config": config_name,
                "n_folds": len(sub),
                "mean_auroc": mean_auroc,
                "std_auroc": std_auroc,
                "mean_auprc": mean_auprc,
                "std_auprc": std_auprc,
                "mean_brier": mean_brier,
                "std_brier": std_brier,
                "perm_p_value": p_perm,
                "random_gene_p_value": p_rand,
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(RESULTS / "nested_loco_summary.csv", index=False)
    print(f"  nested_loco_summary.csv: {len(summary_df)} rows")

    # Print summary
    print(f"\n{'='*70}")
    print("[161] SUMMARY")
    print(f"{'='*70}")
    for _, row in summary_df.iterrows():
        perm_str = f", perm_p={row['perm_p_value']:.4f}" if not np.isnan(row["perm_p_value"]) else ""
        rand_str = f", rand_p={row['random_gene_p_value']:.4f}" if not np.isnan(row["random_gene_p_value"]) else ""
        print(
            f"  {row['target']:10s} | {row['config']:15s} | "
            f"AUROC={row['mean_auroc']:.3f}±{row['std_auroc']:.3f} | "
            f"AUPRC={row['mean_auprc']:.3f}±{row['std_auprc']:.3f} | "
            f"Brier={row['mean_brier']:.3f}±{row['std_brier']:.3f}"
            f"{perm_str}{rand_str}"
        )

    elapsed = time.time() - t0
    print(f"\n[161] Done in {elapsed/60:.1f} min")


if __name__ == "__main__":
    main()
