#!/usr/bin/env python3
"""
Script 172: Cross-Ancestry Disease Progression Model
=====================================================
Classifier using only genes validated across European AND East Asian GWAS
via colocalization (COLOC PP.H4).

Hypothesis: Cross-ancestry validated genes generalize better across cohorts
(lower LOCO variance) even if they have slightly lower mean AUROC, because
shared causal architecture implies robust biology rather than population-
specific LD artifacts.

Three feature sets compared (all LOCO-CV, elastic net):
  1. All-feature baseline:  100 div_* divergence genes from Script 161
  2. Cross-ancestry-only:   Genes with PP.H4 > 0.3 in >= 1 European AND >= 1 East Asian GWAS
  3. European-only COLOC:   Genes with PP.H4 > 0.3 in >= 1 European GWAS but NOT East Asian

Target: fib_ge3 (binary fibrosis >= F3)

Outputs to results/novel_ml/cross_ancestry/:
  - cross_ancestry_genes.csv         — genes validated in both EUR and EAS with PP.H4 values
  - cross_ancestry_classifier_results.csv — LOCO-CV per-fold AUROC for 3 feature sets
  - cross_ancestry_comparison.csv    — summary: mean AUROC, std, CV across feature sets
  - cross_ancestry_feature_importance.csv — per-gene importance across folds
"""

import os
import sys
import time
import gzip
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
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
PROGNOSIS = INTEGRATION / "results/prognosis_v2"
STAGING = INTEGRATION / "results/staging_classifier"
CAUSAL = BASE / "RNA-seq/results/causal_inference"
OUTDIR = INTEGRATION / "results/novel_ml/cross_ancestry"
OUTDIR.mkdir(parents=True, exist_ok=True)

FEATURE_FILE = PROGNOSIS / "feature_matrix_v2.csv"
LABEL_FILE = PROGNOSIS / "labels_v2.csv"
LOGCPM_CACHE = STAGING / "_wenda_logcpm_cache.csv.gz"

# Cross-ancestry comparison (pre-computed by Script 48)
CROSS_ANCESTRY_FILE = CAUSAL / "cross_ancestry/cross_ancestry_comparison.csv"

# Individual COLOC files for building gene symbol -> Ensembl mapping
COLOC_FILES = [
    CAUSAL / "finngen_nafld/coloc_results.csv",
    CAUSAL / "finngen_nash/coloc_results.csv",
    CAUSAL / "finngen_hcc/coloc_results.csv",
    CAUSAL / "bbj_alt/coloc_results.csv",
    CAUSAL / "bbj_ast/coloc_results.csv",
    CAUSAL / "bbj_ggt/coloc_results.csv",
]

PP4_THRESHOLD = 0.3  # Minimum PP.H4 for colocalization evidence
N_INNER_FOLDS = 5
RANDOM_STATE = 42
MAX_ITER = 5000
N_GENES_GRID = [10, 25, 50, 100]
L1_RATIO_GRID = [0.1, 0.5, 0.9]

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
print(f"[172] Using {N_JOBS} parallel jobs")


# ---------------------------------------------------------------------------
# Helper functions (matching Script 161 framework)
# ---------------------------------------------------------------------------
def compute_univariate_auroc(X, y, feature_cols):
    """Compute AUROC for each feature vs binary target. Returns sorted Series."""
    aurocs = {}
    for col in feature_cols:
        vals = X[col].values
        if np.std(vals) == 0:
            aurocs[col] = 0.5
            continue
        try:
            auc = roc_auc_score(y, vals)
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


def fit_expression_model_inner_cv(X_train, y_train, X_test, expression_cols):
    """
    Nested inner CV for expression-based models.
    Feature selection and hyperparameter tuning happen ONLY on training data.
    Returns predictions, coefficients, selected features, best params.
    """
    n_available = len(expression_cols)

    # If very few features, skip grid search for n_genes
    gene_grid = [n for n in N_GENES_GRID if n <= n_available]
    if not gene_grid:
        gene_grid = [n_available]
    elif n_available not in gene_grid:
        gene_grid.append(n_available)

    best_inner_auroc = -1
    best_params = (min(gene_grid), 0.5)

    inner_cv = StratifiedKFold(
        n_splits=N_INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE
    )

    for n_genes in gene_grid:
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
                selected_genes = gene_aurocs.head(
                    min(n_genes, len(gene_aurocs))
                ).index.tolist()

                if len(selected_genes) == 0:
                    inner_aurocs.append(0.5)
                    continue

                scaler = StandardScaler().fit(X_inner_train[selected_genes])
                X_tr = scaler.transform(X_inner_train[selected_genes])
                X_val = scaler.transform(X_inner_val[selected_genes])

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
    selected_genes = gene_aurocs.head(
        min(best_n, len(gene_aurocs))
    ).index.tolist()

    if len(selected_genes) == 0:
        # Fallback: return chance-level predictions
        y_prob = np.full(len(X_test), 0.5)
        return y_prob, {}, [], {"n_genes": 0, "l1_ratio": best_l1, "inner_auroc": 0.5}

    scaler = StandardScaler().fit(X_train[selected_genes])
    X_tr = scaler.transform(X_train[selected_genes])
    X_te = scaler.transform(X_test[selected_genes])

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

    coefs = dict(zip(selected_genes, model.coef_[0]))
    params = {
        "n_genes": best_n,
        "l1_ratio": best_l1,
        "inner_auroc": best_inner_auroc,
    }
    return y_prob, coefs, selected_genes, params


# ---------------------------------------------------------------------------
# Step 1: Build gene symbol -> Ensembl mapping from COLOC results
# ---------------------------------------------------------------------------
def build_symbol_ensembl_map():
    """Build symbol -> ensembl mapping from all COLOC result files."""
    mappings = []
    for f in COLOC_FILES:
        if f.exists():
            df = pd.read_csv(f, usecols=["gene", "ensembl"])
            mappings.append(df)
        else:
            print(f"  WARNING: COLOC file not found: {f}")
    if not mappings:
        raise FileNotFoundError("No COLOC result files found")
    mapping = pd.concat(mappings).drop_duplicates(subset="gene")
    return dict(zip(mapping["gene"], mapping["ensembl"]))


# ---------------------------------------------------------------------------
# Step 2: Load cross-ancestry COLOC data and classify genes
# ---------------------------------------------------------------------------
def load_cross_ancestry_genes(symbol_to_ensembl):
    """
    Load cross-ancestry comparison and classify genes into:
    - cross_ancestry: PP.H4 > threshold in EUR AND EAS
    - eur_only: PP.H4 > threshold in EUR but NOT EAS
    """
    ca = pd.read_csv(CROSS_ANCESTRY_FILE)
    print(f"[172] Cross-ancestry comparison: {len(ca)} genes")

    cross_mask = (ca["best_pp4_eur"] > PP4_THRESHOLD) & (ca["best_pp4_eas"] > PP4_THRESHOLD)
    eur_only_mask = (ca["best_pp4_eur"] > PP4_THRESHOLD) & (ca["best_pp4_eas"] <= PP4_THRESHOLD)

    cross_genes = ca[cross_mask].copy()
    eur_only_genes = ca[eur_only_mask].copy()

    print(f"[172] Cross-ancestry validated (PP.H4 > {PP4_THRESHOLD} EUR+EAS): {len(cross_genes)}")
    print(f"[172] EUR-only (PP.H4 > {PP4_THRESHOLD} EUR only): {len(eur_only_genes)}")

    # Add ensembl IDs
    cross_genes["ensembl"] = cross_genes["gene"].map(symbol_to_ensembl)
    eur_only_genes["ensembl"] = eur_only_genes["gene"].map(symbol_to_ensembl)

    return cross_genes, eur_only_genes, ca


# ---------------------------------------------------------------------------
# Step 3: Load expression data for specific genes from logcpm cache
# ---------------------------------------------------------------------------
def load_expression_for_genes(gene_symbols, symbol_to_ensembl):
    """
    Load log-CPM expression for specified genes from the wenda cache.
    Returns DataFrame: samples x genes (symbol names as columns).
    """
    # Build ensembl -> symbol map (including versioned IDs)
    target_ensembl = {}
    for sym in gene_symbols:
        ens = symbol_to_ensembl.get(sym)
        if ens:
            target_ensembl[ens] = sym
            target_ensembl[ens.split(".")[0]] = sym

    print(f"[172] Loading expression for {len(gene_symbols)} genes "
          f"({len(target_ensembl)} Ensembl IDs to match)...")

    # Read the compressed logcpm cache, extracting only target genes
    rows = {}
    with gzip.open(str(LOGCPM_CACHE), "rt") as fh:
        header = fh.readline().strip().split(",")
        sample_ids = header[1:]  # First column is gene ID

        for line in fh:
            parts = line.strip().split(",", 1)
            gene_id = parts[0]
            gene_base = gene_id.split(".")[0]

            sym = target_ensembl.get(gene_id) or target_ensembl.get(gene_base)
            if sym and sym not in rows:
                vals = [float(x) for x in parts[1].split(",")]
                rows[sym] = vals

    print(f"[172] Found expression for {len(rows)}/{len(gene_symbols)} genes")

    if len(rows) == 0:
        return pd.DataFrame()

    expr_df = pd.DataFrame(rows, index=sample_ids)
    expr_df.index.name = "sample_id"
    expr_df = expr_df.reset_index()

    return expr_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    print(f"[172] Cross-Ancestry Disease Progression Model")
    print(f"[172] Output: {OUTDIR}")
    print(f"{'='*70}\n")

    # -----------------------------------------------------------------------
    # Step 1: Build gene mapping
    # -----------------------------------------------------------------------
    print("[172] Step 1: Building gene symbol -> Ensembl mapping...")
    symbol_to_ensembl = build_symbol_ensembl_map()
    print(f"[172] Mapped {len(symbol_to_ensembl)} gene symbols to Ensembl IDs")

    # -----------------------------------------------------------------------
    # Step 2: Identify cross-ancestry and EUR-only genes
    # -----------------------------------------------------------------------
    print("\n[172] Step 2: Classifying COLOC genes by ancestry...")
    cross_genes_df, eur_only_genes_df, ca_full = load_cross_ancestry_genes(symbol_to_ensembl)

    # Save cross-ancestry gene list
    cross_out = cross_genes_df[["gene", "best_pp4_eur", "n_eur_sources_sig",
                                 "best_eur_source", "best_pp4_eas",
                                 "n_eas_sources_sig", "best_eas_source",
                                 "n_ancestry_sig", "class"]].copy()
    cross_out.to_csv(OUTDIR / "cross_ancestry_genes.csv", index=False)
    print(f"[172] Saved cross_ancestry_genes.csv: {len(cross_out)} genes")

    # -----------------------------------------------------------------------
    # Step 3: Load expression data for COLOC genes
    # -----------------------------------------------------------------------
    print("\n[172] Step 3: Loading expression data for COLOC genes...")
    all_coloc_symbols = list(
        set(cross_genes_df["gene"].tolist() + eur_only_genes_df["gene"].tolist())
    )
    expr_df = load_expression_for_genes(all_coloc_symbols, symbol_to_ensembl)

    if len(expr_df) == 0:
        print("[172] ERROR: No expression data found for COLOC genes. Exiting.")
        sys.exit(1)

    # Identify which genes actually have expression
    expr_gene_cols = [c for c in expr_df.columns if c != "sample_id"]
    cross_expr_genes = [g for g in cross_genes_df["gene"] if g in expr_gene_cols]
    eur_only_expr_genes = [g for g in eur_only_genes_df["gene"] if g in expr_gene_cols]

    print(f"[172] Cross-ancestry genes with expression: {len(cross_expr_genes)}")
    print(f"[172] EUR-only genes with expression: {len(eur_only_expr_genes)}")

    # -----------------------------------------------------------------------
    # Step 4: Load baseline features + labels
    # -----------------------------------------------------------------------
    print("\n[172] Step 4: Loading baseline features and labels...")
    features = pd.read_csv(FEATURE_FILE)
    labels = pd.read_csv(LABEL_FILE)

    # Merge everything
    df = features.merge(labels, on="sample_id", how="inner")
    df = df.merge(expr_df, on="sample_id", how="left")
    print(f"[172] Merged dataset: {len(df)} samples, {df.shape[1]} columns")

    # Identify div_ columns (baseline)
    div_cols = [c for c in features.columns if c.startswith("div_")]
    print(f"[172] Baseline div_ features: {len(div_cols)}")

    # Prefix COLOC expression columns to avoid name collisions
    # (some COLOC genes might share names with other columns)
    coloc_cross_cols = []
    coloc_eur_cols = []

    for g in cross_expr_genes:
        new_name = f"coloc_{g}"
        if g in df.columns and new_name not in df.columns:
            df[new_name] = df[g]
        coloc_cross_cols.append(new_name)

    for g in eur_only_expr_genes:
        new_name = f"coloc_{g}"
        if g in df.columns and new_name not in df.columns:
            df[new_name] = df[g]
        coloc_eur_cols.append(new_name)

    # Verify columns exist
    coloc_cross_cols = [c for c in coloc_cross_cols if c in df.columns]
    coloc_eur_cols = [c for c in coloc_eur_cols if c in df.columns]

    # Check for NaN in COLOC expression features, fill with 0 (lowly expressed)
    for cols in [coloc_cross_cols, coloc_eur_cols]:
        for c in cols:
            n_na = df[c].isna().sum()
            if n_na > 0:
                df[c] = df[c].fillna(0.0)

    print(f"[172] COLOC cross-ancestry feature columns: {len(coloc_cross_cols)}")
    print(f"[172] COLOC EUR-only feature columns: {len(coloc_eur_cols)}")

    # -----------------------------------------------------------------------
    # Step 5: Filter to valid LOCO samples
    # -----------------------------------------------------------------------
    valid_folds = df[
        (df["loco_fold_fibrosis"] != "excluded") & df["loco_fold_fibrosis"].notna()
    ].copy()
    fold_names = sorted(valid_folds["loco_fold_fibrosis"].unique())
    print(f"\n[172] Valid LOCO samples: {len(valid_folds)}, folds: {fold_names}")

    # Target
    target_col = "fib_ge3"

    # -----------------------------------------------------------------------
    # Step 6: Define feature set configurations
    # -----------------------------------------------------------------------
    configs = {
        "M3_baseline_div100": {
            "description": "All 100 div_* divergence genes (Script 161 baseline)",
            "cols": div_cols,
        },
        "cross_ancestry_coloc": {
            "description": f"Cross-ancestry COLOC genes (PP.H4>{PP4_THRESHOLD} EUR+EAS)",
            "cols": coloc_cross_cols,
        },
        "eur_only_coloc": {
            "description": f"EUR-only COLOC genes (PP.H4>{PP4_THRESHOLD} EUR, not EAS)",
            "cols": coloc_eur_cols,
        },
        "all_coloc": {
            "description": f"All COLOC genes (PP.H4>{PP4_THRESHOLD} in any ancestry)",
            "cols": coloc_cross_cols + coloc_eur_cols,
        },
    }

    # Report feature counts
    for cfg_name, cfg in configs.items():
        print(f"  {cfg_name}: {len(cfg['cols'])} features — {cfg['description']}")

    # Warn if cross-ancestry set is very small
    if len(coloc_cross_cols) < 10:
        print(f"\n  WARNING: Only {len(coloc_cross_cols)} cross-ancestry genes available. "
              f"Classifier may be underpowered.")

    # -----------------------------------------------------------------------
    # Step 7: Run LOCO-CV for each configuration
    # -----------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("[172] Step 7: Running LOCO-CV classifiers...")
    print(f"{'='*70}")

    all_results = []
    all_predictions = []
    all_importance = []

    for config_name, config in configs.items():
        expr_cols = config["cols"]

        if len(expr_cols) == 0:
            print(f"\n  --- {config_name}: SKIP (no features) ---")
            continue

        print(f"\n  --- {config_name}: {len(expr_cols)} features ---")
        config_importance = {}

        for fold in fold_names:
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

            if len(y_train.unique()) < 2 or len(y_test.unique()) < 2:
                print(f"    Fold {fold}: SKIP (single class)")
                continue

            # Drop any constant features within this fold's training data
            usable_cols = [c for c in expr_cols
                           if c in train_df.columns and train_df[c].std() > 0]
            if len(usable_cols) == 0:
                print(f"    Fold {fold}: SKIP (all features constant)")
                continue

            # Run nested inner CV
            y_prob, coefs, selected_genes, params = fit_expression_model_inner_cv(
                train_df, y_train, test_df, usable_cols
            )

            # Compute metrics
            auroc, auprc, brier = compute_metrics(y_test.values, y_prob)

            print(
                f"    Fold {fold}: n_train={len(train_df)}, n_test={len(test_df)}, "
                f"AUROC={auroc:.3f}, AUPRC={auprc:.3f}, Brier={brier:.3f}, "
                f"n_sel={params.get('n_genes', len(selected_genes))}"
            )

            all_results.append({
                "config": config_name,
                "fold": fold,
                "n_train": len(train_df),
                "n_test": len(test_df),
                "n_pos_train": int(y_train.sum()),
                "n_pos_test": int(y_test.sum()),
                "n_features_available": len(usable_cols),
                "n_features_selected": params.get("n_genes", len(selected_genes)),
                "auroc": auroc,
                "auprc": auprc,
                "brier": brier,
                "best_l1_ratio": params.get("l1_ratio", np.nan),
                "inner_auroc": params.get("inner_auroc", np.nan),
            })

            for sid, prob, true_label in zip(
                test_df["sample_id"].values, y_prob, y_test.values
            ):
                all_predictions.append({
                    "sample_id": sid,
                    "config": config_name,
                    "predicted_prob": prob,
                    "true_label": int(true_label),
                    "fold": fold,
                })

            for feat, coef in coefs.items():
                if feat not in config_importance:
                    config_importance[feat] = []
                config_importance[feat].append(coef)

        # Summarize importance for this config
        for feat, coef_list in config_importance.items():
            # Map back to gene symbol
            gene_sym = feat.replace("coloc_", "") if feat.startswith("coloc_") else feat
            gene_sym = gene_sym.replace("div_", "") if gene_sym.startswith("div_") else gene_sym
            all_importance.append({
                "feature": feat,
                "gene_symbol": gene_sym,
                "config": config_name,
                "n_folds_selected": len(coef_list),
                "mean_coef": np.mean(coef_list),
                "mean_abs_coef": np.mean(np.abs(coef_list)),
                "std_coef": np.std(coef_list) if len(coef_list) > 1 else 0,
            })

    # -----------------------------------------------------------------------
    # Step 8: Save outputs
    # -----------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("[172] Step 8: Saving outputs...")
    print(f"{'='*70}")

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUTDIR / "cross_ancestry_classifier_results.csv", index=False)
    print(f"  cross_ancestry_classifier_results.csv: {len(results_df)} rows")

    pred_df = pd.DataFrame(all_predictions)
    pred_df.to_csv(OUTDIR / "cross_ancestry_predictions.csv", index=False)
    print(f"  cross_ancestry_predictions.csv: {len(pred_df)} rows")

    imp_df = pd.DataFrame(all_importance)
    imp_df.to_csv(OUTDIR / "cross_ancestry_feature_importance.csv", index=False)
    print(f"  cross_ancestry_feature_importance.csv: {len(imp_df)} rows")

    # -----------------------------------------------------------------------
    # Step 9: Summary comparison
    # -----------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("[172] Step 9: Comparison summary...")
    print(f"{'='*70}")

    comparison_rows = []
    for config_name in configs:
        sub = results_df[results_df["config"] == config_name]
        if len(sub) == 0:
            comparison_rows.append({
                "config": config_name,
                "description": configs[config_name]["description"],
                "n_features_total": len(configs[config_name]["cols"]),
                "n_folds": 0,
                "mean_auroc": np.nan,
                "std_auroc": np.nan,
                "cv_auroc": np.nan,
                "mean_auprc": np.nan,
                "std_auprc": np.nan,
                "mean_brier": np.nan,
                "std_brier": np.nan,
                "min_auroc": np.nan,
                "max_auroc": np.nan,
                "range_auroc": np.nan,
                "median_n_selected": np.nan,
            })
            continue

        mean_auroc = sub["auroc"].mean()
        std_auroc = sub["auroc"].std()
        cv_auroc = std_auroc / mean_auroc if mean_auroc > 0 else np.nan

        comparison_rows.append({
            "config": config_name,
            "description": configs[config_name]["description"],
            "n_features_total": len(configs[config_name]["cols"]),
            "n_folds": len(sub),
            "mean_auroc": mean_auroc,
            "std_auroc": std_auroc,
            "cv_auroc": cv_auroc,
            "mean_auprc": sub["auprc"].mean(),
            "std_auprc": sub["auprc"].std(),
            "mean_brier": sub["brier"].mean(),
            "std_brier": sub["brier"].std(),
            "min_auroc": sub["auroc"].min(),
            "max_auroc": sub["auroc"].max(),
            "range_auroc": sub["auroc"].max() - sub["auroc"].min(),
            "median_n_selected": sub["n_features_selected"].median(),
        })

    comparison_df = pd.DataFrame(comparison_rows)
    comparison_df.to_csv(OUTDIR / "cross_ancestry_comparison.csv", index=False)
    print(f"  cross_ancestry_comparison.csv: {len(comparison_df)} rows")

    # Print summary table
    print(f"\n{'='*70}")
    print("[172] RESULTS SUMMARY")
    print(f"{'='*70}")
    print(f"  Target: fib_ge3 (fibrosis >= F3)")
    print(f"  PP.H4 threshold: {PP4_THRESHOLD}")
    print(f"  Cross-ancestry genes: {len(cross_expr_genes)}")
    print(f"  EUR-only genes: {len(eur_only_expr_genes)}")
    print()

    for _, row in comparison_df.iterrows():
        if row["n_folds"] == 0:
            print(f"  {row['config']:30s} | NO FOLDS (0 features)")
        else:
            print(
                f"  {row['config']:30s} | "
                f"n_feat={row['n_features_total']:4.0f} | "
                f"AUROC={row['mean_auroc']:.3f} +/- {row['std_auroc']:.3f} | "
                f"CV={row['cv_auroc']:.3f} | "
                f"range=[{row['min_auroc']:.3f}, {row['max_auroc']:.3f}] | "
                f"AUPRC={row['mean_auprc']:.3f} | "
                f"sel={row['median_n_selected']:.0f}"
            )

    # Hypothesis test: is cross-ancestry variance lower than baseline?
    baseline_sub = results_df[results_df["config"] == "M3_baseline_div100"]
    cross_sub = results_df[results_df["config"] == "cross_ancestry_coloc"]

    if len(baseline_sub) > 2 and len(cross_sub) > 2:
        print(f"\n  Hypothesis: Cross-ancestry features have lower LOCO variance")
        print(f"  Baseline std:      {baseline_sub['auroc'].std():.4f}")
        print(f"  Cross-ancestry std: {cross_sub['auroc'].std():.4f}")

        # Levene's test for equality of variances
        from scipy import stats
        stat, p = stats.levene(
            baseline_sub["auroc"].dropna().values,
            cross_sub["auroc"].dropna().values
        )
        print(f"  Levene's test: statistic={stat:.3f}, p={p:.4f}")

        # Also report paired difference (fold-matched)
        merged_folds = baseline_sub[["fold", "auroc"]].merge(
            cross_sub[["fold", "auroc"]],
            on="fold",
            suffixes=("_baseline", "_cross"),
        )
        if len(merged_folds) > 0:
            diff = merged_folds["auroc_baseline"] - merged_folds["auroc_cross"]
            print(f"  Paired AUROC difference (baseline - cross): "
                  f"mean={diff.mean():.4f}, median={diff.median():.4f}")
            if len(merged_folds) > 2:
                t_stat, t_p = stats.ttest_rel(
                    merged_folds["auroc_baseline"].values,
                    merged_folds["auroc_cross"].values,
                )
                print(f"  Paired t-test: t={t_stat:.3f}, p={t_p:.4f}")

    # Feature stability analysis
    print(f"\n  Feature stability (genes selected in N/6 folds):")
    for config_name in configs:
        sub_imp = imp_df[imp_df["config"] == config_name]
        if len(sub_imp) == 0:
            continue
        n_total = len(sub_imp)
        n_all_folds = (sub_imp["n_folds_selected"] == sub_imp["n_folds_selected"].max()).sum()
        n_half = (sub_imp["n_folds_selected"] >= 3).sum()
        print(f"    {config_name}: {n_total} total features used, "
              f"{n_all_folds} in all folds, {n_half} in >=3 folds")

    elapsed = time.time() - t0
    print(f"\n[172] Done in {elapsed/60:.1f} min")


if __name__ == "__main__":
    main()
