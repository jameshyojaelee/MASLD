#!/usr/bin/env python3
"""
109_wenda_adaptation.py
WENDA-inspired domain adaptation for cross-cohort MASLD staging.

WENDA (Weighted Elastic Net for Domain Adaptation) reweights features for
cross-cohort consistency by downweighting cohort-specific genes and upweighting
genes with consistent signal across cohorts (Greene Lab, greenelab/wenda_gpu).

This implements the core WENDA idea without the GPU package:
  1. Load merged_dge.rds -> TMM logCPM (1,444 x ~34K genes)
  2. Load modeling_metadata.csv for fib_ge3 and loco_fold_fibrosis
  3. For each gene, compute domain stability = within-cohort variance / total variance
     High stability = gene behaves consistently across cohorts
  4. Weight genes by domain stability score
  5. Run elastic net for F>=3 binary, 6-fold fibrosis LOCO with within-fold
     feature selection (top 3K by stability-weighted variance)
  6. Compare AUROC: weighted vs unweighted (V3 baseline is 0.881)
  7. Also run for NAS>=5 binary, 5-fold NAS LOCO

Output (all to results/staging_classifier/):
  - wenda_results.csv          (per-sample predictions + AUROC per fold)
  - wenda_feature_weights.csv  (per-gene domain stability scores)

SLURM: cpu partition, 8 CPUs, 64G RAM, 48h
Env:   micromamba activate rapids_singlecell

Usage:
  sbatch --job-name=stg109_wenda \
         --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00 \
         --output=logs/109_wenda_%j.out \
         --error=logs/109_wenda_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate rapids_singlecell && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 109_wenda_adaptation.py'"
"""

# Must be set before ANY other imports to prevent CUDA crash on cpu nodes
import os
os.environ["RAPIDS_NO_INITIALIZE"] = "1"
os.environ["CUDF_PANDAS"] = "0"

import sys
import time
import warnings
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
SEED = 42
np.random.seed(SEED)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
OUTDIR = os.path.join(RDIR, "staging_classifier")
LOGDIR = os.path.join(OUTDIR, "logs")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(LOGDIR, exist_ok=True)

RSCRIPT = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
META_PATH = os.path.join(OUTDIR, "modeling_metadata.csv")
DGE_PATH = os.path.join(RDIR, "integration/merged_dge.rds")

# Hyperparameters
K_FEATURES = 3000  # Number of features to select per fold (match V3)
ALPHA = 0.5        # Elastic net mixing (0.5 = balanced L1/L2, same as V3)
STABILITY_FLOOR = 0.01  # Minimum stability weight (avoid zero-weight genes)

# LOCO fold definitions (must match 62_prepare_python_data.py)
FIB_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478",
                "GSE193066", "GSE240729"]
NAS_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066"]

print("=" * 70)
print("109: WENDA-Inspired Domain Adaptation for MASLD Staging")
print("=" * 70)
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Output dir: {OUTDIR}")

# ============================================================
# Step 1: Load merged_dge.rds -> TMM logCPM via Rscript
# ============================================================
print("\n--- Step 1: Load expression data (merged_dge.rds -> TMM logCPM) ---")

if not os.path.exists(DGE_PATH):
    sys.exit(f"ERROR: {DGE_PATH} not found. Run integration pipeline first.")

# Convert DGE object to logCPM CSV via R subprocess
logcpm_cache = os.path.join(OUTDIR, "_wenda_logcpm_cache.csv.gz")

if os.path.exists(logcpm_cache):
    print(f"  Loading cached logCPM: {logcpm_cache}")
    logcpm_df = pd.read_csv(logcpm_cache, index_col=0)
    print(f"  logCPM matrix: {logcpm_df.shape[0]} genes x {logcpm_df.shape[1]} samples")
else:
    print(f"  Converting DGE -> TMM logCPM via Rscript...")
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "logcpm.csv")
        r_code = f"""
        suppressPackageStartupMessages(library(edgeR))
        dge <- readRDS("{DGE_PATH}")
        dge <- calcNormFactors(dge, method = "TMM")
        logcpm <- cpm(dge, log = TRUE, prior.count = 1)
        cat("logCPM:", nrow(logcpm), "genes x", ncol(logcpm), "samples\\n")
        write.csv(logcpm, "{csv_path}", quote = FALSE)
        """
        result = subprocess.run(
            [RSCRIPT, "-e", r_code],
            capture_output=True, text=True, timeout=600
        )
        if result.returncode != 0:
            print(f"R stderr: {result.stderr}")
            sys.exit("ERROR: Failed to convert merged_dge.rds to logCPM")
        if result.stdout:
            print(f"  R output: {result.stdout.strip()}")

        logcpm_df = pd.read_csv(csv_path, index_col=0)
        print(f"  logCPM matrix: {logcpm_df.shape[0]} genes x {logcpm_df.shape[1]} samples")

        # Cache for re-runs
        logcpm_df.to_csv(logcpm_cache, compression="gzip")
        print(f"  Cached to {logcpm_cache}")

# ============================================================
# Step 2: Load metadata
# ============================================================
print("\n--- Step 2: Load modeling metadata ---")
meta = pd.read_csv(META_PATH)
meta["sample_id"] = meta["sample_id"].astype(str).str.strip()
print(f"  Metadata: {len(meta)} samples")

# Align samples between expression and metadata
expr_samples = set(logcpm_df.columns)
meta_samples = set(meta["sample_id"])
common_samples = sorted(expr_samples & meta_samples)
print(f"  Common samples: {len(common_samples)}")

meta = meta[meta["sample_id"].isin(common_samples)].copy()
meta = meta.set_index("sample_id").loc[common_samples].reset_index()
logcpm_df = logcpm_df[common_samples]

# ============================================================
# Step 3: Compute domain stability scores (WENDA core idea)
# ============================================================
print("\n--- Step 3: Compute domain stability scores ---")

# Expression matrix as numpy: genes x samples
X_all = logcpm_df.values  # (n_genes, n_samples)
gene_names = logcpm_df.index.tolist()
n_genes, n_samples = X_all.shape

# Map each sample to its cohort
sample_to_cohort = dict(zip(meta["sample_id"], meta["dataset"]))
cohort_labels = np.array([sample_to_cohort[s] for s in common_samples])
unique_cohorts = np.unique(cohort_labels)
n_cohorts = len(unique_cohorts)
print(f"  Cohorts: {n_cohorts} ({', '.join(unique_cohorts)})")

# For each gene: total variance, within-cohort variance, between-cohort variance
# Within-cohort variance = weighted average of per-cohort variances
# Between-cohort variance = variance of cohort means
# domain_stability = within_var / total_var (high = consistent across cohorts)

total_var = np.var(X_all, axis=1, ddof=1)  # (n_genes,)

# Compute within-cohort variance (pooled)
within_var = np.zeros(n_genes)
cohort_means = np.zeros((n_genes, n_cohorts))

for ci, cohort in enumerate(unique_cohorts):
    mask = cohort_labels == cohort
    n_c = mask.sum()
    if n_c < 2:
        continue
    X_c = X_all[:, mask]
    cohort_means[:, ci] = X_c.mean(axis=1)
    within_var += (n_c - 1) * np.var(X_c, axis=1, ddof=1)

within_var /= (n_samples - n_cohorts)  # Pooled within-cohort variance

# Between-cohort variance: variance of cohort means (weighted by cohort size)
cohort_sizes = np.array([np.sum(cohort_labels == c) for c in unique_cohorts])
grand_mean = X_all.mean(axis=1, keepdims=True)
between_var = np.zeros(n_genes)
for ci, cohort in enumerate(unique_cohorts):
    n_c = cohort_sizes[ci]
    between_var += n_c * (cohort_means[:, ci] - grand_mean.ravel()) ** 2
between_var /= (n_cohorts - 1)

# Domain stability = within / total (1 = all variance is within-cohort = consistent)
# Genes with high between-cohort variance -> low stability -> downweight
# Add floor to avoid zero weights
stability = np.where(total_var > 0, within_var / (within_var + between_var), 0.5)
stability = np.clip(stability, STABILITY_FLOOR, 1.0)

# Create stability DataFrame
stability_df = pd.DataFrame({
    "gene": gene_names,
    "total_var": total_var,
    "within_var": within_var,
    "between_var": between_var,
    "domain_stability": stability,
})
stability_df = stability_df.sort_values("domain_stability", ascending=False)

print(f"  Domain stability: median={np.median(stability):.3f}, "
      f"mean={np.mean(stability):.3f}, "
      f"min={np.min(stability):.4f}, max={np.max(stability):.4f}")
print(f"  Genes with stability > 0.9: {(stability > 0.9).sum()}")
print(f"  Genes with stability < 0.5: {(stability < 0.5).sum()}")

# Save feature weights
weights_path = os.path.join(OUTDIR, "wenda_feature_weights.csv")
stability_df.to_csv(weights_path, index=False)
print(f"  Saved: {weights_path}")


# ============================================================
# Step 4: Helper functions for LOCO elastic net
# ============================================================

def select_features_weighted(logcpm_train, stability_scores, k=K_FEATURES):
    """Select top K features by stability-weighted variance (within-fold only).

    Args:
        logcpm_train: (n_genes, n_train_samples) training expression
        stability_scores: (n_genes,) domain stability weights per gene
        k: number of features to select

    Returns:
        indices of top K genes
    """
    gene_var = np.var(logcpm_train, axis=1, ddof=1)
    weighted_var = gene_var * stability_scores
    top_k = np.argsort(weighted_var)[::-1][:k]
    return top_k


def select_features_unweighted(logcpm_train, k=K_FEATURES):
    """Select top K features by plain variance (within-fold, no weighting).

    This replicates V3 feature selection for the unweighted baseline.
    """
    gene_var = np.var(logcpm_train, axis=1, ddof=1)
    top_k = np.argsort(gene_var)[::-1][:k]
    return top_k


def rank_transform(X):
    """Rank transform columns (samples). Returns (n_features, n_samples)."""
    from scipy.stats import rankdata
    n_feat, n_samp = X.shape
    ranked = np.zeros_like(X, dtype=np.float64)
    for j in range(n_samp):
        r = rankdata(X[:, j], method="average")
        ranked[:, j] = r / len(r)
    return ranked


def run_elastic_net_fold(X_train, y_train, X_test, alpha=ALPHA, seed=SEED):
    """Train elastic net (logistic regression with elastic net penalty).

    Uses sklearn LogisticRegression with saga solver for elastic net support.
    Sweeps C (inverse regularization) via internal 5-fold CV on training data.

    Returns: predicted probabilities for test samples.
    """
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    # Sweep C values (matching glmnet lambda sweep behavior)
    best_auroc = -1
    best_C = 1.0
    best_model = None

    from sklearn.model_selection import StratifiedKFold

    C_values = [0.001, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0]
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)

    for C in C_values:
        inner_probs = np.zeros(len(y_train))
        inner_valid = np.zeros(len(y_train), dtype=bool)

        for inner_train, inner_val in skf.split(X_train_s, y_train):
            if len(np.unique(y_train[inner_train])) < 2:
                continue
            lr = LogisticRegression(
                penalty="elasticnet", solver="saga", l1_ratio=alpha,
                C=C, max_iter=10000, class_weight="balanced",
                random_state=seed, n_jobs=1
            )
            lr.fit(X_train_s[inner_train], y_train[inner_train])
            inner_probs[inner_val] = lr.predict_proba(X_train_s[inner_val])[:, 1]
            inner_valid[inner_val] = True

        if inner_valid.sum() > 0 and len(np.unique(y_train[inner_valid])) >= 2:
            try:
                auroc = roc_auc_score(y_train[inner_valid], inner_probs[inner_valid])
                if auroc > best_auroc:
                    best_auroc = auroc
                    best_C = C
            except ValueError:
                pass

    # Refit on full training data with best C
    final_model = LogisticRegression(
        penalty="elasticnet", solver="saga", l1_ratio=alpha,
        C=best_C, max_iter=10000, class_weight="balanced",
        random_state=seed, n_jobs=1
    )
    final_model.fit(X_train_s, y_train)
    probs = final_model.predict_proba(X_test_s)[:, 1]

    return probs, best_C, final_model


# ============================================================
# Step 5: LOCO evaluation for F>=3 binary (fibrosis)
# ============================================================
print("\n--- Step 5: F>=3 Binary LOCO (6 folds) ---")

# Filter to fibrosis-eligible samples
fib_meta = meta[
    (meta["loco_fold_fibrosis"] != "excluded") &
    (meta["fib_ge3"] >= 0)
].copy()
print(f"  Fibrosis-eligible: {len(fib_meta)} samples")
fib_folds = sorted(fib_meta["loco_fold_fibrosis"].unique())
print(f"  Folds: {fib_folds}")

fib_results_weighted = []
fib_results_unweighted = []

for fold in fib_folds:
    test_mask = fib_meta["loco_fold_fibrosis"] == fold
    train_mask = fib_meta["loco_fold_fibrosis"] != fold

    test_ids = fib_meta[test_mask]["sample_id"].values
    train_ids = fib_meta[train_mask]["sample_id"].values
    y_test = fib_meta[test_mask]["fib_ge3"].values.astype(int)
    y_train = fib_meta[train_mask]["fib_ge3"].values.astype(int)

    if len(test_ids) < 5 or len(train_ids) < 20:
        print(f"  {fold}: SKIP (test={len(test_ids)}, train={len(train_ids)})")
        continue
    if len(np.unique(y_test)) < 2:
        print(f"  {fold}: SKIP (single class in test)")
        continue

    # Get expression for train/test
    train_col_idx = [common_samples.index(s) for s in train_ids]
    test_col_idx = [common_samples.index(s) for s in test_ids]
    X_train_raw = X_all[:, train_col_idx]
    X_test_raw = X_all[:, test_col_idx]

    # --- WEIGHTED (WENDA) ---
    sel_w = select_features_weighted(X_train_raw, stability, k=K_FEATURES)
    X_tr_w = rank_transform(X_train_raw[sel_w, :]).T  # (n_train, k)
    X_te_w = rank_transform(X_test_raw[sel_w, :]).T   # (n_test, k)

    probs_w, best_C_w, _ = run_elastic_net_fold(X_tr_w, y_train, X_te_w)
    auroc_w = roc_auc_score(y_test, probs_w)

    # --- UNWEIGHTED (V3-style baseline) ---
    sel_u = select_features_unweighted(X_train_raw, k=K_FEATURES)
    X_tr_u = rank_transform(X_train_raw[sel_u, :]).T
    X_te_u = rank_transform(X_test_raw[sel_u, :]).T

    probs_u, best_C_u, _ = run_elastic_net_fold(X_tr_u, y_train, X_te_u)
    auroc_u = roc_auc_score(y_test, probs_u)

    # Feature overlap
    overlap = len(set(sel_w) & set(sel_u))
    jaccard = overlap / len(set(sel_w) | set(sel_u))

    print(f"  {fold}: WENDA AUROC={auroc_w:.3f} (C={best_C_w})  |  "
          f"Unweighted AUROC={auroc_u:.3f} (C={best_C_u})  |  "
          f"Feature overlap Jaccard={jaccard:.3f}  "
          f"(n_test={len(y_test)}, n_train={len(y_train)})")

    # Save per-sample predictions
    for i, sid in enumerate(test_ids):
        fib_results_weighted.append({
            "sample_id": sid,
            "fold": fold,
            "target": "fib_ge3",
            "method": "wenda_weighted",
            "true_label": int(y_test[i]),
            "probability": float(probs_w[i]),
            "predicted": int(probs_w[i] > 0.5),
            "fold_auroc": float(auroc_w),
            "best_C": float(best_C_w),
        })
        fib_results_unweighted.append({
            "sample_id": sid,
            "fold": fold,
            "target": "fib_ge3",
            "method": "unweighted_baseline",
            "true_label": int(y_test[i]),
            "probability": float(probs_u[i]),
            "predicted": int(probs_u[i] > 0.5),
            "fold_auroc": float(auroc_u),
            "best_C": float(best_C_u),
        })

# Aggregate fibrosis results
if fib_results_weighted:
    df_w = pd.DataFrame(fib_results_weighted)
    df_u = pd.DataFrame(fib_results_unweighted)

    overall_w = roc_auc_score(df_w["true_label"], df_w["probability"])
    overall_u = roc_auc_score(df_u["true_label"], df_u["probability"])

    print(f"\n  *** F>=3 Overall WENDA AUROC:      {overall_w:.4f} ***")
    print(f"  *** F>=3 Overall Unweighted AUROC:  {overall_u:.4f} ***")
    print(f"  *** V3 Baseline (reference):        0.8810 ***")
    print(f"  *** WENDA delta vs unweighted:      {overall_w - overall_u:+.4f} ***")
else:
    print("  WARNING: No fibrosis fold results produced")
    overall_w = overall_u = None

# ============================================================
# Step 6: LOCO evaluation for NAS>=5 binary
# ============================================================
print("\n--- Step 6: NAS>=5 Binary LOCO (5 folds) ---")

nas_meta = meta[
    (meta["loco_fold_nas"] != "excluded") &
    (meta["nas_ge5"] >= 0)
].copy()
print(f"  NAS-eligible: {len(nas_meta)} samples")
nas_folds = sorted(nas_meta["loco_fold_nas"].unique())
print(f"  Folds: {nas_folds}")

nas_results_weighted = []
nas_results_unweighted = []

for fold in nas_folds:
    test_mask = nas_meta["loco_fold_nas"] == fold
    train_mask = nas_meta["loco_fold_nas"] != fold

    test_ids = nas_meta[test_mask]["sample_id"].values
    train_ids = nas_meta[train_mask]["sample_id"].values
    y_test = nas_meta[test_mask]["nas_ge5"].values.astype(int)
    y_train = nas_meta[train_mask]["nas_ge5"].values.astype(int)

    if len(test_ids) < 5 or len(train_ids) < 20:
        print(f"  {fold}: SKIP (test={len(test_ids)}, train={len(train_ids)})")
        continue
    if len(np.unique(y_test)) < 2:
        print(f"  {fold}: SKIP (single class in test)")
        continue

    train_col_idx = [common_samples.index(s) for s in train_ids]
    test_col_idx = [common_samples.index(s) for s in test_ids]
    X_train_raw = X_all[:, train_col_idx]
    X_test_raw = X_all[:, test_col_idx]

    # --- WEIGHTED (WENDA) ---
    sel_w = select_features_weighted(X_train_raw, stability, k=K_FEATURES)
    X_tr_w = rank_transform(X_train_raw[sel_w, :]).T
    X_te_w = rank_transform(X_test_raw[sel_w, :]).T

    probs_w, best_C_w, _ = run_elastic_net_fold(X_tr_w, y_train, X_te_w)
    auroc_w = roc_auc_score(y_test, probs_w)

    # --- UNWEIGHTED ---
    sel_u = select_features_unweighted(X_train_raw, k=K_FEATURES)
    X_tr_u = rank_transform(X_train_raw[sel_u, :]).T
    X_te_u = rank_transform(X_test_raw[sel_u, :]).T

    probs_u, best_C_u, _ = run_elastic_net_fold(X_tr_u, y_train, X_te_u)
    auroc_u = roc_auc_score(y_test, probs_u)

    overlap = len(set(sel_w) & set(sel_u))
    jaccard = overlap / len(set(sel_w) | set(sel_u))

    print(f"  {fold}: WENDA AUROC={auroc_w:.3f} (C={best_C_w})  |  "
          f"Unweighted AUROC={auroc_u:.3f} (C={best_C_u})  |  "
          f"Feature overlap Jaccard={jaccard:.3f}  "
          f"(n_test={len(y_test)}, n_train={len(y_train)})")

    for i, sid in enumerate(test_ids):
        nas_results_weighted.append({
            "sample_id": sid,
            "fold": fold,
            "target": "nas_ge5",
            "method": "wenda_weighted",
            "true_label": int(y_test[i]),
            "probability": float(probs_w[i]),
            "predicted": int(probs_w[i] > 0.5),
            "fold_auroc": float(auroc_w),
            "best_C": float(best_C_w),
        })
        nas_results_unweighted.append({
            "sample_id": sid,
            "fold": fold,
            "target": "nas_ge5",
            "method": "unweighted_baseline",
            "true_label": int(y_test[i]),
            "probability": float(probs_u[i]),
            "predicted": int(probs_u[i] > 0.5),
            "fold_auroc": float(auroc_u),
            "best_C": float(best_C_u),
        })

# Aggregate NAS results
if nas_results_weighted:
    df_nw = pd.DataFrame(nas_results_weighted)
    df_nu = pd.DataFrame(nas_results_unweighted)

    nas_overall_w = roc_auc_score(df_nw["true_label"], df_nw["probability"])
    nas_overall_u = roc_auc_score(df_nu["true_label"], df_nu["probability"])

    print(f"\n  *** NAS>=5 Overall WENDA AUROC:      {nas_overall_w:.4f} ***")
    print(f"  *** NAS>=5 Overall Unweighted AUROC:  {nas_overall_u:.4f} ***")
    print(f"  *** WENDA delta vs unweighted:        {nas_overall_w - nas_overall_u:+.4f} ***")
else:
    print("  WARNING: No NAS fold results produced")
    nas_overall_w = nas_overall_u = None

# ============================================================
# Step 7: Save all results
# ============================================================
print("\n--- Step 7: Save results ---")

all_results = (fib_results_weighted + fib_results_unweighted +
               nas_results_weighted + nas_results_unweighted)
results_df = pd.DataFrame(all_results)
results_path = os.path.join(OUTDIR, "wenda_results.csv")
results_df.to_csv(results_path, index=False)
print(f"  Saved: {results_path} ({len(results_df)} rows)")

# Summary table
summary_rows = []
if overall_w is not None:
    summary_rows.append({
        "target": "fib_ge3",
        "method": "wenda_weighted",
        "overall_auroc": overall_w,
        "n_samples": len(fib_results_weighted),
        "n_folds": len(fib_folds),
    })
    summary_rows.append({
        "target": "fib_ge3",
        "method": "unweighted_baseline",
        "overall_auroc": overall_u,
        "n_samples": len(fib_results_unweighted),
        "n_folds": len(fib_folds),
    })
    summary_rows.append({
        "target": "fib_ge3",
        "method": "v3_reference",
        "overall_auroc": 0.881,
        "n_samples": None,
        "n_folds": 6,
    })
if nas_overall_w is not None:
    summary_rows.append({
        "target": "nas_ge5",
        "method": "wenda_weighted",
        "overall_auroc": nas_overall_w,
        "n_samples": len(nas_results_weighted),
        "n_folds": len(nas_folds),
    })
    summary_rows.append({
        "target": "nas_ge5",
        "method": "unweighted_baseline",
        "overall_auroc": nas_overall_u,
        "n_samples": len(nas_results_unweighted),
        "n_folds": len(nas_folds),
    })

if summary_rows:
    summary_df = pd.DataFrame(summary_rows)
    print(f"\n  === Summary ===")
    print(summary_df.to_string(index=False))

print(f"\nCompleted: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("=== 109_wenda_adaptation.py finished ===")
