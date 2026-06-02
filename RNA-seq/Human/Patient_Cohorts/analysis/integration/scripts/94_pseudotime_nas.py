#!/usr/bin/env python3
"""
94_pseudotime_nas.py
Compute pseudotime on NAS-VAE embeddings, calibrate to NAS score.

Pipeline:
  1. Load 64-dim embeddings for 1,444 samples
  2. Compute PCA->t-SNE (2D) for visualization
  3. Compute diffusion pseudotime:
     - kNN graph (k=30) via sklearn
     - Diffusion map (eigendecomposition of transition matrix)
     - Root at centroid of lowest-NAS samples (NAS 0-1)
  4. Calibrate pseudotime -> NAS via isotonic regression on 660 labeled samples
  5. 5-fold NAS LOCO: calibrate on 4 folds, predict held-out
  6. Compute NAS QWK and MAE from calibrated predictions
  7. Save UMAP + pseudotime visualization data

Input:
  - results/staging_classifier/nas_embeddings_all_samples.csv (fallback to embeddings_all_samples.csv)
  - results/staging_classifier/modeling_metadata.csv (NAS labels)

Output:
  - results/staging_classifier/pseudotime_nas_results.csv
  - results/staging_classifier/pseudotime_calibration.csv
  - results/staging_classifier/pseudotime_umap_data.csv

SLURM: cpu partition, 8 CPUs, 32GB RAM, 48h
Env:    micromamba activate rapids_singlecell  (or rnaseq)

NOTE: This script avoids scipy entirely to prevent CUDA initialization errors
      in conda envs where scipy is compiled against cupy/CUDA.
      All sparse matrix and statistical operations use numpy directly.
      The matrix is at most ~1444x1444 so dense numpy is perfectly adequate.
"""

# ---------------------------------------------------------------------------
# CRITICAL: Hide all GPUs BEFORE any imports to prevent CUDA initialization.
# scipy in rapids_singlecell env links against cupy, which auto-initializes
# CUDA on import. Setting CUDA_VISIBLE_DEVICES="" prevents this entirely.
# ---------------------------------------------------------------------------
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["RAPIDS_NO_INITIALIZE"] = "1"
os.environ["CUDF_PANDAS"] = "0"

import sys
import time
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

from sklearn.neighbors import NearestNeighbors
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import (
    cohen_kappa_score, mean_absolute_error, balanced_accuracy_score,
)

# ---------------------------------------------------------------------------
# Pure-numpy replacements for scipy functions
# ---------------------------------------------------------------------------

def _spearmanr(x, y):
    """Spearman rank correlation (pure numpy, no scipy)."""
    from numpy import argsort, sqrt
    n = len(x)
    rx = np.empty(n, dtype=np.float64)
    ry = np.empty(n, dtype=np.float64)
    # Rank (average for ties)
    for arr, rank_out in [(x, rx), (y, ry)]:
        order = argsort(argsort(arr).astype(np.float64))
        # Simple ranking: no tie correction needed for continuous pseudotime
        # For integer NAS, ties are rare enough that plain rank is fine
        sorted_idx = argsort(arr)
        ranks = np.empty(n, dtype=np.float64)
        ranks[sorted_idx] = np.arange(1, n + 1, dtype=np.float64)
        # Average ties
        unique_vals = np.unique(arr)
        for v in unique_vals:
            mask = arr == v
            if mask.sum() > 1:
                ranks[mask] = ranks[mask].mean()
        rank_out[:] = ranks

    d = rx - ry
    d2_sum = (d ** 2).sum()
    rho = 1.0 - (6.0 * d2_sum) / (n * (n ** 2 - 1))
    # Approximate p-value via t-distribution (large n approximation)
    if abs(rho) < 1.0:
        t_stat = rho * sqrt((n - 2) / (1 - rho ** 2))
        # Two-tailed p-value from t(n-2): use normal approx for n >> 30
        p = 2.0 * _norm_cdf(-abs(t_stat))
    else:
        p = 0.0
    return rho, p


def _norm_cdf(x):
    """Standard normal CDF approximation (Abramowitz & Stegun)."""
    # Good to ~1e-7 accuracy
    import math
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911
    sign = 1 if x >= 0 else -1
    x = abs(x) / math.sqrt(2.0)
    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x * x)
    return 0.5 * (1.0 + sign * y)


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

NAS_EMB_PATH = os.path.join(OUTDIR, "nas_embeddings_all_samples.csv")
FALLBACK_EMB_PATH = os.path.join(OUTDIR, "embeddings_all_samples.csv")
META_PATH = os.path.join(OUTDIR, "modeling_metadata.csv")

t0 = time.time()

print("=" * 60)
print("94: Pseudotime on VAE Embeddings -> NAS Calibration")
print("=" * 60)
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading data...")

# Embeddings (prefer NAS-VAE, fallback to fibrosis-VAE)
if os.path.exists(NAS_EMB_PATH):
    emb_path_used = NAS_EMB_PATH
    print(f"  Using NAS-VAE embeddings: {NAS_EMB_PATH}")
else:
    emb_path_used = FALLBACK_EMB_PATH
    print(f"  NAS-VAE not found, using fibrosis-VAE: {FALLBACK_EMB_PATH}")

emb_df = pd.read_csv(emb_path_used, index_col="sample_id")
emb_cols = [c for c in emb_df.columns if c.startswith("z")]
X_emb = emb_df[emb_cols].values.astype(np.float64)
sample_ids = emb_df.index.values
N_SAMPLES = len(sample_ids)
print(f"  Embeddings: {N_SAMPLES} samples x {len(emb_cols)} dims")

# Metadata
meta = pd.read_csv(META_PATH, index_col="sample_id")
meta = meta.loc[sample_ids]  # Align to embedding order
print(f"  Metadata: {len(meta)} samples")

# NAS labels
nas_score = meta["nas_group"].values.astype(int)       # raw NAS 0-8, -1=missing
nas_group4 = meta["nas_group4"].values.astype(int)     # 4-group, -1=missing
loco_nas = meta["loco_fold_nas"].values
datasets = meta["dataset"].values

unique_nas_folds = sorted(set(f for f in loco_nas if f != "excluded" and f != "NA"))
has_nas = nas_score >= 0
n_nas_labeled = has_nas.sum()
print(f"  NAS-labeled: {n_nas_labeled} / {N_SAMPLES}")
print(f"  NAS LOCO folds: {len(unique_nas_folds)} ({', '.join(unique_nas_folds)})")
print(f"  NAS range: {nas_score[has_nas].min()}-{nas_score[has_nas].max()}")


# ---------------------------------------------------------------------------
# STEP 1: Dimensionality reduction for visualization (PCA -> t-SNE)
# ---------------------------------------------------------------------------
print("\n=== STEP 1: Dimensionality Reduction (PCA -> t-SNE, CPU) ===")

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

n_pca = min(20, X_emb.shape[1])
print(f"  PCA: {X_emb.shape[1]} dims -> {n_pca} dims")
pca = PCA(n_components=n_pca, random_state=42)
pca_coords = pca.fit_transform(X_emb)

print(f"  t-SNE: {n_pca} dims -> 2 dims (perplexity=30)")
tsne = TSNE(n_components=2, perplexity=30, random_state=42, n_jobs=-1)
umap_coords = tsne.fit_transform(pca_coords)

print(f"  Embedding shape: {umap_coords.shape}")


# ---------------------------------------------------------------------------
# STEP 2: Diffusion pseudotime (pure numpy, no scipy)
# ---------------------------------------------------------------------------
print("\n=== STEP 2: Diffusion Pseudotime ===")

K_NEIGHBORS = 30
N_DIFFUSION_COMPONENTS = 10

# Build kNN graph
print(f"  Building kNN graph (k={K_NEIGHBORS})...")
nn = NearestNeighbors(n_neighbors=K_NEIGHBORS, metric="euclidean", n_jobs=-1)
nn.fit(X_emb)
distances, indices = nn.kneighbors(X_emb)

# Adaptive bandwidth: sigma_i = distance to k-th neighbor
sigma = distances[:, -1].copy()
sigma[sigma < 1e-10] = 1e-10  # Prevent division by zero

# Build dense affinity matrix with adaptive Gaussian kernel
# Dense is fine for N=1,444 (matrix is 1444x1444 ~ 16MB float64)
print("  Computing affinity matrix (dense, numpy)...")
W = np.zeros((N_SAMPLES, N_SAMPLES), dtype=np.float64)
for i in range(N_SAMPLES):
    for j_idx in range(K_NEIGHBORS):
        j = indices[i, j_idx]
        if i == j:
            continue
        d = distances[i, j_idx]
        # Adaptive Gaussian kernel: K(i,j) = exp(-d^2 / (sigma_i * sigma_j))
        aff = np.exp(-(d ** 2) / (sigma[i] * sigma[j]))
        W[i, j] = aff

# Symmetrize
W = (W + W.T) / 2.0

# Row-normalize to transition matrix
print("  Computing transition matrix...")
row_sums = W.sum(axis=1)
row_sums[row_sums < 1e-10] = 1e-10
T = W / row_sums[:, np.newaxis]

# Eigendecomposition (dense numpy -- perfectly tractable for 1444x1444)
print(f"  Eigendecomposition ({N_DIFFUSION_COMPONENTS} components, numpy dense)...")
eigenvalues_all, eigenvectors_all = np.linalg.eigh(T)

# eigh returns ascending order; take top N_DIFFUSION_COMPONENTS (largest)
eigenvalues = eigenvalues_all[-N_DIFFUSION_COMPONENTS:][::-1]
eigenvectors = eigenvectors_all[:, -N_DIFFUSION_COMPONENTS:][:, ::-1]

print(f"  Top eigenvalues: {eigenvalues[:5]}")

# Diffusion map coordinates (skip first trivial eigenvector)
# DC_i = lambda_i * psi_i
diffusion_coords = eigenvectors[:, 1:] * eigenvalues[1:]
print(f"  Diffusion coordinates: {diffusion_coords.shape}")

# Root pseudotime at centroid of lowest-NAS samples
# Find samples with NAS 0-1 as the root population
low_nas_mask = (nas_score >= 0) & (nas_score <= 1)
n_low = low_nas_mask.sum()
print(f"  Root population (NAS 0-1): {n_low} samples")

if n_low >= 3:
    root_centroid = diffusion_coords[low_nas_mask].mean(axis=0)
else:
    # Fallback: use all controls (is_disease == 0 from metadata)
    is_disease = meta["is_disease"].values.astype(int)
    control_mask = is_disease == 0
    root_centroid = diffusion_coords[control_mask].mean(axis=0)
    print(f"  (Fallback) Using {control_mask.sum()} control samples as root")

# Pseudotime = Euclidean distance from root in diffusion space
pseudotime_raw = np.sqrt(((diffusion_coords - root_centroid) ** 2).sum(axis=1))

# Normalize to [0, 1]
pt_min, pt_max = pseudotime_raw.min(), pseudotime_raw.max()
if pt_max > pt_min:
    pseudotime = (pseudotime_raw - pt_min) / (pt_max - pt_min)
else:
    pseudotime = np.zeros_like(pseudotime_raw)

print(f"  Pseudotime range: [{pseudotime.min():.4f}, {pseudotime.max():.4f}]")

# Check orientation: pseudotime should correlate with NAS
rho, pval = _spearmanr(pseudotime[has_nas], nas_score[has_nas].astype(np.float64))
print(f"  Spearman(pseudotime, NAS): rho={rho:.3f}, p={pval:.2e}")

if rho < 0:
    print("  Flipping pseudotime orientation (negative correlation detected)")
    pseudotime = 1.0 - pseudotime
    rho = -rho


# ---------------------------------------------------------------------------
# STEP 3: Calibrate pseudotime -> NAS via isotonic regression
# ---------------------------------------------------------------------------
print("\n=== STEP 3: Isotonic Calibration (full dataset) ===")

# Fit isotonic regression on all labeled samples
iso_full = IsotonicRegression(y_min=0, y_max=8, increasing=True, out_of_bounds="clip")
iso_full.fit(pseudotime[has_nas], nas_score[has_nas])

# Predict NAS from pseudotime for all samples
nas_calibrated_full = iso_full.predict(pseudotime)
# Round to nearest integer for ordinal NAS
nas_pred_full = np.round(nas_calibrated_full).astype(int)
nas_pred_full = np.clip(nas_pred_full, 0, 8)

# Metrics on labeled samples (in-sample, for reference only)
y_true_all = nas_score[has_nas]
y_pred_all = nas_pred_full[has_nas]

qwk_insample = cohen_kappa_score(y_true_all, y_pred_all, weights="quadratic")
mae_insample = mean_absolute_error(y_true_all, y_pred_all)
print(f"  In-sample: QWK={qwk_insample:.3f}, MAE={mae_insample:.3f}")


# ---------------------------------------------------------------------------
# STEP 4: 5-fold NAS LOCO cross-validation
# ---------------------------------------------------------------------------
print("\n=== STEP 4: 5-Fold NAS LOCO Cross-Validation ===")

fold_results = []
all_loo_preds = np.full(N_SAMPLES, np.nan)

for fold_name in unique_nas_folds:
    test_mask = (loco_nas == fold_name) & has_nas
    train_mask = (loco_nas != fold_name) & (loco_nas != "excluded") & \
                 (loco_nas != "NA") & has_nas

    n_train = train_mask.sum()
    n_test = test_mask.sum()

    if n_test < 5 or n_train < 10:
        print(f"  Fold {fold_name}: skipped (train={n_train}, test={n_test})")
        continue

    # Fit isotonic on training folds
    iso_fold = IsotonicRegression(
        y_min=0, y_max=8, increasing=True, out_of_bounds="clip"
    )
    iso_fold.fit(pseudotime[train_mask], nas_score[train_mask])

    # Predict on test fold
    nas_cal_test = iso_fold.predict(pseudotime[test_mask])
    nas_pred_test = np.round(nas_cal_test).astype(int)
    nas_pred_test = np.clip(nas_pred_test, 0, 8)

    y_test = nas_score[test_mask]

    qwk = cohen_kappa_score(y_test, nas_pred_test, weights="quadratic")
    mae = mean_absolute_error(y_test, nas_pred_test)
    ba = balanced_accuracy_score(y_test, nas_pred_test)

    # Also compute NAS 4-group metrics
    # Map NAS -> 4-group: 0-2=0, 3-4=1, 5-6=2, 7-8=3
    def nas_to_4group(nas):
        g = np.full_like(nas, -1)
        g[(nas >= 0) & (nas <= 2)] = 0
        g[(nas >= 3) & (nas <= 4)] = 1
        g[(nas >= 5) & (nas <= 6)] = 2
        g[(nas >= 7)] = 3
        return g

    y_test_4g = nas_to_4group(y_test)
    y_pred_4g = nas_to_4group(nas_pred_test)
    valid_4g = (y_test_4g >= 0) & (y_pred_4g >= 0)

    qwk_4g = np.nan
    if valid_4g.sum() >= 2:
        qwk_4g = cohen_kappa_score(y_test_4g[valid_4g], y_pred_4g[valid_4g],
                                    weights="quadratic")

    fold_results.append({
        "fold": fold_name,
        "n_train": int(n_train),
        "n_test": int(n_test),
        "qwk_ordinal": qwk,
        "mae": mae,
        "balanced_accuracy": ba,
        "qwk_4group": qwk_4g,
    })

    # Store predictions for held-out samples
    all_loo_preds[test_mask] = nas_pred_test

    print(f"  Fold {fold_name}: QWK={qwk:.3f}, MAE={mae:.3f}, "
          f"BA={ba:.3f}, QWK_4g={qwk_4g:.3f} (n={n_test})")

# Overall LOCO metrics (pool all held-out predictions)
loo_valid = ~np.isnan(all_loo_preds) & has_nas
if loo_valid.sum() > 0:
    y_all_loo = nas_score[loo_valid]
    p_all_loo = all_loo_preds[loo_valid].astype(int)

    overall_qwk = cohen_kappa_score(y_all_loo, p_all_loo, weights="quadratic")
    overall_mae = mean_absolute_error(y_all_loo, p_all_loo)
    overall_ba = balanced_accuracy_score(y_all_loo, p_all_loo)

    y_4g = nas_to_4group(y_all_loo)
    p_4g = nas_to_4group(p_all_loo)
    valid_4g = (y_4g >= 0) & (p_4g >= 0)
    overall_qwk_4g = cohen_kappa_score(y_4g[valid_4g], p_4g[valid_4g],
                                        weights="quadratic")

    print(f"\n  Overall LOCO: QWK={overall_qwk:.3f}, MAE={overall_mae:.3f}, "
          f"BA={overall_ba:.3f}, QWK_4g={overall_qwk_4g:.3f} "
          f"(n={loo_valid.sum()})")


# ---------------------------------------------------------------------------
# STEP 5: Save results
# ---------------------------------------------------------------------------
print("\n=== STEP 5: Saving Results ===")

# Per-fold results
fold_df = pd.DataFrame(fold_results)
fold_out = os.path.join(OUTDIR, "pseudotime_nas_results.csv")
fold_df.to_csv(fold_out, index=False)
print(f"  Per-fold results: {fold_out} ({len(fold_df)} rows)")

# Full calibration table (all samples)
cal_df = pd.DataFrame({
    "sample_id": sample_ids,
    "dataset": datasets,
    "pseudotime": pseudotime,
    "nas_score_true": nas_score,
    "nas_calibrated": nas_calibrated_full,
    "nas_predicted": nas_pred_full,
    "loco_fold_nas": loco_nas,
    "has_nas_label": has_nas.astype(int),
})
cal_out = os.path.join(OUTDIR, "pseudotime_calibration.csv")
cal_df.to_csv(cal_out, index=False)
print(f"  Calibration table: {cal_out} ({len(cal_df)} rows)")

# UMAP + pseudotime data for visualization
umap_df = pd.DataFrame({
    "sample_id": sample_ids,
    "dataset": datasets,
    "umap_1": umap_coords[:, 0],
    "umap_2": umap_coords[:, 1],
    "pseudotime": pseudotime,
    "nas_score": nas_score,
    "nas_4group": nas_group4,
    "nas_calibrated": nas_calibrated_full,
})
umap_out = os.path.join(OUTDIR, "pseudotime_umap_data.csv")
umap_df.to_csv(umap_out, index=False)
print(f"  UMAP+pseudotime: {umap_out} ({len(umap_df)} rows)")


# ---------------------------------------------------------------------------
# STEP 6: Visualization
# ---------------------------------------------------------------------------
print("\n=== STEP 6: Visualization ===")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Panel A: UMAP colored by pseudotime
    ax = axes[0, 0]
    sc = ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
                    c=pseudotime, cmap="viridis", s=8, alpha=0.6)
    plt.colorbar(sc, ax=ax, label="Pseudotime")
    ax.set_title("t-SNE — Pseudotime")
    ax.set_xlabel("Dim 1")
    ax.set_ylabel("Dim 2")

    # Panel B: UMAP colored by NAS score (labeled samples only)
    ax = axes[0, 1]
    unlabeled = ~has_nas
    if unlabeled.any():
        ax.scatter(umap_coords[unlabeled, 0], umap_coords[unlabeled, 1],
                   c="lightgrey", s=5, alpha=0.2, label="Unlabeled")
    if has_nas.any():
        sc2 = ax.scatter(umap_coords[has_nas, 0], umap_coords[has_nas, 1],
                         c=nas_score[has_nas], cmap="RdYlBu_r", s=12, alpha=0.7,
                         vmin=0, vmax=8)
        plt.colorbar(sc2, ax=ax, label="NAS Score")
    ax.set_title("t-SNE — NAS Score")
    ax.set_xlabel("Dim 1")
    ax.set_ylabel("Dim 2")

    # Panel C: Pseudotime vs NAS (labeled samples, with isotonic fit)
    ax = axes[1, 0]
    jitter = np.random.RandomState(42).normal(0, 0.15, size=has_nas.sum())
    ax.scatter(pseudotime[has_nas], nas_score[has_nas] + jitter,
               c="steelblue", s=12, alpha=0.4)
    # Isotonic curve
    pt_sorted = np.linspace(pseudotime[has_nas].min(),
                            pseudotime[has_nas].max(), 200)
    nas_curve = iso_full.predict(pt_sorted)
    ax.plot(pt_sorted, nas_curve, "r-", linewidth=2, label="Isotonic fit")
    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("NAS Score (jittered)")
    ax.set_title(f"Pseudotime vs NAS (rho={rho:.3f})")
    ax.legend()

    # Panel D: Calibrated prediction vs true (LOCO held-out)
    ax = axes[1, 1]
    if loo_valid.sum() > 0:
        jitter2 = np.random.RandomState(42).normal(0, 0.15, size=loo_valid.sum())
        jitter3 = np.random.RandomState(43).normal(0, 0.15, size=loo_valid.sum())
        ax.scatter(y_all_loo + jitter2, p_all_loo + jitter3,
                   c="steelblue", s=12, alpha=0.4)
        ax.plot([0, 8], [0, 8], "k--", alpha=0.5, linewidth=1)
        ax.set_xlabel("True NAS")
        ax.set_ylabel("Predicted NAS (LOCO)")
        ax.set_title(f"LOCO Predictions: QWK={overall_qwk:.3f}, MAE={overall_mae:.2f}")
    else:
        ax.text(0.5, 0.5, "No LOCO predictions", transform=ax.transAxes,
                ha="center", va="center")
        ax.set_title("LOCO Predictions")

    plt.tight_layout()
    fig_path = os.path.join(OUTDIR, "pseudotime_nas_visualization.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {fig_path}")

except Exception as e:
    print(f"  Visualization failed (non-fatal): {e}")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
dt = time.time() - t0
print(f"\n{'=' * 60}")
print(f"94_pseudotime_nas.py completed in {dt / 60:.1f} min")
print(f"  Spearman(pseudotime, NAS) = {rho:.3f}")
if loo_valid.sum() > 0:
    print(f"  LOCO QWK (ordinal)  = {overall_qwk:.3f}")
    print(f"  LOCO QWK (4-group)  = {overall_qwk_4g:.3f}")
    print(f"  LOCO MAE            = {overall_mae:.3f}")
print(f"\n  Per-fold results:   {fold_out}")
print(f"  Calibration table:  {cal_out}")
print(f"  UMAP+pseudotime:    {umap_out}")
print(f"{'=' * 60}")
