#!/usr/bin/env python3
"""
162_loco_nmf_subtypes.py
LOCO-NMF Subtypes: recomputes NMF on training folds only to break label
circularity in disease subtype assignments.

For each of the 6 LOCO fibrosis folds:
  1. Train NMF (k=2, n_init=30) on Disease samples from training folds only
  2. Identify S1/S2 by mean fibrosis stage (S2 = higher fibrosis)
  3. Project held-out test samples via NNLS into the NMF basis
  4. Assign test samples to S1/S2 via nearest centroid in H-space
  5. Compute stability metrics (silhouette, centroid correlation, global concordance)

Input:
  - results/integration/merged_dge.rds                (full DGEList, 1,444 samples)
  - results/staging_classifier/modeling_metadata.csv   (fold assignments)
  - results/subtypes/nmf_assignments.csv               (global NMF labels)
  - metadata/unified_metadata.csv                       (for QC pass_technical)
  - qc/sample_qc_report.csv                            (pass_technical flag)

Output (all to results/prognosis_v2/loco_nmf_labels/):
  - fold_{cohort}_labels.csv          per-fold test sample labels
  - loco_nmf_combined.csv             all samples with LOCO-NMF + global S1/S2
  - loco_nmf_stability.csv            stability metrics per fold

SLURM: io partition, 8 CPUs, 32G RAM, 48h
Env:   micromamba activate spatial
"""

import os
import sys
import time
import subprocess
import tempfile
import warnings
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.spatial.distance import cdist
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SEED = 42
np.random.seed(SEED)

# --- Paths ---
BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
INTEGRATION = BASE / "RNA-seq/Human/Patient_Cohorts/analysis/integration"

DGE_PATH       = INTEGRATION / "results/integration/merged_dge.rds"
META_PATH      = INTEGRATION / "results/staging_classifier/modeling_metadata.csv"
GLOBAL_NMF     = BASE / "RNA-seq/results/subtypes/nmf_assignments.csv"
UNIFIED_META   = INTEGRATION / "metadata/unified_metadata.csv"
QC_PATH        = INTEGRATION / "qc/sample_qc_report.csv"

OUT_DIR = INTEGRATION / "results/prognosis_v2/loco_nmf_labels"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# NMF parameters (match original Script 44)
N_TOP_GENES = 5000
K = 2
N_INIT = 30
MAX_ITER = 500

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S"
)
log = logging.getLogger(__name__)


# =============================================================================
# Step 0: Extract batch-corrected logCPM from merged_dge.rds via R subprocess
# =============================================================================
def extract_logcpm_via_r(dge_path, unified_meta_path, qc_path, out_h5):
    """Call R to load DGEList, apply TMM+voom+batch correction, save to HDF5."""
    log.info("Extracting batch-corrected logCPM from DGEList via R...")

    r_script = f"""
suppressPackageStartupMessages({{
  library(edgeR)
  library(limma)
  library(rhdf5)
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
}})

dge <- readRDS("{dge_path}")
cat(sprintf("DGEList: %d genes x %d samples\\n", nrow(dge), ncol(dge)))

# Load metadata + QC to identify Disease + QC-passing samples
meta <- fread("{unified_meta_path}")
qc   <- fread("{qc_path}")
meta  <- meta %>%
  left_join(qc %>% dplyr::select(sample_id, pass_technical), by = "sample_id")

disease_samples <- meta %>%
  dplyr::filter(group_binary == "Disease" & pass_technical == TRUE) %>%
  dplyr::pull(sample_id)

available <- intersect(disease_samples, colnames(dge))
cat(sprintf("Disease + QC-passing + in DGEList: %d samples\\n", length(available)))

# Subset
dge_sub <- dge[, available]
meta_sub <- meta %>% dplyr::filter(sample_id %in% available)

# TMM + voom
dge_sub <- calcNormFactors(dge_sub, method = "TMM")
v <- voom(dge_sub, design = NULL, plot = FALSE)
logcpm <- v$E

# Batch correction (residualize dataset, keeping overall mean)
logcpm_corrected <- removeBatchEffect(logcpm, batch = meta_sub$dataset)

cat(sprintf("Corrected logCPM: %d genes x %d samples\\n",
            nrow(logcpm_corrected), ncol(logcpm_corrected)))

# Save to HDF5
out_file <- "{out_h5}"
if (file.exists(out_file)) file.remove(out_file)
h5createFile(out_file)

h5createDataset(out_file, "logcpm",
                dims = dim(logcpm_corrected),
                chunk = c(min(1000, nrow(logcpm_corrected)), ncol(logcpm_corrected)),
                level = 4)
h5write(logcpm_corrected, out_file, "logcpm")
h5write(rownames(logcpm_corrected), out_file, "gene_ids")
h5write(colnames(logcpm_corrected), out_file, "sample_ids")
h5write(meta_sub$dataset, out_file, "datasets")

cat(sprintf("Saved to %s\\n", out_file))
"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as f:
        f.write(r_script)
        r_tmp = f.name

    try:
        result = subprocess.run(
            ["micromamba", "run", "-n", "rnaseq", "Rscript", r_tmp],
            capture_output=True, text=True, timeout=600
        )
        if result.returncode != 0:
            log.error("R extraction failed:\n%s", result.stderr)
            sys.exit(1)
        log.info("R output:\n%s", result.stdout)
    finally:
        os.unlink(r_tmp)


def load_logcpm(h5_path):
    """Load the batch-corrected logCPM matrix from HDF5.
    R (rhdf5) stores matrices column-major; h5py reads row-major, so the
    matrix appears transposed. We transpose back to (genes, samples)."""
    import h5py
    with h5py.File(h5_path, "r") as f:
        raw = f["logcpm"][:]  # h5py reads as (samples, genes) due to R col-major
        gene_ids = [g.decode() if isinstance(g, bytes) else g for g in f["gene_ids"][:]]
        sample_ids = [s.decode() if isinstance(s, bytes) else s for s in f["sample_ids"][:]]
    # Transpose to (genes, samples) to match R convention
    logcpm = raw.T
    log.info("Loaded logCPM: %d genes x %d samples", logcpm.shape[0], logcpm.shape[1])
    assert logcpm.shape[0] == len(gene_ids), (
        f"Gene count mismatch: matrix has {logcpm.shape[0]} rows but {len(gene_ids)} gene_ids"
    )
    assert logcpm.shape[1] == len(sample_ids), (
        f"Sample count mismatch: matrix has {logcpm.shape[1]} cols but {len(sample_ids)} sample_ids"
    )
    return logcpm, gene_ids, sample_ids


# =============================================================================
# Step 1: LOCO-NMF
# =============================================================================
def select_variable_genes(logcpm, n_top=5000):
    """Select top n variable genes by IQR (gene axis=0, sample axis=1).
    logcpm shape: (genes, samples)."""
    q75 = np.percentile(logcpm, 75, axis=1)
    q25 = np.percentile(logcpm, 25, axis=1)
    iqr = q75 - q25
    top_idx = np.argsort(iqr)[::-1][:n_top]
    return np.sort(top_idx)


def shift_nonneg(mat):
    """Shift each gene (row) so minimum is 0."""
    row_min = mat.min(axis=1, keepdims=True)
    shifted = mat - row_min
    shifted[shifted < 0] = 0
    return shifted


def run_nmf_train(mat_nn, k=2, n_init=30, max_iter=500, seed=42):
    """Run NMF on non-negative matrix. Returns (W, H, model).
    mat_nn: (genes, samples) -> sklearn expects (samples, features).
    So we transpose: X = mat_nn.T  ->  (samples, genes)
    NMF: X ~ W_sk @ H_sk where W_sk=(samples, k), H_sk=(k, genes)
    We return W_bio=(genes, k), H_bio=(k, samples) to match R NMF convention.
    """
    X = mat_nn.T  # (samples, genes)
    model = NMF(
        n_components=k,
        init="nndsvda",
        max_iter=max_iter,
        random_state=seed,
        solver="mu",
        beta_loss="frobenius",
    )
    # Run with multiple inits, keep best
    best_err = np.inf
    best_W = None
    best_H = None
    best_model = None
    for i in range(n_init):
        m = NMF(
            n_components=k,
            init="nndsvda" if i == 0 else "random",
            max_iter=max_iter,
            random_state=seed + i,
            solver="mu",
            beta_loss="frobenius",
        )
        W_sk = m.fit_transform(X)  # (samples, k)
        H_sk = m.components_       # (k, genes)
        err = m.reconstruction_err_
        if err < best_err:
            best_err = err
            best_W = W_sk
            best_H = H_sk
            best_model = m

    # Convert to bioinformatics convention: W_bio=(genes, k), H_bio=(k, samples)
    W_bio = best_H.T   # (genes, k)
    H_bio = best_W.T   # (k, samples)
    log.info("  NMF converged. Reconstruction error: %.4f", best_err)
    return W_bio, H_bio, best_model


def project_nnls(W_bio, test_mat_nn):
    """Project test samples into NMF basis via NNLS.
    W_bio: (genes, k) -- the gene loadings from training
    test_mat_nn: (genes, n_test) -- non-negative test expression
    Returns H_test: (k, n_test) -- component activations for test samples.
    """
    n_test = test_mat_nn.shape[1]
    k = W_bio.shape[1]
    H_test = np.zeros((k, n_test))
    for j in range(n_test):
        h, _ = nnls(W_bio, test_mat_nn[:, j])
        H_test[:, j] = h
    return H_test


def assign_subtypes(H_bio, fibrosis_stages):
    """Assign S1/S2 labels. S2 = component with higher mean fibrosis.
    H_bio: (k, samples), fibrosis_stages: array of numeric fibrosis values.
    Returns labels (array of 'S1'/'S2'), s2_component index.
    """
    k = H_bio.shape[0]
    # Assign each sample to its dominant component
    dominant = np.argmax(H_bio, axis=0)  # component index per sample

    # Mean fibrosis per component
    mean_fib = np.zeros(k)
    for c in range(k):
        mask = dominant == c
        if mask.sum() > 0:
            mean_fib[c] = np.nanmean(fibrosis_stages[mask])

    s2_comp = np.argmax(mean_fib)  # S2 = higher fibrosis
    s1_comp = 1 - s2_comp          # only works for k=2

    labels = np.where(dominant == s2_comp, "S2", "S1")
    log.info("  S1/S2 assignment: S2=component %d (mean fib %.2f), S1=component %d (mean fib %.2f)",
             s2_comp, mean_fib[s2_comp], s1_comp, mean_fib[s1_comp])
    return labels, s2_comp


def assign_test_nearest_centroid(H_train, train_labels, H_test):
    """Assign test samples to S1/S2 via nearest centroid in H-space.
    H_train: (k, n_train), H_test: (k, n_test).
    Returns test_labels: array of 'S1'/'S2'.
    """
    # Compute centroids
    s1_mask = train_labels == "S1"
    s2_mask = train_labels == "S2"
    centroid_s1 = H_train[:, s1_mask].mean(axis=1)
    centroid_s2 = H_train[:, s2_mask].mean(axis=1)

    # Euclidean distance from each test sample to centroids
    # H_test columns are samples
    test_labels = []
    for j in range(H_test.shape[1]):
        h = H_test[:, j]
        d1 = np.linalg.norm(h - centroid_s1)
        d2 = np.linalg.norm(h - centroid_s2)
        test_labels.append("S1" if d1 < d2 else "S2")
    return np.array(test_labels), centroid_s1, centroid_s2


def compute_silhouette(H_bio, labels):
    """Compute mean silhouette score in H-space."""
    if len(set(labels)) < 2:
        return np.nan
    # H_bio: (k, samples) -> transpose to (samples, k) for sklearn
    return silhouette_score(H_bio.T, labels, metric="euclidean")


# =============================================================================
# Main
# =============================================================================
def main():
    t0 = time.time()
    log.info("=== Script 162: LOCO-NMF Subtypes ===")

    # --- Step 0: Extract full logCPM via R ---
    logcpm_h5 = OUT_DIR / "full_logcpm_corrected.h5"
    if not logcpm_h5.exists():
        extract_logcpm_via_r(DGE_PATH, UNIFIED_META, QC_PATH, logcpm_h5)
    else:
        log.info("Using cached logCPM from %s", logcpm_h5)

    logcpm, gene_ids, all_sample_ids = load_logcpm(logcpm_h5)
    # logcpm: (genes, samples)

    # --- Load metadata ---
    meta = pd.read_csv(META_PATH)
    global_nmf = pd.read_csv(GLOBAL_NMF)

    # Build sample->index lookup for the logCPM matrix
    sample_to_idx = {s: i for i, s in enumerate(all_sample_ids)}

    # Identify disease samples present in the logCPM matrix
    disease_in_logcpm = set(all_sample_ids)  # R extraction already filtered to Disease + QC-passing
    log.info("Disease samples in logCPM: %d", len(disease_in_logcpm))

    # Merge metadata with fold info; restrict to samples in logCPM
    meta_valid = meta[
        (meta["loco_fold_fibrosis"] != "excluded") &
        (meta["sample_id"].isin(disease_in_logcpm)) &
        (meta["group_binary"] == "Disease")
    ].copy()
    log.info("Samples with valid fibrosis fold + Disease + in logCPM: %d", len(meta_valid))

    folds = sorted(meta_valid["loco_fold_fibrosis"].unique())
    log.info("LOCO fibrosis folds (%d): %s", len(folds), folds)

    # --- Run LOCO-NMF per fold ---
    all_fold_results = []
    stability_rows = []
    fold_centroids = {}  # fold -> (centroid_s1, centroid_s2) for cross-fold comparison
    fold_gene_loadings = {}  # fold -> (var_gene_idx, W_bio, s2_comp) for cross-fold W comparison

    for fold_name in folds:
        log.info("\n--- Fold: %s (held-out) ---", fold_name)

        # Split train/test
        test_mask = meta_valid["loco_fold_fibrosis"] == fold_name
        train_meta = meta_valid[~test_mask]
        test_meta  = meta_valid[test_mask]
        log.info("  Train: %d samples, Test: %d samples", len(train_meta), len(test_meta))

        # Get sample indices in logCPM matrix
        train_idx = np.array([sample_to_idx[s] for s in train_meta["sample_id"]])
        test_idx  = np.array([sample_to_idx[s] for s in test_meta["sample_id"]])

        # Extract expression submatrices: (genes, n_train), (genes, n_test)
        train_expr = logcpm[:, train_idx]
        test_expr  = logcpm[:, test_idx]

        # Select top 5000 variable genes on TRAINING only
        var_gene_idx = select_variable_genes(train_expr, n_top=N_TOP_GENES)
        train_mat = train_expr[var_gene_idx, :]
        test_mat  = test_expr[var_gene_idx, :]
        log.info("  Selected %d variable genes from training data", len(var_gene_idx))

        # Shift to non-negative (per-gene min from TRAINING data)
        train_min = train_mat.min(axis=1, keepdims=True)
        train_nn = train_mat - train_min
        train_nn[train_nn < 0] = 0
        test_nn = test_mat - train_min   # use training min for test
        test_nn[test_nn < 0] = 0

        # Run NMF on training data
        W_bio, H_train, model = run_nmf_train(
            train_nn, k=K, n_init=N_INIT, max_iter=MAX_ITER, seed=SEED
        )

        # Assign S1/S2 on training data using fibrosis stage
        train_fib = train_meta["fib_stage"].values.astype(float)
        train_labels, s2_comp = assign_subtypes(H_train, train_fib)

        # Project test samples via NNLS
        H_test = project_nnls(W_bio, test_nn)

        # Assign test samples via nearest centroid
        test_labels, c_s1, c_s2 = assign_test_nearest_centroid(
            H_train, train_labels, H_test
        )
        fold_centroids[fold_name] = (c_s1, c_s2)
        fold_gene_loadings[fold_name] = (var_gene_idx, W_bio, s2_comp)

        # Compute stability: silhouette on training data
        train_sil = compute_silhouette(H_train, train_labels)
        test_sil  = compute_silhouette(H_test, test_labels) if len(set(test_labels)) >= 2 else np.nan

        # Subtype distribution
        train_s1 = (train_labels == "S1").sum()
        train_s2 = (train_labels == "S2").sum()
        test_s1  = (test_labels == "S1").sum()
        test_s2  = (test_labels == "S2").sum()

        log.info("  Train: S1=%d, S2=%d (sil=%.3f)", train_s1, train_s2, train_sil)
        log.info("  Test:  S1=%d, S2=%d (sil=%.3f)", test_s1, test_s2,
                 test_sil if not np.isnan(test_sil) else -1)

        # Save per-fold test labels
        fold_df = pd.DataFrame({
            "sample_id": test_meta["sample_id"].values,
            "loco_fold": fold_name,
            "loco_nmf_subtype": test_labels,
            "H_component_0": H_test[0, :],
            "H_component_1": H_test[1, :],
        })
        fold_path = OUT_DIR / f"fold_{fold_name}_labels.csv"
        fold_df.to_csv(fold_path, index=False)
        log.info("  Saved: %s", fold_path)
        all_fold_results.append(fold_df)

        # Mean fibrosis per test subtype
        test_fib = test_meta["fib_stage"].values.astype(float)
        test_s1_fib = np.nanmean(test_fib[test_labels == "S1"]) if test_s1 > 0 else np.nan
        test_s2_fib = np.nanmean(test_fib[test_labels == "S2"]) if test_s2 > 0 else np.nan

        stability_rows.append({
            "fold": fold_name,
            "n_train": len(train_meta),
            "n_test": len(test_meta),
            "train_S1": train_s1,
            "train_S2": train_s2,
            "test_S1": test_s1,
            "test_S2": test_s2,
            "train_silhouette": round(train_sil, 4),
            "test_silhouette": round(test_sil, 4) if not np.isnan(test_sil) else np.nan,
            "test_S1_mean_fib": round(test_s1_fib, 3) if not np.isnan(test_s1_fib) else np.nan,
            "test_S2_mean_fib": round(test_s2_fib, 3) if not np.isnan(test_s2_fib) else np.nan,
            "reconstruction_err": round(model.reconstruction_err_, 4),
        })

    # --- Combine all folds ---
    combined = pd.concat(all_fold_results, ignore_index=True)
    log.info("\nCombined LOCO-NMF labels: %d samples", len(combined))

    # Merge with global NMF labels
    combined = combined.merge(
        global_nmf[["sample_id", "nmf_subtype"]].rename(
            columns={"nmf_subtype": "global_nmf_subtype"}
        ),
        on="sample_id", how="left"
    )

    # Concordance with global NMF
    both = combined.dropna(subset=["global_nmf_subtype"])
    if len(both) > 0:
        concordance = (both["loco_nmf_subtype"] == both["global_nmf_subtype"]).mean()
        log.info("LOCO vs Global NMF concordance: %.3f (%d/%d)",
                 concordance, (both["loco_nmf_subtype"] == both["global_nmf_subtype"]).sum(), len(both))
    else:
        concordance = np.nan

    combined.to_csv(OUT_DIR / "loco_nmf_combined.csv", index=False)
    log.info("Saved: %s", OUT_DIR / "loco_nmf_combined.csv")

    # --- Cross-fold gene loading similarity ---
    # Compare W (gene loadings) across folds using cosine similarity on
    # overlapping genes. This is more meaningful than H-space centroid
    # correlation which is degenerate at k=2.
    from scipy.spatial.distance import cosine as cosine_dist
    fold_names_sorted = sorted(fold_gene_loadings.keys())
    loading_corrs = []
    for i, f1 in enumerate(fold_names_sorted):
        for f2 in fold_names_sorted[i+1:]:
            idx1, W1, s2c1 = fold_gene_loadings[f1]
            idx2, W2, s2c2 = fold_gene_loadings[f2]
            # Find overlapping gene indices
            common = np.intersect1d(idx1, idx2)
            if len(common) < 100:
                loading_corrs.append({
                    "fold_1": f1, "fold_2": f2,
                    "n_common_genes": len(common),
                    "S2_loading_cosine": np.nan,
                    "S1_loading_cosine": np.nan,
                })
                continue
            # Map common genes to positions in each fold's W matrix
            pos1 = np.searchsorted(idx1, common)
            pos2 = np.searchsorted(idx2, common)
            # S2 loading vectors (W column for the S2 component)
            w1_s2 = W1[pos1, s2c1]
            w2_s2 = W2[pos2, s2c2]
            w1_s1 = W1[pos1, 1 - s2c1]
            w2_s1 = W2[pos2, 1 - s2c2]
            cos_s2 = 1.0 - cosine_dist(w1_s2, w2_s2)
            cos_s1 = 1.0 - cosine_dist(w1_s1, w2_s1)
            loading_corrs.append({
                "fold_1": f1, "fold_2": f2,
                "n_common_genes": len(common),
                "S2_loading_cosine": round(cos_s2, 4),
                "S1_loading_cosine": round(cos_s1, 4),
            })

    if loading_corrs:
        cc_df = pd.DataFrame(loading_corrs)
        mean_s1_cos = cc_df["S1_loading_cosine"].mean()
        mean_s2_cos = cc_df["S2_loading_cosine"].mean()
        log.info("Cross-fold gene loading cosine similarity: S1=%.3f, S2=%.3f",
                 mean_s1_cos, mean_s2_cos)
    else:
        mean_s1_cos = np.nan
        mean_s2_cos = np.nan

    # --- Stability metrics ---
    stability = pd.DataFrame(stability_rows)
    stability["global_concordance"] = concordance
    stability["mean_crossfold_S1_loading_cos"] = round(mean_s1_cos, 4) if not np.isnan(mean_s1_cos) else np.nan
    stability["mean_crossfold_S2_loading_cos"] = round(mean_s2_cos, 4) if not np.isnan(mean_s2_cos) else np.nan

    stability.to_csv(OUT_DIR / "loco_nmf_stability.csv", index=False)
    log.info("Saved: %s", OUT_DIR / "loco_nmf_stability.csv")

    if loading_corrs:
        cc_df.to_csv(OUT_DIR / "crossfold_loading_similarity.csv", index=False)
        log.info("Saved: %s", OUT_DIR / "crossfold_loading_similarity.csv")

    # --- Summary ---
    elapsed = time.time() - t0
    log.info("\n=== Summary ===")
    log.info("Folds processed: %d", len(folds))
    log.info("Total LOCO-NMF samples: %d", len(combined))
    log.info("LOCO S1: %d, S2: %d", (combined["loco_nmf_subtype"] == "S1").sum(),
             (combined["loco_nmf_subtype"] == "S2").sum())
    log.info("Global concordance: %.3f", concordance if not np.isnan(concordance) else -1)
    log.info("Mean train silhouette: %.3f", stability["train_silhouette"].mean())
    log.info("Elapsed: %.1f min", elapsed / 60)
    log.info("Done.")


if __name__ == "__main__":
    main()
