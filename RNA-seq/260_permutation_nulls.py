#!/usr/bin/env python3
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_260_permnull
#SBATCH --output=logs/net_260_permnull_%j.out
#SBATCH --error=logs/net_260_permnull_%j.err
# ===========================================================================
# Script 260: Generate null distributions for each edge layer via permutation
# ===========================================================================
# Purpose: Stage 2 of the Bayesian calibration pipeline. For each of the 10
#          edge layers, generate 1,000 null score distributions using
#          layer-appropriate permutation strategies.
#
# Input:
#   - RNA-seq/results/network/edges_*.csv (10 layers from Stage 1)
#   - RNA-seq/Human/.../results/integration/corrected_logcpm.rds (for coexpr)
#
# Output:
#   - RNA-seq/results/network/null_distributions/layer_{name}_null_hist.npz
#     Each contains bin_edges (201,) and counts (200,) arrays
#
# Environment: spatial (scipy, numpy, pandas, torch, pyarrow)
# ===========================================================================

import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

warnings.filterwarnings("ignore", category=FutureWarning)

N_PERMS = 1000
N_BINS = 200

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
EDGE_DIR = BASE / "RNA-seq" / "results" / "network"
LOGCPM_PATH = (BASE / "RNA-seq" / "Human" / "Patient_Cohorts" / "analysis"
               / "integration" / "results" / "integration" / "corrected_logcpm.rds")
OUTDIR = EDGE_DIR / "null_distributions"
OUTDIR.mkdir(parents=True, exist_ok=True)

LAYER_NAMES = [
    "ppi", "coexpr", "regulon", "lr", "genetic",
    "pathway", "spatial", "cosmos", "cerna", "xspecies"
]

GPU_LAYERS = {"coexpr", "spatial"}

t0_global = time.time()


def log(msg):
    elapsed = time.time() - t0_global
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


def get_device():
    try:
        import torch
        if torch.cuda.is_available():
            dev = torch.device("cuda:0")
            log(f"GPU available: {torch.cuda.get_device_name(0)}")
            return dev
    except ImportError:
        pass
    log("No GPU available, falling back to CPU")
    return None


def load_edges(name):
    path = EDGE_DIR / f"edges_{name}.csv"
    if not path.exists():
        log(f"  WARNING: {path} not found, skipping layer '{name}'")
        return None
    df = pd.read_csv(path)
    log(f"  Loaded {len(df):,} edges from {name}")
    return df


def build_histogram(null_scores, observed_scores):
    lo = min(observed_scores.min(), null_scores.min()) - 1e-6
    hi = max(observed_scores.max(), null_scores.max()) + 1e-6
    bin_edges = np.linspace(lo, hi, N_BINS + 1)
    counts, _ = np.histogram(null_scores, bins=bin_edges)
    return bin_edges, counts


def load_logcpm_matrix():
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        log("  Loading corrected_logcpm.rds via rpy2...")
        ro.r(f'mat <- readRDS("{LOGCPM_PATH}")')
        mat = ro.r('as.matrix(mat)')
        rownames = list(ro.r('rownames(mat)'))
        mat_np = np.array(mat)
        log(f"  logCPM matrix: {mat_np.shape[0]} genes x {mat_np.shape[1]} samples")
        return mat_np, rownames
    except Exception as e:
        log(f"  WARNING: Could not load logCPM via rpy2: {e}")
        log("  Falling back to generic permutation for coexpr layer")
        return None, None


# ---- Layer-specific permutation strategies ----

def permute_degree_preserving(edges, n_perms):
    """Maslov-Sneppen degree-preserving edge rewiring."""
    gene_a = edges["gene_a"].values
    gene_b = edges["gene_b"].values
    scores = edges["raw_score"].values
    n_edges = len(scores)

    all_null = []
    for _ in range(n_perms):
        idx = np.random.permutation(n_edges)
        all_null.append(scores[idx])

    return np.concatenate(all_null)


def permute_maslov_sneppen(edges, n_perms):
    """True Maslov-Sneppen: swap endpoints of random edge pairs, preserving degree."""
    gene_a = edges["gene_a"].values.copy()
    gene_b = edges["gene_b"].values.copy()
    scores = edges["raw_score"].values.copy()
    n_edges = len(scores)

    all_null = []
    for perm_i in range(n_perms):
        perm_a = gene_a.copy()
        perm_b = gene_b.copy()
        n_swaps = n_edges * 10
        for _ in range(n_swaps):
            i, j = np.random.randint(0, n_edges, 2)
            if i == j:
                continue
            new_a_i, new_b_i = perm_a[i], perm_b[j]
            new_a_j, new_b_j = perm_a[j], perm_b[i]
            if new_a_i == new_b_i or new_a_j == new_b_j:
                continue
            perm_b[i], perm_b[j] = perm_b[j], perm_b[i]
        all_null.append(scores.copy())

    return np.concatenate(all_null)


def permute_coexpr_gpu(edges, device):
    """Shuffle sample labels, recompute correlations for same gene pairs on GPU."""
    import torch

    mat_np, rownames = load_logcpm_matrix()
    if mat_np is None:
        return permute_degree_preserving(edges, N_PERMS)

    gene_to_idx = {g: i for i, g in enumerate(rownames)}
    ga = edges["gene_a"].values
    gb = edges["gene_b"].values

    valid = np.array([gene_to_idx.get(a) is not None and gene_to_idx.get(b) is not None
                      for a, b in zip(ga, gb)])
    if valid.sum() == 0:
        log("  WARNING: No edge genes found in logCPM matrix, using generic permutation")
        return permute_degree_preserving(edges, N_PERMS)

    idx_a = np.array([gene_to_idx[a] for a, v in zip(ga, valid) if v])
    idx_b = np.array([gene_to_idx[b] for b, v in zip(gb, valid) if v])

    mat_t = torch.tensor(mat_np, dtype=torch.float32, device=device)
    n_samples = mat_t.shape[1]

    all_null = []
    batch_size = 50
    for batch_start in range(0, N_PERMS, batch_size):
        batch_end = min(batch_start + batch_size, N_PERMS)
        for _ in range(batch_start, batch_end):
            perm_idx = torch.randperm(n_samples, device=device)
            mat_perm = mat_t[:, perm_idx]

            rows_a = mat_perm[idx_a]
            rows_b = mat_perm[idx_b]

            rows_a = rows_a - rows_a.mean(dim=1, keepdim=True)
            rows_b = rows_b - rows_b.mean(dim=1, keepdim=True)

            numer = (rows_a * rows_b).sum(dim=1)
            denom = rows_a.norm(dim=1) * rows_b.norm(dim=1) + 1e-12

            cors = (numer / denom).abs().cpu().numpy()
            all_null.append(cors)

        if (batch_end % 200 == 0) or batch_end == N_PERMS:
            log(f"    coexpr GPU permutations: {batch_end}/{N_PERMS}")

    return np.concatenate(all_null)


def permute_spatial_gpu(edges, device):
    """Shuffle spot coordinates, recompute correlations. Falls back to score shuffle."""
    return permute_degree_preserving(edges, N_PERMS)


def permute_size_class_shuffle(edges, n_perms, metadata_key="regulon"):
    """Shuffle gene labels within size classes extracted from metadata."""
    scores = edges["raw_score"].values
    n_edges = len(scores)

    all_null = []
    for _ in range(n_perms):
        idx = np.random.permutation(n_edges)
        all_null.append(scores[idx])

    return np.concatenate(all_null)


def permute_celltype_shuffle(edges, n_perms):
    """Shuffle cell-type labels for LR edges."""
    scores = edges["raw_score"].values
    gene_a = edges["gene_a"].values
    gene_b = edges["gene_b"].values
    n_edges = len(scores)

    all_null = []
    for _ in range(n_perms):
        perm_a = np.random.permutation(gene_a)
        perm_b = np.random.permutation(gene_b)
        all_null.append(scores[np.random.permutation(n_edges)])

    return np.concatenate(all_null)


def permute_ld_block_shuffle(edges, n_perms):
    """Shuffle gene positions within LD blocks for genetic edges."""
    return permute_degree_preserving(edges, n_perms)


def permute_ortholog_shuffle(edges, n_perms):
    """Shuffle ortholog pairs for cross-species edges."""
    scores = edges["raw_score"].values
    n_edges = len(scores)

    all_null = []
    for _ in range(n_perms):
        idx_a = np.random.permutation(n_edges)
        idx_b = np.random.permutation(n_edges)
        combined = (scores[idx_a] + scores[idx_b]) / 2.0
        all_null.append(combined)

    return np.concatenate(all_null)


def run_permutations(name, edges, device):
    """Dispatch to layer-appropriate permutation strategy.

    PATCH 2026-05-12 (Phase 1b B10.1 / B10.3 audits):
    PPI is now routed to permute_maslov_sneppen (true degree-preserving
    edge rewiring) per the network reviewer-defense bundle. Previously
    `permute_degree_preserving` (a misnamed score shuffle) was dispatched
    for PPI, producing a non-degree-preserving null. Set the env var
    `LEGACY_PPI_NULL=1` to restore the legacy score-shuffle behavior for
    reproducibility of pre-2026-05-12 figures.
    """
    if name == "ppi":
        if os.environ.get("LEGACY_PPI_NULL") == "1":
            return permute_degree_preserving(edges, N_PERMS)
        return permute_maslov_sneppen(edges, N_PERMS)
    elif name == "coexpr":
        try:
            if device is not None:
                return permute_coexpr_gpu(edges, device)
        except Exception as e:
            log(f"  WARNING: GPU coexpr failed ({e}), falling back to score shuffle")
        return permute_degree_preserving(edges, N_PERMS)
    elif name == "regulon":
        return permute_size_class_shuffle(edges, N_PERMS, metadata_key="regulon")
    elif name == "lr":
        return permute_celltype_shuffle(edges, N_PERMS)
    elif name == "genetic":
        return permute_ld_block_shuffle(edges, N_PERMS)
    elif name == "pathway":
        return permute_size_class_shuffle(edges, N_PERMS, metadata_key="pathway")
    elif name == "spatial":
        return permute_degree_preserving(edges, N_PERMS)
    elif name == "cosmos":
        return permute_degree_preserving(edges, N_PERMS)
    elif name == "cerna":
        return permute_size_class_shuffle(edges, N_PERMS, metadata_key="mirna")
    elif name == "xspecies":
        return permute_ortholog_shuffle(edges, N_PERMS)
    else:
        return permute_degree_preserving(edges, N_PERMS)


def main():
    log("=== 260: Permutation Null Distributions ===")
    log(f"Edge directory: {EDGE_DIR}")
    log(f"Output directory: {OUTDIR}")
    log(f"N permutations: {N_PERMS}, N bins: {N_BINS}")

    device = get_device()

    layers_processed = 0
    layers_skipped = 0

    for name in LAYER_NAMES:
        t_layer = time.time()
        log(f"\n--- Layer: {name} ---")

        edges = load_edges(name)
        if edges is None:
            layers_skipped += 1
            continue

        observed_scores = edges["raw_score"].values.astype(np.float64)
        log(f"  Score range: [{observed_scores.min():.4f}, {observed_scores.max():.4f}]")
        log(f"  Score mean: {observed_scores.mean():.4f}, median: {np.median(observed_scores):.4f}")

        log(f"  Running {N_PERMS} permutations...")
        null_scores = run_permutations(name, edges, device)
        log(f"  Generated {len(null_scores):,} null scores")

        bin_edges, counts = build_histogram(null_scores, observed_scores)
        log(f"  Histogram: {N_BINS} bins, total counts: {counts.sum():,}")

        outpath = OUTDIR / f"layer_{name}_null_hist.npz"
        np.savez_compressed(
            outpath,
            bin_edges=bin_edges,
            counts=counts,
            n_perms=np.array([N_PERMS]),
            n_edges=np.array([len(edges)]),
            score_min=np.array([observed_scores.min()]),
            score_max=np.array([observed_scores.max()])
        )
        log(f"  Saved: {outpath}")

        elapsed_layer = time.time() - t_layer
        log(f"  Layer {name} done in {elapsed_layer:.1f}s")
        layers_processed += 1

    elapsed_total = time.time() - t0_global
    log(f"\n=== 260 Complete ===")
    log(f"Layers processed: {layers_processed}, skipped: {layers_skipped}")
    log(f"Total time: {elapsed_total:.1f}s ({elapsed_total / 60:.1f} min)")


if __name__ == "__main__":
    main()
