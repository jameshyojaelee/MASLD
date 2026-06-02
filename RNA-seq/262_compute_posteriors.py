#!/usr/bin/env python3
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_262_posteriors
#SBATCH --output=logs/net_262_posteriors_%j.out
#SBATCH --error=logs/net_262_posteriors_%j.err
# ===========================================================================
# Script 262: Compute Bayesian posterior P_m(i,j) for every edge in every layer
# ===========================================================================
# Purpose: For each layer, compute posterior probability that each edge is
#          a true signal (not null) using the two-group mixture model:
#          P_m = 1 - pi0 * f0(s) / f(s), clipped to [0, 1].
#
# Input:
#   - RNA-seq/results/network/edges_{name}.csv (raw edges from Stage 1)
#   - RNA-seq/results/network/null_distributions/layer_{name}_density.pkl
#
# Output:
#   - RNA-seq/results/network/posterior_edges/edges_{name}_posterior.csv
#     Columns: gene_a, gene_b, raw_score, posterior, metadata
#
# Environment: spatial (scipy, numpy, pandas, torch, pyarrow)
# ===========================================================================

import os
import pickle
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

warnings.filterwarnings("ignore", category=FutureWarning)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
EDGE_DIR = BASE / "RNA-seq" / "results" / "network"
NULL_DIR = EDGE_DIR / "null_distributions"
OUTDIR = EDGE_DIR / "posterior_edges"
OUTDIR.mkdir(parents=True, exist_ok=True)

LAYER_NAMES = [
    "ppi", "coexpr", "regulon", "lr", "genetic",
    "pathway", "spatial", "cosmos", "cerna", "xspecies"
]

GPU_LAYERS = {"coexpr", "spatial"}
GPU_THRESHOLD = 500000

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
    log("No GPU available, using CPU for all layers")
    return None


def compute_posterior_cpu(scores, f0_kde, f_obs_kde, pi0):
    """Compute posteriors on CPU using scipy KDE evaluation."""
    n = len(scores)
    posteriors = np.zeros(n, dtype=np.float64)

    chunk_size = 50000
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        s = scores[start:end]

        f0_vals = f0_kde(s)
        f_vals = f_obs_kde(s)

        f_vals = np.maximum(f_vals, 1e-300)

        null_fraction = pi0 * f0_vals / f_vals
        posteriors[start:end] = 1.0 - null_fraction

    posteriors = np.clip(posteriors, 0.0, 1.0)
    return posteriors


def compute_posterior_gpu(scores, f0_kde, f_obs_kde, pi0, device):
    """Compute posteriors on GPU for large edge sets."""
    import torch

    f0_eval = f0_kde(scores)
    f_obs_eval = f_obs_kde(scores)

    f0_t = torch.tensor(f0_eval, dtype=torch.float64, device=device)
    f_obs_t = torch.tensor(f_obs_eval, dtype=torch.float64, device=device)

    f_obs_t = torch.clamp(f_obs_t, min=1e-300)

    null_fraction = pi0 * f0_t / f_obs_t
    posteriors = 1.0 - null_fraction
    posteriors = torch.clamp(posteriors, 0.0, 1.0)

    return posteriors.cpu().numpy()


def process_layer(name, device):
    """Process one layer: load density model, compute posteriors, save augmented parquet."""
    t_layer = time.time()
    log(f"\n--- Layer: {name} ---")

    edge_path = EDGE_DIR / f"edges_{name}.csv"
    density_path = NULL_DIR / f"layer_{name}_density.pkl"

    if not edge_path.exists():
        log(f"  Edge parquet not found: {edge_path}, skipping")
        return None
    if not density_path.exists():
        log(f"  Density pickle not found: {density_path}, skipping")
        return None

    edges = pd.read_csv(edge_path)
    log(f"  Loaded {len(edges):,} edges")

    with open(density_path, "rb") as f:
        density = pickle.load(f)

    f0_kde = density["f0"]
    f_obs_kde = density["f_obs"]
    pi0 = density["pi0"]
    log(f"  pi0 = {pi0:.4f}")

    if f0_kde is None or f_obs_kde is None:
        log("  WARNING: KDE models unavailable, assigning uniform posterior = 0.5")
        edges["posterior"] = 0.5
    else:
        scores = edges["raw_score"].values.astype(np.float64)

        use_gpu = (device is not None and name in GPU_LAYERS and len(scores) > GPU_THRESHOLD)

        if use_gpu:
            log(f"  Computing posteriors on GPU ({len(scores):,} edges)...")
            posteriors = compute_posterior_gpu(scores, f0_kde, f_obs_kde, pi0, device)
        else:
            log(f"  Computing posteriors on CPU ({len(scores):,} edges)...")
            posteriors = compute_posterior_cpu(scores, f0_kde, f_obs_kde, pi0)

        edges["posterior"] = posteriors

    out_df = edges[["gene_a", "gene_b", "raw_score", "posterior", "metadata"]].copy()

    out_path = OUTDIR / f"edges_{name}_posterior.csv"
    out_df.to_parquet(out_path, index=False)
    log(f"  Saved: {out_path}")

    n_sig_05 = (out_df["posterior"] > 0.5).sum()
    n_sig_09 = (out_df["posterior"] > 0.9).sum()
    mean_post = out_df["posterior"].mean()
    median_post = out_df["posterior"].median()

    stats = {
        "layer": name,
        "n_edges": len(out_df),
        "pi0": pi0,
        "mean_posterior": float(mean_post),
        "median_posterior": float(median_post),
        "n_posterior_gt_05": int(n_sig_05),
        "n_posterior_gt_09": int(n_sig_09),
        "frac_gt_05": float(n_sig_05 / len(out_df)) if len(out_df) > 0 else 0.0,
        "frac_gt_09": float(n_sig_09 / len(out_df)) if len(out_df) > 0 else 0.0,
    }

    elapsed_layer = time.time() - t_layer
    log(f"  Posterior stats: mean={mean_post:.4f}, median={median_post:.4f}")
    log(f"  Edges P>0.5: {n_sig_05:,} ({stats['frac_gt_05']:.3f})")
    log(f"  Edges P>0.9: {n_sig_09:,} ({stats['frac_gt_09']:.3f})")
    log(f"  Layer {name} done in {elapsed_layer:.1f}s")

    return stats


def main():
    log("=== 262: Compute Bayesian Posteriors ===")
    log(f"Edge directory: {EDGE_DIR}")
    log(f"Density directory: {NULL_DIR}")
    log(f"Output directory: {OUTDIR}")

    device = get_device()

    summary_rows = []
    layers_processed = 0
    layers_skipped = 0

    for name in LAYER_NAMES:
        stats = process_layer(name, device)
        if stats is not None:
            summary_rows.append(stats)
            layers_processed += 1
        else:
            layers_skipped += 1

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = OUTDIR / "posterior_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        log(f"\nSummary saved: {summary_path}")
        log("\n" + summary_df.to_string(index=False))

    elapsed_total = time.time() - t0_global
    log(f"\n=== 262 Complete ===")
    log(f"Layers processed: {layers_processed}, skipped: {layers_skipped}")
    log(f"Total time: {elapsed_total:.1f}s ({elapsed_total / 60:.1f} min)")


if __name__ == "__main__":
    main()
