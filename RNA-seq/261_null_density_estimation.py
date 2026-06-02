#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_261_nulldens
#SBATCH --output=logs/net_261_nulldens_%j.out
#SBATCH --error=logs/net_261_nulldens_%j.err
# ===========================================================================
# Script 261: Fit smooth null density f0(s) per layer from permutation histograms
# ===========================================================================
# Purpose: For each edge layer, fit a KDE-based null density f0 and observed
#          density f, then estimate pi0 (null proportion) via Storey's method.
#
# Input:
#   - RNA-seq/results/network/null_distributions/layer_{name}_null_hist.npz
#   - RNA-seq/results/network/edges_{name}.csv
#
# Output:
#   - RNA-seq/results/network/null_distributions/layer_{name}_density.pkl
#     Contains: f0 (null KDE), f_obs (observed KDE), pi0, bin_edges,
#               null_counts, obs_counts, obs_scores
#
# Environment: spatial (scipy, numpy, pandas, pyarrow)
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
from scipy.stats import gaussian_kde

warnings.filterwarnings("ignore", category=FutureWarning)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
EDGE_DIR = BASE / "RNA-seq" / "results" / "network"
NULL_DIR = EDGE_DIR / "null_distributions"
NULL_DIR.mkdir(parents=True, exist_ok=True)

LAYER_NAMES = [
    "ppi", "coexpr", "regulon", "lr", "genetic",
    "pathway", "spatial", "cosmos", "cerna", "xspecies"
]

t0_global = time.time()


def log(msg):
    elapsed = time.time() - t0_global
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


def histogram_to_samples(bin_edges, counts, n_samples=50000):
    """Convert histogram back to weighted samples for KDE fitting."""
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    total = counts.sum()
    if total == 0:
        return bin_centers
    probs = counts / total
    indices = np.random.choice(len(bin_centers), size=n_samples, p=probs)
    bin_width = bin_edges[1] - bin_edges[0]
    samples = bin_centers[indices] + np.random.uniform(-bin_width / 2, bin_width / 2, n_samples)
    return samples


def fit_kde(data, bw_method="scott"):
    """Fit gaussian KDE, handling edge cases."""
    data = data[np.isfinite(data)]
    if len(data) < 10:
        log("  WARNING: Too few data points for KDE, returning uniform")
        return None
    spread = data.max() - data.min()
    if spread < 1e-10:
        log("  WARNING: Near-zero variance, adding jitter")
        data = data + np.random.normal(0, 1e-6, len(data))
    try:
        kde = gaussian_kde(data, bw_method=bw_method)
        return kde
    except Exception as e:
        log(f"  WARNING: KDE fitting failed: {e}")
        return None


def estimate_pi0_storey(observed_scores, f0_kde, f_obs_kde):
    """Estimate pi0 using Storey's method: pi0 = #{scores < median} / expected_under_null."""
    if f0_kde is None or f_obs_kde is None:
        log("  WARNING: Cannot estimate pi0 without both KDEs, defaulting to 0.5")
        return 0.5

    med = np.median(observed_scores)

    n_below_median = np.sum(observed_scores <= med)
    n_total = len(observed_scores)
    frac_below = n_below_median / n_total

    eval_grid = np.linspace(observed_scores.min(), med, 500)
    null_cdf_at_median = np.trapezoid(f0_kde(eval_grid), eval_grid)
    null_cdf_at_median = np.clip(null_cdf_at_median, 0.01, 1.0)

    pi0 = frac_below / null_cdf_at_median
    pi0 = np.clip(pi0, 0.01, 0.99)

    return float(pi0)


def process_layer(name):
    """Process one layer: load null histogram + observed edges, fit KDEs, estimate pi0."""
    t_layer = time.time()
    log(f"\n--- Layer: {name} ---")

    null_path = NULL_DIR / f"layer_{name}_null_hist.npz"
    edge_path = EDGE_DIR / f"edges_{name}.csv"

    if not null_path.exists():
        log(f"  Null histogram not found: {null_path}, skipping")
        return False
    if not edge_path.exists():
        log(f"  Edge parquet not found: {edge_path}, skipping")
        return False

    null_data = np.load(null_path)
    bin_edges = null_data["bin_edges"]
    null_counts = null_data["counts"]
    n_perms = int(null_data["n_perms"][0])
    n_edges = int(null_data["n_edges"][0])
    log(f"  Null histogram: {len(null_counts)} bins, {n_perms} perms, {n_edges} edges")

    edges = pd.read_csv(edge_path)
    obs_scores = edges["raw_score"].values.astype(np.float64)
    log(f"  Observed scores: {len(obs_scores):,} edges, range [{obs_scores.min():.4f}, {obs_scores.max():.4f}]")

    log("  Fitting null KDE (f0)...")
    null_samples = histogram_to_samples(bin_edges, null_counts, n_samples=min(100000, n_perms * n_edges))
    f0_kde = fit_kde(null_samples)

    log("  Fitting observed KDE (f_obs)...")
    obs_subsample = obs_scores if len(obs_scores) <= 100000 else np.random.choice(obs_scores, 100000, replace=False)
    f_obs_kde = fit_kde(obs_subsample)

    log("  Estimating pi0 (Storey's method)...")
    pi0 = estimate_pi0_storey(obs_scores, f0_kde, f_obs_kde)
    log(f"  pi0 = {pi0:.4f}")

    obs_hist_counts, _ = np.histogram(obs_scores, bins=bin_edges)

    out_path = NULL_DIR / f"layer_{name}_density.pkl"
    result = {
        "layer_name": name,
        "f0": f0_kde,
        "f_obs": f_obs_kde,
        "pi0": pi0,
        "bin_edges": bin_edges,
        "null_counts": null_counts,
        "obs_counts": obs_hist_counts,
        "obs_scores": obs_scores,
        "n_perms": n_perms,
        "n_edges_raw": n_edges,
    }
    with open(out_path, "wb") as f:
        pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
    log(f"  Saved: {out_path}")

    elapsed_layer = time.time() - t_layer
    log(f"  Layer {name} done in {elapsed_layer:.1f}s")
    return True


def main():
    log("=== 261: Null Density Estimation ===")
    log(f"Edge directory: {EDGE_DIR}")
    log(f"Null directory: {NULL_DIR}")

    layers_processed = 0
    layers_skipped = 0
    summary_rows = []

    for name in LAYER_NAMES:
        success = process_layer(name)
        if success:
            layers_processed += 1
            pkl_path = NULL_DIR / f"layer_{name}_density.pkl"
            with open(pkl_path, "rb") as f:
                result = pickle.load(f)
            summary_rows.append({
                "layer": name,
                "pi0": result["pi0"],
                "n_edges": result["n_edges_raw"],
                "n_perms": result["n_perms"],
                "score_mean": float(result["obs_scores"].mean()),
                "score_median": float(np.median(result["obs_scores"])),
                "score_min": float(result["obs_scores"].min()),
                "score_max": float(result["obs_scores"].max()),
            })
        else:
            layers_skipped += 1

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = NULL_DIR / "density_estimation_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        log(f"\nSummary saved: {summary_path}")
        log("\n" + summary_df.to_string(index=False))

    elapsed_total = time.time() - t0_global
    log(f"\n=== 261 Complete ===")
    log(f"Layers processed: {layers_processed}, skipped: {layers_skipped}")
    log(f"Total time: {elapsed_total:.1f}s ({elapsed_total / 60:.1f} min)")


if __name__ == "__main__":
    main()
