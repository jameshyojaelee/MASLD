#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_263_composite
#SBATCH --output=logs/net_263_composite_%j.out
#SBATCH --error=logs/net_263_composite_%j.err
# ===========================================================================
# Script 263: Compute composite probability across all layers via noisy-OR
# ===========================================================================
# Purpose: For each unique gene pair across all 10 layers, compute the
#          noisy-OR composite: P_composite = 1 - prod(1 - P_m), plus
#          multiplicity K_ij = count of layers with P_m > 0.5.
#
# Input:
#   - RNA-seq/results/network/posterior_edges/edges_{name}_posterior.csv
#
# Output:
#   - RNA-seq/results/network/composite_edges.csv
#     Columns: gene_a, gene_b, p_composite, k_multiplicity, p_ppi, p_coexpr, ...
#   - RNA-seq/results/network/edge_summary.csv
#     Summary stats: total edges per threshold, mean K, layer contribution
#
# Environment: spatial (numpy, pandas, pyarrow)
# ===========================================================================

import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
EDGE_DIR = BASE / "RNA-seq" / "results" / "network"
POST_DIR = EDGE_DIR / "posterior_edges"
OUTDIR = EDGE_DIR
OUTDIR.mkdir(parents=True, exist_ok=True)

LAYER_NAMES = [
    "ppi", "coexpr", "regulon", "lr", "genetic",
    "pathway", "spatial", "cosmos", "cerna", "xspecies"
]

P_COMPOSITE_THRESHOLDS = [0.5, 0.7, 0.8, 0.9, 0.95]

t0_global = time.time()


def log(msg):
    elapsed = time.time() - t0_global
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


def load_posterior_edges():
    """Load all posterior edge parquets and build a unified gene-pair table."""
    layer_data = {}
    for name in LAYER_NAMES:
        path = POST_DIR / f"edges_{name}_posterior.csv"
        if not path.exists():
            log(f"  Posterior parquet not found: {path}, skipping layer '{name}'")
            continue
        df = pd.read_parquet(path, columns=["gene_a", "gene_b", "posterior"])
        df = df.rename(columns={"posterior": f"p_{name}"})
        layer_data[name] = df
        log(f"  Loaded {len(df):,} edges from {name}")
    return layer_data


def merge_layers(layer_data):
    """Merge all layer edges into a single table keyed by (gene_a, gene_b)."""
    if not layer_data:
        log("ERROR: No layer data loaded")
        sys.exit(1)

    merged = None
    for name, df in layer_data.items():
        if merged is None:
            merged = df
        else:
            merged = merged.merge(df, on=["gene_a", "gene_b"], how="outer")
        log(f"  After merging {name}: {len(merged):,} unique gene pairs")

    for name in LAYER_NAMES:
        col = f"p_{name}"
        if col in merged.columns:
            merged[col] = merged[col].fillna(0.0)
        else:
            merged[col] = 0.0

    return merged


def compute_noisy_or(merged):
    """Compute noisy-OR composite and multiplicity for each gene pair."""
    log("Computing noisy-OR composite probabilities...")

    posterior_cols = [f"p_{name}" for name in LAYER_NAMES if f"p_{name}" in merged.columns]

    posterior_matrix = merged[posterior_cols].values.astype(np.float64)

    log_complement = np.log1p(-np.clip(posterior_matrix, 0.0, 1.0 - 1e-15))
    log_product = log_complement.sum(axis=1)
    p_composite = 1.0 - np.exp(log_product)
    p_composite = np.clip(p_composite, 0.0, 1.0)

    merged["p_composite"] = p_composite

    k_multiplicity = (posterior_matrix > 0.5).sum(axis=1)
    merged["k_multiplicity"] = k_multiplicity.astype(np.int32)

    return merged


def filter_and_save(merged):
    """Apply adaptive thresholding and save outputs."""
    n_total = len(merged)
    log(f"\nTotal unique gene pairs: {n_total:,}")

    for thresh in P_COMPOSITE_THRESHOLDS:
        n_pass = (merged["p_composite"] > thresh).sum()
        pct = 100.0 * n_pass / n_total if n_total > 0 else 0.0
        log(f"  P_composite > {thresh}: {n_pass:,} edges ({pct:.2f}%)")

    target_threshold = 0.9
    n_at_target = (merged["p_composite"] > target_threshold).sum()
    if n_at_target < 100:
        log(f"  WARNING: Only {n_at_target} edges at P>{target_threshold}, relaxing to 0.5")
        target_threshold = 0.5

    out_cols = ["gene_a", "gene_b", "p_composite", "k_multiplicity"]
    for name in LAYER_NAMES:
        col = f"p_{name}"
        if col in merged.columns:
            out_cols.append(col)

    merged_sorted = merged.sort_values("p_composite", ascending=False)

    composite_path = OUTDIR / "composite_edges.csv"
    merged_sorted[out_cols].to_csv(composite_path, index=False)
    log(f"\nSaved composite edges: {composite_path}")
    log(f"  Total edges: {len(merged_sorted):,}")

    return merged_sorted, target_threshold


def build_summary(merged):
    """Build edge summary statistics CSV."""
    log("\nBuilding edge summary...")

    summary_rows = []

    for thresh in P_COMPOSITE_THRESHOLDS:
        subset = merged[merged["p_composite"] > thresh]
        n_edges = len(subset)
        mean_k = subset["k_multiplicity"].mean() if n_edges > 0 else 0.0
        median_k = subset["k_multiplicity"].median() if n_edges > 0 else 0.0

        row = {
            "threshold": thresh,
            "n_edges": n_edges,
            "mean_k_multiplicity": float(mean_k),
            "median_k_multiplicity": float(median_k),
            "mean_p_composite": float(subset["p_composite"].mean()) if n_edges > 0 else 0.0,
        }

        for name in LAYER_NAMES:
            col = f"p_{name}"
            if col in merged.columns:
                n_contributing = (subset[col] > 0.5).sum() if n_edges > 0 else 0
                row[f"n_contributing_{name}"] = int(n_contributing)
                row[f"frac_contributing_{name}"] = float(n_contributing / n_edges) if n_edges > 0 else 0.0

        summary_rows.append(row)

    layer_stats = []
    for name in LAYER_NAMES:
        col = f"p_{name}"
        if col not in merged.columns:
            continue
        vals = merged[col].values
        n_present = (vals > 0).sum()
        n_sig = (vals > 0.5).sum()
        layer_stats.append({
            "layer": name,
            "n_edges_present": int(n_present),
            "n_edges_sig": int(n_sig),
            "mean_posterior": float(vals[vals > 0].mean()) if n_present > 0 else 0.0,
            "median_posterior": float(np.median(vals[vals > 0])) if n_present > 0 else 0.0,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = OUTDIR / "edge_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    log(f"Saved threshold summary: {summary_path}")

    layer_df = pd.DataFrame(layer_stats)
    layer_path = OUTDIR / "layer_contribution_summary.csv"
    layer_df.to_csv(layer_path, index=False)
    log(f"Saved layer contribution: {layer_path}")

    log("\nThreshold summary:")
    log(summary_df.to_string(index=False))
    log("\nLayer contribution:")
    log(layer_df.to_string(index=False))

    return summary_df, layer_df


def print_top_edges(merged, n=20):
    """Print top composite edges for sanity check."""
    log(f"\nTop {n} composite edges:")
    top = merged.head(n)
    for _, row in top.iterrows():
        layers_active = []
        for name in LAYER_NAMES:
            col = f"p_{name}"
            if col in merged.columns and row[col] > 0.5:
                layers_active.append(f"{name}={row[col]:.2f}")
        layer_str = ", ".join(layers_active) if layers_active else "none>0.5"
        log(f"  {row['gene_a']:12s} -- {row['gene_b']:12s}  "
            f"P={row['p_composite']:.4f}  K={row['k_multiplicity']}  [{layer_str}]")


def main():
    log("=== 263: Noisy-OR Composite Probabilities ===")
    log(f"Posterior directory: {POST_DIR}")
    log(f"Output directory: {OUTDIR}")

    layer_data = load_posterior_edges()
    merged = merge_layers(layer_data)
    merged = compute_noisy_or(merged)
    merged_sorted, threshold = filter_and_save(merged)
    build_summary(merged_sorted)
    print_top_edges(merged_sorted)

    n_genes_a = merged_sorted["gene_a"].nunique()
    n_genes_b = merged_sorted["gene_b"].nunique()
    n_genes = len(set(merged_sorted["gene_a"]) | set(merged_sorted["gene_b"]))

    elapsed_total = time.time() - t0_global
    log(f"\n=== 263 Complete ===")
    log(f"Unique genes in network: {n_genes:,}")
    log(f"Total composite edges: {len(merged_sorted):,}")
    log(f"Edges at P>{threshold}: {(merged_sorted['p_composite'] > threshold).sum():,}")
    log(f"Total time: {elapsed_total:.1f}s ({elapsed_total / 60:.1f} min)")


if __name__ == "__main__":
    main()
