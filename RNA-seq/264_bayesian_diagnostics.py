#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=net_264_diagnostics
#SBATCH --output=logs/net_264_diagnostics_%j.out
#SBATCH --error=logs/net_264_diagnostics_%j.err
# ===========================================================================
# Script 264: Diagnostic plots for the Bayesian calibration pipeline
# ===========================================================================
# Purpose: Generate per-layer and cross-layer diagnostics to validate the
#          null model fit, posterior calibration, and layer informativeness.
#
# Input:
#   - RNA-seq/results/network/null_distributions/layer_{name}_density.pkl
#   - RNA-seq/results/network/posterior_edges/edges_{name}_posterior.csv
#   - RNA-seq/results/network/composite_edges.csv
#
# Output:
#   - RNA-seq/results/network/diagnostics/*.pdf (per-layer and cross-layer)
#   - RNA-seq/results/network/bayesian_diagnostics.csv
#
# Environment: spatial (scipy, numpy, pandas, matplotlib, pyarrow)
# ===========================================================================

import os
import pickle
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

warnings.filterwarnings("ignore", category=FutureWarning)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
EDGE_DIR = BASE / "RNA-seq" / "results" / "network"
NULL_DIR = EDGE_DIR / "null_distributions"
POST_DIR = EDGE_DIR / "posterior_edges"
DIAG_DIR = EDGE_DIR / "diagnostics"
DIAG_DIR.mkdir(parents=True, exist_ok=True)

LAYER_NAMES = [
    "ppi", "coexpr", "regulon", "lr", "genetic",
    "pathway", "spatial", "cosmos", "cerna", "xspecies"
]

LAYER_COLORS = {
    "ppi": "#1f77b4",
    "coexpr": "#ff7f0e",
    "regulon": "#2ca02c",
    "lr": "#d62728",
    "genetic": "#9467bd",
    "pathway": "#8c564b",
    "spatial": "#e377c2",
    "cosmos": "#7f7f7f",
    "cerna": "#bcbd22",
    "xspecies": "#17becf",
}

t0_global = time.time()


def log(msg):
    elapsed = time.time() - t0_global
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


def load_density(name):
    path = NULL_DIR / f"layer_{name}_density.pkl"
    if not path.exists():
        return None
    with open(path, "rb") as f:
        return pickle.load(f)


def load_posteriors(name):
    path = POST_DIR / f"edges_{name}_posterior.csv"
    if not path.exists():
        return None
    return pd.read_parquet(path)


def plot_layer_diagnostics(name, density, posteriors):
    """Generate 4-panel diagnostic for one layer: density overlay, posterior hist, Q-Q, summary."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"Bayesian Calibration Diagnostics: {name}", fontsize=14, fontweight="bold")
    color = LAYER_COLORS.get(name, "#333333")

    obs_scores = density["obs_scores"]
    f0_kde = density["f0"]
    f_obs_kde = density["f_obs"]
    pi0 = density["pi0"]
    bin_edges = density["bin_edges"]
    null_counts = density["null_counts"]
    obs_counts = density["obs_counts"]

    ax = axes[0, 0]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    bin_width = bin_edges[1] - bin_edges[0]

    null_total = null_counts.sum()
    obs_total = obs_counts.sum()
    if null_total > 0:
        ax.bar(bin_centers, null_counts / (null_total * bin_width),
               width=bin_width * 0.8, alpha=0.3, color="gray", label="Null (permutation)")
    if obs_total > 0:
        ax.bar(bin_centers, obs_counts / (obs_total * bin_width),
               width=bin_width * 0.8, alpha=0.3, color=color, label="Observed")

    if f0_kde is not None and f_obs_kde is not None:
        x_grid = np.linspace(bin_edges[0], bin_edges[-1], 500)
        ax.plot(x_grid, f0_kde(x_grid), "k--", lw=1.5, label=f"f0 KDE (null)")
        ax.plot(x_grid, f_obs_kde(x_grid), "-", color=color, lw=1.5, label="f KDE (obs)")

    ax.set_xlabel("Edge score")
    ax.set_ylabel("Density")
    ax.set_title(f"Null vs Observed Density (pi0={pi0:.3f})")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    post_vals = posteriors["posterior"].values
    ax.hist(post_vals, bins=100, color=color, alpha=0.7, edgecolor="white", linewidth=0.3)
    ax.axvline(0.5, color="red", ls="--", lw=1, label="P=0.5")
    ax.axvline(0.9, color="darkred", ls="--", lw=1, label="P=0.9")
    n_gt05 = (post_vals > 0.5).sum()
    n_gt09 = (post_vals > 0.9).sum()
    ax.set_xlabel("Posterior P_m")
    ax.set_ylabel("Count")
    ax.set_title(f"Posterior Distribution (>0.5: {n_gt05:,}, >0.9: {n_gt09:,})")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    if f0_kde is not None:
        n_quantiles = min(1000, len(obs_scores))
        quantile_probs = np.linspace(0.01, 0.99, n_quantiles)
        obs_quantiles = np.quantile(obs_scores, quantile_probs)

        null_samples_for_qq = density.get("_null_samples_qq", None)
        if null_samples_for_qq is None:
            null_samples_for_qq = np.random.choice(
                bin_centers,
                size=min(100000, len(obs_scores) * 10),
                p=null_counts / null_counts.sum() if null_counts.sum() > 0 else np.ones(len(null_counts)) / len(null_counts)
            )
        null_quantiles = np.quantile(null_samples_for_qq, quantile_probs)

        ax.scatter(null_quantiles, obs_quantiles, s=2, alpha=0.5, color=color)
        lo = min(null_quantiles.min(), obs_quantiles.min())
        hi = max(null_quantiles.max(), obs_quantiles.max())
        ax.plot([lo, hi], [lo, hi], "k--", lw=1, alpha=0.5)
        ax.set_xlabel("Null quantiles")
        ax.set_ylabel("Observed quantiles")
        ax.set_title("Q-Q Plot (Observed vs Null)")
    else:
        ax.text(0.5, 0.5, "KDE unavailable", ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Q-Q Plot (unavailable)")

    ax = axes[1, 1]
    ax.axis("off")
    summary_text = (
        f"Layer: {name}\n"
        f"Total edges: {len(posteriors):,}\n"
        f"pi0 (null fraction): {pi0:.4f}\n"
        f"Mean score: {obs_scores.mean():.4f}\n"
        f"Median score: {np.median(obs_scores):.4f}\n"
        f"Score range: [{obs_scores.min():.4f}, {obs_scores.max():.4f}]\n"
        f"\nPosterior summary:\n"
        f"  Mean: {post_vals.mean():.4f}\n"
        f"  Median: {np.median(post_vals):.4f}\n"
        f"  P > 0.5: {n_gt05:,} ({100*n_gt05/len(post_vals):.1f}%)\n"
        f"  P > 0.9: {n_gt09:,} ({100*n_gt09/len(post_vals):.1f}%)\n"
    )
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=10,
            verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = DIAG_DIR / f"layer_{name}_diagnostics.pdf"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"  Saved: {out_path}")


def plot_cross_layer_heatmap():
    """Cross-layer informativeness heatmap by gene degree class."""
    log("\nGenerating cross-layer heatmap...")

    composite_path = EDGE_DIR / "composite_edges.csv"
    if not composite_path.exists():
        log("  Composite edges not found, skipping cross-layer heatmap")
        return

    composite = pd.read_csv(composite_path)

    all_genes = set(composite["gene_a"]) | set(composite["gene_b"])
    gene_degree = {}
    for g in all_genes:
        d = ((composite["gene_a"] == g) | (composite["gene_b"] == g)).sum()
        gene_degree[g] = d

    degrees = np.array(list(gene_degree.values()))
    if len(degrees) == 0:
        log("  No genes found, skipping heatmap")
        return

    q25, q50, q75 = np.quantile(degrees, [0.25, 0.5, 0.75])
    degree_class = {}
    for g, d in gene_degree.items():
        if d <= q25:
            degree_class[g] = "Low"
        elif d <= q50:
            degree_class[g] = "Med-Low"
        elif d <= q75:
            degree_class[g] = "Med-High"
        else:
            degree_class[g] = "High"

    classes = ["Low", "Med-Low", "Med-High", "High"]
    posterior_cols = [f"p_{name}" for name in LAYER_NAMES if f"p_{name}" in composite.columns]
    available_layers = [name for name in LAYER_NAMES if f"p_{name}" in composite.columns]

    heatmap_data = np.zeros((len(available_layers), len(classes)))
    for i, name in enumerate(available_layers):
        col = f"p_{name}"
        for j, cls in enumerate(classes):
            cls_genes = {g for g, c in degree_class.items() if c == cls}
            mask = composite["gene_a"].isin(cls_genes) | composite["gene_b"].isin(cls_genes)
            subset = composite.loc[mask, col]
            heatmap_data[i, j] = (subset > 0.5).mean() if len(subset) > 0 else 0.0

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(heatmap_data, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(classes)
    ax.set_yticks(range(len(available_layers)))
    ax.set_yticklabels(available_layers)
    ax.set_xlabel("Gene Degree Class")
    ax.set_ylabel("Edge Layer")
    ax.set_title("Layer Informativeness by Gene Degree Class\n(Fraction of edges with P > 0.5)")

    for i in range(len(available_layers)):
        for j in range(len(classes)):
            val = heatmap_data[i, j]
            text_color = "white" if val > 0.5 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=9, color=text_color)

    plt.colorbar(im, ax=ax, label="Fraction P > 0.5")
    plt.tight_layout()
    out_path = DIAG_DIR / "cross_layer_degree_heatmap.pdf"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"  Saved: {out_path}")


def plot_pi0_comparison():
    """Bar chart of pi0 across all layers."""
    log("\nGenerating pi0 comparison plot...")

    layer_pi0 = {}
    for name in LAYER_NAMES:
        density = load_density(name)
        if density is not None:
            layer_pi0[name] = density["pi0"]

    if not layer_pi0:
        log("  No density data found, skipping pi0 plot")
        return

    names = list(layer_pi0.keys())
    values = [layer_pi0[n] for n in names]
    colors = [LAYER_COLORS.get(n, "#333333") for n in names]

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(range(len(names)), values, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_ylabel("pi0 (null proportion)")
    ax.set_title("Estimated Null Proportion per Layer")
    ax.axhline(0.5, color="red", ls="--", lw=1, alpha=0.5, label="pi0 = 0.5")
    ax.set_ylim(0, 1.05)
    ax.legend()

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f"{val:.3f}", ha="center", va="bottom", fontsize=9)

    plt.tight_layout()
    out_path = DIAG_DIR / "pi0_comparison.pdf"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"  Saved: {out_path}")


def build_summary_csv():
    """Build bayesian_diagnostics.csv with per-layer statistics."""
    log("\nBuilding summary CSV...")

    rows = []
    for name in LAYER_NAMES:
        density = load_density(name)
        posteriors = load_posteriors(name)
        if density is None or posteriors is None:
            continue

        post_vals = posteriors["posterior"].values
        obs_scores = density["obs_scores"]

        rows.append({
            "name": name,
            "pi0": density["pi0"],
            "n_edges_raw": len(posteriors),
            "n_edges_posterior_05": int((post_vals > 0.5).sum()),
            "n_edges_posterior_09": int((post_vals > 0.9).sum()),
            "mean_posterior": float(post_vals.mean()),
            "median_posterior": float(np.median(post_vals)),
            "mean_raw_score": float(obs_scores.mean()),
            "median_raw_score": float(np.median(obs_scores)),
            "score_range_min": float(obs_scores.min()),
            "score_range_max": float(obs_scores.max()),
        })

    if rows:
        df = pd.DataFrame(rows)
        out_path = EDGE_DIR / "bayesian_diagnostics.csv"
        df.to_parquet(out_path, index=False)
        log(f"Saved: {out_path}")
        log("\n" + df.to_string(index=False))
    else:
        log("  WARNING: No layers had both density and posterior data")


def main():
    log("=== 264: Bayesian Calibration Diagnostics ===")
    log(f"Network directory: {EDGE_DIR}")
    log(f"Diagnostics output: {DIAG_DIR}")

    layers_plotted = 0
    for name in LAYER_NAMES:
        log(f"\n--- Layer: {name} ---")
        density = load_density(name)
        posteriors = load_posteriors(name)

        if density is None:
            log(f"  Density not found, skipping")
            continue
        if posteriors is None:
            log(f"  Posteriors not found, skipping")
            continue

        plot_layer_diagnostics(name, density, posteriors)
        layers_plotted += 1

    plot_cross_layer_heatmap()
    plot_pi0_comparison()
    build_summary_csv()

    elapsed_total = time.time() - t0_global
    log(f"\n=== 264 Complete ===")
    log(f"Layers plotted: {layers_plotted}")
    log(f"Total time: {elapsed_total:.1f}s ({elapsed_total / 60:.1f} min)")


if __name__ == "__main__":
    main()
