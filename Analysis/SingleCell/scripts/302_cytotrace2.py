#!/usr/bin/env python3
"""
302: CytoTRACE 2 — Cell Potency Prediction.

Predicts differentiation potential per cell for each of 5 core cell types.
Tests whether MASLD disease states show dedifferentiation (increased potency).

Inputs:
    - Cell-type subsets from 300: pseudotime/{celltype}_subset.h5ad

Outputs (to results_gpu_v2/pseudotime/):
    - cytotrace2_scores_{celltype}.csv   (when CytoTRACE 2 is available)
    - cytotrace_v1_scores_{celltype}.csv (when falling back to CytoTRACE v1)
    - cytotrace2_comparison.csv

Usage:
    sbatch run_pseudotime_pipeline.sh
"""

import os
import sys
import warnings
import logging
import gc

import numpy as np
import pandas as pd
from scipy import stats, sparse

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SC_DIR = os.path.join(BASE, "Analysis/SingleCell")
RESULTS = os.path.join(SC_DIR, "results_gpu_v2")
PT_DIR = os.path.join(RESULTS, "pseudotime")

CELL_TYPES = [
    "Hepatocytes",
    "Macrophages",
    "Fibroblasts",
    "Endothelial_cells",
    "Cholangiocytes",
]

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import scanpy as sc

# Try CytoTRACE 2
try:
    import cytotrace2
    HAS_CYTOTRACE2 = True
    log.info("CytoTRACE 2 available")
except ImportError:
    HAS_CYTOTRACE2 = False
    log.warning("CytoTRACE 2 not installed. Using CytoTRACE v1 approximation "
                "(gene count correlation). Install: pip install cytotrace2")

CONDITION_ORDER = {"Healthy": 0, "NAFLD": 1, "MASLD": 1.5, "NASH": 2, "Cirrhotic": 3}


# ---------------------------------------------------------------------------
# CytoTRACE v1 approximation (if v2 not available)
# ---------------------------------------------------------------------------
def cytotrace_v1(adata):
    """Approximate CytoTRACE v1: gene counts correlate with stemness.

    CytoTRACE v1 principle: number of expressed genes per cell
    correlates with differentiation potential. More genes expressed
    = more stem-like (higher potency).
    """
    log.info("Running CytoTRACE v1 approximation...")

    X = adata.X
    if sparse.issparse(X):
        n_genes_per_cell = np.array((X > 0).sum(axis=1)).flatten()
    else:
        n_genes_per_cell = np.sum(X > 0, axis=1)

    # Rank-normalize to [0, 1]
    ranks = stats.rankdata(n_genes_per_cell)
    potency = ranks / len(ranks)

    # Correlate each gene with gene count (CytoTRACE score)
    gene_corrs = []
    X_dense = X.toarray() if sparse.issparse(X) else X
    for i in range(min(5000, X_dense.shape[1])):
        expr = X_dense[:, i]
        if np.std(expr) > 0:
            rho, _ = stats.spearmanr(n_genes_per_cell, expr)
            gene_corrs.append({"gene": adata.var_names[i], "cytotrace_corr": rho})

    gene_df = pd.DataFrame(gene_corrs).sort_values("cytotrace_corr", ascending=False)

    return potency, gene_df


# ---------------------------------------------------------------------------
# CytoTRACE 2 wrapper
# ---------------------------------------------------------------------------
def run_cytotrace2(adata):
    """Run CytoTRACE 2 deep learning model."""
    log.info("Running CytoTRACE 2...")
    results = cytotrace2.cytotrace2(adata)
    return results


# ---------------------------------------------------------------------------
# Statistical comparison across conditions
# ---------------------------------------------------------------------------
def compare_potency_across_conditions(potency_df, celltype):
    """Test potency differences across disease conditions."""
    log.info("Comparing potency across conditions for %s...", celltype)

    results = []

    # Kruskal-Wallis across all conditions
    groups = []
    group_labels = []
    for cond in sorted(potency_df["condition"].unique()):
        vals = potency_df.loc[potency_df["condition"] == cond, "potency"].values
        if len(vals) > 5:
            groups.append(vals)
            group_labels.append(cond)

    if len(groups) >= 2:
        kw_stat, kw_pval = stats.kruskal(*groups)
        log.info("Kruskal-Wallis: H=%.2f, p=%.2e", kw_stat, kw_pval)
        results.append({
            "cell_type": celltype,
            "test": "kruskal_wallis",
            "stat": kw_stat,
            "pval": kw_pval,
            "comparison": "all_conditions",
        })

    # Pairwise: Healthy vs. each disease state
    if "Healthy" in potency_df["condition"].values:
        healthy = potency_df.loc[
            potency_df["condition"] == "Healthy", "potency"
        ].values

        for cond in ["NAFLD", "MASLD", "NASH", "Cirrhotic"]:
            disease = potency_df.loc[
                potency_df["condition"] == cond, "potency"
            ].values
            if len(disease) > 5 and len(healthy) > 5:
                stat, pval = stats.mannwhitneyu(
                    disease, healthy, alternative="two-sided"
                )
                effect = np.mean(disease) - np.mean(healthy)
                log.info("  %s vs Healthy: U=%.0f, p=%.2e, effect=%.3f",
                         cond, stat, pval, effect)
                results.append({
                    "cell_type": celltype,
                    "test": "mannwhitneyu",
                    "stat": stat,
                    "pval": pval,
                    "comparison": f"{cond}_vs_Healthy",
                    "effect_size": effect,
                    "mean_disease": np.mean(disease),
                    "mean_healthy": np.mean(healthy),
                })

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 60)
    log.info("302: CytoTRACE 2 Potency Prediction")
    log.info("=" * 60)

    all_comparisons = []

    for ct in CELL_TYPES:
        log.info("\n" + "=" * 50)
        log.info("Processing: %s", ct)
        log.info("=" * 50)

        h5ad_path = os.path.join(PT_DIR, f"{ct}_subset.h5ad")
        if not os.path.exists(h5ad_path):
            log.warning("Subset not found: %s — skipping", h5ad_path)
            continue

        adata = sc.read_h5ad(h5ad_path)
        log.info("Loaded %s: %s", ct, adata.shape)

        # Run CytoTRACE
        if HAS_CYTOTRACE2:
            try:
                ct2_results = run_cytotrace2(adata)
                potency = ct2_results.obs["CytoTRACE2_Score"].values
                potency_type = "cytotrace2"
            except Exception as e:
                log.error("CytoTRACE 2 failed: %s — falling back to v1", e)
                potency, gene_df = cytotrace_v1(adata)
                potency_type = "cytotrace_v1"
        else:
            potency, gene_df = cytotrace_v1(adata)
            potency_type = "cytotrace_v1"

        log.info("Potency (%s): median=%.3f, range=[%.3f, %.3f]",
                 potency_type, np.median(potency), np.min(potency), np.max(potency))

        # Build output
        score_df = pd.DataFrame(index=adata.obs_names)
        score_df["potency"] = potency
        score_df["potency_method"] = potency_type
        score_df["condition"] = adata.obs["condition"].values
        score_df["leiden_substate"] = adata.obs["leiden_substate"].values
        score_df["dataset"] = adata.obs["dataset"].values

        # Save scores
        out_suffix = "cytotrace2" if potency_type == "cytotrace2" else "cytotrace_v1"
        out_path = os.path.join(PT_DIR, f"{out_suffix}_scores_{ct}.csv")
        score_df.to_csv(out_path)
        log.info("Saved potency scores to %s", out_path)

        # Correlate with Palantir differentiation potential if available
        palantir_path = os.path.join(PT_DIR, f"palantir_pseudotime_{ct}.csv")
        if os.path.exists(palantir_path):
            pr_df = pd.read_csv(palantir_path, index_col=0)
            if "palantir_entropy" in pr_df.columns:
                common = score_df.index.intersection(pr_df.index)
                if len(common) > 30:
                    rho, pval = stats.spearmanr(
                        score_df.loc[common, "potency"],
                        pr_df.loc[common, "palantir_entropy"],
                    )
                    log.info("CytoTRACE vs Palantir entropy: rho=%.3f, p=%.2e",
                             rho, pval)

        # Statistical comparison across conditions
        comparisons = compare_potency_across_conditions(score_df, ct)
        all_comparisons.extend(comparisons)

        del adata
        gc.collect()

    # Save comparison summary
    if all_comparisons:
        comp_df = pd.DataFrame(all_comparisons)
        out_path = os.path.join(PT_DIR, "cytotrace2_comparison.csv")
        comp_df.to_csv(out_path, index=False)
        log.info("\nComparison summary:\n%s", comp_df.to_string())

    log.info("\n=== 302: CytoTRACE 2 COMPLETE ===")


if __name__ == "__main__":
    main()
