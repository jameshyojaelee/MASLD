#!/usr/bin/env python3
"""
09a_zonation_aware_de.py — Zonation-aware differential expression.

Identifies genes differentially expressed in specific zonation bins
(zone-restricted MASLD signatures). Aggregates to sample-level to
avoid pseudo-replication, then uses Wilcoxon rank-sum between conditions.

SLURM: --partition=cpu --cpus=16 --mem=128G --time=8:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_dream_degs,
    build_ensembl_to_symbol_map, save_csv, print_header, print_step,
)


def aggregate_to_sample(adata, gene, bin_label):
    """Aggregate expression to sample-level means within a zonation bin."""
    mask = adata.obs["zonation_bin"] == bin_label
    if mask.sum() == 0:
        return pd.DataFrame()

    adata_bin = adata[mask]
    expr = np.asarray(adata_bin[:, gene].X.todense()).flatten()

    records = []
    for sample in adata_bin.obs["sample_id"].unique():
        smask = adata_bin.obs["sample_id"] == sample
        if smask.sum() < 5:  # Need minimum spots per sample per bin
            continue
        records.append({
            "sample_id": sample,
            "condition": adata_bin.obs.loc[smask, "condition"].iloc[0],
            "mean_expr": float(expr[smask.values].mean()),
            "n_spots": int(smask.sum()),
        })
    return pd.DataFrame(records)


def zonation_aware_de(adata, genes, bin_labels, condition_healthy, condition_masld):
    """Run zone-specific DE for all genes across all bins."""
    results = []
    n_genes = len(genes)

    for i, gene in enumerate(genes):
        if i % 200 == 0:
            print_step(f"Processing gene {i+1}/{n_genes}", i+1, n_genes)

        gene_results = {}
        for bin_label in bin_labels:
            agg = aggregate_to_sample(adata, gene, bin_label)
            if len(agg) < 2:
                continue

            healthy = agg[agg["condition"] == condition_healthy]["mean_expr"].values
            masld = agg[agg["condition"] == condition_masld]["mean_expr"].values

            if len(healthy) < 1 or len(masld) < 1:
                continue

            try:
                U, pval = mannwhitneyu(masld, healthy, alternative="two-sided")
            except Exception:
                U, pval = np.nan, 1.0

            lfc = np.log2(
                max(masld.mean(), 1e-6) / max(healthy.mean(), 1e-6)
            )

            gene_results[bin_label] = {
                "pval": pval,
                "lfc": lfc,
                "mean_healthy": float(healthy.mean()),
                "mean_masld": float(masld.mean()),
                "n_healthy": len(healthy),
                "n_masld": len(masld),
            }

        if not gene_results:
            continue

        # Count significant bins
        pvals = [v["pval"] for v in gene_results.values()]
        n_sig_bins = sum(1 for p in pvals if p < 0.05)
        max_abs_lfc = max(abs(v["lfc"]) for v in gene_results.values())
        min_abs_lfc = min(abs(v["lfc"]) for v in gene_results.values()) if pvals else 0

        # Classification
        if n_sig_bins == 0:
            classification = "Not_significant"
        elif n_sig_bins <= 2 and n_sig_bins < len(gene_results):
            classification = "Zone_restricted"
        elif n_sig_bins == len(gene_results) and max_abs_lfc > 2 * min_abs_lfc:
            classification = "Zone_amplified"
        elif n_sig_bins == len(gene_results):
            classification = "Pan_lobular"
        else:
            classification = "Partial"

        row = {
            "gene": gene,
            "classification": classification,
            "n_sig_bins": n_sig_bins,
            "n_tested_bins": len(gene_results),
            "max_abs_lfc": max_abs_lfc,
        }
        for bin_label, vals in gene_results.items():
            row[f"{bin_label}_pval"] = vals["pval"]
            row[f"{bin_label}_lfc"] = vals["lfc"]
            row[f"{bin_label}_mean_h"] = vals["mean_healthy"]
            row[f"{bin_label}_mean_m"] = vals["mean_masld"]

        results.append(row)

    return pd.DataFrame(results)


def main():
    print_header("09a: Zonation-Aware Differential Expression")

    config = load_config()
    output_dir = RESULTS_DIR / "zonation"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load zonation-annotated data
    zon_h5ad = RESULTS_DIR / "zonation" / "spatial_with_zonation.h5ad"
    if not zon_h5ad.exists():
        print("  ERROR: Run 04a_define_zonation.py first")
        sys.exit(1)
    adata = sc.read_h5ad(zon_h5ad)
    print(f"  Loaded: {adata.n_obs} spots with zonation")

    # Identify conditions
    conditions = adata.obs["condition"].unique().tolist()
    cond_healthy = [c for c in conditions if "healthy" in c.lower() or c == "Healthy"]
    cond_masld = [c for c in conditions if c not in cond_healthy]
    if not cond_healthy or not cond_masld:
        print("  ERROR: Need both healthy and MASLD conditions")
        sys.exit(1)
    print(f"  Conditions: {cond_healthy[0]} vs {cond_masld[0]}")

    # Load dream DEGs and filter to those in spatial data
    dream_degs = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.0)
    if "symbol" in dream_degs.columns:
        deg_genes = dream_degs["symbol"].dropna().unique().tolist()
    else:
        ensembl_to_symbol = build_ensembl_to_symbol_map()
        deg_genes = [ensembl_to_symbol.get(g) for g in dream_degs.iloc[:, 0]]
        deg_genes = [g for g in deg_genes if g is not None]

    available_genes = [g for g in deg_genes if g in adata.var_names]
    print(f"  DEGs in spatial data: {len(available_genes)}/{len(deg_genes)}")

    # Zonation bins
    bin_labels = ["PP1", "PP2", "Mid", "PC2", "PC1"]

    # Run zone-aware DE
    print("\n  Running zonation-aware DE...")
    results = zonation_aware_de(
        adata, available_genes, bin_labels,
        cond_healthy[0], cond_masld[0],
    )

    if len(results) == 0:
        print("  WARNING: No results produced (insufficient samples per bin)")
        return

    # Apply BH FDR correction within each bin
    for bl in bin_labels:
        pcol = f"{bl}_pval"
        if pcol in results.columns:
            pvals = results[pcol].fillna(1.0).values
            _, padj, _, _ = multipletests(pvals, method="fdr_bh")
            results[f"{bl}_padj_bh"] = padj

    # Recompute n_sig_bins and classification using BH-corrected p-values
    # (the initial classification used raw p < 0.05, which inflates significance)
    padj_cols = [f"{bl}_padj_bh" for bl in bin_labels if f"{bl}_padj_bh" in results.columns]
    lfc_cols = [f"{bl}_lfc" for bl in bin_labels if f"{bl}_lfc" in results.columns]
    for idx in results.index:
        tested_bins = [bl for bl in bin_labels if f"{bl}_padj_bh" in results.columns
                       and pd.notna(results.at[idx, f"{bl}_padj_bh"])]
        n_tested = len(tested_bins)
        if n_tested == 0:
            results.at[idx, "n_sig_bins"] = 0
            results.at[idx, "classification"] = "Not_significant"
            continue
        n_sig = sum(1 for bl in tested_bins if results.at[idx, f"{bl}_padj_bh"] < 0.05)
        lfcs = [abs(results.at[idx, f"{bl}_lfc"]) for bl in tested_bins
                if pd.notna(results.at[idx, f"{bl}_lfc"])]
        max_abs_lfc = max(lfcs) if lfcs else 0
        min_abs_lfc = min(lfcs) if lfcs else 0
        results.at[idx, "n_sig_bins"] = n_sig
        results.at[idx, "max_abs_lfc"] = max_abs_lfc
        if n_sig == 0:
            results.at[idx, "classification"] = "Not_significant"
        elif n_sig <= 2 and n_sig < n_tested:
            results.at[idx, "classification"] = "Zone_restricted"
        elif n_sig == n_tested and max_abs_lfc > 2 * min_abs_lfc:
            results.at[idx, "classification"] = "Zone_amplified"
        elif n_sig == n_tested:
            results.at[idx, "classification"] = "Pan_lobular"
        else:
            results.at[idx, "classification"] = "Partial"

    save_csv(results, "zonation_aware_de.csv", subdir="zonation")

    # Extract zone-restricted genes
    zone_restricted = results[results["classification"] == "Zone_restricted"]
    zone_amplified = results[results["classification"] == "Zone_amplified"]
    pan_lobular = results[results["classification"] == "Pan_lobular"]

    print(f"\n  Results ({len(results)} genes tested):")
    for cls in ["Zone_restricted", "Zone_amplified", "Pan_lobular", "Partial", "Not_significant"]:
        n = (results["classification"] == cls).sum()
        pct = n / len(results) * 100
        print(f"    {cls}: {n} ({pct:.1f}%)")

    if len(zone_restricted) > 0:
        save_csv(zone_restricted, "zone_specific_masld_genes.csv", subdir="zonation")
        print(f"\n  Top zone-restricted genes:")
        for _, row in zone_restricted.nlargest(10, "max_abs_lfc").iterrows():
            sig_bins = [bl for bl in bin_labels if row.get(f"{bl}_pval", 1) < 0.05]
            print(f"    {row['gene']}: sig in {sig_bins}, max|LFC|={row['max_abs_lfc']:.2f}")

    print_header("09a: Complete")


if __name__ == "__main__":
    main()
