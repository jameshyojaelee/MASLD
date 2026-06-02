#!/usr/bin/env python3
"""
11_spatial_bulk_validation.py — Direct spatial validation of bulk RNA-seq findings.

Core question: Do bulk MASLD DEGs show spatial dysregulation in independent
Visium datasets (Guilliams et al. + Vu et al.)?

Approach:
  1. For each Visium dataset, compute per-gene spatial metrics
     (Moran's I, mean expression in disease vs control spots)
  2. Test whether dream DEGs have higher spatial autocorrelation than non-DEGs
  3. Test whether Conserved genes are spatially coherent
  4. Correlate bulk dream logFC/tstat with spatial metrics
  5. Compare Guilliams et al. vs Vu et al. for cross-dataset replication

No consensus thresholds, no per-array splitting — just direct gene-level
comparison between bulk and spatial.

SLURM: --partition=cpu --cpus=16 --mem=128G --time=8:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.stats import spearmanr, mannwhitneyu, pearsonr
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_dream_degs,
    load_conserved, build_ensembl_to_symbol_map,
    save_csv, print_header,
)

OUTPUT_DIR = RESULTS_DIR / "validation_bulk"


def compute_spatial_metrics(adata, dataset_label, n_perms=999, n_jobs=8):
    """Donor-aware Moran's I for all genes in a spatial dataset.

    F176 fix: the atlas spatial_max_I / spatial_sig (assembled by 27a from this
    output) previously came from a single pooled-spot Moran's I that treated
    thousands of spots as independent. This routes the computation through
    spatial_stats.morans_i_by_donor — a per-SLICE spatial graph (so spots from
    physically distinct sections are never connected) computed per donor and
    aggregated. `I` is the donor-mean Moran's I; `padj_bh` is set so the
    downstream `padj_bh < 0.05` test in 27a means "spatially autocorrelated in a
    MAJORITY of donor slices" (frac_donors_sig >= 0.5) instead of a pooled p.
    """
    from spatial_stats import morans_i_by_donor
    print(f"  Computing donor-aware spatial metrics for {dataset_label}...")
    print(f"    {adata.n_obs} spots × {adata.n_vars} genes")

    genes_to_test = adata.var_names.tolist()
    if len(genes_to_test) > 5000:
        if "highly_variable" in adata.var.columns:
            genes_to_test = adata.var_names[adata.var["highly_variable"]].tolist()
        else:
            genes_to_test = genes_to_test[:5000]
    donor_col = "sample_id" if "sample_id" in adata.obs.columns else None

    if donor_col is None:
        # Single-slice fallback (no donor column): pooled Moran's I, clearly the
        # degenerate case (1 unit). Kept so the function never hard-fails.
        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)
        sq.gr.spatial_autocorr(adata, mode="moran", genes=genes_to_test,
                               n_perms=n_perms, n_jobs=n_jobs)
        result = adata.uns["moranI"].copy()
        result["padj_bh"] = multipletests(
            result["pval_norm"].fillna(1.0).values, method="fdr_bh")[1]
        result["dataset"] = dataset_label
        return result

    print(f"    Testing {len(genes_to_test)} genes (donor-aware, per-slice graphs)...")
    md = morans_i_by_donor(adata, genes_to_test, donor_col, n_perms=n_perms)
    if md.empty:
        return pd.DataFrame(columns=["I", "padj_bh", "dataset"])
    result = pd.DataFrame(index=md.index)
    result["I"] = md["morans_i_mean"]
    frac = md["frac_donors_sig"] if "frac_donors_sig" in md.columns else 0.0
    # padj_bh encodes the donor-majority rule for 27a's `padj_bh < 0.05` test:
    # 0.0 (=> sig) when autocorrelated in a majority of donor slices, else 1.0.
    result["padj_bh"] = (frac < 0.5).astype(float)
    result["frac_donors_sig"] = frac
    result["n_donors"] = md["n_donors"]
    result["dataset"] = dataset_label

    n_sig = int((result["padj_bh"] < 0.05).sum())
    n_don = int(md["n_donors"].max()) if "n_donors" in md.columns else 0
    print(f"    {n_sig} genes donor-robust spatial autocorr (majority of {n_don} slices)")
    return result


def load_bulk_results():
    """Load dream DEG results with gene symbols."""
    dream = load_dream_degs(padj_thresh=1.0, lfc_thresh=0.0)  # all genes
    gene_col = "symbol" if "symbol" in dream.columns else dream.columns[0]

    # Build clean table: symbol -> logFC, tstat, padj
    bulk = dream[["gene", "logFC", gene_col]].copy()
    if "t" in dream.columns:
        bulk["tstat"] = dream["t"]
    elif "tstat" in dream.columns:
        bulk["tstat"] = dream["tstat"]
    padj_col = "padj" if "padj" in dream.columns else "adj.P.Val"
    bulk["dream_padj"] = dream[padj_col]
    bulk = bulk.rename(columns={gene_col: "symbol"})
    bulk = bulk.dropna(subset=["symbol"])
    bulk = bulk.drop_duplicates(subset=["symbol"], keep="first")
    print(f"  Bulk: {len(bulk)} genes with symbols")
    return bulk


def merge_bulk_spatial(bulk_df, spatial_df):
    """Merge bulk and spatial results on gene symbol."""
    merged = spatial_df.merge(
        bulk_df[["symbol", "logFC", "tstat", "dream_padj"]],
        left_index=True, right_on="symbol", how="inner"
    )
    merged["is_deg"] = merged["dream_padj"] < 0.1
    merged["is_deg_strong"] = (merged["dream_padj"] < 0.1) & (merged["logFC"].abs() > 0.5)
    merged["is_up"] = (merged["dream_padj"] < 0.1) & (merged["logFC"] > 0)
    merged["is_down"] = (merged["dream_padj"] < 0.1) & (merged["logFC"] < 0)
    print(f"    Merged: {len(merged)} genes ({merged['is_deg'].sum()} DEGs)")
    return merged


def test_deg_spatial_enrichment(merged, dataset_label):
    """Test: do DEGs have higher Moran's I than non-DEGs?"""
    results = []

    for label, mask_col in [("All DEGs (padj<0.1)", "is_deg"),
                             ("Strong DEGs (|LFC|>0.5)", "is_deg_strong"),
                             ("Up-regulated DEGs", "is_up"),
                             ("Down-regulated DEGs", "is_down")]:
        deg_I = merged.loc[merged[mask_col], "I"].dropna()
        nondeg_I = merged.loc[~merged[mask_col], "I"].dropna()

        if len(deg_I) < 10 or len(nondeg_I) < 10:
            continue

        stat, pval = mannwhitneyu(deg_I, nondeg_I, alternative="greater")

        results.append({
            "dataset": dataset_label,
            "comparison": label,
            "n_deg": len(deg_I),
            "n_nondeg": len(nondeg_I),
            "median_I_deg": deg_I.median(),
            "median_I_nondeg": nondeg_I.median(),
            "fold_enrichment": deg_I.median() / max(nondeg_I.median(), 1e-6),
            "mannwhitney_pval": pval,
        })

        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
        print(f"    {label}: median I = {deg_I.median():.3f} vs {nondeg_I.median():.3f} "
              f"(fold={deg_I.median()/max(nondeg_I.median(), 1e-6):.2f}x, p={pval:.2e}) {sig}")

    return pd.DataFrame(results)


def test_conserved(merged, cc_genes, dataset_label):
    """Test: do Conserved genes have higher Moran's I?"""
    merged["is_cc"] = merged["symbol"].isin(cc_genes)
    cc_I = merged.loc[merged["is_cc"], "I"].dropna()
    other_I = merged.loc[~merged["is_cc"], "I"].dropna()

    if len(cc_I) < 5:
        print(f"    Conserved: only {len(cc_I)} genes found, skipping")
        return None

    stat, pval = mannwhitneyu(cc_I, other_I, alternative="greater")
    print(f"    Conserved: median I = {cc_I.median():.3f} vs {other_I.median():.3f} "
          f"(n={len(cc_I)}, p={pval:.2e})")

    return {
        "dataset": dataset_label,
        "comparison": "Conserved",
        "n_cc": len(cc_I),
        "n_other": len(other_I),
        "median_I_cc": cc_I.median(),
        "median_I_other": other_I.median(),
        "fold_enrichment": cc_I.median() / max(other_I.median(), 1e-6),
        "mannwhitney_pval": pval,
    }


def correlate_bulk_spatial(merged, dataset_label):
    """Correlate bulk logFC/tstat with spatial Moran's I."""
    results = []

    for bulk_metric in ["logFC", "tstat"]:
        if bulk_metric not in merged.columns:
            continue
        mask = merged[bulk_metric].notna() & merged["I"].notna()
        x = merged.loc[mask, bulk_metric].values
        y = merged.loc[mask, "I"].values

        if len(x) < 50:
            continue

        # Overall correlation
        rho, rho_p = spearmanr(x, y)
        r, r_p = pearsonr(x, y)

        # |bulk_metric| vs I (absolute effect → spatial structure)
        rho_abs, rho_abs_p = spearmanr(np.abs(x), y)

        results.append({
            "dataset": dataset_label,
            "bulk_metric": bulk_metric,
            "n_genes": len(x),
            "spearman_rho": rho,
            "spearman_pval": rho_p,
            "pearson_r": r,
            "pearson_pval": r_p,
            "abs_spearman_rho": rho_abs,
            "abs_spearman_pval": rho_abs_p,
        })

        print(f"    {bulk_metric} vs Moran's I: rho={rho:.3f} (p={rho_p:.2e}), "
              f"|{bulk_metric}| vs I: rho={rho_abs:.3f} (p={rho_abs_p:.2e})")

    return pd.DataFrame(results)


def main():
    print_header("11: Direct Spatial Validation of Bulk RNA-seq")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    config = load_config()
    svg_config = config.get("svg", {"n_perms": 100, "n_jobs": 8})

    # ── Load bulk results ──
    print("  Loading bulk dream results...")
    bulk = load_bulk_results()

    print("  Loading Conserved...")
    cc_genes = set(load_conserved())
    print(f"    {len(cc_genes)} Conserved genes")

    # ── Load spatial datasets ──
    datasets = {}

    # Guilliams et al. (GSE192741)
    guilliams_path = RESULTS_DIR / "preprocessed" / "merged_spatial.h5ad"
    if guilliams_path.exists():
        adata_g = sc.read_h5ad(guilliams_path)
        # Filter to human only if mixed species
        if "species" in adata_g.obs.columns:
            adata_g = adata_g[adata_g.obs["species"] == "human"].copy()
        datasets["Guilliams et al."] = adata_g

    # Vu et al.
    vu_path = RESULTS_DIR / "preprocessed" / "merged_spatial_vu.h5ad"
    if vu_path.exists():
        datasets["Vu et al."] = sc.read_h5ad(vu_path)

    if not datasets:
        print("  ERROR: No spatial datasets found!")
        sys.exit(1)

    # ── Run analysis per dataset ──
    all_enrichment = []
    all_correlations = []
    all_cc_tests = []
    all_spatial = {}

    for ds_label, adata in datasets.items():
        print_header(f"Dataset: {ds_label}")

        # Compute Moran's I
        spatial_df = compute_spatial_metrics(
            adata, ds_label,
            n_perms=svg_config.get("n_perms", 100),
            n_jobs=svg_config.get("n_jobs", 8)
        )

        # Merge with bulk
        print(f"\n  Merging with bulk results...")
        merged = merge_bulk_spatial(bulk, spatial_df)
        all_spatial[ds_label] = merged

        # Save per-dataset spatial metrics
        save_csv(merged, f"spatial_bulk_merged_{ds_label.replace(' ', '_').replace('.', '')}.csv",
                 subdir="validation_bulk")

        # Test DEG enrichment
        print(f"\n  Testing DEG spatial enrichment...")
        enrich = test_deg_spatial_enrichment(merged, ds_label)
        all_enrichment.append(enrich)

        # Test Conserved
        print(f"\n  Testing Conserved...")
        cc_result = test_conserved(merged, cc_genes, ds_label)
        if cc_result:
            all_cc_tests.append(cc_result)

        # Correlate bulk metrics with spatial
        print(f"\n  Correlating bulk vs spatial...")
        corr = correlate_bulk_spatial(merged, ds_label)
        all_correlations.append(corr)

    # ── Cross-dataset concordance ──
    if len(all_spatial) == 2:
        print_header("Cross-Dataset Concordance")
        ds_labels = list(all_spatial.keys())
        m1 = all_spatial[ds_labels[0]].set_index("symbol")
        m2 = all_spatial[ds_labels[1]].set_index("symbol")
        shared = m1.index.intersection(m2.index)
        print(f"  Shared genes: {len(shared)}")

        I1 = m1.loc[shared, "I"]
        I2 = m2.loc[shared, "I"]
        mask = I1.notna() & I2.notna()
        rho, rho_p = spearmanr(I1[mask], I2[mask])
        print(f"  Moran's I correlation: rho={rho:.3f}, p={rho_p:.2e}")

        cross = pd.DataFrame({
            "gene": shared[mask],
            f"morans_I_{ds_labels[0]}": I1[mask].values,
            f"morans_I_{ds_labels[1]}": I2[mask].values,
        })
        save_csv(cross, "cross_dataset_morans_i.csv", subdir="validation_bulk")

    # ── Save combined results ──
    print_header("Summary")

    enrichment_df = pd.concat(all_enrichment, ignore_index=True)
    save_csv(enrichment_df, "deg_spatial_enrichment.csv", subdir="validation_bulk")

    if all_cc_tests:
        cc_df = pd.DataFrame(all_cc_tests)
        save_csv(cc_df, "conserved_spatial.csv", subdir="validation_bulk")

    corr_df = pd.concat(all_correlations, ignore_index=True)
    save_csv(corr_df, "bulk_spatial_correlation.csv", subdir="validation_bulk")

    # Print summary table
    print("\n  DEG Spatial Enrichment:")
    for _, row in enrichment_df.iterrows():
        sig = "***" if row["mannwhitney_pval"] < 0.001 else "ns"
        print(f"    {row['dataset']:20s} | {row['comparison']:30s} | "
              f"fold={row['fold_enrichment']:.2f}x | p={row['mannwhitney_pval']:.2e} {sig}")

    if all_cc_tests:
        print("\n  Conserved:")
        for row in all_cc_tests:
            print(f"    {row['dataset']:20s} | fold={row['fold_enrichment']:.2f}x | "
                  f"p={row['mannwhitney_pval']:.2e}")

    print("\n  Bulk-Spatial Correlations:")
    for _, row in corr_df.iterrows():
        print(f"    {row['dataset']:20s} | {row['bulk_metric']:8s} | "
              f"rho={row['spearman_rho']:.3f} | |metric| rho={row['abs_spearman_rho']:.3f}")

    print_header("11: Complete")


if __name__ == "__main__":
    main()
