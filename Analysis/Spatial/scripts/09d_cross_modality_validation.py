#!/usr/bin/env python3
"""
09d_cross_modality_validation.py — Spatial vs bulk deconvolution validation.

Systematically compares spatial cell2location against bulk MuSiC/BayesPrism:
1. Cell type proportion comparison (Spearman rho of rank ordering)
2. Per-cell-type gene co-localization validation (not just hepatocytes)
3. Negative hepatocyte correlation investigation with multiple thresholds

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, mannwhitneyu, pearsonr
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, C2L_PREFIX, load_config, load_deconvolved_adata,
    load_deconv_scores, build_ensembl_to_symbol_map,
    save_csv, print_header, print_step,
)
# F-coverage (2026-06-01): the deconvolved object's .X is RAW UMI counts, so
# per-gene Mann-Whitney / fold-change must run on a log1p-CPM layer, not .X.
from spatial_stats import ensure_lognorm


def _gene_expr(adata, gene, layer):
    """Dense 1-D expression vector for one gene from a named layer."""
    col = adata[:, gene].layers[layer]
    return np.asarray(col.todense()).flatten() if hasattr(col, "todense") \
        else np.asarray(col).flatten()


def validate_per_celltype(adata_vis, deconv_scores, ensembl_to_symbol,
                          hep_threshold=0.5, expr_layer="lognorm"):
    """Validate gene co-localization for each cell type, not just hepatocytes.

    For each category in deconv_attribution_scores (e.g., Hepatocyte_intrinsic,
    Macrophage_intrinsic, etc.), tests whether attributed genes show higher
    expression in spots dominated by that cell type.

    Expression is read from ``expr_layer`` (log1p-CPM), NOT raw ``.X``: the
    deconvolved object's .X is raw UMI counts, so a raw Mann-Whitney would
    conflate library depth with expression (F-coverage 2026-06-01).
    """
    if "category" not in deconv_scores.columns:
        return pd.DataFrame()

    abundances = adata_vis.obsm.get("q05_cell_abundance_w_sf")
    if abundances is None or not isinstance(abundances, pd.DataFrame):
        print("  WARNING: No cell type abundances found")
        return pd.DataFrame()

    # Map Ensembl to symbols
    deconv_scores = deconv_scores.copy()
    deconv_scores["symbol"] = deconv_scores["gene"].map(ensembl_to_symbol)
    unmapped = deconv_scores["symbol"].isna()
    if unmapped.any():
        deconv_scores.loc[unmapped, "symbol"] = (
            deconv_scores.loc[unmapped, "gene"]
            .str.split(".").str[0]
            .map(ensembl_to_symbol)
        )
    deconv_scores = deconv_scores.dropna(subset=["symbol"])

    # Available genes in spatial
    available = set(adata_vis.var_names)

    # Get cell type categories
    categories = deconv_scores["category"].unique()
    total_abundances = abundances.sum(axis=1)

    results = []
    for category in categories:
        cat_genes = deconv_scores[deconv_scores["category"] == category]
        cat_genes_in_spatial = cat_genes[cat_genes["symbol"].isin(available)]

        if len(cat_genes_in_spatial) < 5:
            continue

        # Find matching cell type column in abundances
        cell_type = category.replace("_intrinsic", "").replace("_", " ")
        matching_cols = [c for c in abundances.columns
                        if cell_type.lower() in c.lower()]

        if not matching_cols:
            print(f"    {category}: no matching abundance column")
            continue

        ct_col = matching_cols[0]
        ct_frac = abundances[ct_col] / total_abundances.clip(lower=1e-10)

        # Dynamic threshold: use top quartile as "dominated" spots
        threshold = max(ct_frac.quantile(0.75), 0.1)
        dominated = ct_frac > threshold
        n_dom = dominated.sum()
        n_other = (~dominated).sum()

        if n_dom < 50 or n_other < 50:
            # Use simple median split instead
            threshold = ct_frac.median()
            dominated = ct_frac > threshold
            n_dom = dominated.sum()
            n_other = (~dominated).sum()

        if n_dom < 20 or n_other < 20:
            print(f"    {category}: insufficient spots for comparison")
            continue

        # Test each gene
        n_validated = 0
        n_tested = 0
        pvals = []

        for _, row in cat_genes_in_spatial.iterrows():
            gene = row["symbol"]
            # log1p-CPM layer, not raw .X (depth-confounded) — F-coverage fix.
            expr = _gene_expr(adata_vis, gene, expr_layer)
            expr_dom = expr[dominated.values]
            expr_other = expr[~dominated.values]

            try:
                _, pval = mannwhitneyu(expr_dom, expr_other, alternative="greater")
            except Exception:
                pval = 1.0

            pvals.append(pval)
            n_tested += 1

        # BH correction
        if pvals:
            _, padj, _, _ = multipletests(pvals, method="fdr_bh")
            n_validated = (padj < 0.05).sum()

        pct = n_validated / max(n_tested, 1) * 100
        results.append({
            "category": category,
            "cell_type_column": ct_col,
            "threshold": threshold,
            "n_genes_bulk": len(cat_genes),
            "n_genes_spatial": len(cat_genes_in_spatial),
            "n_tested": n_tested,
            "n_validated": n_validated,
            "pct_validated": pct,
            "n_dominated_spots": n_dom,
            "n_other_spots": n_other,
        })
        print(f"    {category}: {n_validated}/{n_tested} validated ({pct:.1f}%)")

    return pd.DataFrame(results)


def investigate_negative_correlation(adata_vis, deconv_scores, ensembl_to_symbol,
                                      thresholds=[0.5, 0.6, 0.7, 0.8],
                                      expr_layer="lognorm"):
    """Investigate negative bulk-spatial hepatocyte correlation at multiple thresholds.

    The hepatocyte spatial fold-change (mean_hep / mean_other) is computed on
    the log1p-CPM ``expr_layer``, NOT raw ``.X``: the prior raw-count fold-change
    conflated per-spot sequencing depth with expression (F-coverage 2026-06-01),
    which could itself produce the apparent negative concordance.
    """
    if "category" not in deconv_scores.columns:
        return pd.DataFrame()

    hep_genes = deconv_scores[deconv_scores["category"] == "Hepatocyte_intrinsic"].copy()
    hep_genes["symbol"] = hep_genes["gene"].map(ensembl_to_symbol)
    unmapped = hep_genes["symbol"].isna()
    if unmapped.any():
        hep_genes.loc[unmapped, "symbol"] = (
            hep_genes.loc[unmapped, "gene"]
            .str.split(".").str[0]
            .map(ensembl_to_symbol)
        )
    hep_genes = hep_genes.dropna(subset=["symbol"])
    hep_genes = hep_genes[hep_genes["symbol"].isin(adata_vis.var_names)]

    if len(hep_genes) < 10:
        return pd.DataFrame()

    abundances = adata_vis.obsm.get("q05_cell_abundance_w_sf")
    if abundances is None:
        return pd.DataFrame()

    hep_col = [c for c in abundances.columns if "hepatocyte" in c.lower()][0]
    total = abundances.sum(axis=1)
    hep_frac = abundances[hep_col] / total.clip(lower=1e-10)

    results = []
    for threshold in thresholds:
        hep_spots = hep_frac > threshold
        n_hep = hep_spots.sum()
        n_other = (~hep_spots).sum()

        if n_hep < 20 or n_other < 20:
            continue

        # Compute spatial_fc for each gene at this threshold
        spatial_fcs = []
        attr_raws = []

        for _, row in hep_genes.iterrows():
            gene = row["symbol"]
            # log1p-CPM layer, not raw .X (depth-confounded) — F-coverage fix.
            expr = _gene_expr(adata_vis, gene, expr_layer)
            mean_hep = expr[hep_spots.values].mean()
            mean_other = expr[~hep_spots.values].mean()
            fc = mean_hep / max(mean_other, 1e-6)
            spatial_fcs.append(fc)
            attr_raws.append(row.get("attribution_raw", np.nan))

        spatial_fcs = np.array(spatial_fcs)
        attr_raws = np.array(attr_raws)
        valid = np.isfinite(spatial_fcs) & np.isfinite(attr_raws)

        if valid.sum() < 5:
            continue

        rho, p_rho = spearmanr(attr_raws[valid], spatial_fcs[valid])
        r, p_r = pearsonr(attr_raws[valid], spatial_fcs[valid])

        # Rank-biserial: is attribution_raw higher for validated genes (FC>1)?
        validated = spatial_fcs[valid] > 1
        if validated.sum() > 0 and (~validated).sum() > 0:
            try:
                _, rb_p = mannwhitneyu(
                    attr_raws[valid][validated],
                    attr_raws[valid][~validated],
                    alternative="two-sided"
                )
            except Exception:
                rb_p = 1.0
        else:
            rb_p = np.nan

        results.append({
            "hep_threshold": threshold,
            "n_hep_spots": n_hep,
            "n_other_spots": n_other,
            "pct_hep_spots": n_hep / (n_hep + n_other) * 100,
            "n_genes": valid.sum(),
            "spearman_rho": rho,
            "spearman_pval": p_rho,
            "pearson_r": r,
            "pearson_pval": p_r,
            "pct_fc_gt1": (spatial_fcs[valid] > 1).mean() * 100,
            "median_spatial_fc": float(np.median(spatial_fcs[valid])),
            "rank_biserial_pval": rb_p,
        })

    return pd.DataFrame(results)


def main():
    print_header("09d: Cross-Modality Deconvolution Validation")

    config = load_config()
    output_dir = RESULTS_DIR / "cell2location"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    adata_vis = load_deconvolved_adata()
    print(f"  Loaded: {adata_vis.n_obs} spots, {adata_vis.n_vars} genes")

    # Build a log1p-CPM layer once: the deconvolved object's .X is raw UMI
    # counts, so all per-gene co-localization / fold-change tests below must run
    # on this normalized layer (F-coverage 2026-06-01).
    expr_layer = ensure_lognorm(adata_vis)
    print(f"  Using normalized expression layer: '{expr_layer}' (log1p-CPM)")

    deconv_scores = load_deconv_scores()
    print(f"  Bulk deconv scores: {len(deconv_scores)} genes")

    ensembl_to_symbol = build_ensembl_to_symbol_map()

    # 1. Per-cell-type validation
    print("\n  --- Per-Cell-Type Validation ---")
    ct_validation = validate_per_celltype(
        adata_vis, deconv_scores, ensembl_to_symbol,
        hep_threshold=config.get("validation", {}).get("hep_abundance_threshold", 0.5),
        expr_layer=expr_layer,
    )
    if len(ct_validation) > 0:
        save_csv(ct_validation, "per_celltype_concordance.csv", subdir="cell2location")

    # 2. Investigate negative hepatocyte correlation
    print("\n  --- Hepatocyte Correlation Investigation ---")
    corr_investigation = investigate_negative_correlation(
        adata_vis, deconv_scores, ensembl_to_symbol,
        thresholds=[0.5, 0.6, 0.7, 0.8, 0.9],
        expr_layer=expr_layer,
    )
    if len(corr_investigation) > 0:
        save_csv(corr_investigation, "hep_correlation_threshold_sensitivity.csv",
                 subdir="cell2location")

        print(f"\n  Correlation at different thresholds:")
        for _, row in corr_investigation.iterrows():
            print(f"    hep_frac>{row['hep_threshold']}: "
                  f"rho={row['spearman_rho']:.3f} (p={row['spearman_pval']:.2e}), "
                  f"hep spots={row['pct_hep_spots']:.0f}%, "
                  f"FC>1: {row['pct_fc_gt1']:.0f}%")

    # 3. Cell type proportion comparison
    print("\n  --- Cell Type Proportion Comparison ---")
    abundances = adata_vis.obsm.get("q05_cell_abundance_w_sf")
    if abundances is not None and isinstance(abundances, pd.DataFrame):
        # Aggregate to overall proportions
        mean_props = abundances.mean()
        total = mean_props.sum()
        spatial_props = (mean_props / total).sort_values(ascending=False)

        print(f"  Spatial cell type proportions (mean across spots):")
        for ct, prop in spatial_props.items():
            ct_clean = ct.replace(C2L_PREFIX, "") if isinstance(ct, str) else ct
            if prop > 0.01:
                print(f"    {ct_clean}: {prop:.3f}")

        # Save for cross-modality comparison
        props_df = pd.DataFrame({
            "cell_type": [ct.replace(C2L_PREFIX, "") if isinstance(ct, str) else ct
                         for ct in spatial_props.index],
            "spatial_proportion": spatial_props.values,
        })
        save_csv(props_df, "cross_modality_validation.csv", subdir="cell2location")

    print_header("09d: Complete")


if __name__ == "__main__":
    main()
