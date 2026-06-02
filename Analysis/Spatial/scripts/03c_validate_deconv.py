#!/usr/bin/env python3
"""
03c_validate_deconv.py — Validate spatial deconvolution against bulk results.

Key validations:
1. Aggregate spatial cell type proportions → compare to MuSiC/BayesPrism
2. Hepatocyte-intrinsic DEG co-localization with hepatocyte spots
3. Bulk-spatial proportion correlation (Spearman)

SLURM: --partition=cpu --cpus=8 --mem=64G --time=2:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_deconvolved_adata, load_deconv_scores,
    build_ensembl_to_symbol_map, save_csv, print_header, print_step,
)
from spatial_stats import ensure_lognorm, compare_two_groups_many_genes


def validate_hepatocyte_intrinsic(adata_vis, deconv_scores, hep_threshold):
    """Validate bulk hepatocyte-intrinsic DEGs against spatial cell2location.

    Uses `category == "Hepatocyte_intrinsic"` from deconv_attribution_scores.csv,
    maps Ensembl IDs to symbols, then tests enrichment in hepatocyte-rich spots
    with Wilcoxon rank-sum + BH FDR correction.
    """
    # Filter to Hepatocyte_intrinsic genes
    if "category" not in deconv_scores.columns:
        print("  WARNING: No 'category' column in deconv scores")
        return pd.DataFrame()

    hep_genes = deconv_scores[deconv_scores["category"] == "Hepatocyte_intrinsic"].copy()
    print(f"  Hepatocyte_intrinsic genes from bulk: {len(hep_genes)}")
    if len(hep_genes) == 0:
        return pd.DataFrame()

    # Map Ensembl IDs to symbols
    ensembl_to_symbol = build_ensembl_to_symbol_map()
    hep_genes["symbol"] = hep_genes["gene"].map(ensembl_to_symbol)
    # Fallback: strip version suffix
    unmapped = hep_genes["symbol"].isna()
    if unmapped.any():
        hep_genes.loc[unmapped, "symbol"] = (
            hep_genes.loc[unmapped, "gene"]
            .str.split(".").str[0]
            .map(ensembl_to_symbol)
        )
    hep_genes = hep_genes.dropna(subset=["symbol"])
    print(f"  After symbol mapping: {len(hep_genes)} genes")

    # Get hepatocyte abundance per spot
    abundances = adata_vis.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundances, pd.DataFrame):
        print("  WARNING: Abundances not in DataFrame format")
        return pd.DataFrame()

    hep_col_name = [c for c in abundances.columns if "hepatocyte" in c.lower()]
    if not hep_col_name:
        print("  WARNING: No 'Hepatocyte' column in cell2location abundances")
        return pd.DataFrame()
    hep_col_name = hep_col_name[0]
    total = abundances.sum(axis=1)
    hep_frac = abundances[hep_col_name] / total.clip(lower=1e-10)
    hep_spots = hep_frac > hep_threshold
    n_hep = hep_spots.sum()
    n_other = (~hep_spots).sum()
    print(f"  Hep-rich spots (frac>{hep_threshold}): {n_hep}, other: {n_other}")

    if n_hep == 0 or n_other == 0:
        print("  WARNING: Cannot partition spots into hep/other groups")
        return pd.DataFrame()

    # Intersect with spatial var_names
    available = set(adata_vis.var_names)
    hep_genes = hep_genes[hep_genes["symbol"].isin(available)]
    print(f"  Genes present in spatial data: {len(hep_genes)}")

    # Donor-aware test (fixes F016 raw-counts + F018 spot pseudoreplication).
    # hep-rich vs other is a WITHIN-donor factor (each donor's slice has both),
    # so expr ~ hep_status + (1|donor) via spatial_stats (MixedLM; donor-mean
    # Mann-Whitney fallback). Expression is read from log1p-CPM, NOT raw counts.
    lognorm = ensure_lognorm(adata_vis)
    adata_vis.obs["hep_status"] = np.where(hep_spots.values, "hep_rich", "other")
    if adata_vis.obs["hep_status"].nunique() < 2:
        print("  WARNING: Cannot partition spots into hep/other groups")
        return pd.DataFrame()

    genes = list(hep_genes["symbol"])
    res = compare_two_groups_many_genes(
        adata_vis, genes, group_col="hep_status", donor_col="sample_id",
        layer=lognorm, fdr=True,
    )
    if len(res) == 0:
        return res

    # `effect` = mean(other) - mean(hep_rich) on the log1p-CPM scale (groups are
    # sorted alphabetically: 'hep_rich' < 'other'). Flip so positive = hep higher.
    res = res.rename(columns={"gene": "symbol"})
    res["hep_minus_other_lognorm"] = -res["effect"]
    res["spatial_fc"] = res["hep_minus_other_lognorm"]  # downstream-compat alias (now a log-scale diff)
    meta = hep_genes.drop_duplicates("symbol").set_index("symbol")
    res["ensembl_id"] = res["symbol"].map(meta["gene"])
    res["category"] = res["symbol"].map(meta["category"])
    if "attribution_raw" in hep_genes.columns:
        res["attribution_raw"] = res["symbol"].map(meta["attribution_raw"])
    if "logFC_adj" in hep_genes.columns:
        res["bulk_logFC_adj"] = res["symbol"].map(meta["logFC_adj"])
    res["validated"] = (res["padj"] < 0.05) & (res["hep_minus_other_lognorm"] > 0)

    n_don = int(res["n_donors"].dropna().iloc[0]) if res["n_donors"].notna().any() else 0
    pct = res["validated"].mean() * 100
    print(f"  Hepatocyte-intrinsic validation (donor-aware, n_donors={n_don}): "
          f"{pct:.1f}% of {len(res)} genes validated (padj<0.05 & hep>other)")
    print(f"  Validated: {int(res['validated'].sum())}/{len(res)} | "
          f"test methods: {dict(res['method'].value_counts())}")
    return res


def compute_bulk_spatial_correlation(hep_validation_df):
    """Compute Spearman correlation between bulk attribution_raw and spatial_fc."""
    if len(hep_validation_df) == 0:
        return
    valid = hep_validation_df[["attribution_raw", "spatial_fc"]].dropna()
    if len(valid) < 5:
        print("  Too few genes for bulk-spatial correlation")
        return
    rho, pval = spearmanr(valid["attribution_raw"], valid["spatial_fc"])
    print(f"\n  Bulk-spatial correlation:")
    print(f"    Spearman rho = {rho:.3f}, p = {pval:.2e} (n={len(valid)})")
    r, p_pearson = pearsonr(valid["attribution_raw"], valid["spatial_fc"])
    print(f"    Pearson r = {r:.3f}, p = {p_pearson:.2e}")


def compute_proportion_correlation(adata_vis):
    """Compute per-sample spatial cell type proportions."""
    abundances = adata_vis.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(abundances, pd.DataFrame):
        return pd.DataFrame()

    # Add sample_id
    props = abundances.copy()
    props["sample_id"] = adata_vis.obs["sample_id"].values
    sample_props = props.groupby("sample_id").mean()
    # Normalize to proportions
    sample_props = sample_props.div(sample_props.sum(axis=1), axis=0)
    return sample_props


def main():
    print_header("03c: Validate Spatial Deconvolution")

    config = load_config()
    output_dir = RESULTS_DIR / "cell2location"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load deconvolved spatial data
    adata_vis = load_deconvolved_adata()
    print(f"  Loaded: {adata_vis.n_obs} spots, {adata_vis.n_vars} genes")

    # 1. Cell type proportion summary
    print("\n  --- Cell type proportions ---")
    props = compute_proportion_correlation(adata_vis)
    if len(props) > 0:
        save_csv(props, "spatial_cell_type_proportions.csv",
                 subdir="cell2location")
        print(f"  Mean proportions across {len(props)} samples:")
        for ct in props.columns:
            print(f"    {ct}: {props[ct].mean():.3f} ± {props[ct].std():.3f}")

    # 2. Hepatocyte-intrinsic validation
    print("\n  --- Hepatocyte-intrinsic DEG validation ---")
    try:
        deconv_scores = load_deconv_scores()
        hep_threshold = config.get("validation", {}).get("hep_abundance_threshold", 0.5)
        hep_validation = validate_hepatocyte_intrinsic(adata_vis, deconv_scores, hep_threshold=hep_threshold)
        if len(hep_validation) > 0:
            save_csv(hep_validation, "hep_intrinsic_validation.csv",
                     subdir="cell2location")
            compute_bulk_spatial_correlation(hep_validation)
    except Exception as e:
        print(f"  WARNING: Could not validate hepatocyte DEGs: {e}")
        import traceback
        traceback.print_exc()

    # 3. Dominant cell type distribution
    print("\n  --- Dominant cell type per spot ---")
    if "cell_type_dominant" in adata_vis.obs.columns:
        ct_counts = adata_vis.obs["cell_type_dominant"].value_counts()
        for ct, n in ct_counts.items():
            pct = n / adata_vis.n_obs * 100
            print(f"    {ct}: {n} ({pct:.1f}%)")

        # Per-condition comparison
        for cond in adata_vis.obs["condition"].unique():
            mask = adata_vis.obs["condition"] == cond
            print(f"\n  Condition: {cond} ({mask.sum()} spots)")
            ct_cond = adata_vis.obs.loc[mask, "cell_type_dominant"].value_counts()
            for ct, n in ct_cond.head(5).items():
                pct = n / mask.sum() * 100
                print(f"    {ct}: {n} ({pct:.1f}%)")

    print_header("03c: Complete")


if __name__ == "__main__":
    main()
