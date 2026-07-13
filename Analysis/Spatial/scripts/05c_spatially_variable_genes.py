#!/usr/bin/env python3
"""
05c_spatially_variable_genes.py — Identify spatially variable genes (SVGs).

Uses Moran's I to detect genes with significant spatial autocorrelation.
Compares SVGs between healthy and MASLD to find disease-emergent spatial patterns.

SLURM: --partition=cpu --cpus=16 --mem=128G --time=8:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_spatial_adata, load_dream_degs,
    save_csv, print_header, print_step,
)
from spatial_stats import morans_i_by_donor

# Donor (INDIVIDUAL) column. GSE192741 has 5 slides (sample_id
# JBO014/015/018/019/022) but only 4 distinct individuals: JBO014 and JBO015 are
# both individual H35 (Steatotic), so slide-level counting double-counts H35 and
# inflates the Steatotic unit n to 3. The true experimental units are 2 Healthy
# (H36=JBO018, H38=JBO022) and 2 Steatotic (H35=JBO014+JBO015, H37=JBO019)
# individuals — i.e. n=2-vs-2, not 2-vs-3 slides. Aggregating Moran's I to the
# INDIVIDUAL (not the slide) fixes the remaining pseudoreplication
# (F009/F044/F183/C021): the per-individual mean I and majority-significance vote
# are now over biological replicates. NOTE: morans_i_by_donor groups by this
# column and builds ONE generic-kNN graph per group; for H35 that pools spots
# from the two physically distinct JBO014/JBO015 slides into a single graph (a
# residual F229-style cross-slide-edge concern). The clean fix — build the
# spatial graph per slide and aggregate to the individual — lives in
# spatial_stats.morans_i_by_donor (outside this file's scope).
DONOR_COL = "individual"


def identify_svgs(adata, condition, config):
    """Identify spatially variable genes for one condition (donor-aware).

    F009/F044/F229/F183/C021: previously this pooled ALL spots of a condition
    into ONE generic-kNN spatial graph (connecting spots across physically
    distinct slides whose coordinates share a frame) and tested thousands of
    non-independent spots as if independent (n=spots, not n=donors). We now
    aggregate Moran's I to the INDIVIDUAL (donor=individual; JBO014+JBO015 = one
    individual H35, so the Steatotic unit n is 2 not 3) via
    spatial_stats.morans_i_by_donor.

    Output column meanings CHANGED but names are preserved for the downstream
    atlas contract (06_integration.py maps 'I'->spatial_morans_i, 'svg'->
    spatial_is_svg): 'I' is now the per-donor MEAN Moran's I (not a single
    pooled statistic); 'svg' now requires autocorrelation significance in a
    MAJORITY of donor slides (frac_donors_sig >= 0.5) plus mean I above
    threshold, rather than a spot-level permutation p-value.
    """
    adata_sub = adata[adata.obs["condition"] == condition].copy()
    if adata_sub.n_obs < 100:
        print(f"    WARNING: Only {adata_sub.n_obs} spots for {condition}")
        return pd.DataFrame()
    if DONOR_COL not in adata_sub.obs.columns:
        print(f"    WARNING: '{DONOR_COL}' missing; cannot donor-block SVGs")
        return pd.DataFrame()

    # Use HVGs to limit computation
    genes_to_test = adata_sub.var_names[adata_sub.var["highly_variable"]].tolist()
    n_donors = adata_sub.obs[DONOR_COL].nunique()
    print(f"    Testing {len(genes_to_test)} HVGs across {n_donors} donor slide(s)...")

    # Per-slide Moran's I, aggregated across donors (per-slice neighbor graphs).
    donor_df = morans_i_by_donor(
        adata_sub, genes=genes_to_test, donor_col=DONOR_COL,
        n_perms=config["n_perms"], n_neighs=6,
    )
    if donor_df.empty:
        print(f"    WARNING: no per-donor Moran's I computed for {condition}")
        return pd.DataFrame()

    svg_df = donor_df.copy()
    # Preserve the downstream column name 'I' but it is now the donor-mean I.
    svg_df["I"] = svg_df["morans_i_mean"]
    if "frac_donors_sig" not in svg_df.columns:
        # Single-donor condition (or no p-values): fall back to mean-I threshold,
        # significance-in-majority cannot be assessed.
        svg_df["frac_donors_sig"] = np.nan
    # SVG = autocorrelation significant in a majority of donor slides AND mean I
    # above the biological-floor threshold. With n<2 donors frac is NaN, so we
    # require only the mean-I floor (descriptive, single-slide).
    maj_sig = (svg_df["frac_donors_sig"].fillna(1.0) >= 0.5)
    svg_df["svg"] = maj_sig & (svg_df["I"] > config["min_morans_i"])
    svg_df = svg_df.sort_values("I", ascending=False)
    return svg_df


def differential_svgs(svg_healthy, svg_masld):
    """Find genes gaining/losing spatial structure in disease."""
    shared = svg_healthy.index.intersection(svg_masld.index)
    if len(shared) == 0:
        return pd.DataFrame()

    diff = pd.DataFrame({
        "morans_I_healthy": svg_healthy.loc[shared, "I"],
        "morans_I_masld": svg_masld.loc[shared, "I"],
        "svg_healthy": svg_healthy.loc[shared, "svg"],
        "svg_masld": svg_masld.loc[shared, "svg"],
    })
    diff["delta_I"] = diff["morans_I_masld"] - diff["morans_I_healthy"]
    diff["category"] = "stable"
    diff.loc[diff["svg_masld"] & ~diff["svg_healthy"], "category"] = "disease_emergent_SVG"
    diff.loc[~diff["svg_masld"] & diff["svg_healthy"], "category"] = "disease_lost_SVG"
    return diff.sort_values("delta_I", key=abs, ascending=False)


def main():
    print_header("05c: Spatially Variable Gene Detection")

    config = load_config()
    svg_config = config["svg"]
    output_dir = RESULTS_DIR / "svg"
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_spatial_adata()
    print(f"  Loaded: {adata.n_obs} spots, {adata.n_vars} genes")

    conditions = adata.obs["condition"].unique().tolist()
    svg_results = {}

    for condition in conditions:
        print(f"\n  Identifying SVGs: {condition}...")
        svg_df = identify_svgs(adata, condition, svg_config)
        if len(svg_df) > 0:
            svg_results[condition] = svg_df
            n_svg = svg_df["svg"].sum()
            save_csv(svg_df, f"svgs_{condition}.csv", subdir="svg")
            print(f"    SVGs: {n_svg} / {len(svg_df)} tested "
                  f"(donor-aware: significant in majority of slides)")
            print(f"    Top 10 by donor-mean Moran's I:")
            for gene, row in svg_df.head(10).iterrows():
                frac = row.get("frac_donors_sig", np.nan)
                print(f"      {gene}: mean_I={row['I']:.3f}, "
                      f"frac_donors_sig={frac:.2f}, n_donors={int(row['n_donors'])}")

    # Differential SVGs
    if len(svg_results) >= 2:
        cond_h = [c for c in conditions if "healthy" in c.lower() or c == "Healthy"]
        cond_m = [c for c in conditions if c not in cond_h]
        if cond_h and cond_m:
            print(f"\n  Differential SVGs: {cond_h[0]} vs {cond_m[0]}...")
            diff = differential_svgs(svg_results[cond_h[0]], svg_results[cond_m[0]])
            save_csv(diff, "differential_svgs.csv", subdir="svg")

            for cat in ["disease_emergent_SVG", "disease_lost_SVG", "stable"]:
                n = (diff["category"] == cat).sum()
                print(f"    {cat}: {n}")

    # Cross-reference with dream DEGs
    print("\n  Cross-referencing with bulk DEGs...")
    try:
        dream_degs = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.8)
        gene_col = "symbol" if "symbol" in dream_degs.columns else dream_degs.columns[0]
        deg_set = set(dream_degs[gene_col])

        for condition, svg_df in svg_results.items():
            svg_genes = set(svg_df[svg_df["svg"]].index)
            overlap = svg_genes & deg_set
            pct = len(overlap) / max(len(svg_genes), 1) * 100
            print(f"    {condition}: {len(overlap)}/{len(svg_genes)} SVGs are dream DEGs ({pct:.1f}%)")
    except Exception as e:
        print(f"    WARNING: Could not cross-reference: {e}")

    print_header("05c: Complete")


if __name__ == "__main__":
    main()
