#!/usr/bin/env python3
"""
04b_deg_zonation_mapping.py — Map MASLD DEGs onto liver zonation gradient.

Classifies each DEG as periportal/pericentral/pan-lobular and computes
zonation disruption scores between healthy and MASLD conditions.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=6:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import kruskal, spearmanr

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_dream_degs, load_conserved,
    build_ensembl_to_symbol_map, save_csv, print_header, print_step,
)
from spatial_stats import ensure_lognorm, spearman_by_donor, _get_expr

# Donor (slide) column. GSE192741 = JBO014/015/018/019/022 (n=5). Zonation
# tests are blocked on this so spots are not treated as independent (F027:
# Visium spots are spatially autocorrelated; effective n = donors, not spots).
DONOR_COL = "sample_id"


def classify_degs_zonation_vectorized(adata, deg_genes):
    """Classify all DEGs' zonation preferences (donor-aware).

    F027: previously the Kruskal-Wallis bin test and the Spearman zonation
    correlation pooled all ~6,546 spots as independent observations, but the
    cohort has only ~5 donors and Visium spots are strongly spatially
    autocorrelated (pseudoreplication). We now block on donor (sample_id):

      * spearman_rho / spearman_pval are the MEAN per-donor Spearman rho and a
        donor-level Wilcoxon p across donors (spatial_stats.spearman_by_donor),
        NOT a spot-level rho/p. Column names are preserved for downstream
        (06_integration.py, 09b/09c, 07_spatial_figures.py).
      * kruskal_stat / kruskal_pval are computed on per-donor-bin MEAN
        expression (one value per donor x bin), so the effective n is donors,
        not spots. The classification thresholds are unchanged but now rest on
        donor-level statistics (descriptive at n~5 donors — see figure caveat).

    Returns DataFrame with one row per gene that exists in adata.
    """
    # Filter to genes present in spatial data
    available = [g for g in deg_genes if g in adata.var_names]
    print(f"  DEGs present in spatial data: {len(available)}/{len(deg_genes)}")
    if not available:
        return pd.DataFrame()

    if DONOR_COL not in adata.obs.columns:
        print(f"  WARNING: '{DONOR_COL}' missing; zonation tests fall back to "
              f"spot-level (pseudoreplicated)")
        donors = np.array(["_all"] * adata.n_obs)
    else:
        donors = adata.obs[DONOR_COL].astype(str).values
    n_donors = len(np.unique(donors))
    print(f"  Donor-blocked zonation tests across {n_donors} donor(s)")

    # Read log1p-CPM, not raw counts, for per-gene correlation tests (depth
    # confound; the deconvolved/zonation h5ad .X may be raw counts).
    layer = ensure_lognorm(adata)

    z_scores = adata.obs["zonation_score"].values
    bins = adata.obs["zonation_bin"].astype(str).values
    bin_labels = ["PP1", "PP2", "Mid", "PC2", "PC1"]

    results = []
    for i, gene in enumerate(available):
        if i % 500 == 0:
            print_step(f"Processing gene {i+1}/{len(available)}")

        expr = _get_expr(adata, gene, layer=layer)

        # Mean expression per bin (pooled, descriptive only).
        bin_means = {}
        for b in bin_labels:
            mask = bins == b
            bin_means[b] = float(expr[mask].mean()) if mask.sum() > 0 else 0.0

        # Donor-aware Kruskal-Wallis across zonation bins: use per-donor-bin
        # mean expression as the observation unit (n = donors x bins, not
        # spots). This removes spot-level pseudoreplication.
        edf = pd.DataFrame({"expr": expr, "bin": bins, "donor": donors})
        cell_means = (edf.groupby(["donor", "bin"])["expr"].mean()
                      .reset_index())
        groups = [cell_means.loc[cell_means["bin"] == b, "expr"].values
                  for b in bin_labels
                  if (cell_means["bin"] == b).sum() > 0]
        groups = [g for g in groups if len(g) > 0]
        try:
            stat, pval = kruskal(*groups) if len(groups) >= 2 else (0.0, 1.0)
        except Exception:
            stat, pval = 0.0, 1.0

        # Donor-blocked Spearman: mean per-donor rho + donor-level Wilcoxon p.
        sb = spearman_by_donor(expr, z_scores, donors, min_per_donor=20)
        rho = sb["mean_rho"]
        rho_p = sb.get("pval", np.nan)
        if not np.isfinite(rho):
            rho = 0.0

        # Classification (thresholds unchanged; now on donor-level statistics).
        if not np.isfinite(pval) or pval > 0.05:
            zonation_class = "Pan-lobular"
        elif rho > 0.2:
            zonation_class = "Pericentral-enriched"
        elif rho < -0.2:
            zonation_class = "Periportal-enriched"
        else:
            zonation_class = "Pan-lobular"

        results.append({
            "gene": gene, "zonation_class": zonation_class,
            "kruskal_stat": stat, "kruskal_pval": pval,
            "spearman_rho": rho, "spearman_pval": rho_p,
            "n_donors": int(sb["n_donors"]),
            **{f"mean_{b}": v for b, v in bin_means.items()},
        })

    return pd.DataFrame(results)


def compute_zonation_disruption(adata, condition_healthy, condition_masld, marker_genes):
    """Compute per-gene zonation disruption between conditions (donor-aware).

    F027: the healthy gradient is pooled over only ~2 donors and the MASLD side
    over ~3, so this is severely underpowered. We now compute the healthy
    reference rho as the MEAN of per-donor (per-slide) Spearman rho rather than
    one pooled spot-level rho, and tag every row with n_healthy_donors /
    underpowered=True so downstream/figures treat it as descriptive, not a
    quantitative disruption test. Expression is read from the log1p-CPM layer,
    not raw counts.
    """
    results = []
    adata_h = adata[adata.obs["condition"] == condition_healthy].copy()
    adata_m = adata[adata.obs["condition"] == condition_masld].copy()
    layer_h = ensure_lognorm(adata_h)
    layer_m = ensure_lognorm(adata_m)

    h_donors = adata_h.obs["sample_id"].astype(str).values
    n_h_donors = int(len(np.unique(h_donors)))

    for gene in marker_genes:
        if gene not in adata.var_names:
            continue

        # Healthy gradient: mean of per-donor (per-slide) Spearman rho so the
        # reference is not a single pooled-spot statistic (F027).
        h_expr = _get_expr(adata_h, gene, layer=layer_h)
        h_z = adata_h.obs["zonation_score"].values
        sb_h = spearman_by_donor(h_expr, h_z, h_donors, min_per_donor=20)
        rho_h = sb_h["mean_rho"]
        if not np.isfinite(rho_h):
            rho_h = 0.0

        # Per MASLD sample (already donor-resolved)
        m_expr_all = _get_expr(adata_m, gene, layer=layer_m)
        for sample in adata_m.obs["sample_id"].unique():
            mask = (adata_m.obs["sample_id"] == sample).values
            d_expr = m_expr_all[mask]
            d_z = adata_m.obs.loc[mask, "zonation_score"].values
            try:
                rho_d = spearmanr(d_expr, d_z).correlation
            except Exception:
                rho_d = 0.0
            if not np.isfinite(rho_d):
                rho_d = 0.0

            # Disruption = absolute difference in zonation correlation
            # (avoids extreme scores when rho_h ≈ 0)
            disruption = abs(rho_h - rho_d)
            results.append({
                "gene": gene, "sample": sample,
                "rho_healthy": rho_h, "rho_masld": rho_d,
                "disruption_score": disruption,
                "n_healthy_donors": n_h_donors,
                "underpowered": True,  # n=2 vs 3 donors — descriptive only
            })

    return pd.DataFrame(results)


def main():
    print_header("04b: DEG Zonation Mapping")

    config = load_config()
    zon_config = config["zonation"]
    output_dir = RESULTS_DIR / "zonation"

    # Load zonation-annotated spatial data
    zon_h5ad = RESULTS_DIR / "zonation" / "spatial_with_zonation.h5ad"
    if not zon_h5ad.exists():
        print("  ERROR: Run 04a_define_zonation.py first")
        sys.exit(1)
    adata = sc.read_h5ad(zon_h5ad)
    print(f"  Loaded: {adata.n_obs} spots with zonation scores")

    # Load ALL significant dream DEGs (padj<0.1, no LFC filter for zonation)
    dream_degs = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.0)
    padj_col = "padj" if "padj" in dream_degs.columns else "adj.P.Val"

    # Use symbol column; map Ensembl IDs if needed
    if "symbol" in dream_degs.columns:
        deg_genes = dream_degs["symbol"].dropna().tolist()
    else:
        ensembl_to_symbol = build_ensembl_to_symbol_map()
        deg_genes = dream_degs[dream_degs.columns[0]].map(ensembl_to_symbol).dropna().tolist()

    deg_genes = list(set(deg_genes))  # deduplicate
    print(f"  Total DEGs to classify: {len(deg_genes)}")

    # Classify all DEGs (vectorized)
    print("\n  Classifying DEG zonation preferences...")
    deg_zon = classify_degs_zonation_vectorized(adata, deg_genes)
    save_csv(deg_zon, "deg_zonation_classification.csv", subdir="zonation")

    # Summary
    if len(deg_zon) > 0:
        print(f"\n  Zonation classification ({len(deg_zon)} genes):")
        for cls, n in deg_zon["zonation_class"].value_counts().items():
            pct = n / len(deg_zon) * 100
            print(f"    {cls}: {n} ({pct:.1f}%)")

    # Cross-reference Conserved
    conserved = load_conserved()
    if conserved and len(deg_zon) > 0:
        core_in_zon = deg_zon[deg_zon["gene"].isin(conserved)]
        print(f"\n  Conserved in classified DEGs: {len(core_in_zon)}")
        if len(core_in_zon) > 0:
            for cls, n in core_in_zon["zonation_class"].value_counts().items():
                print(f"    {cls}: {n}")

    # Zonation disruption analysis
    conditions = adata.obs["condition"].unique()
    if len(conditions) >= 2:
        cond_healthy = [c for c in conditions if "healthy" in c.lower() or c == "Healthy"]
        cond_masld = [c for c in conditions if c not in cond_healthy]
        if cond_healthy and cond_masld:
            print(f"\n  Computing zonation disruption: {cond_healthy[0]} vs {cond_masld[0]}...")
            all_markers = zon_config["periportal_markers"] + zon_config["pericentral_markers"]
            disruption = compute_zonation_disruption(
                adata, cond_healthy[0], cond_masld[0], all_markers
            )
            save_csv(disruption, "zonation_disruption_scores.csv", subdir="zonation")
            if len(disruption) > 0:
                mean_d = disruption["disruption_score"].mean()
                print(f"  Mean disruption score: {mean_d:.3f}")

    print_header("04b: Complete")


if __name__ == "__main__":
    main()
