#!/usr/bin/env python3
"""
14d_ontrac_zonation_overlay.py — Overlay niche trajectory on zonation and domains.

Correlates NT scores with marker-gene zonation from 04a, cross-tabulates
niche clusters with 05d spatial domains (ARI, NMI), and attempts Vu et al.
cross-dataset TAG replication.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_zonation_scores,
    save_csv, print_header, print_step,
)


def load_nt_scores():
    """Load niche trajectory scores from 14b."""
    path = RESULTS_DIR / "ontrac" / "niche_trajectory_scores.csv"
    if not path.exists():
        print("  ERROR: Run 14b_ontrac_run.py first")
        sys.exit(1)
    return pd.read_csv(path, index_col=0)


def load_tags():
    """Load trajectory-associated genes from 14c."""
    path = RESULTS_DIR / "ontrac" / "trajectory_associated_genes.csv"
    if not path.exists():
        print("  WARNING: No TAG results from 14c")
        return pd.DataFrame()
    return pd.read_csv(path, index_col=0)


def correlate_nt_with_zonation(nt_df, zon_df):
    """Correlate NT score with zonation scores per spot (Spearman)."""
    shared = nt_df.index.intersection(zon_df.index)
    if len(shared) < 50:
        print(f"  WARNING: Only {len(shared)} shared spots for NT-zonation correlation")
        return pd.DataFrame()

    nt_vals = nt_df.loc[shared, "nt_score"].values
    results = []

    # Correlate with each zonation column
    zon_cols = [c for c in zon_df.columns if "score" in c.lower() or "zone" in c.lower()
                or c in ["zonation_score", "periportal_score", "pericentral_score"]]
    if not zon_cols:
        # Use all numeric columns
        zon_cols = zon_df.select_dtypes(include=[np.number]).columns.tolist()

    for col in zon_cols:
        vals = zon_df.loc[shared, col].values
        valid = ~(np.isnan(nt_vals) | np.isnan(vals))
        if valid.sum() < 50:
            continue

        rho, pval = spearmanr(nt_vals[valid], vals[valid])
        results.append({
            "zonation_metric": col,
            "n_spots": valid.sum(),
            "spearman_rho": rho,
            "pval": pval,
        })
        sig = "*" if pval < 0.05 else "ns"
        print(f"    NT vs {col}: rho={rho:.3f}, p={pval:.2e} {sig}")

    return pd.DataFrame(results)


def compare_niche_clusters_with_domains(nt_df):
    """Cross-tabulate niche clusters with spatial domains (ARI, NMI)."""
    # Load spatial domains from 05d
    dom_path = RESULTS_DIR / "domains" / "spatial_domains.csv"
    if not dom_path.exists():
        print("  WARNING: No spatial domain data (run 05d first)")
        return pd.DataFrame()

    dom_df = pd.read_csv(dom_path, index_col=0)
    if "spatial_domain" not in dom_df.columns:
        print("  WARNING: No spatial_domain column in domain results")
        return pd.DataFrame()

    shared = nt_df.index.intersection(dom_df.index)
    if len(shared) < 50:
        print(f"  WARNING: Only {len(shared)} shared spots for domain comparison")
        return pd.DataFrame()

    niche_labels = nt_df.loc[shared, "niche_cluster"].values
    domain_labels = dom_df.loc[shared, "spatial_domain"].values

    # Compute agreement metrics
    ari = adjusted_rand_score(domain_labels, niche_labels)
    nmi = normalized_mutual_info_score(domain_labels, niche_labels)

    print(f"    ARI (niche clusters vs spatial domains): {ari:.3f}")
    print(f"    NMI (niche clusters vs spatial domains): {nmi:.3f}")

    # Cross-tabulation
    ct = pd.crosstab(
        pd.Series(niche_labels, name="niche_cluster"),
        pd.Series(domain_labels, name="spatial_domain"),
    )
    ct_norm = ct.div(ct.sum(axis=1), axis=0)

    result = pd.DataFrame([{
        "n_spots": len(shared),
        "n_niche_clusters": len(np.unique(niche_labels)),
        "n_spatial_domains": len(np.unique(domain_labels)),
        "adjusted_rand_index": ari,
        "normalized_mutual_info": nmi,
    }])

    # Append crosstab as separate output
    save_csv(ct_norm, "niche_domain_crosstab.csv", subdir="ontrac")

    return result


def vu_cross_dataset_replication(tag_df):
    """Attempt TAG replication with Vu et al. dataset.

    Checks for Vu ONTraC results or Vu expression data for correlation-based
    replication. Reports Jaccard overlap if available.
    """
    results = []

    # Check for Vu ONTraC TAG results
    vu_tag_path = RESULTS_DIR / "ontrac" / "vu_trajectory_associated_genes.csv"
    if vu_tag_path.exists():
        vu_tags = pd.read_csv(vu_tag_path, index_col=0)
        gse_tags = set(tag_df[tag_df["is_tag"]]["gene"]) if "is_tag" in tag_df.columns else set()
        vu_tag_set = set(vu_tags[vu_tags.get("is_tag", True)].index
                         if "is_tag" in vu_tags.columns
                         else vu_tags.index)

        overlap = gse_tags & vu_tag_set
        union = gse_tags | vu_tag_set
        jaccard = len(overlap) / max(len(union), 1)

        results.append({
            "comparison": "GSE192741_vs_Vu_TAGs",
            "n_gse_tags": len(gse_tags),
            "n_vu_tags": len(vu_tag_set),
            "n_overlap": len(overlap),
            "jaccard": jaccard,
        })
        print(f"    GSE192741 vs Vu TAGs: overlap={len(overlap)}, Jaccard={jaccard:.3f}")
        if overlap:
            print(f"    Replicated TAGs: {sorted(overlap)[:20]}")
        return pd.DataFrame(results)

    # Fallback: try expression-based TAG replication on Vu data
    vu_h5ad = RESULTS_DIR / "preprocessed" / "merged_spatial_vu.h5ad"
    if not vu_h5ad.exists():
        vu_h5ad = RESULTS_DIR / "validation_vu" / "merged_spatial_vu.h5ad"

    if not vu_h5ad.exists():
        print("  WARNING: No Vu et al. data found for cross-dataset replication")
        results.append({
            "comparison": "GSE192741_vs_Vu_TAGs",
            "n_gse_tags": len(tag_df[tag_df.get("is_tag", False) == True])
            if "is_tag" in tag_df.columns else 0,
            "n_vu_tags": 0,
            "n_overlap": 0,
            "jaccard": np.nan,
            "note": "Vu data not available",
        })
        return pd.DataFrame(results)

    # Load Vu data and compute gene-level spatial autocorrelation as proxy
    print("  Loading Vu et al. data for expression-based replication...")
    try:
        import squidpy as sq

        vu_adata = sc.read_h5ad(vu_h5ad)
        print(f"    Vu: {vu_adata.n_obs} spots, {vu_adata.n_vars} genes")

        # Get GSE192741 TAGs
        gse_tags = tag_df[tag_df.get("is_tag", False) == True] if "is_tag" in tag_df.columns else tag_df
        gse_tag_genes = set(gse_tags["gene"]) if "gene" in gse_tags.columns else set(gse_tags.index)

        # Genes present in both
        shared_genes = gse_tag_genes & set(vu_adata.var_names)
        print(f"    GSE192741 TAGs present in Vu: {len(shared_genes)}/{len(gse_tag_genes)}")

        if len(shared_genes) > 0:
            # Compute Moran's I for TAGs in Vu as spatial validation proxy
            sq.gr.spatial_neighbors(vu_adata, coord_type="generic", n_neighs=6)
            sq.gr.spatial_autocorr(
                vu_adata, mode="moran", genes=list(shared_genes)[:500],
                n_perms=100, n_jobs=8,
            )
            vu_morans = vu_adata.uns["moranI"]
            n_spatially_sig = (vu_morans["pval_norm"] < 0.05).sum()

            results.append({
                "comparison": "GSE192741_TAGs_spatial_in_Vu",
                "n_gse_tags": len(gse_tag_genes),
                "n_vu_tags": len(shared_genes),
                "n_overlap": n_spatially_sig,
                "jaccard": n_spatially_sig / max(len(shared_genes), 1),
                "note": "Spatial autocorrelation proxy (Moran I p<0.05 in Vu)",
            })
            print(f"    TAGs spatially significant in Vu: {n_spatially_sig}/{len(shared_genes)}")
    except Exception as e:
        print(f"  WARNING: Vu replication failed: {e}")
        results.append({
            "comparison": "GSE192741_vs_Vu_TAGs",
            "n_gse_tags": len(tag_df[tag_df.get("is_tag", False) == True])
            if "is_tag" in tag_df.columns else 0,
            "n_vu_tags": 0,
            "n_overlap": 0,
            "jaccard": np.nan,
            "note": f"Replication failed: {str(e)[:100]}",
        })

    return pd.DataFrame(results)


def main():
    print_header("14d: Zonation and Domain Overlay")

    config = load_config()
    output_dir = RESULTS_DIR / "ontrac"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load NT scores
    print_step("Loading niche trajectory scores")
    nt_df = load_nt_scores()
    print(f"  NT scores: {len(nt_df)} spots")

    # Load zonation scores from 04a
    print_step("Correlating NT with zonation")
    zon_df = load_zonation_scores()
    if len(zon_df) > 0:
        zon_corr = correlate_nt_with_zonation(nt_df, zon_df)
        if len(zon_corr) > 0:
            save_csv(zon_corr, "zonation_correlation.csv", subdir="ontrac")
    else:
        print("  WARNING: No zonation scores available")

    # Cross-tabulate niche clusters with spatial domains
    print_step("Comparing niche clusters with spatial domains")
    domain_comp = compare_niche_clusters_with_domains(nt_df)
    if len(domain_comp) > 0:
        save_csv(domain_comp, "domain_niche_comparison.csv", subdir="ontrac")

    # Vu cross-dataset TAG replication
    print_step("Cross-dataset TAG replication (Vu et al.)")
    tag_df = load_tags()
    if len(tag_df) > 0:
        replication = vu_cross_dataset_replication(tag_df)
        if len(replication) > 0:
            save_csv(replication, "cross_dataset_tag_replication.csv", subdir="ontrac")
    else:
        print("  WARNING: No TAGs from 14c to replicate")

    print_header("14d: Complete")


if __name__ == "__main__":
    main()
