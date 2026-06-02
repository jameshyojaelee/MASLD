#!/usr/bin/env python3
"""
04a_define_zonation.py — Define liver zonation axis from spatial data.

Computes continuous periportal-pericentral zonation score per spot using
established marker genes. Validates via spatial autocorrelation (Moran's I).

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_deconvolved_adata,
    save_checkpoint, save_csv, check_checkpoint,
    print_header, print_step,
)


def compute_zonation_score(adata, periportal, pericentral, ctrl_multiplier=5):
    """Compute continuous periportal-pericentral zonation score per spot."""
    pp_present = [g for g in periportal if g in adata.var_names]
    pc_present = [g for g in pericentral if g in adata.var_names]
    print(f"  Periportal markers present: {len(pp_present)}/{len(periportal)}: {pp_present}")
    print(f"  Pericentral markers present: {len(pc_present)}/{len(pericentral)}: {pc_present}")

    if len(pp_present) < 2 or len(pc_present) < 2:
        print("  WARNING: Too few zonation markers detected")

    sc.tl.score_genes(adata, gene_list=pp_present, score_name="periportal_score",
                      ctrl_size=len(pp_present) * ctrl_multiplier)
    sc.tl.score_genes(adata, gene_list=pc_present, score_name="pericentral_score",
                      ctrl_size=len(pc_present) * ctrl_multiplier)

    adata.obs["zonation_score"] = (
        adata.obs["pericentral_score"] - adata.obs["periportal_score"]
    )

    # Quantile bins
    adata.obs["zonation_bin"] = pd.qcut(
        adata.obs["zonation_score"], q=5,
        labels=["PP1", "PP2", "Mid", "PC2", "PC1"],
        duplicates="drop",
    )
    return adata


def validate_zonation(adata, condition_label="all"):
    """Validate zonation gradient via spatial autocorrelation."""
    # Ensure spatial neighbors built
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)

    # Per-sample Moran's I on zonation score
    results = []
    for sample in adata.obs["sample_id"].unique():
        mask = adata.obs["sample_id"] == sample
        adata_sub = adata[mask].copy()
        if adata_sub.n_obs < 50:
            continue
        try:
            # Build a minimal AnnData with zonation_score as a single "gene"
            # and recompute spatial neighbors on it (avoids obsp shape mismatch)
            adata_tmp = sc.AnnData(
                X=np.array(adata_sub.obs[["zonation_score"]].values, dtype=np.float32),
                obs=adata_sub.obs.copy(),
            )
            adata_tmp.obsm["spatial"] = adata_sub.obsm["spatial"].copy()
            sq.gr.spatial_neighbors(adata_tmp, coord_type="generic", n_neighs=6)
            sq.gr.spatial_autocorr(
                adata_tmp, mode="moran", n_perms=199,
            )
            mi = adata_tmp.uns["moranI"].iloc[0]
            results.append({
                "sample_id": sample,
                "condition": adata_sub.obs["condition"].iloc[0],
                "morans_i": mi["I"],
                "morans_pval": mi["pval_norm"],
                "n_spots": adata_sub.n_obs,
            })
        except Exception as e:
            print(f"    WARNING: Moran's I failed for {sample}: {e}")
            results.append({
                "sample_id": sample,
                "condition": adata_sub.obs["condition"].iloc[0],
                "morans_i": np.nan,
                "morans_pval": np.nan,
                "n_spots": adata_sub.n_obs,
            })

    return pd.DataFrame(results)


def main():
    print_header("04a: Define Liver Zonation Axis")

    config = load_config()
    zon_config = config["zonation"]
    output_dir = RESULTS_DIR / "zonation"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load deconvolved spatial data (has cell type assignments)
    adata = load_deconvolved_adata()
    print(f"  Loaded: {adata.n_obs} spots")

    # Need log-normalized expression for gene scoring
    # If only raw counts, normalize
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Compute zonation score
    print("\n  Computing zonation scores...")
    adata = compute_zonation_score(
        adata,
        periportal=zon_config["periportal_markers"],
        pericentral=zon_config["pericentral_markers"],
        ctrl_multiplier=zon_config["ctrl_gene_multiplier"],
    )

    print(f"\n  Zonation score range: [{adata.obs['zonation_score'].min():.3f}, "
          f"{adata.obs['zonation_score'].max():.3f}]")
    print(f"  Zonation bin distribution:")
    for b, n in adata.obs["zonation_bin"].value_counts().sort_index().items():
        print(f"    {b}: {n} spots ({n/adata.n_obs*100:.1f}%)")

    # Validate spatial autocorrelation
    print("\n  Validating zonation gradient (Moran's I)...")
    validation = validate_zonation(adata)
    save_csv(validation, "zonation_validation.csv", subdir="zonation")

    valid_samples = validation[validation["morans_i"] > zon_config["min_morans_i"]]
    print(f"\n  Validation: {len(valid_samples)}/{len(validation)} samples have "
          f"Moran's I > {zon_config['min_morans_i']}")
    if len(validation) > 0:
        print(f"  Mean Moran's I: {validation['morans_i'].mean():.3f}")
        for cond in validation["condition"].unique():
            cond_mi = validation[validation["condition"] == cond]["morans_i"].mean()
            print(f"    {cond}: {cond_mi:.3f}")

    # Save zonation scores
    zon_df = adata.obs[["sample_id", "dataset", "condition",
                         "zonation_score", "zonation_bin",
                         "periportal_score", "pericentral_score"]].copy()
    save_csv(zon_df, "zonation_scores.csv", subdir="zonation")

    # Save AnnData with zonation
    save_checkpoint(adata, "spatial_with_zonation.h5ad", subdir="zonation")

    print_header("04a: Complete")


if __name__ == "__main__":
    main()
