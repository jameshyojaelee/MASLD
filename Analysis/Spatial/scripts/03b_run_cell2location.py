#!/usr/bin/env python3
"""
03b_run_cell2location.py — Run cell2location spatial deconvolution.

Maps scRNA-seq cell type signatures onto Visium spatial data.
Produces per-spot cell type abundance estimates.

SLURM: --partition=gpu --gres=gpu:1 --cpus=8 --mem=128G --time=12:00:00
"""

import pathlib
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, init_spatial_gpu, load_spatial_adata,
    save_checkpoint, check_checkpoint, print_header, print_step,
)


def main():
    parser = argparse.ArgumentParser(description="Run cell2location spatial deconvolution")
    parser.add_argument("--force", action="store_true", help="Rerun even if output exists")
    args = parser.parse_args()

    print_header("03b: cell2location Spatial Deconvolution")

    config = load_config()
    c2l_config = config["cell2location"]["spatial"]
    ref_model_dir = RESULTS_DIR / "cell2location" / "reference_model"
    spatial_model_dir = RESULTS_DIR / "cell2location" / "spatial_model"
    spatial_model_dir.mkdir(parents=True, exist_ok=True)

    # Check if already done
    output_h5ad = spatial_model_dir / "spatial_deconvolved.h5ad"
    if not args.force and output_h5ad.exists():
        print("  Deconvolved data exists. Use --force to rerun.")
        return

    # Check reference model exists
    inf_aver_path = ref_model_dir / "inf_aver.csv"
    if not inf_aver_path.exists():
        print(f"  ERROR: Reference signatures not found at {inf_aver_path}")
        print("  Run 03a_prepare_c2l_reference.py first.")
        sys.exit(1)

    # Initialize GPU
    init_spatial_gpu()
    print("  GPU initialized")

    import cell2location

    # Load reference signatures
    print(f"  Loading reference signatures: {inf_aver_path}")
    inf_aver = pd.read_csv(inf_aver_path, index_col=0)
    print(f"  Reference: {inf_aver.shape[0]} genes × {inf_aver.shape[1]} cell types")
    print(f"  Cell types: {list(inf_aver.columns)}")

    # Load spatial data (raw counts)
    print("\n  Loading spatial data...")
    adata_vis = load_spatial_adata("merged_spatial_raw.h5ad")
    # Ensure raw counts
    if "counts" in adata_vis.layers:
        adata_vis.X = adata_vis.layers["counts"].copy()
    print(f"  Spatial data: {adata_vis.n_obs} spots × {adata_vis.n_vars} genes")

    # Intersect genes
    shared_genes = adata_vis.var_names.intersection(inf_aver.index)
    adata_vis = adata_vis[:, shared_genes].copy()
    inf_aver = inf_aver.loc[shared_genes, :]
    print(f"  Shared genes: {len(shared_genes)}")

    # Filter genes with very low expression (cell2location recommendation)
    sc.pp.filter_genes(adata_vis, min_cells=5)
    inf_aver = inf_aver.loc[inf_aver.index.isin(adata_vis.var_names), :]
    print(f"  After filter: {adata_vis.n_vars} genes")

    # Setup spatial model
    print("\n  Setting up cell2location spatial model...")
    cell2location.models.Cell2location.setup_anndata(
        adata_vis,
        batch_key="sample_id",
    )

    mod_spatial = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=c2l_config["N_cells_per_location"],
        detection_alpha=c2l_config["detection_alpha"],
    )
    print(f"  Model parameters: {sum(p.numel() for p in mod_spatial.module.parameters()):,}")

    # Train
    print(f"\n  Training spatial model (max_epochs={c2l_config['max_epochs']})...")
    train_kwargs = dict(
        max_epochs=c2l_config["max_epochs"],
        batch_size=None,  # Auto (full dataset)
        train_size=1.0,
    )
    try:
        mod_spatial.train(accelerator="gpu", **train_kwargs)
    except TypeError:
        mod_spatial.train(use_gpu=True, **train_kwargs)

    # Check convergence
    history = mod_spatial.history["elbo_train"]
    last_10 = history.iloc[-10:].values.flatten()
    print(f"\n  Final ELBO (last 10): mean={last_10.mean():.0f}, std={last_10.std():.0f}")

    # Export posterior
    print("\n  Exporting posterior estimates...")
    sample_kwargs = {
        "num_samples": c2l_config["num_posterior_samples"],
        "batch_size": 2500,
    }
    try:
        adata_vis = mod_spatial.export_posterior(adata_vis, sample_kwargs=sample_kwargs)
    except TypeError:
        sample_kwargs["use_gpu"] = True
        adata_vis = mod_spatial.export_posterior(adata_vis, sample_kwargs=sample_kwargs)

    # Assign dominant cell type per spot
    abundances = adata_vis.obsm["q05_cell_abundance_w_sf"]
    if isinstance(abundances, pd.DataFrame):
        adata_vis.obs["cell_type_dominant"] = abundances.columns[
            abundances.values.argmax(axis=1)
        ]
        total = abundances.values.sum(axis=1)
        adata_vis.obs["cell_type_confidence"] = (
            abundances.values.max(axis=1) / np.maximum(total, 1e-10)
        )
        # Also store as separate obs columns for easy plotting
        for ct in abundances.columns:
            adata_vis.obs[f"c2l_{ct}"] = abundances[ct].values
    else:
        # NumPy array format
        ct_names = inf_aver.columns.tolist()
        adata_vis.obs["cell_type_dominant"] = [
            ct_names[i] for i in abundances.argmax(axis=1)
        ]
        total = abundances.sum(axis=1)
        adata_vis.obs["cell_type_confidence"] = (
            abundances.max(axis=1) / np.maximum(total, 1e-10)
        )

    # Save model and results
    print(f"\n  Saving model to {spatial_model_dir}...")
    mod_spatial.save(str(spatial_model_dir), overwrite=True)

    print("  Saving deconvolved AnnData...")
    adata_vis.write_h5ad(output_h5ad)

    # Summary
    print(f"\n  === Deconvolution Summary ===")
    print(f"  Spots: {adata_vis.n_obs}")
    print(f"  Genes: {adata_vis.n_vars}")
    print(f"  Cell types: {len(inf_aver.columns)}")
    print(f"  Dominant cell type distribution:")
    for ct, count in adata_vis.obs["cell_type_dominant"].value_counts().items():
        pct = count / adata_vis.n_obs * 100
        print(f"    {ct}: {count} spots ({pct:.1f}%)")

    print_header("03b: Complete")


if __name__ == "__main__":
    main()
