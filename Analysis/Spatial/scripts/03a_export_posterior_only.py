#!/usr/bin/env python3
"""
03a_export_posterior_only.py — Resume from saved reference model, export posterior.

Use this after 03a_prepare_c2l_reference.py saved model.pt but crashed during
export_posterior. Skips the 250-epoch training.
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, init_spatial_gpu,
    print_header,
)


def main():
    print_header("03a Resume: Export Posterior from Saved Model")

    config = load_config()
    c2l_config = config["cell2location"]["reference"]
    ref_model_dir = RESULTS_DIR / "cell2location" / "reference_model"

    if not (ref_model_dir / "model.pt").exists():
        print("  ERROR: No saved model found. Run 03a_prepare_c2l_reference.py first.")
        sys.exit(1)

    init_spatial_gpu()
    import cell2location

    # Reload the reference atlas (same preprocessing as 03a)
    atlas_path = PROJECT_ROOT / config["paths"]["scrna_atlas"]
    print(f"  Loading reference atlas: {atlas_path}")
    adata_ref = sc.read_h5ad(atlas_path)
    print(f"  Reference: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes")

    # Identify cell type and batch columns
    ct_col = None
    for candidate in ["cell_type", "majority_voting", "predicted_labels", "leiden"]:
        if candidate in adata_ref.obs.columns:
            ct_col = candidate
            break
    batch_col = None
    for candidate in ["dataset", "batch", "sample", "sample_id"]:
        if candidate in adata_ref.obs.columns:
            batch_col = candidate
            break
    print(f"  Cell type: '{ct_col}', Batch: '{batch_col}'")

    # Subsample (same as training)
    max_cells = c2l_config["max_cells"]
    if adata_ref.n_obs > max_cells:
        sc.pp.subsample(adata_ref, n_obs=max_cells, random_state=42)
        print(f"  Subsampled to {adata_ref.n_obs} cells")

    # Get raw counts from .raw
    if adata_ref.raw is not None:
        adata_raw = adata_ref.raw.to_adata()
        adata_ref = ad.AnnData(
            X=adata_raw.X,
            obs=adata_ref.obs,
            var=adata_raw.var,
        )
        print(f"  Raw counts: {adata_ref.n_obs} × {adata_ref.n_vars}")
    elif "counts" in adata_ref.layers:
        adata_ref.X = adata_ref.layers["counts"].copy()

    # Filter genes (same as training)
    sc.pp.filter_genes(adata_ref, min_cells=10)
    spatial_checkpoint = RESULTS_DIR / "preprocessed" / "merged_spatial.h5ad"
    if spatial_checkpoint.exists():
        adata_sp = sc.read_h5ad(spatial_checkpoint, backed="r")
        shared_genes = adata_ref.var_names.intersection(adata_sp.var_names)
        adata_ref = adata_ref[:, shared_genes].copy()
        del adata_sp
    print(f"  After filtering: {adata_ref.n_vars} genes")

    # Setup and load saved model
    print("\n  Loading saved model...")
    setup_kwargs = {"labels_key": ct_col}
    if batch_col:
        setup_kwargs["batch_key"] = batch_col
    cell2location.models.RegressionModel.setup_anndata(adata_ref, **setup_kwargs)
    mod_ref = cell2location.models.RegressionModel.load(str(ref_model_dir), adata_ref)
    print("  Model loaded successfully")

    # Export posterior (without use_gpu)
    print("  Exporting posterior signatures...")
    sample_kwargs = {
        "num_samples": 1000,
        "batch_size": c2l_config["batch_size"],
    }
    try:
        adata_ref = mod_ref.export_posterior(adata_ref, sample_kwargs=sample_kwargs)
    except TypeError:
        sample_kwargs["use_gpu"] = True
        adata_ref = mod_ref.export_posterior(adata_ref, sample_kwargs=sample_kwargs)

    # Save inf_aver
    if "means_per_cluster_mu_fg" in adata_ref.varm:
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"]
        if isinstance(inf_aver, pd.DataFrame):
            inf_aver.to_csv(ref_model_dir / "inf_aver.csv")
        else:
            pd.DataFrame(
                inf_aver,
                index=adata_ref.var_names,
                columns=adata_ref.obs[ct_col].cat.categories
                if hasattr(adata_ref.obs[ct_col], "cat") else sorted(adata_ref.obs[ct_col].unique()),
            ).to_csv(ref_model_dir / "inf_aver.csv")
        print(f"  Saved inf_aver: {inf_aver.shape}")
    else:
        print("  WARNING: means_per_cluster_mu_fg not found in varm")

    adata_ref.write_h5ad(ref_model_dir / "reference_posterior.h5ad")

    print(f"\n  === Reference Model Summary ===")
    print(f"  Cells: {adata_ref.n_obs}")
    print(f"  Genes: {adata_ref.n_vars}")
    print(f"  Cell types: {adata_ref.obs[ct_col].nunique()}")
    print(f"  Model dir: {ref_model_dir}")

    print_header("03a Resume: Complete")


if __name__ == "__main__":
    main()
