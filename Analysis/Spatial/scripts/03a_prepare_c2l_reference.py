#!/usr/bin/env python3
"""
03a_prepare_c2l_reference.py — Train cell2location reference model from scRNA-seq atlas.

Uses the 524K-cell CellTypist-annotated atlas as reference to estimate
per-gene, per-cell-type expression signatures (NB regression model).

SLURM: --partition=gpu --gres=gpu:1 --cpus=8 --mem=128G --time=6:00:00
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
    PROJECT_ROOT, RESULTS_DIR, load_config, init_spatial_gpu,
    check_checkpoint, print_header, print_step,
)


def main():
    parser = argparse.ArgumentParser(description="Train cell2location reference model")
    parser.add_argument("--force", action="store_true", help="Retrain even if model exists")
    args = parser.parse_args()

    print_header("03a: Train cell2location Reference Model")

    config = load_config()
    c2l_config = config["cell2location"]["reference"]
    ref_model_dir = RESULTS_DIR / "cell2location" / "reference_model"
    ref_model_dir.mkdir(parents=True, exist_ok=True)

    # Check if model already trained
    if not args.force and (ref_model_dir / "model.pt").exists():
        print("  Reference model already exists. Use --force to retrain.")
        return

    # Initialize GPU
    init_spatial_gpu()
    print("  GPU initialized")

    # Import cell2location (requires GPU)
    import cell2location

    # Load scRNA-seq reference atlas
    atlas_path = PROJECT_ROOT / config["paths"]["scrna_atlas"]
    print(f"  Loading reference atlas: {atlas_path}")
    adata_ref = sc.read_h5ad(atlas_path)
    print(f"  Reference: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes")

    # Identify cell type column
    ct_col = None
    for candidate in ["cell_type", "majority_voting", "predicted_labels", "leiden"]:
        if candidate in adata_ref.obs.columns:
            ct_col = candidate
            break
    if ct_col is None:
        print("  ERROR: No cell type column found in reference atlas")
        sys.exit(1)
    print(f"  Cell type column: '{ct_col}' ({adata_ref.obs[ct_col].nunique()} types)")
    print(f"  Cell types: {sorted(adata_ref.obs[ct_col].unique())}")

    # Identify batch column
    batch_col = None
    for candidate in ["dataset", "batch", "sample", "sample_id"]:
        if candidate in adata_ref.obs.columns:
            batch_col = candidate
            break
    print(f"  Batch column: '{batch_col}' ({adata_ref.obs[batch_col].nunique() if batch_col else 'None'} batches)")

    # Subsample if too large (OOM protection for L40S 48GB)
    max_cells = c2l_config["max_cells"]
    if adata_ref.n_obs > max_cells:
        print(f"  Subsampling {adata_ref.n_obs} → {max_cells} cells (stratified by {ct_col})")
        # Stratified sampling: equal allocation per cell type to preserve proportions
        ct_series = adata_ref.obs[ct_col]
        n_types = ct_series.nunique()
        per_type = max_cells // n_types
        indices = []
        for ct_val, ct_idx in adata_ref.obs.groupby(ct_col).groups.items():
            n_sample = min(len(ct_idx), per_type)
            rng = np.random.RandomState(42)
            indices.extend(rng.choice(ct_idx, size=n_sample, replace=False).tolist())
        adata_ref = adata_ref[indices].copy()
        print(f"  After stratified subsample: {adata_ref.n_obs} cells ({n_types} types)")

    # Ensure raw counts in .X (cell2location requires integer counts)
    if "counts" in adata_ref.layers:
        adata_ref.X = adata_ref.layers["counts"].copy()
        print("  Using 'counts' layer as raw counts")
    elif "raw_counts" in adata_ref.layers:
        adata_ref.X = adata_ref.layers["raw_counts"].copy()
        print("  Using 'raw_counts' layer as raw counts")
    elif adata_ref.raw is not None:
        # Standard scanpy convention: .raw stores pre-normalization data
        print("  Extracting raw counts from .raw")
        adata_raw = adata_ref.raw.to_adata()
        # Use raw X (integer counts) but keep current obs/var metadata
        adata_ref = ad.AnnData(
            X=adata_raw.X,
            obs=adata_ref.obs,
            var=adata_raw.var,
        )
        print(f"  Raw counts: {adata_ref.n_obs} cells × {adata_ref.n_vars} genes")
    else:
        # Check if .X already contains integers
        import scipy.sparse as sp
        x_sample = adata_ref.X[:100, :100]
        if sp.issparse(x_sample):
            x_sample = x_sample.toarray()
        non_zero = x_sample[x_sample != 0]
        if len(non_zero) > 0 and np.allclose(non_zero, non_zero.astype(int)):
            print("  .X appears to contain integer counts, using as-is")
        else:
            print("  WARNING: No raw counts found. .X contains non-integer values.")
            print("  Attempting to reverse log1p normalization...")
            import scipy.sparse as sp
            if sp.issparse(adata_ref.X):
                adata_ref.X = np.expm1(adata_ref.X.toarray())
            else:
                adata_ref.X = np.expm1(adata_ref.X)
            adata_ref.X = np.round(adata_ref.X).astype(int)

    # Filter genes: minimum 10 cells expressing
    sc.pp.filter_genes(adata_ref, min_cells=10)
    print(f"  After gene filter: {adata_ref.n_vars} genes")

    # Load spatial data gene list (if available) for intersection
    spatial_checkpoint = RESULTS_DIR / "preprocessed" / "merged_spatial.h5ad"
    if spatial_checkpoint.exists():
        print(f"  Loading spatial gene list for intersection...")
        adata_sp = sc.read_h5ad(spatial_checkpoint, backed="r")
        shared_genes = adata_ref.var_names.intersection(adata_sp.var_names)
        adata_ref = adata_ref[:, shared_genes].copy()
        print(f"  After spatial intersection: {adata_ref.n_vars} genes")
        del adata_sp

    # Setup and train NB regression reference model
    print(f"\n  Setting up cell2location RegressionModel...")
    setup_kwargs = {"labels_key": ct_col}
    if batch_col:
        setup_kwargs["batch_key"] = batch_col
    cell2location.models.RegressionModel.setup_anndata(adata_ref, **setup_kwargs)

    mod_ref = cell2location.models.RegressionModel(adata_ref)
    print(f"  Model parameters: {sum(p.numel() for p in mod_ref.module.parameters()):,}")

    print(f"\n  Training reference model (max_epochs={c2l_config['max_epochs']})...")
    # Note: use_gpu deprecated in newer scvi-tools; GPU auto-detected via accelerator
    train_kwargs = dict(
        max_epochs=c2l_config["max_epochs"],
        batch_size=c2l_config["batch_size"],
        train_size=c2l_config["train_size"],
        lr=c2l_config["lr"],
    )
    # Try with accelerator="gpu" (modern API), fall back to use_gpu=True (legacy)
    try:
        mod_ref.train(accelerator="gpu", **train_kwargs)
    except TypeError:
        mod_ref.train(use_gpu=True, **train_kwargs)

    # Check convergence
    print("\n  Checking convergence...")
    history = mod_ref.history["elbo_train"]
    last_10 = history.iloc[-10:].values.flatten()
    first_10 = history.iloc[:10].values.flatten()
    improvement = (first_10.mean() - last_10.mean()) / abs(first_10.mean()) * 100
    print(f"  ELBO improvement: {improvement:.1f}% (first 10 vs last 10 epochs)")
    if abs(last_10[-1] - last_10[0]) / abs(last_10[0]) > 0.05:
        print("  WARNING: Model may not have converged (>5% change in last 10 epochs)")

    # Save model
    print(f"\n  Saving model to {ref_model_dir}...")
    mod_ref.save(str(ref_model_dir), overwrite=True)

    # Export reference signatures
    print("  Exporting posterior signatures...")
    sample_kwargs = {
        "num_samples": 1000,
        "batch_size": c2l_config["batch_size"],
    }
    # Try modern API first (no use_gpu), fall back to legacy
    try:
        adata_ref = mod_ref.export_posterior(adata_ref, sample_kwargs=sample_kwargs)
    except TypeError:
        sample_kwargs["use_gpu"] = True
        adata_ref = mod_ref.export_posterior(adata_ref, sample_kwargs=sample_kwargs)

    # Save inf_aver (per-gene, per-cell-type mean expression)
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

    # Save reference AnnData with posterior
    adata_ref.write_h5ad(ref_model_dir / "reference_posterior.h5ad")

    print(f"\n  === Reference Model Summary ===")
    print(f"  Cells: {adata_ref.n_obs}")
    print(f"  Genes: {adata_ref.n_vars}")
    print(f"  Cell types: {adata_ref.obs[ct_col].nunique()}")
    print(f"  Model dir: {ref_model_dir}")

    print_header("03a: Complete")


if __name__ == "__main__":
    main()
