#!/usr/bin/env python3
"""
15c_prepare_spatial_for_gsmap.py — Prepare spatial h5ad files for gsMap.

Splits GSE192741 and Vu et al. spatial data into per-sample h5ad files
with raw counts and spatial coordinates, as required by gsMap.

gsMap requires:
  - adata.X or adata.layers['counts'] with raw integer counts
  - adata.obsm['spatial'] with spatial coordinates
  - Gene names as human gene symbols (var_names)

SLURM: --partition=cpu --cpus=8 --mem=64G --time=2:00:00
"""

import pathlib
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp_sparse

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_dataset_config,
    print_header, print_step,
)

OUTPUT_DIR = PROJECT_ROOT / "Analysis/Spatial/data/gsmap_input"


def ensure_raw_counts(adata, name):
    """Ensure adata.X contains raw integer counts.

    gsMap needs raw counts (not normalized/log-transformed).
    Check layers['counts'], layers['raw_counts'], then X.
    """
    for layer_name in ["counts", "raw_counts", "count"]:
        if layer_name in adata.layers:
            print(f"    Using layer '{layer_name}' as raw counts")
            adata.X = adata.layers[layer_name].copy()
            return adata

    # Check if X looks like raw counts (integers, max > 10)
    if sp_sparse.issparse(adata.X):
        sample_data = adata.X[:100, :100].toarray()
    else:
        sample_data = adata.X[:100, :100]

    is_integer = np.allclose(sample_data, np.round(sample_data), atol=0.01, equal_nan=True)
    max_val = np.nanmax(sample_data)

    if is_integer and max_val > 5:
        print(f"    X appears to contain raw counts (integer, max={max_val:.0f})")
        return adata

    print(f"    WARNING: X may not be raw counts for {name} "
          f"(integer={is_integer}, max={max_val:.1f}). Using as-is.")
    return adata


def ensure_spatial_coords(adata, name):
    """Ensure adata.obsm['spatial'] exists."""
    if "spatial" in adata.obsm:
        return True
    # Try X_spatial
    if "X_spatial" in adata.obsm:
        adata.obsm["spatial"] = adata.obsm["X_spatial"]
        return True
    print(f"    ERROR: No spatial coordinates in {name}")
    return False


def split_and_save_dataset(adata, sample_col, dataset_prefix, dataset_label):
    """Split merged AnnData by sample and save per-sample h5ad files.

    Parameters
    ----------
    adata : AnnData
        Merged spatial data with raw counts.
    sample_col : str
        Column in adata.obs containing sample IDs.
    dataset_prefix : str
        Prefix for output filenames (e.g., 'gse192741', 'vu').
    dataset_label : str
        Label for logging.

    Returns
    -------
    list of (sample_name, output_path) tuples.
    """
    samples = adata.obs[sample_col].unique().tolist()
    results = []

    for i, sample_id in enumerate(sorted(samples), 1):
        print_step(f"{dataset_label}: {sample_id}", i, len(samples))

        mask = adata.obs[sample_col] == sample_id
        sub = adata[mask].copy()

        # Ensure raw counts
        sub = ensure_raw_counts(sub, sample_id)

        # Ensure spatial coords
        if not ensure_spatial_coords(sub, sample_id):
            continue

        # Minimal h5ad: X (counts), obs, var, obsm['spatial']
        out_adata = ad.AnnData(
            X=sub.X.copy(),
            obs=sub.obs[[]].copy(),  # minimal obs (just index)
            var=sub.var[[]].copy(),  # minimal var (just gene names)
            obsm={"spatial": sub.obsm["spatial"].copy()},
        )
        # Ensure gene names are var_names
        out_adata.var_names_make_unique()

        # Add condition annotation for Cauchy combination
        if "condition" in sub.obs.columns:
            out_adata.obs["condition"] = sub.obs["condition"].values.copy()

        # Store counts in layers['counts'] as gsMap expects
        out_adata.layers["counts"] = out_adata.X.copy()

        # Save
        safe_name = str(sample_id).replace("/", "_").replace(" ", "_")
        out_path = OUTPUT_DIR / f"{dataset_prefix}_{safe_name}.h5ad"
        out_adata.write_h5ad(out_path)
        size_mb = out_path.stat().st_size / 1e6
        print(f"      {out_adata.n_obs} spots x {out_adata.n_vars} genes "
              f"({size_mb:.1f} MB)")

        results.append((f"{dataset_prefix}_{safe_name}", str(out_path)))

    return results


def process_gse192741():
    """Load and split GSE192741 human samples."""
    # Use raw merged checkpoint
    raw_path = RESULTS_DIR / "preprocessed" / "merged_spatial_raw.h5ad"
    if not raw_path.exists():
        # Fall back to processed version
        raw_path = RESULTS_DIR / "preprocessed" / "merged_spatial.h5ad"

    if not raw_path.exists():
        print("  ERROR: GSE192741 spatial data not found")
        return []

    print(f"  Loading GSE192741: {raw_path}")
    adata = sc.read_h5ad(raw_path)
    print(f"    {adata.n_obs} spots x {adata.n_vars} genes")

    # Only human samples
    if "species" in adata.obs.columns:
        adata = adata[adata.obs["species"] == "human"].copy()
        print(f"    Human only: {adata.n_obs} spots")

    # Identify sample column
    sample_col = None
    for col in ["sample_id", "sample", "library_id"]:
        if col in adata.obs.columns:
            sample_col = col
            break
    if sample_col is None:
        print("  ERROR: Cannot find sample column in GSE192741")
        return []

    print(f"    Sample column: '{sample_col}', "
          f"samples: {adata.obs[sample_col].nunique()}")

    return split_and_save_dataset(adata, sample_col, "gse192741", "GSE192741")


def process_vu():
    """Load and split Vu et al. arrays."""
    raw_path = RESULTS_DIR / "preprocessed" / "merged_spatial_vu_raw.h5ad"
    if not raw_path.exists():
        raw_path = RESULTS_DIR / "preprocessed" / "merged_spatial_vu.h5ad"

    if not raw_path.exists():
        print("  ERROR: Vu et al. spatial data not found")
        return []

    print(f"\n  Loading Vu et al.: {raw_path}")
    adata = sc.read_h5ad(raw_path)
    print(f"    {adata.n_obs} spots x {adata.n_vars} genes")

    # Identify sample column
    sample_col = None
    for col in ["sample_id", "sample", "array_id", "library_id"]:
        if col in adata.obs.columns:
            sample_col = col
            break
    if sample_col is None:
        print("  ERROR: Cannot find sample column in Vu et al.")
        return []

    print(f"    Sample column: '{sample_col}', "
          f"arrays: {adata.obs[sample_col].nunique()}")

    return split_and_save_dataset(adata, sample_col, "vu", "Vu et al.")


def write_h5ad_yaml(all_samples, output_dir):
    """Write YAML mapping sample_name -> h5ad path for gsMap."""
    import yaml

    # GSE192741 YAML
    gse_samples = {name: path for name, path in all_samples
                   if name.startswith("gse192741")}
    if gse_samples:
        yaml_path = output_dir / "gse192741_h5ad.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(gse_samples, f, default_flow_style=False, sort_keys=False)
        print(f"  GSE192741 YAML: {yaml_path} ({len(gse_samples)} samples)")

    # Vu YAML
    vu_samples = {name: path for name, path in all_samples
                  if name.startswith("vu")}
    if vu_samples:
        yaml_path = output_dir / "vu_h5ad.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(vu_samples, f, default_flow_style=False, sort_keys=False)
        print(f"  Vu YAML: {yaml_path} ({len(vu_samples)} samples)")


def main():
    print_header("15c: Prepare Spatial Data for gsMap")

    config = load_config()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_samples = []

    # Process GSE192741
    gse_results = process_gse192741()
    all_samples.extend(gse_results)

    # Process Vu et al.
    vu_results = process_vu()
    all_samples.extend(vu_results)

    # Write YAML configs
    if all_samples:
        write_h5ad_yaml(all_samples, OUTPUT_DIR)

    print(f"\n  Total: {len(all_samples)} h5ad files written to {OUTPUT_DIR}")

    print_header("15c: Complete")


if __name__ == "__main__":
    main()
