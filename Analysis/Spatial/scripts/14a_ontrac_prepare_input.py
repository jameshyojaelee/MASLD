#!/usr/bin/env python3
"""
14a_ontrac_prepare_input.py — Prepare ONTraC input from deconvolved spatial data.

Extracts per-spot metadata (barcode, sample, dominant cell type, coordinates)
and cell-type composition matrix from cell2location results for ONTraC v2 mode.

SLURM: --partition=cpu --cpus=4 --mem=32G --time=1:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, C2L_PREFIX, load_config, load_deconvolved_adata,
    strip_c2l_prefix, save_csv, print_header, print_step,
)


def extract_cell_type_composition(adata):
    """Extract cell2location proportions per spot, strip prefix from names.

    Returns DataFrame: spots x cell types (proportions summing to 1 per spot).
    """
    # Find cell2location columns in obsm
    c2l_key = None
    for key in adata.obsm.keys():
        if "q05" in key and "cell_abundance" in key:
            c2l_key = key
            break

    if c2l_key is not None:
        obsm_data = adata.obsm[c2l_key]
        # Use .values to extract numpy array (avoids column-name mismatch when
        # constructing DataFrame with renamed columns from an existing DataFrame)
        if hasattr(obsm_data, "values"):
            raw_values = obsm_data.values
            col_names = [strip_c2l_prefix(c) for c in obsm_data.columns]
        else:
            raw_values = obsm_data
            col_names = [f"CellType_{i}" for i in range(obsm_data.shape[1])]
        comp = pd.DataFrame(
            raw_values,
            index=adata.obs_names,
            columns=col_names,
        )
    else:
        # Fallback: look for columns in obs starting with C2L_PREFIX
        c2l_cols = [c for c in adata.obs.columns if c.startswith(C2L_PREFIX)]
        if not c2l_cols:
            print("  ERROR: No cell2location proportions found in obsm or obs")
            sys.exit(1)
        comp = adata.obs[c2l_cols].copy()
        comp.columns = [strip_c2l_prefix(c) for c in comp.columns]

    # Normalize rows to sum to 1 (proportions)
    row_sums = comp.sum(axis=1)
    row_sums = row_sums.replace(0, 1)  # avoid division by zero
    comp = comp.div(row_sums, axis=0)

    return comp


def determine_dominant_cell_type(comp_df):
    """Assign dominant cell type per spot (highest proportion)."""
    return comp_df.idxmax(axis=1)


def extract_spatial_coords(adata):
    """Extract spatial coordinates from adata.obsm['spatial']."""
    if "spatial" not in adata.obsm:
        print("  ERROR: No spatial coordinates in adata.obsm['spatial']")
        sys.exit(1)

    coords = pd.DataFrame(
        adata.obsm["spatial"][:, :2],
        index=adata.obs_names,
        columns=["x", "y"],
    )
    return coords


def main():
    print_header("14a: Prepare ONTraC Input")

    config = load_config()
    ontrac_config = config["ontrac"]
    output_dir = RESULTS_DIR / "ontrac" / "input"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load deconvolved spatial data
    print_step("Loading deconvolved spatial data")
    adata = load_deconvolved_adata()
    print(f"  Loaded: {adata.n_obs} spots, {adata.n_vars} genes")

    # Extract cell-type composition
    print_step("Extracting cell-type composition")
    comp = extract_cell_type_composition(adata)
    print(f"  Composition matrix: {comp.shape[0]} spots x {comp.shape[1]} cell types")
    print(f"  Cell types: {list(comp.columns)}")

    # Determine dominant cell type
    print_step("Determining dominant cell type per spot")
    dominant_ct = determine_dominant_cell_type(comp)
    ct_counts = dominant_ct.value_counts()
    for ct, n in ct_counts.items():
        print(f"    {ct}: {n} spots ({n / len(dominant_ct) * 100:.1f}%)")

    # Extract spatial coordinates
    print_step("Extracting spatial coordinates")
    coords = extract_spatial_coords(adata)

    # Extract sample IDs
    sample_col = "sample_id" if "sample_id" in adata.obs.columns else "dataset"
    sample_ids = adata.obs[sample_col]

    # Build ONTraC metadata CSV
    print_step("Building ONTraC metadata CSV")
    metadata = pd.DataFrame({
        "Cell_ID": adata.obs_names,
        "Sample": sample_ids.values,
        "Cell_Type": dominant_ct.values,
        "x": coords["x"].values,
        "y": coords["y"].values,
    })
    metadata.index = adata.obs_names

    # Add condition for downstream analysis
    if "condition" in adata.obs.columns:
        metadata["condition"] = adata.obs["condition"].values

    save_csv(metadata, "ontrac_metadata.csv", subdir="ontrac/input")

    # Save composition matrix for ONTraC v2 mode
    print_step("Saving cell-type composition matrix")
    save_csv(comp, "ontrac_cell_type_composition.csv", subdir="ontrac/input")

    # Summary
    print(f"\n  Summary:")
    print(f"    Total spots: {len(metadata)}")
    print(f"    Samples: {metadata['Sample'].nunique()}")
    print(f"    Cell types: {metadata['Cell_Type'].nunique()}")
    print(f"    Coordinate range: x=[{coords['x'].min():.0f}, {coords['x'].max():.0f}], "
          f"y=[{coords['y'].min():.0f}, {coords['y'].max():.0f}]")

    print_header("14a: Complete")


if __name__ == "__main__":
    main()
