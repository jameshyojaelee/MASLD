#!/usr/bin/env python3
"""Regenerate only panels 1, 5, 10 with updated font sizes (PDF only)."""

import sys
sys.path.insert(0, "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures")

from spatial_presentation_panels import (
    panel_01_deconvolution, panel_05_svg_heatmap, panel_10_communication_waterfall,
    RESULTS_DIR, OUT_DIR
)
import scanpy as sc

print("Regenerating panels 1, 5, 10")
print(f"Output: {OUT_DIR}\n")

# Panel 1 needs zonation h5ad (has c2l columns)
zon_path = RESULTS_DIR / "zonation" / "spatial_with_zonation.h5ad"
print("Loading zonation h5ad...")
adata_z = sc.read_h5ad(zon_path)
print(f"  Loaded: {adata_z.shape}\n")

panel_01_deconvolution(adata_z)
panel_05_svg_heatmap()
panel_10_communication_waterfall()

print("\nDone!")
