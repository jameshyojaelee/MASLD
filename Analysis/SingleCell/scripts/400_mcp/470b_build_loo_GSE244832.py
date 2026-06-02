#!/usr/bin/env python
"""Build the critical LOO h5ad (GSE244832 = 62.5% of cells)."""
from __future__ import annotations
import os
from pathlib import Path
import anndata as ad

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
inp = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/atlas_cnmf_global.h5ad"
out_dir = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs"

atlas = ad.read_h5ad(inp)
for d in ["GSE244832", "GSE202379", "GSE185477"]:
    sub = atlas[atlas.obs["dataset"] != d].copy()
    out = out_dir / f"atlas_cnmf_global_loo_{d}.h5ad"
    if not out.exists():
        print(f"writing {out.name} ({sub.shape[0]:,} cells)")
        sub.write_h5ad(out, compression="gzip")
    else:
        print(f"exists: {out.name}")
