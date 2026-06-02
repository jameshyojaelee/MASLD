#!/usr/bin/env python
"""Export scATAC UMAP triad data for Fig 1 (cell_type / F_stage / condition).

Reads label-transferred snapATAC2 h5ad (anndata only — no scanpy), joins per-donor
metadata (F_stage_augmented, condition), and writes a flat TSV for the R figure
script `scripts/figures/fig1_scatac_umap_triad.R`.

Output columns: x, y, cell_type, donor, F_stage_augmented, condition
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
H5AD_PATH = (
    PROJECT_ROOT
    / "Analysis/ATAC/Human_Multiome/results/label_transfer"
    / "snapatac2_label_transferred.h5ad"
)
DONOR_META_PATH = (
    PROJECT_ROOT
    / "Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv"
)
OUT_PATH = (
    PROJECT_ROOT
    / "Analysis/ATAC/Human_Multiome/results/snapatac2/umap_export.tsv.gz"
)


def main() -> int:
    print(f"[load] {H5AD_PATH}", flush=True)
    # Backed='r' avoids loading the giant peak matrix into memory; we only
    # touch obs + obsm.
    adata = ad.read_h5ad(H5AD_PATH, backed="r")
    print(f"[load] shape={adata.shape}", flush=True)

    if "X_umap" not in adata.obsm:
        print("[ERROR] X_umap missing from obsm", file=sys.stderr)
        return 1

    umap = np.asarray(adata.obsm["X_umap"])
    if umap.shape[1] < 2:
        print(f"[ERROR] X_umap shape {umap.shape} not 2-D", file=sys.stderr)
        return 1

    # Use canonical label-transferred annotation (Hepatocytes, Fibroblasts,
    # Macrophages, ... — 10 categories matching snRNA atlas), NOT the legacy
    # gene-activity `cell_type` column.
    ct_col = "cell_type_transferred" if "cell_type_transferred" in adata.obs.columns else "cell_type"
    obs = adata.obs[["donor_id", ct_col]].copy()
    obs = obs.rename(columns={ct_col: "cell_type"})
    obs["x"] = umap[:, 0]
    obs["y"] = umap[:, 1]
    obs["donor"] = obs["donor_id"].astype(str)
    obs = obs.drop(columns=["donor_id"])
    obs["cell_type"] = obs["cell_type"].astype(str)

    # Close the backed handle now that we have what we need.
    adata.file.close()

    print(f"[load] donor metadata {DONOR_META_PATH}", flush=True)
    meta = pd.read_csv(DONOR_META_PATH, sep="\t", dtype=str)
    keep = ["donor_id", "F_stage_augmented", "condition"]
    missing = [c for c in keep if c not in meta.columns]
    if missing:
        print(f"[ERROR] donor metadata missing cols: {missing}", file=sys.stderr)
        return 1
    meta = meta[keep].rename(columns={"donor_id": "donor"})

    print(f"[merge] cells={len(obs):,}  donors_in_meta={meta['donor'].nunique()}",
          flush=True)
    out = obs.merge(meta, on="donor", how="left")

    n_missing = int(out["condition"].isna().sum())
    if n_missing:
        print(f"[warn] {n_missing} cells have no donor metadata match",
              flush=True)

    # Column order — match the figure script expectations.
    out = out[["x", "y", "cell_type", "donor", "F_stage_augmented", "condition"]]

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    print(f"[write] {OUT_PATH}  rows={len(out):,}", flush=True)
    out.to_csv(OUT_PATH, sep="\t", index=False, compression="gzip")
    print("[done]", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
