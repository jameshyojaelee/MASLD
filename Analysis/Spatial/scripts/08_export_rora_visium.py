#!/usr/bin/env python
"""
08_export_rora_visium.py — Extract per-spot RORA expression from the human
Visium atlas for the Fig 3c RORA case-study panel (Panel F).

Replaces an ad-hoc `/tmp/rora_visium.csv` step. Concatenates the two human
Visium sources that the figure was built against (15 sections total):
  - merged_spatial.h5ad      : GSE192741, 5 sections (Healthy + Steatotic;
                               JBO014/015/018/019/022)
  - merged_spatial_vu.h5ad   : Vu et al. 2025, 10 sections (MASLD_spectrum;
                               VLP115/116/119/120/121 x {A,D})
and writes a stable, persistent CSV:

    Analysis/Spatial/results/rora_case_study/rora_visium.csv

Output columns consumed by scripts/figures/fig_rora_case_study.R (Panel F):
    sample_id, condition, x, y, RORA
One row per Visium spot, all sections. Both .X layers are normalized log1p
from 02_build_anndata.py, so they are directly comparable.
"""
import pathlib
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

BASE = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
H5ADS = [
    BASE / "Analysis/Spatial/results/preprocessed/merged_spatial.h5ad",     # GSE192741
    BASE / "Analysis/Spatial/results/preprocessed/merged_spatial_vu.h5ad",  # Vu 2025
]
OUT_DIR = BASE / "Analysis/Spatial/results/rora_case_study"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT_DIR / "rora_visium.csv"

GENE = "RORA"


def extract_one(h5ad: pathlib.Path) -> pd.DataFrame:
    print(f"Loading {h5ad}")
    adata = ad.read_h5ad(h5ad)
    print(f"  shape: {adata.shape}")

    # Keep human spots only
    if "species" in adata.obs.columns:
        before = adata.n_obs
        adata = adata[adata.obs["species"].astype(str).str.lower().isin(
            ["human", "homo_sapiens", "hsapiens", "unknown", "nan"])].copy()
        if before != adata.n_obs:
            print(f"  species filter: {before} -> {adata.n_obs}")

    if GENE not in adata.var_names:
        sys.exit(f"ERROR: {GENE} not in var_names of {h5ad.name} ({adata.n_vars} genes)")

    expr = adata[:, GENE].X
    expr = (np.asarray(expr.todense()).ravel() if sp.issparse(expr)
            else np.asarray(expr).ravel())

    if "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"])
        x, y = coords[:, 0], coords[:, 1]
    else:
        x = adata.obs.get("array_col", pd.Series(np.arange(adata.n_obs))).to_numpy()
        y = adata.obs.get("array_row", pd.Series(np.arange(adata.n_obs))).to_numpy()

    cond_col = "condition" if "condition" in adata.obs.columns else None
    df = pd.DataFrame({
        "sample_id": adata.obs["sample_id"].astype(str).to_numpy(),
        "condition": (adata.obs[cond_col].astype(str).to_numpy()
                      if cond_col else np.full(adata.n_obs, "unknown")),
        "x": x.astype(float),
        "y": y.astype(float),
        "RORA": expr.astype(float),
    })
    return df.dropna(subset=["x", "y"]).reset_index(drop=True)


def main():
    parts = [extract_one(h) for h in H5ADS if h.exists()]
    if not parts:
        sys.exit("ERROR: no input h5ad found")
    df = pd.concat(parts, ignore_index=True)

    df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {len(df):,} spots x {df.shape[1]} cols -> {OUT_CSV}")
    print("\n=== sections (sample_id) ===")
    summ = (df.groupby(["sample_id", "condition"])
              .agg(n_spots=("RORA", "size"), mean_RORA=("RORA", "mean"))
              .reset_index())
    print(summ.to_string(index=False))
    print(f"\nn sections: {df['sample_id'].nunique()}")
    print(f"condition values: {sorted(df['condition'].unique().tolist())}")


if __name__ == "__main__":
    main()
