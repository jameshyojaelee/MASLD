#!/usr/bin/env python3
"""Downsampled UMAP embedding for the /single-cell page (data contract §8 item 6).

Source: the integrated human scRNA atlas
  Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad
  (1,232,318 cells / 7 datasets / 275 samples). It ALREADY carries a precomputed
  UMAP in .obsm['X_umap'] and cell_type / dataset / condition_harmonized in .obs,
  so nothing is recomputed. Opened backed='r' — only .obs + .obsm are read, never
  the expression matrix .X.

Downsample: ~60,000 cells, stratified by cell_type x condition_harmonized, with a
floor of ~200 cells per cell type (where available) so rare types survive. Fixed
seed for reproducibility.

Output (-> --output-dir):
  sc_umap_downsampled.parquet  {cell_id, umap_x, umap_y, cell_type, dataset, condition_harmonized}

Run (compute node):
  srun --partition=cpu --qos=interactive --mem=64G --time=48:00:00 --job-name=portal \
    micromamba run -n spatial python scripts/portal/generate_sc_umap_downsampled.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from _portal_io import project_root, write_parquet

SEED = 42
TARGET_TOTAL = 60_000
FLOOR_PER_CELLTYPE = 200
H5AD_REL = ("Analysis/SingleCell/integration/output/human/"
            "scalesc_human_annotated_celltypist.h5ad")


def _decode(a):
    """bytes -> str for object/bytes arrays; passthrough otherwise."""
    if len(a) and isinstance(a[0], (bytes, bytes)):
        return np.array([x.decode() if isinstance(x, bytes) else x for x in a])
    return a


def read_obs_col(obs: h5py.Group, col: str) -> np.ndarray:
    """Read an obs column from an h5ad `obs` group, decoding anndata categoricals
    (group of categories+codes; code -1 -> None) or plain datasets. Reads directly
    from HDF5 so it is NOT filtered by anndata's `column-order` attr (which omits
    `condition_harmonized` in this file)."""
    o = obs[col]
    if isinstance(o, h5py.Group):
        cats = _decode(o["categories"][:])
        codes = o["codes"][:]
        out = np.array([cats[c] if c >= 0 else None for c in codes], dtype=object)
        return out
    return _decode(o[:])


def largest_remainder(total: int, sizes: dict) -> dict:
    """Distribute `total` across groups proportional to `sizes` (Hamilton method).
    Never allocates more than a group's own size."""
    s = sum(sizes.values())
    if s == 0 or total <= 0:
        return {k: 0 for k in sizes}
    raw = {k: total * v / s for k, v in sizes.items()}
    floor = {k: int(np.floor(x)) for k, x in raw.items()}
    alloc = {k: min(floor[k], sizes[k]) for k in sizes}
    remaining = total - sum(alloc.values())
    # hand out leftover by largest fractional remainder, respecting group caps
    order = sorted(sizes, key=lambda k: raw[k] - floor[k], reverse=True)
    i = 0
    guard = 0
    while remaining > 0 and guard < 10 * len(order) + 10:
        k = order[i % len(order)]
        if alloc[k] < sizes[k]:
            alloc[k] += 1
            remaining -= 1
        i += 1
        guard += 1
    return alloc


def main():
    root = project_root()
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default=str(root / "masld-atlas-v2/public/data"))
    ap.add_argument("--h5ad", default=str(root / H5AD_REL))
    args = ap.parse_args()
    out = Path(args.output_dir)

    print(f"Reading (h5py, obs + X_umap only) {args.h5ad}")
    with h5py.File(args.h5ad, "r") as f:
        obs = f["obs"]
        assert "X_umap" in f["obsm"], "no precomputed X_umap in .obsm"
        umap = f["obsm/X_umap"][:, :2]
        idx_key = obs.attrs.get("_index", "_index")
        idx_key = idx_key.decode() if isinstance(idx_key, bytes) else idx_key
        df = pd.DataFrame({
            "cell_id": _decode(obs[idx_key][:]).astype(str),
            "umap_x": umap[:, 0],
            "umap_y": umap[:, 1],
            "cell_type": read_obs_col(obs, "cell_type"),
            "dataset": read_obs_col(obs, "dataset"),
            "condition_harmonized": read_obs_col(obs, "condition_harmonized"),
        })
    n_total = len(df)
    print(f"  n_obs={n_total}")
    # drop cells lacking a UMAP coord or a cell-type label
    df = df[np.isfinite(df["umap_x"]) & np.isfinite(df["umap_y"])]
    df = df[df["cell_type"].notna()
            & ~df["cell_type"].astype(str).isin(["nan", "NA", "None", ""])]
    df = df.reset_index(drop=True)
    n_clean = len(df)
    print(f"  usable cells (finite UMAP + labeled): {n_clean}")

    rng_base = SEED
    picks = []
    for ct, g in df.groupby("cell_type", sort=True):
        ct_size = len(g)
        prop = round(TARGET_TOTAL * ct_size / n_clean)
        ct_target = min(max(prop, min(FLOOR_PER_CELLTYPE, ct_size)), ct_size)
        cond_sizes = g["condition_harmonized"].value_counts().to_dict()
        alloc = largest_remainder(ct_target, cond_sizes)
        for cond, sub in g.groupby("condition_harmonized", sort=True):
            k = min(alloc.get(cond, 0), len(sub))
            if k > 0:
                picks.append(sub.sample(n=k, random_state=rng_base))
        rng_base += 1

    sampled = pd.concat(picks, ignore_index=True)
    sampled = sampled.sort_values(["cell_type", "condition_harmonized"]).reset_index(drop=True)
    sampled["umap_x"] = sampled["umap_x"].round(3)
    sampled["umap_y"] = sampled["umap_y"].round(3)

    # sort_col=None: not symbol-keyed, keep the cell_type-grouped order
    write_parquet(sampled, out / "sc_umap_downsampled.parquet", sort_col=None)

    ct_counts = sampled["cell_type"].value_counts().sort_index()
    print(f"  -> sc_umap_downsampled.parquet: {len(sampled)} rows, "
          f"{sampled['cell_type'].nunique()} cell types")
    print("  per-cell-type counts:")
    for ct, n in ct_counts.items():
        print(f"      {ct}: {n}")


if __name__ == "__main__":
    main()
