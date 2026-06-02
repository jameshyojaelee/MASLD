#!/usr/bin/env python3
"""Aggregate per-cell consensus stage signature scores to (cell_type x signature)
mean matrix for fig2 cell-type transition program panels.

Inputs:
  - scored_atlas.h5ad  (obs/cell_type + obs/_index = 1.23M cell barcodes)
  - scores_consensus.csv.gz  (cell barcode index x stage / ct signatures)

Output:
  - celltype_stage_mean_consensus.csv  (rows = cell_type, cols = signature)
  - celltype_stage_n_cells.csv         (cells per cell_type used in mean)
"""

import os
import h5py
import numpy as np
import pandas as pd

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SIG_DIR = os.path.join(
    BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures"
)
H5AD = os.path.join(SIG_DIR, "scored_atlas.h5ad")
SCORES = os.path.join(SIG_DIR, "scores_consensus.csv.gz")
OUT_MEAN = os.path.join(SIG_DIR, "celltype_stage_mean_consensus.csv")
OUT_N = os.path.join(SIG_DIR, "celltype_stage_n_cells.csv")


def read_obs_categorical(f, key):
    obj = f[f"obs/{key}"]
    cats = obj["categories"][:]
    if cats.dtype.kind == "O":
        cats = np.array([c.decode("utf-8") if isinstance(c, bytes) else str(c)
                         for c in cats])
    elif cats.dtype.kind == "S":
        cats = cats.astype(str)
    codes = obj["codes"][:]
    return cats[codes]


def main():
    print(f"Reading cell_type + barcodes from {H5AD}")
    with h5py.File(H5AD, "r") as f:
        cell_types = read_obs_categorical(f, "cell_type")
        index_obj = f["obs/_index"]
        barcodes = index_obj[:]
        if barcodes.dtype.kind in ("S", "O"):
            barcodes = np.array([b.decode("utf-8") if isinstance(b, bytes) else str(b)
                                 for b in barcodes])
    print(f"  {len(barcodes):,} cells, {len(np.unique(cell_types))} cell types")

    ct_lookup = pd.Series(cell_types, index=pd.Index(barcodes, name="cell"))

    print(f"Reading {SCORES}")
    scores = pd.read_csv(SCORES, index_col=0)
    print(f"  shape: {scores.shape[0]:,} cells x {scores.shape[1]} signatures")

    common = scores.index.intersection(ct_lookup.index)
    if len(common) < len(scores):
        print(f"  joining: {len(common):,} of {len(scores):,} cells matched")
    scores = scores.loc[common]
    scores["cell_type"] = ct_lookup.loc[common].values

    mean_mat = scores.groupby("cell_type").mean(numeric_only=True)
    n_cells = scores.groupby("cell_type").size().rename("n_cells")

    mean_mat.to_csv(OUT_MEAN)
    n_cells.to_csv(OUT_N)
    print(f"  wrote {OUT_MEAN}  ({mean_mat.shape[0]} x {mean_mat.shape[1]})")
    print(f"  wrote {OUT_N}     ({len(n_cells)} cell types)")


if __name__ == "__main__":
    main()
