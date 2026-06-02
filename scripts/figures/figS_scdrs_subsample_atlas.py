#!/usr/bin/env python
"""
figS_scdrs_subsample_atlas.py
=============================
Create a 25% deterministic subsample of the preprocessed scDRS atlas so
that scoring runs fit on the cpu partition (210GB cap; full atlas needs
~470GB peak on bigmem). Used by the parallel cpu pipeline for the
hep-expression-matched null bias check.

Subsampling is stratified by (cell_type, disease_stage_coarse) so the rare
populations (basophils, plasma cells) retain enough cells for group-analysis
z-scores. SCDRS_PARAM (per-gene control bin assignments) is computed on the
FULL atlas and copied unchanged into the subsample — population-level
normalization is preserved.

Outputs:
  Analysis/SingleCell/results_gpu_v2/disease_signatures/
    scdrs_preprocessed_atlas_25pct.h5ad
    scdrs_preprocessed_params_25pct.pkl   (identical to full; mirrored for path symmetry)
    subsample_25pct_cell_ids.txt          (audit log)
"""
from __future__ import annotations

import os
import pickle
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SIG_DIR        = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
FULL_ATLAS     = SIG_DIR / "scdrs_preprocessed_atlas.h5ad"
FULL_PARAMS    = SIG_DIR / "scdrs_preprocessed_params.pkl"
SUB_ATLAS      = SIG_DIR / "scdrs_preprocessed_atlas_25pct.h5ad"
SUB_PARAMS     = SIG_DIR / "scdrs_preprocessed_params_25pct.pkl"
SUB_AUDIT      = SIG_DIR / "subsample_25pct_cell_ids.txt"

FRAC = 0.25
RANDOM_SEED = 42
STRATA_COLS = ["cell_type", "disease_stage_coarse"]


def log(msg: str) -> None:
    print(f"[subsample] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def main() -> None:
    log(f"loading full atlas: {FULL_ATLAS}")
    t0 = time.time()
    adata = ad.read_h5ad(FULL_ATLAS)
    log(f"  full atlas: {adata.n_obs:,} x {adata.n_vars:,};  loaded in {time.time()-t0:.1f}s")

    # Stratified random sample at FRAC of each (CT × stage) bin
    rng = np.random.default_rng(RANDOM_SEED)
    obs = adata.obs.copy()
    log(f"  stratifying by {STRATA_COLS}")
    obs["_strata"] = obs[STRATA_COLS].astype(str).agg("|".join, axis=1)
    pick_idx = []
    for _, group_idx in obs.groupby("_strata", observed=True).groups.items():
        n_in = len(group_idx)
        n_pick = int(round(n_in * FRAC))
        if n_pick == 0 and n_in > 0:
            n_pick = 1
        picked = rng.choice(list(group_idx), n_pick, replace=False)
        pick_idx.extend(picked.tolist())
    pick_idx = sorted(set(pick_idx))
    log(f"  total picked: {len(pick_idx):,}  ({len(pick_idx)/adata.n_obs*100:.1f}% of full)")

    sub = adata[pick_idx, :].copy()
    log(f"  subsample: {sub.n_obs:,} x {sub.n_vars:,}")

    # Load SCDRS_PARAM from the full atlas's params pickle.
    # IMPORTANT: do NOT stash it in sub.uns — it contains pandas Series objects
    # that anndata's h5ad writer can't serialize (IORegistryError on COV_GENE_MEAN).
    # The scoring script reads SCDRS_PARAM from the separate pickle path and
    # attaches it back to the loaded h5ad at scoring time, so write_h5ad must
    # NOT see it.
    with open(FULL_PARAMS, "rb") as f:
        params = pickle.load(f)
    if "SCDRS_PARAM" in sub.uns:
        del sub.uns["SCDRS_PARAM"]
    log(f"  SCDRS_PARAM kept out of uns; will be loaded from pickle at scoring time")
    log(f"  SCDRS_PARAM keys: {list(params.keys())}")

    log(f"writing subsample atlas: {SUB_ATLAS}")
    t1 = time.time()
    sub.write_h5ad(SUB_ATLAS, compression="gzip")
    log(f"  wrote in {time.time()-t1:.1f}s")

    with open(SUB_PARAMS, "wb") as f:
        pickle.dump(params, f)
    log(f"  wrote {SUB_PARAMS}")

    # Audit log of which cells were picked
    with open(SUB_AUDIT, "w") as f:
        for c in sub.obs_names:
            f.write(f"{c}\n")
    log(f"  wrote {SUB_AUDIT}")

    # Strata composition summary
    log("strata sizes in subsample (top 20 by count):")
    sub_obs = sub.obs[STRATA_COLS].copy()
    sub_obs["_strata"] = sub_obs.astype(str).agg("|".join, axis=1)
    counts = sub_obs["_strata"].value_counts().head(20)
    for s, n in counts.items():
        log(f"  {s}: {n:,}")

    log("done.")


if __name__ == "__main__":
    main()
