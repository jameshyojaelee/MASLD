"""352a2_add_harmony.py — add X_harmony to the v2 annotated atlas.

The earlier exploration agent reported X_harmony was in the v2 atlas obsm.
Direct inspection shows obsm_keys=['X_pca','X_scVI','X_umap'] only -- the
01_integrate_atlas_v2 Harmony step did not persist its output. Compute it
now in-place so the three-way method test has a defensible Harmony latent.

Strategy: run scanpy.external.pp.harmony_integrate on X_pca with batch_key
'dataset' (same as scVI). Adds obsm['X_harmony'] (50-dim).
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import harmonypy as hm

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS = ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_annotated.h5ad"


def main():
    log.info(f"loading: {ATLAS}")
    a = ad.read_h5ad(ATLAS)
    log.info(f"  shape: {a.shape}; obsm: {list(a.obsm)}")
    if "X_harmony" in a.obsm:
        log.info("X_harmony already present; skipping")
        return 0
    if "X_pca" not in a.obsm:
        log.error("X_pca not in obsm; cannot harmonize")
        return 1
    if "dataset" not in a.obs.columns:
        log.error("'dataset' obs column missing")
        return 1

    log.info("running harmonypy.run_harmony on X_pca (batch_key='dataset')")
    pca_mat = np.asarray(a.obsm["X_pca"], dtype=np.float64)
    meta_df = a.obs[["dataset"]].copy()
    ho = hm.run_harmony(pca_mat, meta_df, "dataset",
                        max_iter_harmony=10)
    # ho.Z_corr is (n_pcs, n_cells); transpose to (n_cells, n_pcs)
    z_corr = np.asarray(ho.Z_corr)
    if z_corr.shape[0] == pca_mat.shape[1] and z_corr.shape[1] == pca_mat.shape[0]:
        z_corr = z_corr.T
    elif z_corr.shape != pca_mat.shape:
        log.error(f"unexpected Z_corr shape {z_corr.shape}; pca was {pca_mat.shape}")
        return 2
    a.obsm["X_harmony"] = z_corr.astype(np.float32)
    log.info(f"  X_harmony shape: {a.obsm['X_harmony'].shape}")

    log.info(f"writing back -> {ATLAS}")
    a.write_h5ad(ATLAS, compression="gzip")
    log.info("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
