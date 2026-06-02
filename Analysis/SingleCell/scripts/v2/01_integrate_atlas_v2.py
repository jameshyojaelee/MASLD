#!/usr/bin/env python3
"""
Phase 0.5 — Step 2 / v2 atlas
============================

Concatenate the 7 per-dataset filtered h5ads (from Step 1) and re-run
the integration pipeline (no scVI in this script — scVI retraining is
Agent 3's job). Produces a CPU-friendly atlas based on harmony or
PCA-only reductions.

Pipeline:
    raw counts -> normalize_total(1e4) -> log1p -> HVG (seurat_v3,
    batch_key='dataset', n_top_genes=3000) -> PCA (50) -> harmony (if
    available; otherwise PCA-only is the fallback) -> neighbors -> UMAP
    -> Leiden.

Inputs:
    Analysis/SingleCell/results_gpu_v2_phase05/atlas/per_dataset/{dataset_id}_v2.h5ad

Outputs:
    Analysis/SingleCell/results_gpu_v2_phase05/atlas/scalesc_human_annotated_celltypist_v2.h5ad
    Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2.h5ad
    Analysis/SingleCell/results_gpu_v2_phase05/atlas/integration_summary_v2.json

Env: rapids_singlecell.
SBATCH (bigmem / qos=interactive / 16 cpus / 500G / 24h).
"""

from __future__ import annotations

import gc
import glob
import json
import logging
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# Determinism
import random
os.environ["PYTHONHASHSEED"] = "42"
random.seed(42)
np.random.seed(42)

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SC_DIR = BASE / "Analysis" / "SingleCell"
PER_DATASET_DIR = SC_DIR / "results_gpu_v2_phase05" / "atlas" / "per_dataset"

OUT_FULL = SC_DIR / "results_gpu_v2_phase05" / "atlas" / "scalesc_human_annotated_celltypist_v2.h5ad"
OUT_HEP = SC_DIR / "results_gpu_v2_phase05" / "atlas" / "hepatocyte_atlas_v2.h5ad"
OUT_SUMMARY = SC_DIR / "results_gpu_v2_phase05" / "atlas" / "integration_summary_v2.json"

N_TOP_GENES = 3000
N_PCS = 50
LEIDEN_RES = 1.0


def _normalize_var_index(adata: ad.AnnData) -> ad.AnnData:
    adata.var_names = pd.Index(adata.var_names.astype(str))
    adata.var_names_make_unique()
    return adata


def _load_one(path: Path) -> ad.AnnData:
    a = sc.read_h5ad(path)
    a = _normalize_var_index(a)
    return a


def _intersect_genes(adatas: list[ad.AnnData]) -> list[str]:
    if not adatas:
        return []
    common = set(adatas[0].var_names)
    for a in adatas[1:]:
        common &= set(a.var_names)
    log.info("Common gene set across datasets: %d", len(common))
    return sorted(common)


def main() -> int:
    OUT_FULL.parent.mkdir(parents=True, exist_ok=True)

    # 1. Load all per-dataset filtered h5ads
    paths = sorted(glob.glob(str(PER_DATASET_DIR / "*_v2.h5ad")))
    log.info("Per-dataset filtered files: %d", len(paths))
    if not paths:
        log.error("No per-dataset h5ads found under %s", PER_DATASET_DIR)
        return 1
    for p in paths:
        log.info("  %s", p)

    log.info("Reading datasets ...")
    adatas = []
    for p in paths:
        a = _load_one(Path(p))
        log.info("  loaded %s: %s", Path(p).name, a.shape)
        adatas.append(a)

    # 2. Intersect var to a common gene panel
    common_genes = _intersect_genes(adatas)
    if len(common_genes) < 1000:
        log.error("Too few common genes (%d) — abort", len(common_genes))
        return 2

    for i, a in enumerate(adatas):
        adatas[i] = a[:, common_genes].copy()
        del a
    gc.collect()

    log.info("Concatenating ...")
    atlas = ad.concat(
        adatas,
        join="inner",
        merge="same",
        label="dataset",  # safety: relabel even though obs has 'dataset'
        index_unique=None,
    )
    # restore obs.dataset from the per-dataset file (concat label may stomp it)
    # The per-dataset h5ad already has obs['dataset']; ad.concat does NOT overwrite obs
    # unless label='dataset' AND obs is missing. Force it explicitly.
    if "dataset" not in atlas.obs.columns or atlas.obs["dataset"].isna().any():
        # Reconstruct from per-dataset file order
        labels = []
        for a in adatas:
            n = a.n_obs
            ds_vals = a.obs.get("dataset", pd.Series(["unknown"] * n)).astype(str)
            labels.append(ds_vals.values)
        atlas.obs["dataset"] = np.concatenate(labels)

    log.info("Concatenated atlas: %s", atlas.shape)
    del adatas
    gc.collect()

    # 3. Save raw counts in a layer + .raw
    atlas.layers["counts"] = atlas.X.copy()
    atlas.raw = atlas.copy()

    # 4. Normalize + log1p
    log.info("normalize_total(target_sum=1e4) + log1p ...")
    sc.pp.normalize_total(atlas, target_sum=1e4)
    sc.pp.log1p(atlas)

    # 5. HVG seurat_v3 with batch_key='dataset' on RAW counts
    log.info("HVG (seurat_v3, batch_key='dataset', n_top_genes=%d) ...", N_TOP_GENES)
    try:
        sc.pp.highly_variable_genes(
            atlas, n_top_genes=N_TOP_GENES, flavor="seurat_v3",
            batch_key="dataset", layer="counts",
        )
    except Exception as e:
        log.warning("seurat_v3 HVG failed (%s) — falling back to seurat flavor on log-data", e)
        sc.pp.highly_variable_genes(
            atlas, n_top_genes=N_TOP_GENES, flavor="seurat",
            batch_key="dataset",
        )

    # exclude MT/RP* from HVG list
    mt_mask = atlas.var_names.str.startswith("MT-")
    ribo_mask = atlas.var_names.str.startswith("RPS") | atlas.var_names.str.startswith("RPL")
    excl = mt_mask | ribo_mask
    atlas.var.loc[excl, "highly_variable"] = False
    n_hvg = int(atlas.var["highly_variable"].sum())
    log.info("HVG count after MT/RP exclusion: %d", n_hvg)

    # 6. Scale + PCA on HVG only
    log.info("Subsetting to HVG and scaling ...")
    atlas_hvg = atlas[:, atlas.var["highly_variable"]].copy()
    sc.pp.scale(atlas_hvg, max_value=10)

    log.info("PCA (n_comps=%d) ...", N_PCS)
    sc.pp.pca(atlas_hvg, n_comps=N_PCS, zero_center=True, svd_solver="arpack",
              random_state=42)

    # 7. Harmony (best-effort; pure-CPU fallback to PCA only)
    use_rep = "X_pca"
    try:
        import harmonypy
        log.info("Running Harmony on PCA (batch='dataset') ...")
        import numpy as _np
        hpy = harmonypy.run_harmony(
            atlas_hvg.obsm["X_pca"], atlas_hvg.obs, ["dataset"],
            random_state=42, max_iter_harmony=20,
        )
        atlas_hvg.obsm["X_harmony"] = _np.array(hpy.Z_corr).T
        use_rep = "X_harmony"
        log.info("Harmony done: %s", atlas_hvg.obsm["X_harmony"].shape)
    except Exception as exc:
        log.warning("Harmony unavailable / failed (%s) — falling back to PCA (no batch correction)", exc)

    # 8. Carry the chosen rep back to the full atlas
    atlas.obsm["X_pca"] = atlas_hvg.obsm["X_pca"]
    if "X_harmony" in atlas_hvg.obsm:
        atlas.obsm["X_harmony"] = atlas_hvg.obsm["X_harmony"]

    # 9. Neighbors + UMAP + Leiden (CPU, deterministic)
    log.info("Neighbors (use_rep=%s, n_neighbors=15) ...", use_rep)
    sc.pp.neighbors(atlas, use_rep=use_rep, n_neighbors=15, random_state=42)
    log.info("UMAP ...")
    sc.tl.umap(atlas, random_state=42)
    log.info("Leiden (resolution=%.2f) ...", LEIDEN_RES)
    sc.tl.leiden(atlas, resolution=LEIDEN_RES, random_state=42, key_added="leiden_v2")

    del atlas_hvg
    gc.collect()

    # 10. Persist
    log.info("Writing full atlas (%s) ...", atlas.shape)
    atlas.write_h5ad(OUT_FULL, compression="gzip")
    log.info("Wrote %s", OUT_FULL)

    # 11. Hepatocyte subset
    if "cell_type" not in atlas.obs.columns:
        log.warning("cell_type missing — skipping hepatocyte subset")
        hep_n = 0
    else:
        hep_mask = atlas.obs["cell_type"].astype(str) == "Hepatocytes"
        hep_n = int(hep_mask.sum())
        log.info("Hepatocyte subset: %d cells", hep_n)
        if hep_n > 0:
            atlas_hep = atlas[hep_mask].copy()
            atlas_hep.write_h5ad(OUT_HEP, compression="gzip")
            log.info("Wrote %s (%s)", OUT_HEP, atlas_hep.shape)
            del atlas_hep
            gc.collect()

    summary = {
        "v2_atlas_path": str(OUT_FULL),
        "v2_hepatocyte_path": str(OUT_HEP),
        "n_cells_v2": int(atlas.n_obs),
        "n_genes_v2": int(atlas.n_vars),
        "n_hvg": n_hvg,
        "use_rep": use_rep,
        "leiden_resolution": LEIDEN_RES,
        "n_hepatocytes_v2": hep_n,
        "dataset_breakdown": atlas.obs["dataset"].astype(str).value_counts().to_dict(),
    }
    with open(OUT_SUMMARY, "w") as f:
        json.dump(summary, f, indent=2)
    log.info("Wrote summary %s", OUT_SUMMARY)
    log.info("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
