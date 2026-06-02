#!/usr/bin/env python
"""
run_scdrs_stage1_preprocess.py
==============================
Stage 1 of parallelized scDRS pipeline. Builds the .gs file, applies the
prep-method filter, runs scDRS preprocess (the slow ~60min step), and saves
the preprocessed AnnData + .gs file to disk so Stage 2 (SLURM array of
trait-scoring jobs) can load it once per task.

Refactored from run_scdrs_F_transitions.py. Identical inputs & filters.

Outputs:
  Analysis/SingleCell/results_gpu_v2/disease_signatures/
    F_transitions_scdrs.gs
    scdrs_gs_summary.csv
    scdrs_preprocessed_atlas.h5ad   (909237 x ~37k, .uns['SCDRS_PARAM'] populated)
    scdrs_obs_metadata.csv          (cell_id, cell_type, sample, dataset for stage 3)
"""
from __future__ import annotations

import gc
import os
import pickle
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# NumPy 2.0 shim
if not hasattr(np, "float_"):
    np.float_ = np.float64  # type: ignore[attr-defined]
if not hasattr(np, "int_"):
    np.int_ = np.int64  # type: ignore[attr-defined]

import scdrs

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
DEG_CSV   = BASE / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_consecutive_dream.csv"
ATLAS_H5  = BASE / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
UMAP_META = BASE / "Analysis/SingleCell/results_gpu_v2/atlas_umap_for_fig2.csv.gz"

OUT_DIR = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
GS_PATH         = OUT_DIR / "F_transitions_scdrs.gs"
GS_SUMMARY_PATH = OUT_DIR / "scdrs_gs_summary.csv"
PP_ATLAS_PATH   = OUT_DIR / "scdrs_preprocessed_atlas.h5ad"
PP_PARAM_PATH   = OUT_DIR / "scdrs_preprocessed_params.pkl"  # SCDRS_PARAM dict (pd.Series-safe)
OBS_META_PATH   = OUT_DIR / "scdrs_obs_metadata.csv"

PADJ_CUT = 0.05
LFC_CUT  = 0.5
PREP_KEEP = ["unsorted", "nuclei"]
TRANSITIONS = ["F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3"]
RANDOM_SEED = 42


def log(msg: str) -> None:
    print(f"[scdrs stage1] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def strip_ensembl_version(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace(r"\.[0-9]+$", "", regex=True)


def build_gs_lines(deg_df: pd.DataFrame, atlas_var: pd.DataFrame):
    var_map = (
        atlas_var.reset_index().rename(columns={atlas_var.index.name or "index": "var_name"})
        [["gene_ids", "var_name"]]
        .drop_duplicates("gene_ids")
    )
    var_map["gene_ids"] = var_map["gene_ids"].astype(str)
    deg = deg_df.copy()
    deg["gene_ensg"] = strip_ensembl_version(deg["gene"])
    deg = deg.merge(var_map, left_on="gene_ensg", right_on="gene_ids", how="left")
    log(f"DEG rows: {len(deg):,};  mapped: {deg['var_name'].notna().sum():,}")
    deg["sig"] = (deg["padj"] < PADJ_CUT) & (deg["logFC"].abs() > LFC_CUT)
    deg_sig = deg[deg["sig"] & deg["var_name"].notna()].copy()
    deg_sig["abs_t"] = deg_sig["t"].abs()

    lines = ["TRAIT\tGENESET"]
    summary = []
    for tr in TRANSITIONS:
        sub = deg_sig[deg_sig["contrast"] == tr]
        signed = sub[["var_name", "abs_t"]].drop_duplicates("var_name")
        up = sub[sub["t"] > 0][["var_name", "abs_t"]].drop_duplicates("var_name")
        dn = sub[sub["t"] < 0][["var_name", "abs_t"]].drop_duplicates("var_name")
        for mode, gs_df in [("signed", signed), ("up", up), ("down", dn)]:
            tag = f"{tr}__{mode}"
            if len(gs_df) == 0:
                summary.append({"trait": tag, "n_genes": 0})
                continue
            payload = ",".join(f"{g}:{w:.6g}" for g, w in zip(gs_df["var_name"], gs_df["abs_t"]))
            lines.append(f"{tag}\t{payload}")
            summary.append({"trait": tag, "n_genes": int(len(gs_df))})
    return lines, pd.DataFrame(summary)


def main() -> None:
    np.random.seed(RANDOM_SEED)

    log("STEP 1: load bulk DEG table")
    deg_df = pd.read_csv(DEG_CSV, usecols=["gene", "logFC", "padj", "t", "contrast"])
    log(f"contrasts: {sorted(deg_df['contrast'].unique())}")

    log("STEP 2: load atlas .h5ad")
    adata = sc.read_h5ad(ATLAS_H5)
    log(f"atlas: {adata.n_obs:,} x {adata.n_vars:,};  raw: {adata.raw is not None}")

    log("STEP 3: merge prep-method via positional alignment")
    umap_df = pd.read_csv(UMAP_META)
    if len(umap_df) != adata.n_obs:
        raise RuntimeError(f"UMAP meta {len(umap_df):,} != n_obs {adata.n_obs:,}")
    if "cell_type" in umap_df.columns:
        match = (umap_df["cell_type"].astype(str).values == adata.obs["cell_type"].astype(str).values).mean()
        log(f"  positional cell_type match: {match:.4f}")
        if match < 0.99:
            raise RuntimeError("positional misalignment")
    adata.obs["preparation_method"] = umap_df["preparation_method"].astype(str).values
    if "disease_stage_coarse" in umap_df.columns:
        adata.obs["disease_stage_coarse"] = umap_df["disease_stage_coarse"].astype(str).values

    log("STEP 4: filter prep-method")
    keep = adata.obs["preparation_method"].isin(PREP_KEEP).values
    log(f"  keep: {int(keep.sum()):,} / {adata.n_obs:,}")
    adata = adata[keep].copy()
    gc.collect()

    log("STEP 5: build .gs")
    lines, gs_summary = build_gs_lines(deg_df, adata.var)
    gs_summary.to_csv(GS_SUMMARY_PATH, index=False)
    GS_PATH.write_text("\n".join(lines) + "\n")
    log(f"  wrote {GS_PATH} ({len(lines)-1} traits)")
    log(f"  GS summary:\n{gs_summary.to_string(index=False)}")

    log("STEP 6: replace .X with raw counts")
    if adata.raw is None:
        raise RuntimeError("no .raw")
    raw_X = adata.raw.X
    if raw_X.shape[0] != adata.n_obs:
        raw_X = adata.raw[adata.obs_names].X
    if not sp.issparse(raw_X):
        raw_X = sp.csr_matrix(raw_X)
    if not (adata.var_names == adata.raw.var_names).all():
        common_v = adata.var_names.intersection(adata.raw.var_names)
        log(f"  var mismatch; intersect = {len(common_v):,}")
        raw_full = adata.raw.to_adata()
        raw_aligned = raw_full[adata.obs_names, list(common_v)]
        raw_X = raw_aligned.X
        if not sp.issparse(raw_X):
            raw_X = sp.csr_matrix(raw_X)
        adata = adata[:, list(common_v)].copy()
    adata.X = raw_X.astype(np.float32)
    log(f"  X: sparse={sp.issparse(adata.X)}  dtype={adata.X.dtype}  shape={adata.X.shape}")

    log("STEP 7: compute n_genes / pct_mito covariates")
    sc.pp.calculate_qc_metrics(adata, qc_vars=[], percent_top=None, log1p=False, inplace=True)
    mt_mask = np.asarray(adata.var_names.astype(str).str.upper().str.startswith("MT-"))
    if mt_mask.sum() == 0:
        mt_mask = np.asarray(adata.var["gene_ids"].astype(str).str.upper().str.startswith("MT-"))
    mt_idx = np.where(mt_mask)[0]
    if len(mt_idx) > 0:
        mt_counts = np.asarray(adata.X[:, mt_idx].sum(axis=1)).ravel()
        total_counts = np.asarray(adata.X.sum(axis=1)).ravel()
        adata.obs["pct_mito"] = (mt_counts / np.maximum(total_counts, 1) * 100.0).astype(np.float32)
    else:
        adata.obs["pct_mito"] = np.float32(0)
    adata.obs["n_genes"] = adata.obs["n_genes_by_counts"].astype(np.float32)
    log(f"  pct_mito median={float(adata.obs['pct_mito'].median()):.3f}")

    log("STEP 8: assemble covariates (dataset one-hot + n_genes + pct_mito)")
    dataset_dummies = pd.get_dummies(
        adata.obs["dataset"].astype("category"),
        prefix="ds", drop_first=True, dtype=np.float32,
    )
    cov = pd.concat([dataset_dummies, adata.obs[["n_genes", "pct_mito"]].astype(np.float32)], axis=1)
    cov.index = adata.obs_names
    log(f"  cov shape: {cov.shape}")

    log("STEP 9: scdrs.preprocess (the ~60min step)")
    t0 = time.time()
    scdrs.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20)
    log(f"  preprocess done in {time.time()-t0:.1f}s; SCDRS_PARAM keys: {list(adata.uns['SCDRS_PARAM'].keys())}")

    log("STEP 10: write preprocessed atlas + obs metadata")
    # Save thin obs metadata for stage 3 aggregation
    obs_keep = ["cell_type", "sample", "dataset", "preparation_method", "disease_stage_coarse"]
    obs_keep = [c for c in obs_keep if c in adata.obs.columns]
    adata.obs[obs_keep].rename_axis("cell_id").reset_index().to_csv(OBS_META_PATH, index=False)
    log(f"  wrote {OBS_META_PATH}  shape={adata.obs[obs_keep].shape}")

    # Pop SCDRS_PARAM out of adata.uns -- contains pd.Series (e.g., COV_GENE_MEAN)
    # which anndata h5ad cannot serialize. Save as pickle alongside the h5ad and
    # reattach in Stage 2.
    scdrs_param = adata.uns.pop("SCDRS_PARAM")
    log(f"  popped SCDRS_PARAM (keys: {list(scdrs_param.keys())}) for separate pickle write")

    t0 = time.time()
    adata.write_h5ad(PP_ATLAS_PATH, compression="gzip", compression_opts=4)
    log(f"  wrote {PP_ATLAS_PATH}  in {time.time()-t0:.1f}s  ({PP_ATLAS_PATH.stat().st_size/1e9:.2f} GB)")

    t0 = time.time()
    with open(PP_PARAM_PATH, "wb") as f:
        pickle.dump(scdrs_param, f, protocol=pickle.HIGHEST_PROTOCOL)
    log(f"  wrote {PP_PARAM_PATH}  in {time.time()-t0:.1f}s  ({PP_PARAM_PATH.stat().st_size/1e6:.1f} MB)")

    log("STAGE 1 DONE")


if __name__ == "__main__":
    main()
