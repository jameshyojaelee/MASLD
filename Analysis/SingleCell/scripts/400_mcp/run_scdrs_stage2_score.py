#!/usr/bin/env python
"""
run_scdrs_stage2_score.py
=========================
Stage 2 of parallelized scDRS pipeline. Loads the preprocessed atlas from
Stage 1, scores ONE trait specified via SCDRS_TRAIT env var (set by SLURM
array task index), and writes per-cell + per-celltype outputs for that trait.

Submitted as a SLURM array (--array=0-11), one task per trait.

Inputs:
  SCDRS_TRAIT_INDEX env var (0-11, maps to .gs file row order)
  scdrs_preprocessed_atlas.h5ad
  F_transitions_scdrs.gs

Outputs (per trait):
  scdrs_trait_outputs/{trait}_cell_scores.parquet
  scdrs_trait_outputs/{trait}_celltype_enrichment.csv  (downstream_group_analysis)
"""
from __future__ import annotations

import gc
import os
import pickle
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# NumPy 2.0 shim
if not hasattr(np, "float_"):
    np.float_ = np.float64  # type: ignore[attr-defined]
if not hasattr(np, "int_"):
    np.int_ = np.int64  # type: ignore[attr-defined]

import scdrs
from scdrs.method import downstream_group_analysis

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SIG_DIR = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
PP_ATLAS_PATH = SIG_DIR / "scdrs_preprocessed_atlas.h5ad"
PP_PARAM_PATH = SIG_DIR / "scdrs_preprocessed_params.pkl"
GS_PATH = SIG_DIR / "F_transitions_scdrs.gs"
OUT_DIR = SIG_DIR / "scdrs_trait_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_CTRL = 1000
WEIGHT_OPT = "vs"
RANDOM_SEED = 42


def log(msg: str) -> None:
    print(f"[scdrs stage2 t={TRAIT_INDEX}] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


# Read trait index from SLURM array task index or env var
TRAIT_INDEX = int(os.environ.get("SCDRS_TRAIT_INDEX",
                                  os.environ.get("SLURM_ARRAY_TASK_ID", "0")))


def main() -> None:
    np.random.seed(RANDOM_SEED + TRAIT_INDEX)

    log("STEP 1: load .gs file and select trait")
    gs_dict = scdrs.util.load_gs(str(GS_PATH))
    trait_names = list(gs_dict.keys())
    if TRAIT_INDEX < 0 or TRAIT_INDEX >= len(trait_names):
        raise RuntimeError(f"trait index {TRAIT_INDEX} out of range (0..{len(trait_names)-1})")
    trait = trait_names[TRAIT_INDEX]
    genes, weights = gs_dict[trait]
    log(f"  trait: {trait};  n_genes={len(genes)}")

    out_score = OUT_DIR / f"{trait}_cell_scores.parquet"
    out_enrich = OUT_DIR / f"{trait}_celltype_enrichment.csv"
    if out_score.exists() and out_enrich.exists():
        log(f"  outputs already exist; skipping. To re-run, delete {out_score} and {out_enrich}")
        return

    log("STEP 2: load preprocessed atlas + reattach SCDRS_PARAM from pickle")
    t0 = time.time()
    adata = ad.read_h5ad(PP_ATLAS_PATH)
    log(f"  atlas: {adata.n_obs:,} x {adata.n_vars:,};  loaded in {time.time()-t0:.1f}s")
    with open(PP_PARAM_PATH, "rb") as f:
        adata.uns["SCDRS_PARAM"] = pickle.load(f)
    log(f"  SCDRS_PARAM reattached;  keys: {list(adata.uns['SCDRS_PARAM'].keys())}")

    log("STEP 3: restrict trait genes to atlas vars")
    gene_set = pd.Series(genes)
    weight_set = pd.Series(weights)
    in_atlas = gene_set.isin(adata.var_names)
    if int(in_atlas.sum()) == 0:
        log(f"  WARNING: 0 trait genes in atlas; writing empty outputs")
        pd.DataFrame().to_parquet(out_score)
        pd.DataFrame().to_csv(out_enrich)
        return
    gene_use = gene_set[in_atlas].tolist()
    weight_use = weight_set[in_atlas].tolist()
    log(f"  trait genes in atlas: {len(gene_use)} / {len(gene_set)}")

    log("STEP 4: score_cell")
    t1 = time.time()
    df_score = scdrs.score_cell(
        data=adata,
        gene_list=gene_use,
        gene_weight=weight_use,
        ctrl_match_key="mean_var",
        n_ctrl=N_CTRL,
        weight_opt=WEIGHT_OPT,
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        random_seed=RANDOM_SEED + TRAIT_INDEX,
        verbose=False,
    )
    log(f"  done in {time.time()-t1:.1f}s; rows={len(df_score)}")
    log(f"  norm_score: median={df_score['norm_score'].median():.3f}  "
        f"  zscore: median={df_score['zscore'].median():.3f}  "
        f"  mc_pval<0.05: {(df_score['mc_pval']<0.05).sum()}/{len(df_score)}")

    log("STEP 5: write per-cell scores parquet")
    rec = df_score[["raw_score", "norm_score", "mc_pval", "pval", "nlog10_pval", "zscore"]].copy()
    rec["trait"] = trait
    rec.index.name = "cell_id"
    rec.reset_index().to_parquet(out_score, index=False, compression="zstd")
    log(f"  wrote {out_score}  shape={rec.shape}")

    log("STEP 6: per-celltype group analysis")
    try:
        res = downstream_group_analysis(
            adata=adata,
            df_full_score=df_score,
            group_cols=["cell_type"],
            fdr_thresholds=[0.05, 0.1, 0.2],
        )
        ct_df = res.get("cell_type")
        if ct_df is not None:
            ct_df = ct_df.copy()
            ct_df["trait"] = trait
            ct_df["cell_type"] = ct_df.index
            ct_df = ct_df.reset_index(drop=True)
            ct_df.to_csv(out_enrich, index=False)
            log(f"  wrote {out_enrich}  shape={ct_df.shape}  assoc_mcp_min={ct_df['assoc_mcp'].min():.3g}")
    except Exception as exc:
        log(f"  group analysis FAILED: {exc}")
        pd.DataFrame().to_csv(out_enrich)

    log(f"STAGE 2 trait {TRAIT_INDEX} ({trait}) DONE")


if __name__ == "__main__":
    main()
