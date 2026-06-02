#!/usr/bin/env python
"""
run_scdrs_score_25pct.py
========================
Score one trait on the 25% subsample atlas, for the parallel cpu pipeline
of the hep-expression-matched null bias check. Variant of
Analysis/SingleCell/scripts/400_mcp/run_scdrs_stage2_score_gwas.py.

Env vars:
  SCDRS_GS         "real" | "null"      — which .gs file to read
  SCDRS_TRAIT_INDEX  integer            — trait row index in selected .gs file
                                          (set by SLURM array task ID)

Outputs (per trait, under scdrs_*_25pct/):
  {trait}_celltype_enrichment.csv
"""
from __future__ import annotations

import os
import pickle
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

if not hasattr(np, "float_"):
    np.float_ = np.float64  # type: ignore[attr-defined]
if not hasattr(np, "int_"):
    np.int_ = np.int64      # type: ignore[attr-defined]

import scdrs
from scdrs.method import downstream_group_analysis

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SIG_DIR = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
PP_ATLAS_PATH = SIG_DIR / "scdrs_preprocessed_atlas_25pct.h5ad"
PP_PARAM_PATH = SIG_DIR / "scdrs_preprocessed_params_25pct.pkl"

GS_KIND = os.environ.get("SCDRS_GS", "null").lower()
TRAIT_INDEX = int(os.environ.get(
    "SCDRS_TRAIT_INDEX",
    os.environ.get("SLURM_ARRAY_TASK_ID", "0"),
))

if GS_KIND == "real":
    GS_PATH = SIG_DIR / "GWAS_anchored_scdrs.gs"
    OUT_DIR = SIG_DIR / "scdrs_gwas_trait_outputs_25pct"
    # Score only the 4 polyfun phenotypes from the 17-trait real file
    REAL_TARGETS = [
        "gwas_coloc_nafld_specific_polyfun",
        "gwas_coloc_liver_enzymes_polyfun",
        "gwas_coloc_pdff_polyfun",
        "gwas_coloc_cirrhosis_hcc_polyfun",
    ]
elif GS_KIND == "null":
    GS_PATH = SIG_DIR / "GWAS_anchored_hep_matched_nulls.gs"
    OUT_DIR = SIG_DIR / "scdrs_gwas_hep_null_outputs_25pct"
    REAL_TARGETS = None
else:
    raise ValueError(f"SCDRS_GS must be 'real' or 'null'; got {GS_KIND!r}")

OUT_DIR.mkdir(parents=True, exist_ok=True)
N_CTRL = 1000
WEIGHT_OPT = "vs"
RANDOM_SEED = 42


def log(msg: str) -> None:
    tag = f"{GS_KIND}_25pct t={TRAIT_INDEX}"
    print(f"[scdrs {tag}] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def main() -> None:
    np.random.seed(RANDOM_SEED + TRAIT_INDEX)

    log(f"STEP 1: load .gs ({GS_PATH.name}) and select trait")
    gs_dict = scdrs.util.load_gs(str(GS_PATH))
    trait_names = list(gs_dict.keys())
    if REAL_TARGETS is not None:
        # Filter to polyfun phenotypes only
        trait_names = [t for t in trait_names if t in REAL_TARGETS]
        if not trait_names:
            raise RuntimeError("No real polyfun targets found in .gs")
    if TRAIT_INDEX < 0 or TRAIT_INDEX >= len(trait_names):
        raise RuntimeError(f"trait index {TRAIT_INDEX} out of range (0..{len(trait_names)-1})")
    trait = trait_names[TRAIT_INDEX]
    genes, weights = gs_dict[trait]
    log(f"  trait: {trait};  n_genes={len(genes)}")

    out_enrich = OUT_DIR / f"{trait}_celltype_enrichment.csv"
    out_sentinel = OUT_DIR / f"{trait}_DONE.txt"
    if out_enrich.exists() and out_sentinel.exists():
        log("  outputs already exist; skipping")
        return

    log(f"STEP 2: load subsample atlas + SCDRS_PARAM (path: {PP_ATLAS_PATH})")
    t0 = time.time()
    adata = ad.read_h5ad(PP_ATLAS_PATH)
    log(f"  atlas: {adata.n_obs:,} x {adata.n_vars:,};  loaded in {time.time()-t0:.1f}s")
    with open(PP_PARAM_PATH, "rb") as f:
        adata.uns["SCDRS_PARAM"] = pickle.load(f)

    log("STEP 3: restrict trait genes to atlas vars")
    gs_genes = pd.Series(genes)
    gs_w     = pd.Series(weights)
    in_atlas = gs_genes.isin(adata.var_names)
    if int(in_atlas.sum()) == 0:
        log("  WARNING: 0 trait genes in atlas; writing empty outputs")
        pd.DataFrame().to_csv(out_enrich)
        out_sentinel.write_text("empty\n")
        return
    gene_use   = gs_genes[in_atlas].tolist()
    weight_use = gs_w[in_atlas].tolist()
    log(f"  trait genes in atlas: {len(gene_use)} / {len(gs_genes)}")

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
    log(f"  scored in {time.time()-t1:.1f}s;  rows={len(df_score)}")

    log("STEP 5: per-celltype group analysis")
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
            out_sentinel.write_text("done\n")
            log(f"  wrote {out_enrich}  shape={ct_df.shape}")
    except Exception as exc:
        log(f"  group analysis FAILED: {exc}")
        pd.DataFrame().to_csv(out_enrich)

    log(f"DONE  trait_index={TRAIT_INDEX}  trait={trait}  ({GS_KIND}_25pct)")


if __name__ == "__main__":
    main()
