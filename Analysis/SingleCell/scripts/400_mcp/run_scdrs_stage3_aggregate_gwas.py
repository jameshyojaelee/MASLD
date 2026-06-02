#!/usr/bin/env python
"""
run_scdrs_stage3_aggregate.py
=============================
Stage 3 of parallelized scDRS pipeline. Aggregates per-trait outputs from
Stage 2 array into final figure-ready tables.

Inputs:
  scdrs_trait_outputs/{trait}_cell_scores.parquet (12 files)
  scdrs_trait_outputs/{trait}_celltype_enrichment.csv (12 files)
  scdrs_obs_metadata.csv (from Stage 1)

Outputs:
  scdrs_cell_scores.parquet                       (combined per-cell × 12 traits)
  scdrs_cell_scores_macrophages_signed.csv        (R-friendly slice)
  scdrs_celltype_norm_signed_median.csv           (per CT × transition)
  scdrs_celltype_enrichment.csv                   (combined group analysis)
  scdrs_signed_vs_split_concordance.csv           (signed vs up-down QC)
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SIG_DIR = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
TRAIT_DIR = SIG_DIR / "scdrs_gwas_trait_outputs"
OBS_META_PATH = SIG_DIR / "scdrs_obs_metadata.csv"

CELL_SCORE_PATH = SIG_DIR / "scdrs_gwas_cell_scores.parquet"
MAC_CSV_PATH    = SIG_DIR / "scdrs_gwas_cell_scores_macrophages.csv"
CT_MEDIAN_PATH  = SIG_DIR / "scdrs_gwas_celltype_score_median.csv"
CT_ENRICH_PATH  = SIG_DIR / "scdrs_gwas_celltype_enrichment.csv"
CONCORD_PATH    = SIG_DIR / "scdrs_gwas_signed_vs_split_concordance.csv"


def log(msg: str) -> None:
    print(f"[scdrs stage3] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def main() -> None:
    log("STEP 1: gather per-trait outputs")
    score_files = sorted(TRAIT_DIR.glob("*_cell_scores.parquet"))
    enrich_files = sorted(TRAIT_DIR.glob("*_celltype_enrichment.csv"))
    log(f"  found {len(score_files)} cell-score files, {len(enrich_files)} enrichment files")

    log("STEP 2: load obs metadata")
    obs_meta = pd.read_csv(OBS_META_PATH)
    log(f"  obs: {obs_meta.shape}")

    log("STEP 3: combine per-cell scores")
    cs_dfs = []
    for f in score_files:
        df = pd.read_parquet(f)
        if df.empty:
            log(f"  skipping empty {f.name}")
            continue
        cs_dfs.append(df)
    if not cs_dfs:
        raise RuntimeError("no per-trait cell scores")
    cell_scores = pd.concat(cs_dfs, ignore_index=True)
    # GWAS traits are flat-named like `gwas_coloc_pooled_pp05_1kg` or `sceqtl_hep`.
    # Parse panel suffix (_1kg/_polyfun) and trait_class (panel-free name).
    def parse_trait(t: str):
        if t.endswith("_1kg"):
            return t[:-4], "1kg"
        if t.endswith("_polyfun"):
            return t[:-8], "polyfun"
        return t, "sceqtl"  # sceQTL traits have no panel suffix
    parsed = cell_scores["trait"].map(parse_trait)
    cell_scores["trait_class"] = parsed.map(lambda x: x[0])
    cell_scores["panel"] = parsed.map(lambda x: x[1])
    cell_scores.to_parquet(CELL_SCORE_PATH, index=False, compression="zstd")
    log(f"  wrote {CELL_SCORE_PATH}  shape={cell_scores.shape}")

    log("STEP 4: macrophage-only score CSV slice")
    mac_ids = obs_meta.loc[obs_meta["cell_type"] == "Macrophages", "cell_id"]
    mac_slice = cell_scores[cell_scores["cell_id"].isin(mac_ids)][
        ["cell_id", "trait", "trait_class", "panel", "norm_score", "zscore"]
    ]
    mac_slice.to_csv(MAC_CSV_PATH, index=False)
    log(f"  wrote {MAC_CSV_PATH}  shape={mac_slice.shape}")

    log("STEP 5: per-celltype median per (cell_type x trait_class x panel)")
    cs_ct = cell_scores.merge(
        obs_meta[["cell_id", "cell_type"]], on="cell_id", how="left"
    )
    ct_median = (
        cs_ct.groupby(["cell_type", "trait_class", "panel"])
        .agg(scdrs_norm_median=("norm_score", "median"),
             scdrs_zscore_median=("zscore", "median"),
             n_cells=("norm_score", "size"))
        .reset_index()
    )
    ct_median.to_csv(CT_MEDIAN_PATH, index=False)
    log(f"  wrote {CT_MEDIAN_PATH}  shape={ct_median.shape}")

    log("STEP 6: combine enrichment files")
    enrich_dfs = []
    for f in enrich_files:
        try:
            df = pd.read_csv(f)
            if df.empty:
                continue
            enrich_dfs.append(df)
        except Exception as exc:
            log(f"  skipping {f.name}: {exc}")
    if enrich_dfs:
        enrich_df = pd.concat(enrich_dfs, ignore_index=True)
        parsed_e = enrich_df["trait"].map(parse_trait)
        enrich_df["trait_class"] = parsed_e.map(lambda x: x[0])
        enrich_df["panel"] = parsed_e.map(lambda x: x[1])
        enrich_df.to_csv(CT_ENRICH_PATH, index=False)
        log(f"  wrote {CT_ENRICH_PATH}  shape={enrich_df.shape}")
    else:
        log("  WARNING: no enrichment records")
        enrich_df = pd.DataFrame()

    log("STEP 7: 1KG vs PolyFun cross-panel concordance per (cell_type x trait_class)")
    paired = ct_median[ct_median["panel"].isin(["1kg", "polyfun"])]
    if not paired.empty:
        wide = paired.pivot_table(
            index=["cell_type", "trait_class"],
            columns="panel",
            values="scdrs_zscore_median",
        ).reset_index()
        wide.columns.name = None
        if {"1kg", "polyfun"}.issubset(wide.columns):
            wide["delta_zscore"] = wide["polyfun"] - wide["1kg"]
            wide["abs_delta"] = wide["delta_zscore"].abs()
            wide["sign_concordant"] = (
                np.sign(wide["1kg"].fillna(0)) == np.sign(wide["polyfun"].fillna(0))
            )
            wide.to_csv(CONCORD_PATH, index=False)
            try:
                rho = wide[["1kg", "polyfun"]].corr(method="spearman").iloc[0, 1]
            except Exception:
                rho = float("nan")
            log(f"  wrote {CONCORD_PATH};  Spearman(1kg, polyfun) zscore = {rho:.3f}; "
                f"median |Δ| = {wide['abs_delta'].median():.3f}; "
                f"sign-concordant = {int(wide['sign_concordant'].sum())}/{len(wide)}")
        else:
            log("  WARNING: only one panel present; cross-panel concordance skipped")
    else:
        log("  no panel-tagged rows; cross-panel concordance skipped")

    log("STAGE 3 DONE")


if __name__ == "__main__":
    main()
