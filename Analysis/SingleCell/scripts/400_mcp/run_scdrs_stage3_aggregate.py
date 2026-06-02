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
TRAIT_DIR = SIG_DIR / "scdrs_trait_outputs"
OBS_META_PATH = SIG_DIR / "scdrs_obs_metadata.csv"

CELL_SCORE_PATH = SIG_DIR / "scdrs_cell_scores.parquet"
MAC_CSV_PATH    = SIG_DIR / "scdrs_cell_scores_macrophages_signed.csv"
CT_MEDIAN_PATH  = SIG_DIR / "scdrs_celltype_norm_signed_median.csv"
CT_ENRICH_PATH  = SIG_DIR / "scdrs_celltype_enrichment.csv"
CONCORD_PATH    = SIG_DIR / "scdrs_signed_vs_split_concordance.csv"


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
    parts = cell_scores["trait"].str.split("__", n=1, expand=True)
    cell_scores["transition"] = parts[0]
    cell_scores["mode"] = parts[1]
    cell_scores.to_parquet(CELL_SCORE_PATH, index=False, compression="zstd")
    log(f"  wrote {CELL_SCORE_PATH}  shape={cell_scores.shape}")

    log("STEP 4: macrophage-only signed CSV slice")
    mac_ids = obs_meta.loc[obs_meta["cell_type"] == "Macrophages", "cell_id"]
    mac_slice = cell_scores[
        (cell_scores["mode"] == "signed")
        & cell_scores["cell_id"].isin(mac_ids)
    ][["cell_id", "transition", "norm_score", "zscore"]]
    mac_slice.to_csv(MAC_CSV_PATH, index=False)
    log(f"  wrote {MAC_CSV_PATH}  shape={mac_slice.shape}")

    log("STEP 5: per-celltype median of signed norm_score")
    cs_signed = cell_scores[cell_scores["mode"] == "signed"].merge(
        obs_meta[["cell_id", "cell_type"]], on="cell_id", how="left"
    )
    ct_median = (
        cs_signed.groupby(["cell_type", "transition"])
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
        parts = enrich_df["trait"].str.split("__", n=1, expand=True)
        enrich_df["transition"] = parts[0]
        enrich_df["mode"] = parts[1]
        enrich_df.to_csv(CT_ENRICH_PATH, index=False)
        log(f"  wrote {CT_ENRICH_PATH}  shape={enrich_df.shape}")
    else:
        log("  WARNING: no enrichment records")
        enrich_df = pd.DataFrame()

    log("STEP 7: signed vs (up - down) concordance")
    if not enrich_df.empty and {"signed", "up", "down"}.issubset(set(enrich_df["mode"].unique())):
        ct_z = (
            cs_signed.groupby(["cell_type", "transition"])["zscore"]
            .median().rename("z_signed").reset_index()
        )
        # up + down medians per cell_type x transition
        cs_split = cell_scores[cell_scores["mode"].isin(["up", "down"])].merge(
            obs_meta[["cell_id", "cell_type"]], on="cell_id", how="left"
        )
        ct_split = (
            cs_split.groupby(["cell_type", "transition", "mode"])["zscore"]
            .median().reset_index()
        )
        up = ct_split[ct_split["mode"] == "up"].rename(columns={"zscore": "z_up"}).drop(columns="mode")
        dn = ct_split[ct_split["mode"] == "down"].rename(columns={"zscore": "z_down"}).drop(columns="mode")
        merged = ct_z.merge(up, on=["cell_type", "transition"], how="outer") \
                     .merge(dn, on=["cell_type", "transition"], how="outer")
        merged["z_split"] = merged["z_up"].fillna(0) - merged["z_down"].fillna(0)
        merged["abs_diff"] = (merged["z_signed"] - merged["z_split"]).abs()
        merged["discordant_sign"] = (
            (np.sign(merged["z_signed"].fillna(0)) != np.sign(merged["z_split"].fillna(0)))
            & (merged["z_signed"].abs() > 0.1)
            & (merged["z_split"].abs() > 0.1)
        )
        merged.to_csv(CONCORD_PATH, index=False)
        log(f"  wrote {CONCORD_PATH};  median |signed - split|={merged['abs_diff'].median():.3f}")
    else:
        log("  cannot build concordance (modes missing)")

    log("STAGE 3 DONE")


if __name__ == "__main__":
    main()
