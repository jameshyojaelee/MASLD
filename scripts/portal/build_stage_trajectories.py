#!/usr/bin/env python3
"""build_stage_trajectories.py

Aggregate per-cohort stage-vs-F0 DE (from per_cohort_contrasts.R) into a single
weighted-meta-analysis trajectory file. Also emit an NAS trajectory from the
existing one-vs-rest Dream results.

Outputs go under:
  RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier/
    - fibrosis_stage_trajectory_vs_F0.csv       (new; only if per-cohort files exist)
    - fibrosis_stage_trajectory_vs_rest.csv     (shim over one_vs_rest_fibrosis_dream.csv)
    - fibrosis_stage_trajectory_vs_healthy.csv  (shim over fibrosis_vs_healthy_dream.csv, if present)
    - nas_trajectory.csv                        (from one_vs_rest_nas_dream.csv)
    - nas_trajectory_vs_nas0.csv                (shim over nas_vs_nas0_dream.csv, if present)
    - nas_trajectory_vs_healthy.csv             (shim over nas_vs_healthy_dream.csv, if present)

Usage:
  python scripts/portal/build_stage_trajectories.py
"""
from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
STAGING_DIR = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier"
)
PER_STUDY_CONTRASTS_DIR = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study_contrasts"
)

F_STAGES = ["F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0"]


def _meta_combine(group: pd.DataFrame) -> pd.Series:
    """Weighted meta-analysis across cohorts: inverse-variance weights from t stats.

    Fallback: if t/P.Value yield no SE, use simple mean logFC + Fisher on P.
    """
    g = group.dropna(subset=["logFC", "P.Value"])
    if g.empty:
        return pd.Series({"logFC": np.nan, "padj": np.nan, "n_cohorts": 0})
    # Derive SE from t: SE = logFC / t (guard t==0)
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.where(g["t"].abs() > 1e-9, g["logFC"] / g["t"], np.nan)
    se = np.abs(se)
    valid = np.isfinite(se) & (se > 0)
    if valid.sum() >= 1:
        w = 1.0 / (se[valid] ** 2)
        lfc = float(np.sum(w * g["logFC"].values[valid]) / np.sum(w))
        # Combined P via Stouffer using z = lfc / sqrt(1/sum(w))
        z = lfc / math.sqrt(1.0 / np.sum(w))
        from math import erfc, sqrt

        p = erfc(abs(z) / sqrt(2))
    else:
        lfc = float(g["logFC"].mean())
        # Fisher combined p
        chi = -2.0 * np.log(np.clip(g["P.Value"].values, 1e-300, 1.0)).sum()
        from scipy.stats import chi2  # type: ignore

        p = float(chi2.sf(chi, df=2 * len(g)))
    return pd.Series({"logFC": lfc, "P": p, "n_cohorts": int(len(g))})


def build_vs_control() -> Path | None:
    """Aggregate F{k}_vs_F0 per-cohort CSVs, if present."""
    any_present = any(
        (PER_STUDY_CONTRASTS_DIR / s).is_dir() and
        any((PER_STUDY_CONTRASTS_DIR / s).glob("*_de.csv"))
        for s in F_STAGES
    )
    if not any_present:
        print("[vs_F0] no per-cohort fibrosis contrast outputs yet — skipping.")
        return None

    rows = []
    for stage in F_STAGES:
        sdir = PER_STUDY_CONTRASTS_DIR / stage
        if not sdir.is_dir():
            continue
        for fp in sorted(sdir.glob("*_de.csv")):
            cohort = fp.name.replace("_de.csv", "")
            df = pd.read_csv(fp, usecols=["gene", "logFC", "t", "P.Value", "padj"])
            df["stage"] = stage
            df["cohort"] = cohort
            rows.append(df)
    if not rows:
        print("[vs_F0] no rows aggregated.")
        return None

    big = pd.concat(rows, ignore_index=True)
    print(f"[vs_F0] {len(big):,} per-cohort rows across {big['cohort'].nunique()} cohorts")

    meta = (
        big.groupby(["gene", "stage"], sort=False)
        .apply(_meta_combine)
        .reset_index()
    )
    # BH-adjust within stage
    from statsmodels.stats.multitest import multipletests  # type: ignore

    meta["padj"] = np.nan
    for s in meta["stage"].unique():
        m = meta["stage"] == s
        pv = meta.loc[m, "P"].fillna(1.0).values
        _, padj, *_ = multipletests(pv, method="fdr_bh")
        meta.loc[m, "padj"] = padj
    meta = meta.rename(columns={"logFC": "logFC"})[
        ["gene", "stage", "logFC", "padj", "n_cohorts"]
    ]
    out = STAGING_DIR / "fibrosis_stage_trajectory_vs_F0.csv"
    meta.to_csv(out, index=False)
    print(f"[vs_F0] wrote {len(meta):,} rows -> {out}")
    return out


def build_vs_rest_shim() -> Path:
    src = STAGING_DIR / "one_vs_rest_fibrosis_dream.csv"
    out = STAGING_DIR / "fibrosis_stage_trajectory_vs_rest.csv"
    df = pd.read_csv(src, usecols=["gene", "stage_label", "logFC", "padj"])
    df = df.rename(columns={"stage_label": "stage"})
    df["n_cohorts"] = None
    df.to_csv(out, index=False)
    print(f"[vs_rest] shimmed {len(df):,} rows -> {out}")
    return out


def build_nas() -> Path:
    src = STAGING_DIR / "one_vs_rest_nas_dream.csv"
    out = STAGING_DIR / "nas_trajectory.csv"
    df = pd.read_csv(src, usecols=["gene", "stage_label", "logFC", "padj"])
    df = df.rename(columns={"stage_label": "nas_score"})
    df["n_cohorts"] = None
    df.to_csv(out, index=False)
    print(f"[nas] wrote {len(df):,} rows -> {out}")
    return out


def build_fibrosis_vs_healthy_shim() -> Path | None:
    """Shim fibrosis_vs_healthy_dream.csv into the stage-trajectory schema.

    Schema: gene, stage, logFC, padj, n_cohorts
    If the source is not present yet (dream job still running), skip.
    """
    src = STAGING_DIR / "fibrosis_vs_healthy_dream.csv"
    out = STAGING_DIR / "fibrosis_stage_trajectory_vs_healthy.csv"
    if not src.exists():
        print("[fib_vs_healthy] source not present yet — skipping")
        return None
    df = pd.read_csv(src, usecols=["gene", "stage_label", "logFC", "padj"])
    df = df.rename(columns={"stage_label": "stage"})
    df["n_cohorts"] = None
    df.to_csv(out, index=False)
    print(f"[fib_vs_healthy] shimmed {len(df):,} rows -> {out}")
    return out


def build_nas_baseline_shim(src_name: str, out_name: str, tag: str) -> Path | None:
    """Shim NAS-vs-fixed-baseline dream CSV into the nas_trajectory schema.

    Schema: gene, nas_score, logFC, padj, n_cohorts
    If the source file is not present yet (dream job still running), skip.
    """
    src = STAGING_DIR / src_name
    if not src.exists():
        print(f"[{tag}] source not present yet ({src_name}) -- skipping")
        return None
    df = pd.read_csv(src, usecols=["gene", "stage_label", "logFC", "padj"])
    df = df.rename(columns={"stage_label": "nas_score"})
    df["n_cohorts"] = None
    out = STAGING_DIR / out_name
    df.to_csv(out, index=False)
    print(f"[{tag}] shimmed {len(df):,} rows -> {out}")
    return out


def main() -> None:
    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    build_vs_rest_shim()
    build_nas()
    build_vs_control()
    build_fibrosis_vs_healthy_shim()
    build_nas_baseline_shim(
        "nas_vs_nas0_dream.csv", "nas_trajectory_vs_nas0.csv", "nas_vs_nas0"
    )
    build_nas_baseline_shim(
        "nas_vs_healthy_dream.csv", "nas_trajectory_vs_healthy.csv", "nas_vs_healthy"
    )


if __name__ == "__main__":
    main()
