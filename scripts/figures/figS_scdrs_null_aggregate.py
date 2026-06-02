#!/usr/bin/env python
"""
figS_scdrs_null_aggregate.py
============================
Aggregate the 200 null scDRS group-analysis outputs (50 nulls × 4 phenotypes,
hepatocyte-expression-matched) with the 4 real panel outputs to build the
bias-controlled comparison.

For each (real_phenotype × cell_type):
  - z_observed = real panel scDRS group z-score
  - z_null_dist = z-scores from the 50 hep-expression-matched nulls
  - empirical_p = (#{null_z >= z_observed} + 1) / 51   [one-sided]
  - delta_z = z_observed − mean(z_null_dist)
  - z_null_mean, z_null_sd, q025, q975 = summary stats for the null distribution

Outputs:
  Analysis/SingleCell/results_gpu_v2/disease_signatures/
    scdrs_gwas_hep_null_summary.csv         (one row per CT × phenotype)
    scdrs_gwas_hep_null_z_long.csv          (long: trait, cell_type, z, is_real)
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
SIG_DIR = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
REAL_DIR = SIG_DIR / "scdrs_gwas_trait_outputs"
NULL_DIR = SIG_DIR / "scdrs_gwas_hep_null_outputs"
OUT_SUMMARY = SIG_DIR / "scdrs_gwas_hep_null_summary.csv"
OUT_LONG    = SIG_DIR / "scdrs_gwas_hep_null_z_long.csv"

TARGET_PHENOS = [
    "gwas_coloc_nafld_specific_polyfun",
    "gwas_coloc_liver_enzymes_polyfun",
    "gwas_coloc_pdff_polyfun",
    "gwas_coloc_cirrhosis_hcc_polyfun",
]


def log(msg: str) -> None:
    print(f"[null_aggregate] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def read_enrich(path: Path, trait: str, is_real: bool, real_pheno: str) -> pd.DataFrame:
    if not path.exists():
        log(f"  WARNING: missing {path}")
        return pd.DataFrame()
    df = pd.read_csv(path)
    if "cell_type" not in df.columns or "assoc_mcz" not in df.columns:
        log(f"  WARNING: {path} missing expected columns")
        return pd.DataFrame()
    df = df[["cell_type", "assoc_mcz", "assoc_mcp"]].copy()
    df["trait"] = trait
    df["is_real"] = is_real
    df["real_pheno"] = real_pheno
    return df


def main() -> None:
    log(f"reading real panels from {REAL_DIR}")
    real_parts = []
    for pheno in TARGET_PHENOS:
        f = REAL_DIR / f"{pheno}_celltype_enrichment.csv"
        d = read_enrich(f, trait=pheno, is_real=True, real_pheno=pheno)
        log(f"  {pheno}: {len(d)} CT rows")
        real_parts.append(d)
    real_df = pd.concat(real_parts, ignore_index=True)
    log(f"  combined real: {len(real_df)} rows")

    log(f"reading null panels from {NULL_DIR}")
    null_parts = []
    for pheno in TARGET_PHENOS:
        for i in range(50):
            null_name = f"{pheno}__null{i:03d}"
            f = NULL_DIR / f"{null_name}_celltype_enrichment.csv"
            d = read_enrich(f, trait=null_name, is_real=False, real_pheno=pheno)
            null_parts.append(d)
    null_df = pd.concat(null_parts, ignore_index=True)
    log(f"  combined null: {len(null_df)} rows")
    log(f"  null traits collected: {null_df['trait'].nunique()}/200")
    log(f"  null CT rows per pheno:")
    log(null_df.groupby("real_pheno").size().to_string())

    long_df = pd.concat([real_df, null_df], ignore_index=True)
    long_df.to_csv(OUT_LONG, index=False)
    log(f"wrote {OUT_LONG}  shape={long_df.shape}")

    log("computing observed vs null summary per (CT × real_pheno)")
    rows = []
    for (pheno, ct), g in long_df.groupby(["real_pheno", "cell_type"]):
        real = g[g.is_real]
        null = g[~g.is_real]
        if len(real) == 0 or len(null) == 0:
            continue
        z_obs = float(real["assoc_mcz"].iloc[0])
        z_null = null["assoc_mcz"].dropna().values
        if len(z_null) == 0:
            continue
        empirical_p = (np.sum(z_null >= z_obs) + 1) / (len(z_null) + 1)
        delta_z = z_obs - float(np.mean(z_null))
        rows.append({
            "real_pheno": pheno,
            "cell_type": ct,
            "z_observed": z_obs,
            "n_nulls": len(z_null),
            "z_null_mean": float(np.mean(z_null)),
            "z_null_sd": float(np.std(z_null, ddof=1)) if len(z_null) > 1 else np.nan,
            "z_null_q025": float(np.quantile(z_null, 0.025)),
            "z_null_q975": float(np.quantile(z_null, 0.975)),
            "delta_z": delta_z,
            "empirical_p": empirical_p,
            "bias_controlled_sig": bool(empirical_p < 0.05),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_SUMMARY, index=False)
    log(f"wrote {OUT_SUMMARY}  shape={summary.shape}")

    log("Hepatocyte rows (highlight for the figure):")
    log(summary[summary["cell_type"] == "Hepatocytes"].to_string(index=False))


if __name__ == "__main__":
    main()
