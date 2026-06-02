#!/usr/bin/env python
"""
figS_scdrs_null_aggregate_25pct.py
==================================
Aggregate the 25%-subsample variant of the hep-matched null bias check.
Same logic as figS_scdrs_null_aggregate.py but reads from the *_25pct
output dirs (real and null both scored on the same 225K-cell subsample, so
the observed-vs-null comparison is paired).

Outputs:
  Analysis/SingleCell/results_gpu_v2/disease_signatures/
    scdrs_gwas_hep_null_summary_25pct.csv
    scdrs_gwas_hep_null_z_long_25pct.csv
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
SIG_DIR  = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
REAL_DIR = SIG_DIR / "scdrs_gwas_trait_outputs_25pct"
NULL_DIR = SIG_DIR / "scdrs_gwas_hep_null_outputs_25pct"
OUT_SUMMARY = SIG_DIR / "scdrs_gwas_hep_null_summary_25pct.csv"
OUT_LONG    = SIG_DIR / "scdrs_gwas_hep_null_z_long_25pct.csv"

TARGET_PHENOS = [
    "gwas_coloc_nafld_specific_polyfun",
    "gwas_coloc_liver_enzymes_polyfun",
    "gwas_coloc_pdff_polyfun",
    "gwas_coloc_cirrhosis_hcc_polyfun",
]


def log(msg: str) -> None:
    print(f"[null_aggregate_25pct] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


def read_enrich(path: Path, trait: str, is_real: bool, real_pheno: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    if "cell_type" not in df.columns or "assoc_mcz" not in df.columns:
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

    log(f"reading null panels from {NULL_DIR}")
    null_parts = []
    for pheno in TARGET_PHENOS:
        for i in range(50):
            null_name = f"{pheno}__null{i:03d}"
            f = NULL_DIR / f"{null_name}_celltype_enrichment.csv"
            d = read_enrich(f, trait=null_name, is_real=False, real_pheno=pheno)
            null_parts.append(d)
    null_df = pd.concat(null_parts, ignore_index=True)
    n_null_traits = null_df["trait"].nunique()
    log(f"  combined null: {len(null_df)} rows  ({n_null_traits}/200 traits collected)")

    long_df = pd.concat([real_df, null_df], ignore_index=True)
    long_df.to_csv(OUT_LONG, index=False)
    log(f"wrote {OUT_LONG}  shape={long_df.shape}")

    log("summarizing observed vs null per (CT × real_pheno)")
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

    log("Hepatocyte rows:")
    log(summary[summary["cell_type"] == "Hepatocytes"].to_string(index=False))


if __name__ == "__main__":
    main()
