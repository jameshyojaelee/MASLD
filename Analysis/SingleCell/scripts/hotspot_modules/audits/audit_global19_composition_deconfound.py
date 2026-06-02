#!/usr/bin/env python
"""Audit: does global__19 disease β survive after controlling for cholangiocyte
cell fraction? Tests the trans-differentiation interpretation vs the
cholangiocyte-expansion-only alternative.

Models compared:
  M1 unadjusted:   score ~ disease_stage_numeric + (1|dataset)
  M2 adjusted:     score ~ disease_stage_numeric + frac_Cholangiocytes + (1|dataset)
  M3 hep_proxy:    score ~ disease_stage_numeric + frac_Hepatocytes + (1|dataset)

Output: TSV row with β / p / retention for each model.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
HS = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
META = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
OUT = HS.parent / "hotspot_modules/audits" / "global19_deconfound.tsv"
OUT.parent.mkdir(parents=True, exist_ok=True)

ds = pd.read_csv(HS / "donor_scores_all.tsv", sep="\t")
ds = ds[ds.module == "global__19"].copy()
meta = pd.read_csv(META, sep="\t")
df = ds.merge(meta, on="sample", how="inner")
df = df.rename(columns={"frac_Cholangiocytes": "frac_chol",
                        "frac_Hepatocytes": "frac_hep"})
df = df.dropna(subset=["disease_stage_numeric", "frac_chol", "frac_hep", "score"])

def fit(formula, label):
    m = smf.mixedlm(formula, data=df, groups=df["dataset"]).fit()
    return {
        "model": label,
        "formula": formula,
        "beta_stage": float(m.params["disease_stage_numeric"]),
        "p_stage": float(m.pvalues["disease_stage_numeric"]),
        "beta_frac_chol": float(m.params.get("frac_chol", np.nan)),
        "p_frac_chol": float(m.pvalues.get("frac_chol", np.nan)),
        "beta_frac_hep": float(m.params.get("frac_hep", np.nan)),
        "p_frac_hep": float(m.pvalues.get("frac_hep", np.nan)),
        "n_donors": int(len(df)),
    }

rows = [
    fit("score ~ disease_stage_numeric", "M1_unadjusted"),
    fit("score ~ disease_stage_numeric + frac_chol", "M2_chol_adjusted"),
    fit("score ~ disease_stage_numeric + frac_hep", "M3_hep_adjusted"),
    fit("score ~ disease_stage_numeric + frac_chol + frac_hep", "M4_both_adjusted"),
]
result = pd.DataFrame(rows)
unadj = result.loc[0, "beta_stage"]
result["beta_stage_retained_pct"] = 100 * result["beta_stage"] / unadj
print(result.to_string(index=False))
result.to_csv(OUT, sep="\t", index=False)
print(f"\nWrote {OUT}")
