#!/usr/bin/env python
"""
496_dataset_confound.py — Test whether cNMF programs are dataset-driven artifacts.

Approach: for each program, decompose donor-score variance into between-dataset vs
within-dataset. Programs where >70% of variance is between-dataset are likely batch artifacts.
Also fit ANOVA: program_score ~ dataset + disease_stage to see which contributes more.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
OUT = MCP / "reviewer_defense"


def main():
    donor_meta = pd.read_csv(MCP / "inputs/donor_metadata.tsv", sep="\t").set_index("sample")
    scores = pd.read_csv(MCP / "integration/program_donor_scores_wide.tsv.gz", sep="\t", index_col=0)
    common = donor_meta.index.intersection(scores.index)
    donor_meta = donor_meta.loc[common]
    scores = scores.loc[common]

    rows = []
    for prog in scores.columns:
        y = pd.to_numeric(scores[prog], errors="coerce")
        ds = donor_meta["dataset"]
        stage = pd.to_numeric(donor_meta["disease_stage_numeric"], errors="coerce")
        # ANOVA decomposition
        valid = y.notna() & ds.notna()
        y0 = y[valid]; ds0 = ds[valid]
        # Between-dataset SS / total SS
        overall_mean = y0.mean()
        tss = ((y0 - overall_mean) ** 2).sum()
        ds_means = y0.groupby(ds0, observed=True).mean()
        ds_sizes = y0.groupby(ds0, observed=True).size()
        bss = (ds_sizes * (ds_means - overall_mean) ** 2).sum()
        dataset_frac_var = float(bss / tss) if tss > 0 else np.nan

        # F-test: dataset as categorical vs null
        try:
            ds_dummies = pd.get_dummies(ds0, drop_first=True)
            X_ds = sm.add_constant(ds_dummies.astype(float))
            m_ds = sm.OLS(y0.astype(float).values, X_ds.values).fit()
            f_ds_p = float(m_ds.f_pvalue)
            r2_ds = float(m_ds.rsquared)
        except Exception:
            f_ds_p = np.nan; r2_ds = np.nan

        # Stage as continuous (after dataset adjustment)
        valid_s = valid & stage.notna()
        if valid_s.sum() > 30:
            y1 = y[valid_s].astype(float); s1 = stage[valid_s].astype(float); d1 = ds[valid_s]
            try:
                dsd = pd.get_dummies(d1, drop_first=True).astype(float)
                X_full = sm.add_constant(pd.concat([pd.Series(s1.values, name="stage"), dsd.reset_index(drop=True)], axis=1))
                m_full = sm.OLS(y1.values, X_full.values).fit()
                r2_full = float(m_full.rsquared)
                # Dataset-only fit on stage-valid subset
                X_dsonly = sm.add_constant(dsd.astype(float).values)
                m_dsonly = sm.OLS(y1.values, X_dsonly).fit()
                r2_dsonly = float(m_dsonly.rsquared)
                stage_added_var = r2_full - r2_dsonly
                # partial stage effect
                stage_beta = float(m_full.params[1])
                stage_p = float(m_full.pvalues[1])
            except Exception as e:
                r2_full = np.nan; stage_added_var = np.nan; stage_beta = np.nan; stage_p = np.nan
        else:
            r2_full = stage_added_var = stage_beta = stage_p = np.nan

        rows.append({
            "program": prog,
            "n_donors": int(valid.sum()),
            "dataset_frac_var": dataset_frac_var,
            "dataset_f_p": f_ds_p,
            "r2_dataset_only": r2_ds,
            "r2_with_stage": r2_full,
            "stage_added_var": stage_added_var,
            "stage_beta_adj_ds": stage_beta,
            "stage_p_adj_ds": stage_p,
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "dataset_confound.tsv", sep="\t", index=False)
    n_batch = int((df["dataset_frac_var"] > 0.7).sum())
    n_stage = int((df["stage_p_adj_ds"] < 0.05).sum())
    print(f"[496] Programs with >70% variance from dataset: {n_batch}/{len(df)}")
    print(f"[496] Programs with stage effect (after dataset adjustment) p<0.05: {n_stage}/{len(df)}")
    print(df.sort_values("dataset_frac_var", ascending=False).to_string())


if __name__ == "__main__":
    main()
