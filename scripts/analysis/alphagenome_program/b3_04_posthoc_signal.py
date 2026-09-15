#!/usr/bin/env python3
"""B3 post hoc (NOT prespecified): does the skill association survive measured signal strength?

Reason it is run: the producing deposit
GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables/region_covariate_associations.tsv
reports that `skill_shipped_mean` correlates 0.286 with `signal_mean` and 0.226 with `signal_sd`,
far more strongly than with promoter (-0.058) or log10_dist_tss (+0.053). The hashed B3
prespecification adjusts only for promoter, log10_dist_tss and gc, so the prespecified adjusted slope
cannot rule out that the skill association is measured H3K27ac signal strength. This file tests that
and is labelled post hoc everywhere it is reported.

Same loader, same 1-Mb blocks, same 10,000 draws, same seed as b3_02_analysis.py.
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

SEED, B = 20260914, 10000
STATS = ["mean_abs", "p95_abs", "conc_top5"]
PRIMARY_GROUP = {"ATAC": "adult_liver", "DNASE": "adult_liver", "CHIP_HISTONE": "adult_liver",
                 "CHIP_TF": "adult_liver", "CAGE": "adult_liver", "AVI_SCORE": "avi_all"}
MODELS = {
    "M0_unadjusted": [],
    "M1_prespecified_promoter_tss_gc": ["promoter", "log10_dist_tss", "gc"],
    "M2_posthoc_plus_signal": ["promoter", "log10_dist_tss", "gc", "signal_mean", "signal_sd"],
    "M3_posthoc_signal_only": ["signal_mean", "signal_sd"],
}


def z(v):
    s = v.std()
    return (v - v.mean()) / s if s > 0 else np.zeros_like(v)


def slopes(Y, cov, names):
    """Standardised partial slope of each Y column on skill, adjusting for `names`."""
    Z = (Y - Y.mean(0)) / np.where(Y.std(0) > 0, Y.std(0), np.nan)
    D = [np.ones(len(Y)), z(cov["skill_shipped_mean"])]
    for n in names:
        D.append(cov[n].astype(float) if n == "promoter" else z(cov[n]))
    D = np.column_stack(D)
    beta, *_ = np.linalg.lstsq(D, Z, rcond=None)
    return beta[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tab = os.path.join(a.out, "tables")
    sm = pd.read_csv(os.path.join(tab, "region_summaries.tsv.gz"), sep="\t")
    sm = sm[sm.group != "ERROR"]
    sm = sm[[PRIMARY_GROUP[s] == g for s, g in zip(sm.scorer, sm.group)]]
    samp = pd.read_csv(os.path.join(a.out, "prespec", "sample_regions.tsv"), sep="\t")

    long = sm.melt(id_vars=["region_key", "scorer"], value_vars=STATS, var_name="stat", value_name="v")
    long["col"] = long.scorer + "|" + long.stat
    wide = long.pivot_table(index="region_key", columns="col", values="v", aggfunc="first")
    dat = samp.set_index("region_key").loc[wide.index]
    cols, M = list(wide.columns), wide.to_numpy(float)
    if np.isnan(M).any():
        raise SystemExit("missing values in the post hoc matrix")
    cov = {c: dat[c].to_numpy(float) for c in
           ["skill_shipped_mean", "promoter", "log10_dist_tss", "gc", "signal_mean", "signal_sd"]}
    ub, binv = np.unique(dat.mb_block.to_numpy(), return_inverse=True)
    by_block = [np.flatnonzero(binv == i) for i in range(len(ub))]

    obs = {m: slopes(M, cov, v) for m, v in MODELS.items()}
    rng = np.random.default_rng(SEED)
    draws = {m: np.empty((B, len(cols)), np.float32) for m in MODELS}
    for b in range(B):
        idx = np.concatenate([by_block[i] for i in rng.integers(0, len(ub), len(ub))])
        cb = {k: v[idx] for k, v in cov.items()}
        for m, v in MODELS.items():
            draws[m][b] = slopes(M[idx], cb, v)

    rows = []
    for m in MODELS:
        for ci, c in enumerate(cols):
            d = draws[m][:, ci].astype(float)
            d = d[np.isfinite(d)]
            lo, hi = np.percentile(d, [2.5, 97.5])
            p = min(1.0, 2 * min((d <= 0).sum() + 1, (d >= 0).sum() + 1) / (len(d) + 1))
            rows.append({"model": m, "adjusted_for": ",".join(MODELS[m]) or "none",
                         "scorer": c.split("|")[0], "stat": c.split("|")[1],
                         "skill_slope_std": float(obs[m][ci]), "ci_lo": float(lo), "ci_hi": float(hi),
                         "p_boot": float(p), "prespecified": m in ("M0_unadjusted",
                                                                  "M1_prespecified_promoter_tss_gc")})
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(tab, "b3_posthoc_signal_adjustment.tsv"), sep="\t",
               index=False, float_format="%.6g")
    pd.set_option("display.width", 250)
    for stat in STATS:
        print(f"\n=== skill slope on {stat}, by adjustment set ===")
        print(res[res.stat == stat].pivot_table(index="scorer", columns="model",
                                                values="skill_slope_std").to_string())
    print("\nfull table -> tables/b3_posthoc_signal_adjustment.tsv")


if __name__ == "__main__":
    main()
