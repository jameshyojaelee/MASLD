#!/usr/bin/env python3
"""B3 step 2: skill, promoter and TSS-distance associations of predicted sequence sensitivity.

Every statistic is recomputed end to end inside each 1-Mb-block bootstrap draw through the same code
path that produced the observed value. Blocks are the resampling unit. Standardisation happens inside
the draw, never once on the observed sample.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

SEED = 20260914
B = 10000
STATS = ["mean_abs", "p95_abs", "conc_top5"]
CONTRASTS = ["T1_skill_Q5_minus_Q1", "T2_skill_slope_std", "T3_promoter_std",
             "T4_log10_dist_tss_std", "T5_skill_slope_std_adjusted"]
PRIMARY_GROUP = {"ATAC": "adult_liver", "DNASE": "adult_liver", "CHIP_HISTONE": "adult_liver",
                 "CHIP_TF": "adult_liver", "CAGE": "adult_liver", "AVI_SCORE": "avi_all"}


def z(v):
    s = v.std()
    return (v - v.mean()) / s if s > 0 else np.zeros_like(v)


def contrasts(Y, skill, prom, ldist, gc, q1, q5):
    """Y: (m, K) native statistic columns. Returns (5, K) contrast estimates."""
    Z = (Y - Y.mean(0)) / np.where(Y.std(0) > 0, Y.std(0), np.nan)
    x1, x3, x4 = z(skill), z(ldist), z(gc)
    out = np.empty((5, Y.shape[1]), float)
    out[0] = Y[q5].mean(0) - Y[q1].mean(0)                       # native units
    out[1] = (Z * x1[:, None]).mean(0)
    out[2] = Z[prom == 1].mean(0) - Z[prom == 0].mean(0)
    out[3] = (Z * x3[:, None]).mean(0)
    D = np.column_stack([np.ones(len(x1)), x1, prom.astype(float), x3, x4])
    bad = ~np.isfinite(Z).all(0)          # a column with any missing region is not fit jointly here
    beta, *_ = np.linalg.lstsq(D, np.nan_to_num(Z), rcond=None)
    out[4] = np.where(bad, np.nan, beta[1])
    return out


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    idx = np.flatnonzero(ok)
    o = idx[np.argsort(p[idx], kind="mergesort")]
    n = len(o)
    adj = p[o] * n / np.arange(1, n + 1)
    q[o] = np.minimum.accumulate(adj[::-1])[::-1].clip(max=1.0)
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out, tab = args.out, os.path.join(args.out, "tables")

    sm = pd.read_csv(os.path.join(tab, "region_summaries.tsv.gz"), sep="\t")
    samp = pd.read_csv(os.path.join(out, "prespec", "sample_regions.tsv"), sep="\t")
    if (sm.group == "ERROR").any():
        sm.loc[sm.group == "ERROR"].to_csv(os.path.join(tab, "region_read_errors.tsv"), sep="\t", index=False)
        sm = sm[sm.group != "ERROR"]
    print("[load] summaries", sm.shape, "sample", samp.shape, flush=True)

    # one wide region x (scorer,group,stat,scale) matrix
    long = sm.melt(id_vars=["region_key", "scorer", "group"],
                   value_vars=STATS + [s + "_q" for s in STATS],
                   var_name="stat", value_name="value")
    long["scale"] = np.where(long.stat.str.endswith("_q"), "quantile", "native")
    long["stat"] = long.stat.str.replace("_q$", "", regex=True)
    long["col"] = long.scorer + "|" + long.group + "|" + long.stat + "|" + long.scale
    wide = long.pivot_table(index="region_key", columns="col", values="value", aggfunc="first")
    dat = samp.set_index("region_key").loc[wide.index]
    cols = list(wide.columns)
    M = wide.to_numpy(float)
    nan_by_col = pd.Series(np.isnan(M).sum(0), index=cols)
    nan_by_col[nan_by_col > 0].to_csv(os.path.join(tab, "columns_with_missing_regions.tsv"), sep="\t")
    print("[load] columns with any missing region:", int((nan_by_col > 0).sum()), "of", len(cols), flush=True)
    # every contrast is a complete-case statistic; a channel missing any sampled region is dropped
    # whole and named in columns_with_missing_regions.tsv rather than partially imputed.
    keep = nan_by_col.to_numpy() == 0
    if not keep.all():
        print("[load] dropping channels with missing regions:",
              [c for c, k in zip(cols, keep) if not k], flush=True)
        M = M[:, keep]
        cols = [c for c, k in zip(cols, keep) if k]

    skill = dat.skill_shipped_mean.to_numpy(float)
    prom = dat.promoter.to_numpy(int)
    ldist = dat.log10_dist_tss.to_numpy(float)
    gc = dat.gc.to_numpy(float)
    quint = dat.skill_quintile.to_numpy(int)
    blocks = dat.mb_block.to_numpy()
    if not np.isfinite(ldist).all() or not np.isfinite(gc).all():
        raise SystemExit("[load] non-finite covariate in the sample; refusing to proceed")

    ub, binv = np.unique(blocks, return_inverse=True)
    by_block = [np.flatnonzero(binv == i) for i in range(len(ub))]
    print(f"[boot] {len(ub)} 1-Mb blocks over {len(dat)} regions "
          f"(median {np.median([len(b) for b in by_block]):.0f} regions/block)", flush=True)

    obs = contrasts(M, skill, prom, ldist, gc, quint == 1, quint == 5)

    rng = np.random.default_rng(SEED)
    draws = np.empty((B, 5, len(cols)), np.float32)
    t0 = time.time()
    for b in range(B):
        pick = rng.integers(0, len(ub), len(ub))
        idx = np.concatenate([by_block[i] for i in pick])
        draws[b] = contrasts(M[idx], skill[idx], prom[idx], ldist[idx], gc[idx],
                             quint[idx] == 1, quint[idx] == 5)
        if (b + 1) % 1000 == 0:
            print(f"[boot] {b+1}/{B} draws, {time.time()-t0:.0f}s", flush=True)

    rows = []
    for ci, c in enumerate(cols):
        scorer, group, stat, scale = c.split("|")
        for ti, tname in enumerate(CONTRASTS):
            d = draws[:, ti, ci].astype(float)
            d = d[np.isfinite(d)]
            if len(d) < B * 0.9 or not np.isfinite(obs[ti, ci]):
                rows.append({"scorer": scorer, "group": group, "stat": stat, "scale": scale,
                             "contrast": tname, "observed": obs[ti, ci], "ci_lo": np.nan,
                             "ci_hi": np.nan, "p_boot": np.nan, "n_draws": len(d)})
                continue
            lo, hi = np.percentile(d, [2.5, 97.5])
            p = min(1.0, 2 * min((d <= 0).sum() + 1, (d >= 0).sum() + 1) / (len(d) + 1))
            rows.append({"scorer": scorer, "group": group, "stat": stat, "scale": scale,
                         "contrast": tname, "observed": float(obs[ti, ci]), "ci_lo": float(lo),
                         "ci_hi": float(hi), "p_boot": float(p), "n_draws": len(d)})
    res = pd.DataFrame(rows)
    res["is_primary_group"] = [PRIMARY_GROUP[s] == g for s, g in zip(res.scorer, res.group)]
    res["family"] = np.where(res.is_primary_group, "primary_" + res.scale, "hepg2_" + res.scale)
    res["q_bh"] = np.nan
    for fam, sub in res.groupby("family"):
        res.loc[sub.index, "q_bh"] = bh(sub.p_boot.to_numpy())
    res["p_floor"] = 2.0 / (B + 1)
    res = res.sort_values(["family", "scorer", "stat", "contrast"])
    res.to_csv(os.path.join(tab, "b3_contrasts.tsv"), sep="\t", index=False, float_format="%.6g")
    print("[boot] wrote b3_contrasts.tsv:", res.shape, flush=True)

    # descriptive quintile profile, native scale, primary groups
    desc = []
    for c in cols:
        scorer, group, stat, scale = c.split("|")
        if scale != "native" or PRIMARY_GROUP[scorer] != group:
            continue
        v = M[:, cols.index(c)]
        for q in range(1, 6):
            m = quint == q
            desc.append({"scorer": scorer, "group": group, "stat": stat, "quintile": q,
                         "n": int(m.sum()), "mean": float(np.nanmean(v[m])),
                         "sd": float(np.nanstd(v[m])), "median": float(np.nanmedian(v[m]))})
    pd.DataFrame(desc).to_csv(os.path.join(tab, "b3_quintile_profile.tsv"), sep="\t",
                              index=False, float_format="%.6g")

    with open(os.path.join(tab, "b3_boot_meta.json"), "w") as fh:
        json.dump({"B": B, "seed": SEED, "n_regions": int(len(dat)), "n_blocks": int(len(ub)),
                   "p_floor": 2.0 / (B + 1), "resampling_unit": "1-Mb block (mb_block)",
                   "standardisation": "recomputed inside every draw",
                   "families": sorted(res.family.unique().tolist())}, fh, indent=1)
    print("[done]", flush=True)


if __name__ == "__main__":
    sys.exit(main())
