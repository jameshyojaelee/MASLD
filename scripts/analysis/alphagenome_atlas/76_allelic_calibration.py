#!/usr/bin/env python3
"""Step 76 (P3): is the Atlas ATAC quantile CALIBRATED against measured allelic direction?

Sign concordance pooled over all sites blends strata the model is confident about with strata it is not.
The question that matters for using these scores is whether the confident predictions are the right ones:
does agreement with the measurement rise with the predicted effect size? A model whose agreement is flat in
|quantile| is right on average but carries no usable per-variant information.

Reported per cohort: equal-count quantile bins of |predicted quantile| with concordance in each, a Spearman
trend across sites, and a 1-Mb block bootstrap of that trend (one site per block per draw, because sites at
one locus are in LD). Agreement is also tabulated against the number of het donors, which separates "the
model is uncertain" from "the measurement is noisy".

Outputs (tables/): allelic_calibration_bins.tsv, allelic_calibration.json
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
COHORTS = (("gse281367", "allelic_sites.tsv", ""), ("gse244832", "gse244832_allelic_sites.tsv", "gse244832_"))
N_BINS = 4
DRAWS = 2000
SEED = 76076


def quantile_bins(absq: np.ndarray, n_bins: int) -> np.ndarray:
    """Equal-count bins of |quantile|, ordered by magnitude; degenerate input collapses rather than raising."""
    a = np.asarray(absq, float)
    try:
        b = pd.qcut(a, q=n_bins, labels=False, duplicates="drop")
    except (ValueError, IndexError):
        return np.zeros(a.shape, int)
    return np.asarray(pd.Series(b).fillna(0), int)


def calibration_trend(absq: np.ndarray, agree: np.ndarray) -> float:
    """Spearman between |predicted quantile| and per-site agreement. Positive means calibrated."""
    a, g = np.asarray(absq, float), np.asarray(agree, float)
    ok = ~np.isnan(a) & ~np.isnan(g)
    # 3 is the smallest n for which a rank correlation is not degenerate; on real data (hundreds of sites)
    # this guard never binds, it only keeps tiny bootstrap draws from returning nonsense.
    if ok.sum() < 3 or len(np.unique(g[ok])) < 2 or len(np.unique(a[ok])) < 2:
        return float("nan")
    return float(spearmanr(a[ok], g[ok]).correlation)


def block_bootstrap_trend(absq: np.ndarray, agree: np.ndarray, blocks: np.ndarray,
                          draws: int = DRAWS, seed: int = SEED) -> dict:
    """Block bootstrap of the trend: one site per sampled 1-Mb block, so LD cannot inflate the interval."""
    rng = np.random.default_rng(seed)
    uniq = np.unique(blocks)
    members = {b: np.flatnonzero(blocks == b) for b in uniq}
    vals = []
    for _ in range(draws):
        picked = rng.choice(uniq, size=len(uniq), replace=True)
        take = np.array([rng.choice(members[b]) for b in picked])
        t = calibration_trend(np.asarray(absq)[take], np.asarray(agree)[take])
        if t == t:
            vals.append(t)
    if not vals:
        return {"mean": float("nan"), "ci95": [float("nan"), float("nan")], "n_blocks": int(len(uniq))}
    v = np.asarray(vals)
    return {"mean": float(v.mean()), "ci95": [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))],
            "n_blocks": int(len(uniq)), "n_draws": int(v.size)}


def main() -> None:
    out, rows = {}, []
    for cohort, fname, prefix in COHORTS:
        path = TABLES / fname
        if not path.exists():
            continue
        d = pd.read_csv(path, sep="\t")
        d = d[d.atac_quantile.notna() & (d.atac_quantile != 0)].copy()
        if len(d) < 20:
            out[cohort] = {"verdict": "cannot conclude: fewer than 20 scored sites"}
            continue
        d["agree"] = (np.sign(d.atac_quantile) == np.sign(d.mean_log2_alt_over_ref)).astype(int)
        d["absq"] = d.atac_quantile.abs()
        d["bin"] = quantile_bins(d.absq.to_numpy(), N_BINS)
        for b, g in d.groupby("bin"):
            rows.append({"cohort": cohort, "bin": int(b) + 1, "n": int(len(g)),
                         "median_abs_quantile": float(g.absq.median()),
                         "concordance": float(g.agree.mean()),
                         "marginal_expected": la.marginal_expected_concordance(
                             np.sign(g.atac_quantile), np.sign(g.mean_log2_alt_over_ref))})
        trend = calibration_trend(d.absq.to_numpy(), d.agree.to_numpy())
        boot = block_bootstrap_trend(d.absq.to_numpy(), d.agree.to_numpy(), d.block.to_numpy())
        by_donor = {str(k): {"n": int(len(g)), "concordance": float(g.agree.mean())}
                    for k, g in d.groupby(pd.cut(d.n_het_donors, [2, 3, 4, 100], labels=["3", "4", "5+"]), observed=True)}
        first = d[d.bin == d.bin.min()].agree.mean()
        last = d[d.bin == d.bin.max()].agree.mean()
        out[cohort] = {"n_sites": int(len(d)), "overall_concordance": float(d.agree.mean()),
                       "trend_spearman": trend, "trend_block_bootstrap": boot,
                       "weakest_bin_concordance": float(first), "strongest_bin_concordance": float(last),
                       "concordance_by_n_het_donors": by_donor,
                       "verdict": ("calibrated: agreement rises with predicted effect size and the block-bootstrap "
                                   "interval excludes zero" if boot["ci95"][0] > 0 else
                                   "not established: the block-bootstrap interval for the trend includes zero")}
    la.write_tsv_once(TABLES / "allelic_calibration_bins.tsv", rows, sorted({k for r in rows for k in r}))
    json.dump(out, (TABLES / "allelic_calibration.json").open("w"), indent=1, default=float)
    for c, v in out.items():
        if "trend_spearman" in v:
            la.log(f"{c}: trend {v['trend_spearman']:+.3f} boot {v['trend_block_bootstrap']['mean']:+.3f} "
                   f"CI{[round(x, 3) for x in v['trend_block_bootstrap']['ci95']]}; "
                   f"weakest {v['weakest_bin_concordance']:.3f} -> strongest {v['strongest_bin_concordance']:.3f}")


if __name__ == "__main__":
    main()
