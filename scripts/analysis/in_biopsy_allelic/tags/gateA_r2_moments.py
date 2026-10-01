#!/usr/bin/env python3
"""Model A-final: covariate moments for the conditional lambda (spec v1.2 section 3 items 6-8).

For each development F cohort, the design X is the spec z columns in order (b, d, FFPE, sex, age/10,
axis_EUR_AFR, axis_EUR_EAS; v1.4 drops c), mean-filled within cohort, plus a missing indicator for every
column that is partly missing within the cohort (appended in base-column order), then centred
within cohort; S is centred the same way. Written (aggregate only; no per-person values):
  r2_moments_observed.json   per cohort: n, retained column names, effective rank p
                             (numpy.linalg.matrix_rank of centred X), XtX, XtS, StS, raw R2, adjusted
                             R2 = 1 - (1 - R2)(n - 1)/(n - p - 1), floored at 0
  p2_resample_development.tsv.gz  the fixed P2 resamples of F (spec 5.2): draw b = 0..B-1,
                             rng = numpy default_rng(SeedSequence([20260923, 2, b])); strata = cohorts
                             in sorted order; within a stratum, clusters sorted by cluster_id and
                             drawn with rng.integers(0, n_clusters, n_clusters); columns draw, stratum,
                             slot, cluster_id (pseudonymous individual IDs)
  r2_by_resample.tsv.gz      per draw x cohort: n (with multiplicity), p, raw R2, adjusted R2, with the
                             mean fill, indicators and centring recomputed on the resample using
                             multiplicities
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

Z = ["z_b_log_e", "z_d_dup", "z_ffpe", "z_sex_female", "z_age10", "axis_EUR_AFR", "axis_EUR_EAS"]  # spec v1.4


def design(g):
    """Centred X (with indicators) and centred S for one cohort sample (rows may repeat)."""
    cols, names = [], []
    miss = []
    for c in Z:
        v = g[c].to_numpy(float)
        m = np.isnan(v)
        if m.all():
            continue
        v = np.where(m, np.nanmean(v), v)
        cols.append(v - v.mean()); names.append(c)
        if m.any():
            miss.append((c, m.astype(float)))
    for c, m in miss:
        cols.append(m - m.mean()); names.append(f"{c}_missing")
    X = np.column_stack(cols) if cols else np.zeros((len(g), 0))
    S = g["S"].to_numpy(float)
    return X, S - S.mean(), names


def r2_stats(X, S):
    n = len(S)
    p = int(np.linalg.matrix_rank(X)) if X.shape[1] else 0
    sts = float(S @ S)
    if sts == 0 or p == 0:
        return n, p, float("nan"), float("nan")
    beta, *_ = np.linalg.lstsq(X, S, rcond=None)
    raw = float(1 - ((S - X @ beta) @ (S - X @ beta)) / sts)
    adj = max(0.0, 1 - (1 - raw) * (n - 1) / (n - p - 1)) if n - p - 1 > 0 else float("nan")
    return n, p, raw, adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--individuals", required=True, help="Gate A design_individuals.tsv")
    ap.add_argument("--crosswalk", required=True)
    ap.add_argument("--ancestry", required=True, help="restricted library_ancestry_pcs.tsv (axes read, never written)")
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    ind = pd.read_csv(a.individuals, sep="\t")
    ind = ind[ind["set"] == "development"].copy()
    xw = pd.read_csv(a.crosswalk, sep="\t", dtype=str)
    xw = xw[(xw["unit_library"] == "True") & (xw["run_role"] == "development")][["run", "individual_id"]]
    anc = pd.read_csv(a.ancestry, sep="\t", usecols=["run", "axis_EUR_AFR", "axis_EUR_EAS"])
    ind = ind.merge(xw, on="individual_id", how="left").merge(anc, on="run", how="left").drop(columns=["run"])

    obs = {}
    for c, g in ind.groupby("cohort"):
        X, S, names = design(g)
        n, p, raw, adj = r2_stats(X, S)
        obs[c] = {"n": n, "columns": names, "rank_p": p, "XtX": (X.T @ X).tolist(), "XtS": (X.T @ S).tolist(),
                  "StS": float(S @ S), "r2_raw": raw, "r2_adjusted": adj}
    (out / "r2_moments_observed.json").write_text(json.dumps(obs, indent=1))

    strata = {c: sorted(g["cluster_id"].unique()) for c, g in ind.groupby("cohort")}
    members = ind.groupby("cluster_id")
    draws, stats = [], []
    for b in range(a.draws):
        rng = np.random.default_rng(np.random.SeedSequence([20260923, 2, b]))
        for c in sorted(strata):
            cl = strata[c]
            pick = rng.integers(0, len(cl), len(cl))
            chosen = [cl[j] for j in pick]
            draws.append(pd.DataFrame({"draw": b, "stratum": c, "slot": np.arange(len(chosen)), "cluster_id": chosen}))
            g = pd.concat([members.get_group(k) for k in chosen], ignore_index=True)
            g = g[g["cohort"] == c]
            n, p, raw, adj = r2_stats(*design(g)[:2])
            stats.append({"draw": b, "cohort": c, "n": n, "rank_p": p, "r2_raw": raw, "r2_adjusted": adj})
    pd.concat(draws).to_csv(out / "p2_resample_development.tsv.gz", sep="\t", index=False, compression="gzip")
    st = pd.DataFrame(stats)
    st.to_csv(out / "r2_by_resample.tsv.gz", sep="\t", index=False, compression="gzip")
    summary = {c: {"n": v["n"], "rank_p": v["rank_p"], "r2_raw": v["r2_raw"], "r2_adjusted": v["r2_adjusted"],
                   "resample_r2_adjusted_p5_p50_p95": st.loc[st["cohort"] == c, "r2_adjusted"].quantile([0.05, 0.5, 0.95]).round(4).tolist()}
               for c, v in obs.items()}
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
