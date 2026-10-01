#!/usr/bin/env python3
"""C5: Model C sealed test (PRESPEC_C section 8, Amendments 3-5). Run ONCE.

Fits B0 and M1 exactly as in c4_ridge.py on all 179 development donors, predicts
the 45 sealed donors, and reports per primary target the sealed R^2(M1) - R^2(B0).
Null: 1,000 shuffles of the image block among sealed donors within Hardy x
ischemic-tertile strata, with the development fit fixed. BH within the primary
family. Two versions side by side: primary, and residualized on the two RNA
negative-control scores with OLS coefficients estimated on development donors only.
Refuses to run if the sealed output directory already exists.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from statsmodels.stats.multitest import multipletests

import c4_ridge as c4

SEED = 20260923


def fit_predict(Xb_dev, E_dev, Y_dev, Xb_s, E_s, n_pc):
    mu, sd = Xb_dev.mean(0), np.where(Xb_dev.std(0) > 0, Xb_dev.std(0), 1.0)
    xd, xs = (Xb_dev - mu) / sd, (Xb_s - mu) / sd
    m0 = RidgeCV(alphas=c4.ALPHAS, alpha_per_target=True).fit(xd, Y_dev)
    b0_dev = m0.predict(xd).reshape(len(xd), -1)
    b0_s = m0.predict(xs).reshape(len(xs), -1)
    pca = PCA(n_components=min(n_pc, len(E_dev) - 1), random_state=SEED).fit(E_dev)
    m1 = RidgeCV(alphas=c4.ALPHAS, alpha_per_target=True).fit(pca.transform(E_dev), Y_dev - b0_dev)
    return b0_s, (m1, pca)


def r2(y, p):
    return 1 - ((y - p) ** 2).sum(0) / ((y - y.mean(0)) ** 2).sum(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--c2", required=True)
    ap.add_argument("--n-null", type=int, default=1000)
    a = ap.parse_args()
    here = Path(__file__).resolve().parent
    c4.check_prespec(here)
    out = Path(a.c2) / "c5_sealed_test"
    if out.exists():
        raise SystemExit(f"{out} exists: the sealed test runs once")
    out.mkdir()
    c2 = Path(a.c2)
    ceil = pd.read_csv(c2 / "targets/program_ceilings.tsv", sep="\t").set_index("program_uid")["ceiling"]
    dev = pd.read_csv(c2 / "targets/targets_development.tsv", sep="\t").merge(
        pd.read_csv(c2 / "embed/slide_features.tsv.gz", sep="\t"), on="SUBJID")
    sea = pd.read_csv(c2 / "targets/sealed/targets_sealed.tsv", sep="\t").merge(
        pd.read_csv(c2 / "embed/sealed/slide_features.tsv.gz", sep="\t"), on="SUBJID")
    prog = [c for c in dev.columns if c.startswith("hotspot_")]
    primary = [p for p in prog if ceil.get(p, np.nan) >= 0.6 and dev[p].notna().all() and sea[p].notna().all()]
    neg = ["neg_hypoxia", "neg_immediate_early"]
    emb = [c for c in dev.columns if c.startswith("emb_")]
    both = pd.concat([dev.assign(_s=0), sea.assign(_s=1)], ignore_index=True)
    Xb = c4.baseline(both).to_numpy(float)
    Xb_dev, Xb_s = Xb[both["_s"] == 0], Xb[both["_s"] == 1]
    E_dev, E_s = dev[emb].to_numpy(float), sea[emb].to_numpy(float)
    strata = sea["hardy"].fillna(-1).astype(str) + "_" + pd.qcut(sea["ischemic_min"].rank(method="first"), 3, labels=False).astype(str)
    rng = np.random.default_rng(SEED)
    results = {}
    for version in ("primary", "residualized"):
        Yd, Ys = dev[primary].to_numpy(float), sea[primary].to_numpy(float)
        if version == "residualized":
            Zd = np.column_stack([np.ones(len(dev)), dev[neg].to_numpy(float)])
            Zs = np.column_stack([np.ones(len(sea)), sea[neg].to_numpy(float)])
            beta = np.linalg.lstsq(Zd, Yd, rcond=None)[0]
            Yd, Ys = Yd - Zd @ beta, Ys - Zs @ beta
        b0_s, (m1, pca) = fit_predict(Xb_dev, E_dev, Yd, Xb_s, E_s, 256)
        m1_s = b0_s + m1.predict(pca.transform(E_s)).reshape(len(sea), -1)
        gain = r2(Ys, m1_s) - r2(Ys, b0_s)
        null = np.zeros((a.n_null, len(primary)))
        for i in range(a.n_null):
            perm = np.arange(len(sea))
            for _, idx in sea.groupby(strata).groups.items():
                idx = np.asarray(list(idx)); perm[idx] = rng.permutation(idx)
            null[i] = r2(Ys, b0_s + m1.predict(pca.transform(E_s[perm])).reshape(len(sea), -1)) - r2(Ys, b0_s)
        p = (1 + (null >= gain).sum(0)) / (1 + a.n_null)
        q = multipletests(p, method="fdr_bh")[1]
        tab = pd.DataFrame({"target": primary, "r2_B0": r2(Ys, b0_s), "r2_M1": r2(Ys, m1_s), "gain": gain,
                            "p_perm": p, "q_bh": q, "readable": (q < 0.05) & (gain >= 0.05)})
        tab.to_csv(out / f"sealed_gains_{version}.tsv", sep="\t", index=False)
        results[version] = {"n_readable": int(tab["readable"].sum()), "median_gain": float(np.median(gain))}
    a1 = pd.read_csv(out / "sealed_gains_primary.tsv", sep="\t").set_index("target")["readable"]
    a2 = pd.read_csv(out / "sealed_gains_residualized.tsv", sep="\t").set_index("target")["readable"]
    results["readable_in_both"] = int((a1 & a2).sum())
    results["n_sealed_donors"] = int(len(sea)); results["n_primary_targets"] = len(primary)
    (out / "c5_summary.json").write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
