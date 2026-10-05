#!/usr/bin/env python3
"""C6: Model C external test in TCGA-LIHC non-tumor liver (PRESPEC_C Amendment 6). Run ONCE.

Family: the 16 programs that met the readability rule in development and sealed GTEx
donors, primary and residualized (c5_sealed_test/readability_four_criteria.tsv).
Model: the image part of M1, fitted on the 179 GTEx development donors exactly as
in c5_sealed_test.py and not refitted. GTEx B0 covariates (ischemic time, Hardy
class, autolysis, RIN) do not exist in TCGA, so only the image part is applied.
Statistic per program: partial Spearman between the image prediction and the
measured score given fat-vacuole fraction, age and sex; donor is the unit.
Null: 1,000 shuffles of the image prediction among donors within tissue source
site (covariates left paired). Primary endpoint: median partial Spearman over the
family. Per program: one-sided p, BH over the family. Two versions side by side:
primary, and residualized on the two RNA negative-control scores.
Before any TCGA value is read, the script recomputes the GTEx sealed gains and
stops unless they match the saved table. Refuses to run if the output exists.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

import c4_ridge as c4
import c5_sealed_test as c5

SEED = 20261003
NEG = ["neg_hypoxia", "neg_immediate_early"]
MIN_RHO = 0.30


def partial_spearman(pred, Y, C):
    """Pearson correlation of rank residuals after regressing ranks on ranked covariates."""
    Z = np.column_stack([np.ones(len(C))] + [rankdata(C[:, j]) for j in range(C.shape[1])])
    H = np.eye(len(C)) - Z @ np.linalg.pinv(Z)
    rp = H @ np.apply_along_axis(rankdata, 0, pred)
    ry = H @ np.apply_along_axis(rankdata, 0, Y)
    return (rp * ry).sum(0) / np.sqrt((rp ** 2).sum(0) * (ry ** 2).sum(0))


def resid_on_neg(Y, N, beta=None):
    Z = np.column_stack([np.ones(len(N)), N])
    if beta is None:
        beta = np.linalg.lstsq(Z, Y, rcond=None)[0]
    return Y - Z @ beta, beta


def reproduce_sealed(c2, dev, emb):
    """Recompute the GTEx sealed primary gains with the c5 code path and compare."""
    saved = pd.read_csv(c2 / "c5_sealed_test/sealed_gains_primary.tsv", sep="\t")
    sea = pd.read_csv(c2 / "targets/sealed/targets_sealed.tsv", sep="\t").merge(
        pd.read_csv(c2 / "embed/sealed/slide_features.tsv.gz", sep="\t"), on="SUBJID")
    primary = list(saved["target"])
    both = pd.concat([dev.assign(_s=0), sea.assign(_s=1)], ignore_index=True)
    Xb = c4.baseline(both).to_numpy(float)
    Xb_dev, Xb_s = Xb[both["_s"] == 0], Xb[both["_s"] == 1]
    Ys = sea[primary].to_numpy(float)
    b0_s, (m1, pca) = c5.fit_predict(Xb_dev, dev[emb].to_numpy(float), dev[primary].to_numpy(float),
                                     Xb_s, sea[emb].to_numpy(float), 256)
    m1_s = b0_s + m1.predict(pca.transform(sea[emb].to_numpy(float))).reshape(len(sea), -1)
    gain = c5.r2(Ys, m1_s) - c5.r2(Ys, b0_s)
    diff = float(np.abs(gain - saved["gain"].to_numpy()).max())
    if diff > 1e-8:
        raise SystemExit(f"GTEx sealed gains do not reproduce (max abs difference {diff}); TCGA not read")
    return diff, Xb_dev  # B0 design as built in c5 (median fill and Hardy levels over all 224 donors)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--c2", required=True, help="GTEx execution directory (he-molecular-c2-*)")
    ap.add_argument("--tcga-targets", required=True)
    ap.add_argument("--tcga-features", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-null", type=int, default=1000)
    a = ap.parse_args()
    here = Path(__file__).resolve().parent
    c4.check_prespec(here)
    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"{out} exists: the external test runs once")
    c2 = Path(a.c2)

    dev = pd.read_csv(c2 / "targets/targets_development.tsv", sep="\t").merge(
        pd.read_csv(c2 / "embed/slide_features.tsv.gz", sep="\t"), on="SUBJID")
    emb = [c for c in dev.columns if c.startswith("emb_")]
    crit = pd.read_csv(c2 / "c5_sealed_test/readability_four_criteria.tsv", sep="\t")
    family = list(crit.loc[crit["all_four"], "target"])
    if len(family) != 16:
        raise SystemExit(f"expected 16 programs in the family, found {len(family)}")
    repro_diff, Xb_dev = reproduce_sealed(c2, dev, emb)
    out.mkdir(parents=True)

    tc = pd.read_csv(a.tcga_targets, sep="\t").merge(pd.read_csv(a.tcga_features, sep="\t"), on="SUBJID")
    tc = tc[tc["n_tiles"] > 0].reset_index(drop=True)
    tested = [p for p in family if p in tc.columns and tc[p].notna().all()]
    C = np.column_stack([tc["fat_fraction"], tc["age"].fillna(tc["age"].median()), (tc["gender"] == "male").astype(float)])
    E_dev, E_t = dev[emb].to_numpy(float), tc[emb].to_numpy(float)
    strata = tc.groupby("tss").indices
    rng = np.random.default_rng(SEED)
    perms = np.tile(np.arange(len(tc)), (a.n_null, 1))
    for i in range(a.n_null):
        for idx in strata.values():
            perms[i, idx] = rng.permutation(idx)

    results = {"n_donors": int(len(tc)), "family": len(family), "tested": len(tested),
               "untestable": sorted(set(family) - set(tested)),
               "donors_fixed_under_null": int(sum(len(v) == 1 for v in strata.values())),
               "gtex_sealed_reproduction_max_abs_diff": repro_diff}
    tabs = {}
    for version in ("primary", "residualized"):
        Yd, Yt = dev[tested].to_numpy(float), tc[tested].to_numpy(float)
        if version == "residualized":
            Yd, _ = resid_on_neg(Yd, dev[NEG].to_numpy(float))
            Yt, _ = resid_on_neg(Yt, tc[NEG].to_numpy(float))  # within TCGA; no image value involved
        _, (m1, pca) = c5.fit_predict(Xb_dev, E_dev, Yd, Xb_dev[:1], E_dev[:1], 256)
        pred = m1.predict(pca.transform(E_t)).reshape(len(tc), -1)
        rho = partial_spearman(pred, Yt, C)
        null = np.stack([partial_spearman(pred[p], Yt, C) for p in perms])
        p = (1 + (null >= rho).sum(0)) / (1 + a.n_null)
        q = multipletests(p, method="fdr_bh")[1]
        tab = pd.DataFrame({"target": tested, "partial_spearman": rho, "p_perm": p, "q_bh": q,
                            "transfers": (q < 0.05) & (rho >= MIN_RHO)})
        tab.to_csv(out / f"external_{version}.tsv", sep="\t", index=False)
        tabs[version] = tab.set_index("target")["transfers"]
        med, med_null = float(np.median(rho)), np.median(null, axis=1)
        results[version] = {"median_partial_spearman": med,
                            "p_median": float((1 + (med_null >= med).sum()) / (1 + a.n_null)),
                            "n_transfer": int(tab["transfers"].sum())}
    results["transfer_in_both"] = int((tabs["primary"] & tabs["residualized"]).sum())

    # negative controls: same fit and statistic on the two RNA negative-control scores
    _, (m1, pca) = c5.fit_predict(Xb_dev, E_dev, dev[NEG].to_numpy(float), Xb_dev[:1], E_dev[:1], 256)
    pred = m1.predict(pca.transform(E_t)).reshape(len(tc), -1)
    Yn = tc[NEG].to_numpy(float)
    rho = partial_spearman(pred, Yn, C)
    null = np.stack([partial_spearman(pred[p], Yn, C) for p in perms])
    pd.DataFrame({"target": NEG, "partial_spearman": rho,
                  "p_perm": (1 + (null >= rho).sum(0)) / (1 + a.n_null)}).to_csv(
        out / "external_negative_controls.tsv", sep="\t", index=False)
    results["negative_controls"] = dict(zip(NEG, map(float, rho)))
    (out / "c6_summary.json").write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
