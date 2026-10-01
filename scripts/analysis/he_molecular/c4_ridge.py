#!/usr/bin/env python3
"""C4: Model C primary analysis on development donors (PRESPEC_C section 6).

B0 = fat-vacuole fraction + ischemic time + Hardy class + autolysis score + RIN +
     age bracket + sex.  M1 = B0 prediction + ridge on embedding PCs (fit on training
     folds, not rescaled) fitted to B0's training residuals (Amendment 3).
Ridge, penalty chosen per target by efficient leave-one-out (GCV) on the training
fold; outer donor 5-fold x 3 repeats (seed 20260923). Metric per target: held-out
R^2(M1) - R^2(B0), R^2 from pooled out-of-fold predictions, averaged over repeats.
Null: 1,000 refits with the embedding block shuffled across donors within strata of
Hardy class x ischemic-time tertile (B0 left paired). BH within the primary family
(programs with ceiling >= 0.6). "Readable from H&E": BH q < 0.05 and gain >= 0.05.
Negative controls reported beside: hypoxia and immediate-early scores, and how well
embeddings alone predict ischemic time, Hardy class and autolysis.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import RepeatedKFold
from statsmodels.stats.multitest import multipletests

SEED = 20260923
ALPHAS = np.logspace(-2, 8, 40)
AGE = {"20-29": 0, "30-39": 1, "40-49": 2, "50-59": 3, "60-69": 4, "70-79": 5}


def check_prespec(here):
    for line in (here / "PRESPEC.sha256").read_text().splitlines():
        want, name = line.split()
        if hashlib.sha256((here / name).read_bytes()).hexdigest() != want:
            raise SystemExit(f"{name} does not match PRESPEC.sha256")


def baseline(d):
    b = pd.DataFrame(index=d.index)
    b["fat_fraction"] = d["fat_fraction"]
    b["ischemic_min"] = d["ischemic_min"]
    b["rin"] = d["rin"]
    b["age"] = d["age_bracket"].map(AGE)
    b["male"] = (d["sex"] == 1).astype(float)
    b["autolysis"] = d["autolysis"]
    b["autolysis_missing"] = d["autolysis"].isna().astype(float)
    hardy = d["hardy"].fillna(-1).astype(int)
    for h in sorted(hardy.unique())[1:]:
        b[f"hardy_{h}"] = (hardy == h).astype(float)
    return b.fillna(b.median())


def standardize(train, test):
    mu, sd = train.mean(0), train.std(0)
    sd = np.where(sd > 0, sd, 1.0)
    return (train - mu) / sd, (test - mu) / sd


def oof_r2(X_b0, E, Y, n_pc, splits):
    """Held-out R^2 for B0 and M1, pooled out-of-fold predictions per repeat, averaged.
    M1 = B0 prediction + ridge on (unscaled) embedding PC scores fitted to the B0
    training residuals, so the gain is what the image adds beyond B0 (Amendment 3)."""
    r2 = {"B0": [], "M1": []}
    for rep_folds in splits:
        pred = {"B0": np.zeros_like(Y), "M1": np.zeros_like(Y)}
        for tr, te in rep_folds:
            xb_tr, xb_te = standardize(X_b0[tr], X_b0[te])
            m0 = RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(xb_tr, Y[tr])
            fit_tr = m0.predict(xb_tr).reshape(len(tr), -1)
            pred["B0"][te] = m0.predict(xb_te).reshape(len(te), -1)
            pca = PCA(n_components=min(n_pc, len(tr) - 1), random_state=SEED).fit(E[tr])
            e_tr, e_te = pca.transform(E[tr]), pca.transform(E[te])  # not rescaled: trailing PCs stay small
            m1 = RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(e_tr, Y[tr] - fit_tr)
            pred["M1"][te] = pred["B0"][te] + m1.predict(e_te).reshape(len(te), -1)
        sst = ((Y - Y.mean(0)) ** 2).sum(0)
        for name in r2:
            r2[name].append(1 - ((Y - pred[name]) ** 2).sum(0) / sst)
    return {k: np.mean(v, 0) for k, v in r2.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True, help="targets_development.tsv")
    ap.add_argument("--ceilings", required=True)
    ap.add_argument("--features", required=True, help="development slide_features.tsv.gz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-null", type=int, default=1000)
    ap.add_argument("--n-pc", type=int, default=256)
    ap.add_argument("--residualize-negatives", action="store_true",
                    help="Amendment 4: remove from each program score its OLS fit on the two RNA negative-control scores")
    a = ap.parse_args()
    here = Path(__file__).resolve().parent
    check_prespec(here)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    t = pd.read_csv(a.targets, sep="\t")
    f = pd.read_csv(a.features, sep="\t")
    d = t.merge(f, on="SUBJID").reset_index(drop=True)
    if d["sealed"].any():
        raise SystemExit("sealed donor in development input")
    ceil = pd.read_csv(a.ceilings, sep="\t").set_index("program_uid")["ceiling"]
    prog = [c for c in t.columns if c.startswith("hotspot_")]
    primary = [p for p in prog if ceil.get(p, np.nan) >= 0.6 and d[p].notna().all()]
    neg = ["neg_hypoxia", "neg_immediate_early"]
    targets = primary + neg
    Y = d[targets].to_numpy(float)
    if a.residualize_negatives:
        Z = np.column_stack([np.ones(len(d)), d[neg].to_numpy(float)])
        beta = np.linalg.lstsq(Z, Y[:, :len(primary)], rcond=None)[0]
        Y[:, :len(primary)] = Y[:, :len(primary)] - Z @ beta
    E = d[[c for c in d.columns if c.startswith("emb_")]].to_numpy(float)
    X_b0 = baseline(d).to_numpy(float)
    rkf = RepeatedKFold(n_splits=5, n_repeats=3, random_state=SEED)
    folds = list(rkf.split(d))
    splits = [folds[i * 5:(i + 1) * 5] for i in range(3)]

    obs = oof_r2(X_b0, E, Y, a.n_pc, splits)
    gain = obs["M1"] - obs["B0"]

    strata = d["hardy"].fillna(-1).astype(str) + "_" + pd.qcut(d["ischemic_min"].rank(method="first"), 3, labels=False).astype(str)
    rng = np.random.default_rng(SEED)
    null = np.zeros((a.n_null, len(targets)))
    for i in range(a.n_null):
        perm = np.arange(len(d))
        for _, idx in d.groupby(strata).groups.items():
            idx = np.asarray(list(idx))
            perm[idx] = rng.permutation(idx)
        r = oof_r2(X_b0, E[perm], Y, a.n_pc, splits)
        null[i] = r["M1"] - r["B0"]
    p = (1 + (null >= gain).sum(0)) / (1 + a.n_null)
    res = pd.DataFrame({"target": targets, "ceiling": [ceil.get(x, np.nan) for x in targets],
                        "r2_B0": obs["B0"], "r2_M1": obs["M1"], "gain": gain, "p_perm": p,
                        "null_q95": np.quantile(null, 0.95, axis=0), "family": ["primary"] * len(primary) + ["negative_control"] * len(neg)})
    fam = res["family"] == "primary"
    res.loc[fam, "q_bh"] = multipletests(res.loc[fam, "p_perm"], method="fdr_bh")[1]
    res["readable_from_he"] = fam & (res["q_bh"] < 0.05) & (res["gain"] >= 0.05)
    res.to_csv(out / "c4_program_gains.tsv", sep="\t", index=False)

    # embeddings alone -> postmortem covariates (negative-control check)
    nuis = {}
    for col in ("ischemic_min", "autolysis", "hardy"):
        ok = d[col].notna().to_numpy()
        sub_folds = list(RepeatedKFold(n_splits=5, n_repeats=3, random_state=SEED).split(np.where(ok)[0]))
        sub_splits = [sub_folds[i * 5:(i + 1) * 5] for i in range(3)]
        r = oof_r2(np.zeros((ok.sum(), 1)), E[ok], d.loc[ok, col].to_numpy(float)[:, None], a.n_pc, sub_splits)
        nuis[col] = float(r["M1"][0])
    summary = {"n_development_donors": int(len(d)), "n_primary_targets": len(primary),
               "n_readable_from_he": int(res["readable_from_he"].sum()),
               "median_gain_primary": float(res.loc[fam, "gain"].median()),
               "negative_control_gains": res.loc[~fam, ["target", "gain", "p_perm"]].to_dict("records"),
               "embeddings_predict_postmortem_r2": nuis, "n_null": a.n_null,
               "exposure": "Phikon-v2 pretraining included GTEx slides (images only): likely image exposure, no target exposure"}
    (out / "c4_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
