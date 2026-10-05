#!/usr/bin/env python3
"""F2: Model F statistics (PRESPEC_F sections 4-5). Run ONCE.

Reads the program scores written by f1_score_programs.R. Human reference: per-cohort
OLS of the program score on fibrosis stage and sex, fixed-effect meta-analysis over the
four cohorts, and the leave-one-cohort-out ceiling (V1). Mouse: per-contrast OLS of the
program score on case arm and covariates; fidelity = Spearman over programs against
the human stage coefficients, with a within-strata label-shuffle null and a bootstrap
interval over animals. V2: same-diet-class versus different-class similarity of mouse
program profiles, class labels permuted over contrasts. V3: leave-one-series-out.
Refuses to run if the output directory exists.
"""
import argparse
import hashlib
import json
import subprocess
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

SEED = 20261003
META = {"human": ["sample_id", "dataset", "fibrosis_stage", "inferred_sex"],
        "mouse": ["contrast", "series", "diet_class", "sample", "arm", "timepoint", "sex"]}


def ols(Y, X):
    """Coefficient, standard error and two-sided p for column 1 of X, per column of Y."""
    XtXi = np.linalg.pinv(X.T @ X)
    B = XtXi @ X.T @ Y
    df = len(X) - np.linalg.matrix_rank(X)
    s2 = ((Y - X @ B) ** 2).sum(0) / df
    se = np.sqrt(s2 * XtXi[1, 1])
    return B[1], se, 2 * stats.t.sf(np.abs(B[1] / se), df)


def standardize(Y):
    return (Y - Y.mean(0)) / Y.std(0, ddof=1)


def spearman(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return float(stats.spearmanr(a[ok], b[ok])[0]), int(ok.sum())


def bh(p):
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def mouse_design(d):
    cols = [np.ones(len(d)), (d["arm"] == "case").to_numpy(float)]
    strata = pd.Series("all", index=d.index)
    for v in ("timepoint", "sex"):
        if d[v].notna().all() and d[v].nunique() > 1:
            cols += [x.to_numpy(float) for _, x in pd.get_dummies(d[v], drop_first=True).items()]
            strata = strata + "|" + d[v].astype(str)
    return np.column_stack(cols), strata.to_numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True, help="directory written by f1_score_programs.R")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-null", type=int, default=1000)
    ap.add_argument("--n-class-perm", type=int, default=10000)
    a = ap.parse_args()
    here = Path(__file__).resolve().parent
    subprocess.run(["sha256sum", "-c", "--quiet", "PRESPEC_F.sha256"], cwd=here, check=True)
    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"{out} exists: the fidelity analysis runs once")
    sc = Path(a.scores)
    hum = pd.read_csv(sc / "human_program_scores.tsv", sep="\t")
    mou = pd.read_csv(sc / "mouse_program_scores.tsv", sep="\t", dtype={"sample": str})
    progs = [c for c in hum.columns if c not in META["human"]]
    assert progs == [c for c in mou.columns if c not in META["mouse"]] and len(progs) == 117
    out.mkdir(parents=True)
    rng = np.random.default_rng(SEED)

    # ---- human reference and ceiling (V1) ----
    cohorts = sorted(hum["dataset"].unique())
    hb, hs = {}, {}
    for ds in cohorts:
        d = hum[hum["dataset"] == ds]
        X = np.column_stack([np.ones(len(d)), d["fibrosis_stage"].to_numpy(float), (d["inferred_sex"] == "M").to_numpy(float)])
        hb[ds], hs[ds], _ = ols(standardize(d[progs].to_numpy(float)), X)
    hb, hs = pd.DataFrame(hb, index=progs), pd.DataFrame(hs, index=progs)

    def meta(cols):
        w = 1 / hs[cols] ** 2
        return (w * hb[cols]).sum(axis=1, skipna=False) / w.sum(axis=1, skipna=False), np.sqrt(1 / w.sum(axis=1, skipna=False))

    h_beta, h_se = meta(cohorts)
    h_p = pd.Series(2 * stats.norm.sf(np.abs(h_beta / h_se)), index=progs)
    human = pd.DataFrame({"program_uid": progs, "beta_stage": h_beta.to_numpy(), "se": h_se.to_numpy(), "p": h_p.to_numpy(),
                          "q_bh": bh(h_p.to_numpy())})
    human["testable"] = human["beta_stage"].notna()
    human["supported"] = human["q_bh"] < 0.05
    human = pd.concat([human, hb.add_prefix("beta_").reset_index(drop=True)], axis=1)
    human.to_csv(out / "human_reference.tsv", sep="\t", index=False)
    loco = {ds: spearman(hb[ds].to_numpy(), meta([c for c in cohorts if c != ds])[0].to_numpy())[0] for ds in cohorts}
    ceiling = float(np.median(list(loco.values())))
    hvec, hsup = human["beta_stage"].to_numpy(), human["supported"].to_numpy()

    # ---- mouse effects, fidelity, program states ----
    rows, states, prof = [], [], {}
    for cid, d in mou.groupby("contrast", sort=True):
        d = d.reset_index(drop=True)
        X, strata = mouse_design(d)
        Y = standardize(d[progs].to_numpy(float))
        beta, se, p = ols(Y, X)
        prof[cid] = beta
        rho, n_prog = spearman(beta, hvec)
        null = np.empty(a.n_null)
        boot = np.empty(a.n_null)
        arm = d["arm"].to_numpy()
        for i in range(a.n_null):
            Xp = X.copy()
            for s in np.unique(strata):
                idx = np.flatnonzero(strata == s)
                Xp[idx, 1] = X[rng.permutation(idx), 1]
            null[i] = spearman(ols(Y, Xp)[0], hvec)[0]
            b = np.concatenate([rng.choice(np.flatnonzero(arm == g), (arm == g).sum()) for g in ("case", "control")])
            boot[i] = spearman(ols(standardize(d[progs].to_numpy(float)[b]), X[b])[0], hvec)[0]
        q = bh(p)
        state = np.where(~np.isfinite(beta) | ~np.isfinite(hvec), "untestable",
                 np.where(~hsup, "human_unsupported",
                 np.where(q >= 0.05, "mouse_unchanged",
                 np.where(np.sign(beta) == np.sign(hvec), "concordant", "discordant"))))
        states.append(pd.DataFrame({"contrast": cid, "program_uid": progs, "mouse_beta": beta, "mouse_se": se, "mouse_p": p,
                                    "mouse_q_bh": q, "human_beta_stage": hvec, "state": state}))
        rows.append({"contrast": cid, "series": d["series"].iloc[0], "diet_class": d["diet_class"].iloc[0],
                     "n_case": int((arm == "case").sum()), "n_control": int((arm == "control").sum()), "n_programs": n_prog,
                     "fidelity": rho, "ci_low": float(np.nanpercentile(boot, 2.5)), "ci_high": float(np.nanpercentile(boot, 97.5)),
                     "p_perm": float((1 + (null >= rho).sum()) / (1 + a.n_null)) if np.isfinite(rho) else np.nan, "fidelity_over_ceiling": rho / ceiling,
                     **{f"n_{k}": int((state == k).sum()) for k in ("concordant", "discordant", "mouse_unchanged", "human_unsupported", "untestable")}})
    fid = pd.DataFrame(rows)
    fid["q_bh"] = bh(fid["p_perm"].to_numpy())
    fid.to_csv(out / "contrast_fidelity.tsv", sep="\t", index=False)
    pd.concat(states).to_csv(out / "program_states.tsv", sep="\t", index=False)

    # ---- V2: diet-class specificity of mouse program profiles ----
    ids = list(fid.loc[fid["fidelity"].notna(), "contrast"])  # a contrast with no testable program has no profile
    cls = fid.set_index("contrast")["diet_class"]
    ser = fid.set_index("contrast")["series"]
    pairs = [(i, j) for i, j in combinations(ids, 2) if ser[i] != ser[j]]
    r = np.array([spearman(prof[i], prof[j])[0] for i, j in pairs])
    pd.DataFrame({"a": [i for i, _ in pairs], "b": [j for _, j in pairs], "spearman": r,
                  "same_class": [cls[i] == cls[j] for i, j in pairs]}).to_csv(out / "mouse_pairwise_similarity.tsv", sep="\t", index=False)

    def v2(labels):
        same = np.array([labels[i] == labels[j] for i, j in pairs])
        if not same.any() or same.all():
            return np.nan, np.nan
        per = [r[same & np.array([labels[i] == c for i, _ in pairs])].mean() for c in sorted(set(labels.values()))
               if (same & np.array([labels[i] == c for i, _ in pairs])).any()]
        return r[same].mean() - r[~same].mean(), np.mean(per) - r[~same].mean()

    obs = v2(cls[ids].to_dict())
    nullv = np.array([v2(dict(zip(ids, rng.permutation(cls[ids].to_numpy())))) for _ in range(a.n_class_perm)])
    v2_res = {"n_pairs_same": int(sum(cls[i] == cls[j] for i, j in pairs)), "n_pairs_different": int(sum(cls[i] != cls[j] for i, j in pairs)),
              "mean_same": float(r[[cls[i] == cls[j] for i, j in pairs]].mean()), "mean_different": float(r[[cls[i] != cls[j] for i, j in pairs]].mean()),
              "difference": float(obs[0]), "p_perm": float((1 + np.nansum(nullv[:, 0] >= obs[0])) / (1 + np.isfinite(nullv[:, 0]).sum())),
              "difference_class_balanced": float(obs[1]),
              "p_perm_class_balanced": float((1 + np.nansum(nullv[:, 1] >= obs[1])) / (1 + np.isfinite(nullv[:, 1]).sum()))}

    # ---- V3: leave-one-series-out ----
    v3 = []
    for i in ids:
        same = [j for j in ids if cls[j] == cls[i] and ser[j] != ser[i]]
        other = [j for j in ids if cls[j] != cls[i] and ser[j] != ser[i]]
        if same:
            v3.append({"contrast": i, "diet_class": cls[i], "n_same_class": len(same),
                       "rho_same_class_mean": spearman(prof[i], np.nanmean([prof[j] for j in same], 0))[0],
                       "rho_other_class_mean": spearman(prof[i], np.nanmean([prof[j] for j in other], 0))[0]})
    pd.DataFrame(v3).to_csv(out / "leave_one_series_out.tsv", sep="\t", index=False)

    summary = {"seed": SEED, "n_null": a.n_null, "n_class_perm": a.n_class_perm,
               "human": {"participants": {ds: int((hum["dataset"] == ds).sum()) for ds in cohorts},
                         "programs_testable": int(human["testable"].sum()), "programs_supported": int(human["supported"].sum()),
                         "leave_one_cohort_out_spearman": loco, "ceiling_median": ceiling,
                         "interpretable": bool(ceiling >= 0.5)},
               "mouse": {"contrasts": len(fid), "contrasts_testable": len(ids), "fidelity_q_lt_0.05": int((fid["q_bh"] < 0.05).sum()),
                         "fidelity_median": float(fid["fidelity"].median())},
               "v2_class_specificity": v2_res, "v2_endpoint_met": bool(v2_res["p_perm"] < 0.05),
               "inputs_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(sc.glob("*_program_scores.tsv"))}}
    (out / "f2_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(fid.to_string(index=False))


if __name__ == "__main__":
    main()
