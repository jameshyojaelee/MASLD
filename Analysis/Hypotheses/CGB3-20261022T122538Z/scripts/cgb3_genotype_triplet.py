#!/usr/bin/env python
"""CGB3: total vs stage-adjusted genotype effects on the continuum, with a
no-direct-effect simulation. Follows PRESPEC_CGB3.md."""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CG0-20261022T122300Z/scripts")
import cg_latent as cl  # noqa: E402

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
GENO = f"{ROOT}/Analysis/Hypotheses/H3-20261009T180331Z/results/calls-20261009T181453Z/person_genotypes.tsv"
VARIANTS = {"rs738409": "PNPLA3", "rs58542926": "TM6SF2", "rs1260326": "GCKR", "rs2642438": "MTARC1"}
COHORTS = ["GSE130970", "GSE135251", "GSE162694", "GSE213621"]
SIGMAS = [0.4655, 0.5092, 0.6433]
SEED = 20261022
NSIM = 2000


def load(units_path):
    u = pd.read_csv(units_path, sep="\t")
    g = pd.read_csv(GENO, sep="\t")
    lib = g.assign(lib=g.libraries.str.split(";|,")).explode("lib")
    lib["lib"] = lib.lib.str.strip()
    d = u.merge(lib, left_on="sample_id", right_on="lib", how="left", suffixes=("", "_g"))
    d = d[d.dataset.isin(COHORTS) & d.fibrosis_stage.notna()]
    d = d[~((d.dataset == "GSE213621") & (d.fibrosis_stage == 0))].copy()
    for rs in VARIANTS:
        d[rs] = d[f"{rs}_call"].map({"RR": 0, "RA": 1, "AA": 2})
    d["group"] = d.fibrosis_stage.astype(int)
    d["sex"] = (d.inferred_sex == "M").astype(int)
    sys_cols = [c for c in u.columns if c.startswith("sys_")]
    for ds, idx in d.groupby("dataset").groups.items():
        for c in ["continuum_z"] + sys_cols:
            v = d.loc[idx, c]
            d.loc[idx, c] = (v - v.mean()) / v.std(ddof=1)
        s = d.loc[idx, "group"]
        d.loc[idx, "stage_z"] = (s - s.mean()) / s.std(ddof=1)
    X = d[sys_cols].fillna(0).values
    X = X - X.mean(0)
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    for i in range(3):
        d[f"sysPC{i+1}"] = X @ Vt[i]
    return d


def ols_coef(df, formula, term):
    f = smf.ols(formula, data=df).fit()
    return f.params[term], f.bse[term]


def meta(est):
    b = np.array([e[0] for e in est])
    se = np.array([e[1] for e in est])
    w = 1 / se ** 2
    bm = (w * b).sum() / w.sum()
    sm = np.sqrt(1 / w.sum())
    return bm, sm, 2 * stats.norm.sf(abs(bm / sm))


def triplet(d, rs, y):
    per = []
    for ds, s in d.dropna(subset=[rs, y]).groupby("dataset"):
        if s[rs].nunique() < 2 or len(s) < 30:
            continue
        gs = ols_coef(s, f"stage_z ~ {rs} + sex", rs)
        tot = ols_coef(s, f"{y} ~ {rs} + sex", rs)
        adj = ols_coef(s, f"{y} ~ {rs} + C(group) + sex", rs)
        rho = stats.pearsonr(smf.ols("stage_z ~ sex", data=s).fit().resid,
                             smf.ols(f"{y} ~ sex", data=s).fit().resid)[0]
        per.append(dict(dataset=ds, n=len(s), maf=s[rs].mean() / 2, b_gs=gs[0], se_gs=gs[1], b_tot=tot[0],
                        se_tot=tot[1], b_adj=adj[0], se_adj=adj[1], rho=rho))
    p = pd.DataFrame(per)
    out = {}
    for k in ["gs", "tot", "adj"]:
        bm, sm, pv = meta(list(zip(p[f"b_{k}"], p[f"se_{k}"])))
        out[f"b_{k}"], out[f"se_{k}"], out[f"p_{k}"] = bm, sm, pv
    out["rho_mean"] = float(np.average(p.rho, weights=p.n))
    out["bound_B"] = -out["b_gs"] * out["rho_mean"]
    return p, out


def simulate_null(d, rs, y, per, sigma, rng, cache):
    """b_adj meta distribution when G acts on Y only through latent severity."""
    sims = np.empty(NSIM)
    setups = []
    for _, r in per.iterrows():
        s = d[(d.dataset == r.dataset)].dropna(subset=[rs, y])
        gmap = cl.group_map(r.dataset)
        key = (r.dataset, sigma)
        if key not in cache:
            par = cl.fit_latent(s.group.value_counts().to_dict(), gmap, sigma, SEED)
            cache[key] = (par, cl.group_moments(par, gmap, SEED))
        par, mom = cache[key]
        fit = smf.ols(f"{y} ~ C(group) + sex", data=s).fit()
        ms = s.sex.mean()
        means = {k: float(fit.predict(pd.DataFrame({"group": [k], "sex": [ms]}))[0]) for k in sorted(s.group.unique())}
        c0, a = cl.fit_a(means, s.group.value_counts().to_dict(), mom)
        var_u = max(float(np.var(s[y], ddof=1) - a ** 2 * par["s"] ** 2), 1e-3)
        G = s[rs].values.astype(float)
        sex = s.sex.values
        sd_stage = s.group.std(ddof=1)
        # delta in latent units reproducing observed b_gs (stage_z per allele), two-pass rescale
        delta = r.b_gs * sd_stage
        for _ in range(2):
            bs = []
            for _ in range(200):
                L = rng.normal(par["mu"], par["s"], len(G)) + delta * (G - G.mean())
                K = cl.to_group(cl.read_stage(L, par["sd_e"], rng), gmap).astype(float)
                bs.append(np.polyfit(G, (K - K.mean()) / max(K.std(ddof=1), 1e-9), 1)[0])
            m = np.mean(bs)
            if abs(m) > 1e-6:
                delta *= r.b_gs / m
        setups.append((par, gmap, G, sex, delta, c0, a, var_u))
    for i in range(NSIM):
        est = []
        for par, gmap, G, sex, delta, c0, a, var_u in setups:
            L = rng.normal(par["mu"], par["s"], len(G)) + delta * (G - G.mean())
            K = cl.to_group(cl.read_stage(L, par["sd_e"], rng), gmap)
            Y = c0 + a * L + rng.normal(0, np.sqrt(var_u), len(G))
            D = np.column_stack([np.ones(len(G)), G, sex] + [(K == k).astype(float) for k in np.unique(K)[1:]])
            beta, res, rank, _ = np.linalg.lstsq(D, Y, rcond=None)
            resid = Y - D @ beta
            s2 = resid @ resid / (len(Y) - rank)
            cov = s2 * np.linalg.pinv(D.T @ D)
            est.append((beta[1], np.sqrt(cov[1, 1])))
        sims[i] = meta(est)[0]
    return sims


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--units", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=False)
    rng = np.random.default_rng(SEED)
    d = load(a.units)
    res = {"n_units_with_stage": int(len(d)),
           "n_by_cohort": d.dataset.value_counts().to_dict(),
           "genotyped": {rs: int(d[rs].notna().sum()) for rs in VARIANTS}}
    # selection check
    res["sex_check"] = {rs: float(stats.chi2_contingency(pd.crosstab(d[rs], d.sex))[1]) for rs in VARIANTS}
    cache, rows, percohort = {}, [], []
    for rs, gene in VARIANTS.items():
        for y in ["continuum_z", "sysPC1", "sysPC2", "sysPC3"]:
            per, out = triplet(d, rs, y)
            per["variant"], per["outcome"] = rs, y
            percohort.append(per)
            row = dict(variant=rs, gene=gene, outcome=y, **out)
            if y == "continuum_z":
                for sg in SIGMAS:
                    sims = simulate_null(d, rs, y, per, sg, rng, cache)
                    lo, hi = np.quantile(sims, [0.025, 0.975])
                    row[f"null_mean_{sg}"] = float(sims.mean())
                    row[f"null_lo_{sg}"], row[f"null_hi_{sg}"] = float(lo), float(hi)
                    outside = (out["b_adj"] < lo) or (out["b_adj"] > hi)
                    opposite = np.sign(out["b_adj"] - sims.mean()) != np.sign(sims.mean()) if abs(sims.mean()) > 0 else True
                    row[f"direct_{sg}"] = bool(outside and opposite)
                row["direct_effect_all_sigma"] = all(row[f"direct_{sg}"] for sg in SIGMAS)
            rows.append(row)
    t = pd.DataFrame(rows)
    prim = t.outcome == "continuum_z"
    t.loc[prim, "holm_p_tot"] = multipletests(t.loc[prim, "p_tot"], method="holm")[1]
    t.to_csv(f"{a.out}/triplets.tsv", sep="\t", index=False)
    pd.concat(percohort).to_csv(f"{a.out}/triplets_per_cohort.tsv", sep="\t", index=False)
    res["primary"] = t[prim].to_dict(orient="records")
    json.dump(res, open(f"{a.out}/results.json", "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
