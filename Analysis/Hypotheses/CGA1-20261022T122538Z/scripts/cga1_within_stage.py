#!/usr/bin/env python
"""CGA1: within-stage system variation, severity share vs orthogonal structure.
Follows PRESPEC_CGA1.md (sha256 in PRESPEC_CGA1.md.sha256)."""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from sklearn.mixture import GaussianMixture

sys.path.insert(0, "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CG0-20261022T122300Z/scripts")
import cg_latent as cl  # noqa: E402

SEED = 20261022
SIGMAS = [0.0, 0.4655, 0.5092, 0.6433]
PRIMARY = ["GSE162694", "GSE213621"]
SENS = ["GSE135251", "GSE130970"]
NBOOT = 1000
NMIX = 1000


def cohort_frame(u, ds, groups=None):
    d = u[(u.dataset == ds) & u.fibrosis_stage.notna() & u.continuum_z.notna()].copy()
    if ds == "GSE213621":
        d = d[d.fibrosis_stage > 0]
    d["group"] = d.fibrosis_stage.astype(int)
    if groups is not None:
        d = d[d.group.isin(groups)]
    d["sex"] = (d.inferred_sex == "M").astype(int)
    sys_cols = [c for c in u.columns if c.startswith("sys_")]
    for c in sys_cols + ["continuum_z"]:
        d[c] = (d[c] - d[c].mean()) / d[c].std(ddof=1)
    return d.reset_index(drop=True), sys_cols


def resid(d, col):
    return smf.ols(f"{col} ~ C(group) + sex", data=d).fit().resid.values


def observed_c(d):
    fit = smf.ols("continuum_z ~ C(group) + sex", data=d).fit()
    ms = d.sex.mean()
    means = {k: float(fit.predict(pd.DataFrame({"group": [k], "sex": [ms]}))[0]) for k in sorted(d.group.unique())}
    n = d.group.value_counts().to_dict()
    return means, n, float(fit.mse_resid)


def part1(d, ds, latent_cache):
    gmap = cl.group_map(ds)
    counts = d.group.value_counts().to_dict()
    rows = []
    rng = np.random.default_rng(SEED)
    for sigma in SIGMAS:
        key = (ds, sigma)
        if key not in latent_cache:
            par = cl.fit_latent(counts, gmap, sigma, SEED)
            latent_cache[key] = (par, cl.group_moments(par, gmap, SEED))
        par, mom = latent_cache[key]
        means, n, mse = observed_c(d)
        c0, a = cl.fit_a(means, n, mom)
        wvar = sum(n[k] * mom[k][1] for k in n if np.isfinite(mom[k][1])) / sum(n.values())
        pi = a ** 2 * wvar / mse
        boots = []
        for _ in range(NBOOT):
            b = d.iloc[rng.integers(0, len(d), len(d))]
            if b.group.nunique() < len(n):
                continue
            mb, nb, mseb = observed_c(b)
            _, ab = cl.fit_a(mb, nb, mom)
            wb = sum(nb[k] * mom[k][1] for k in nb if np.isfinite(mom[k][1])) / sum(nb.values())
            boots.append(ab ** 2 * wb / mseb)
        var_u = float(np.var(d.continuum_z, ddof=1) - a ** 2 * par["s"] ** 2)
        rows.append(dict(dataset=ds, sigma=sigma, n=len(d), mu=par["mu"], s=par["s"], sd_e=par["sd_e"],
                         sigma_achieved=par["sigma_achieved"], capped=par["capped"], c0=c0, a=a,
                         within_var_L=wvar, within_var_C=mse, pi_sev=pi,
                         pi_lo=np.quantile(boots, 0.025), pi_hi=np.quantile(boots, 0.975), var_u=var_u))
    return pd.DataFrame(rows)


def part2_3(d, sys_cols):
    R = np.column_stack([resid(d, c) for c in sys_cols])
    c = resid(d, "continuum_z")
    gamma = R.T @ c / (c @ c)
    Fc = float((gamma ** 2).sum() * c.var() / R.var(axis=0).sum())
    Rp = R - np.outer(c, gamma)
    Rp_c = Rp - Rp.mean(axis=0)
    U, S, Vt = np.linalg.svd(Rp_c, full_matrices=False)
    var_exp = S ** 2 / (S ** 2).sum()
    return dict(Fc=Fc, gamma=gamma, Vt=Vt[:3], var_exp=var_exp[:3], scores=Rp_c @ Vt[:3].T, Rp=Rp)


def congruence(Va, Vb):
    M = np.abs(Va @ Vb.T)
    return M


def comp_r2(d, score):
    bp = [c for c in d.columns if c.startswith("bp_") and c != "bp_Hepatocytes"]
    X = d[bp].values
    X = np.column_stack([np.ones(len(X)), X])
    beta, *_ = np.linalg.lstsq(X, score, rcond=None)
    pred = X @ beta
    r2 = 1 - ((score - pred) ** 2).sum() / ((score - score.mean()) ** 2).sum()
    p = X.shape[1] - 1
    return float(r2), float(1 - (1 - r2) * (len(score) - 1) / (len(score) - p - 1))


def f2_mixture(d, ds, latent_cache, p1, rng):
    x = d.loc[d.group == 2, "continuum_z"].values.reshape(-1, 1)
    if len(x) < 15:
        return []

    def dbic(v):
        b1 = GaussianMixture(1, random_state=SEED).fit(v).bic(v)
        b2 = GaussianMixture(2, random_state=SEED, n_init=3).fit(v).bic(v)
        return b1 - b2

    obs = dbic(x)
    rows = []
    gmap = cl.group_map(ds)
    for sigma in SIGMAS:
        par, _ = latent_cache[(ds, sigma)]
        pr = p1[p1.sigma == sigma].iloc[0]
        sd_u = np.sqrt(max(pr.var_u, 1e-4))
        null = []
        for _ in range(NMIX):
            L = rng.normal(par["mu"], par["s"], 20 * len(x))
            g = cl.to_group(cl.read_stage(L, par["sd_e"], rng), gmap)
            Lk = L[g == 2][: len(x)]
            if len(Lk) < len(x):
                continue
            cs = pr.c0 + pr.a * Lk + rng.normal(0, sd_u, len(Lk))
            null.append(dbic(cs.reshape(-1, 1)))
        null = np.array(null)
        rows.append(dict(dataset=ds, sigma=sigma, n_f2=len(x), dbic_obs=obs, null_q95=np.quantile(null, 0.95),
                         p_emp=(1 + np.sum(null >= obs)) / (len(null) + 1)))
    return rows


def windows(d, scores, axis_names):
    rows = []
    for k, s in d.assign(**{n: scores[:, i] for i, n in enumerate(axis_names)}).groupby("group"):
        s = s.sort_values("continuum_z")
        for start in range(0, max(len(s) - 30, 0) + 1, 10):
            w = s.iloc[start:start + 30]
            rows.append(dict(group=k, start=start, n=len(w), continuum_mid=w.continuum_z.median(),
                             **{n: w[n].mean() for n in axis_names}))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--units", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=False)
    u = pd.read_csv(a.units, sep="\t")
    rng = np.random.default_rng(SEED)
    cache, res = {}, {}

    p1_all, p23 = [], {}
    frames = {}
    for ds in PRIMARY + SENS:
        d, sys_cols = cohort_frame(u, ds)
        frames[ds] = d
        p1_all.append(part1(d, ds, cache))
        p23[ds] = part2_3(d, sys_cols)
    p1 = pd.concat(p1_all)
    p1.to_csv(f"{a.out}/part1_severity_share.tsv", sep="\t", index=False)
    res["part1"] = p1[["dataset", "sigma", "n", "a", "pi_sev", "pi_lo", "pi_hi", "capped"]].to_dict(orient="records")
    res["part2_Fc"] = {ds: p23[ds]["Fc"] for ds in p23}

    A, B = p23[PRIMARY[0]], p23[PRIMARY[1]]
    M = congruence(A["Vt"], B["Vt"])
    res["part3_congruence_matrix_162694_rows_213621_cols"] = M.round(3).tolist()
    res["part3_var_explained"] = {PRIMARY[0]: A["var_exp"].tolist(), PRIMARY[1]: B["var_exp"].tolist()}
    sys_cols = [c for c in u.columns if c.startswith("sys_")]
    load = pd.DataFrame({"system": sys_cols, **{f"{PRIMARY[0]}_PC{i+1}": A["Vt"][i] for i in range(3)},
                         **{f"{PRIMARY[1]}_PC{i+1}": B["Vt"][i] for i in range(3)},
                         "gamma_" + PRIMARY[0]: A["gamma"], "gamma_" + PRIMARY[1]: B["gamma"]})
    load.to_csv(f"{a.out}/part3_orthogonal_loadings.tsv", sep="\t", index=False)
    axes = []
    for i in range(3):
        j = int(np.argmax(M[i]))
        rep = bool(M[i, j] >= 0.85)
        r2a = comp_r2(frames[PRIMARY[0]], A["scores"][:, i])
        r2b = comp_r2(frames[PRIMARY[1]], B["scores"][:, j])
        sexp = [float(smf.ols("y ~ sex", data=frames[ds].assign(y=sc)).fit().pvalues["sex"])
                for ds, sc in [(PRIMARY[0], A["scores"][:, i]), (PRIMARY[1], B["scores"][:, j])]]
        axes.append(dict(pc_162694=i + 1, best_pc_213621=j + 1, congruence=float(M[i, j]), replicated=rep,
                         comp_r2_adj_162694=r2a[1], comp_r2_adj_213621=r2b[1], comp_r2_raw=[r2a[0], r2b[0]],
                         sex_p=sexp))
    res["part3_axes"] = axes

    mix = []
    for ds in PRIMARY:
        mix += f2_mixture(frames[ds], ds, cache, p1[p1.dataset == ds], rng)
    pd.DataFrame(mix).to_csv(f"{a.out}/part4_f2_mixture.tsv", sep="\t", index=False)
    res["part4"] = mix

    # Part 5: F0-F2 restriction
    r = {}
    rp = {}
    for ds, groups in [("GSE162694", [0, 1, 2]), ("GSE213621", [1, 2])]:
        d, sc = cohort_frame(u, ds, groups)
        rp[ds] = part2_3(d, sc)
        sub = part1(d, ds, cache)
        r[ds] = dict(n=len(d), Fc=rp[ds]["Fc"], pi_sev=dict(zip(sub.sigma, sub.pi_sev)))
    rep_restr = []
    for ax in axes:
        i, j = ax["pc_162694"] - 1, ax["best_pc_213621"] - 1
        rep_restr.append(dict(pc_162694=i + 1, restricted_congruence_162694=float(np.max(np.abs(rp["GSE162694"]["Vt"] @ A["Vt"][i]))),
                              restricted_congruence_213621=float(np.max(np.abs(rp["GSE213621"]["Vt"] @ B["Vt"][j])))))
    r["axis_persistence"] = rep_restr
    res["part5_F0_F2"] = r

    # display windows (descriptive)
    pd.DataFrame(windows(frames[PRIMARY[0]], A["scores"], ["PC1", "PC2", "PC3"])).assign(dataset=PRIMARY[0]).to_csv(
        f"{a.out}/display_windows_{PRIMARY[0]}.tsv", sep="\t", index=False)
    pd.DataFrame(windows(frames[PRIMARY[1]], B["scores"], ["PC1", "PC2", "PC3"])).assign(dataset=PRIMARY[1]).to_csv(
        f"{a.out}/display_windows_{PRIMARY[1]}.tsv", sep="\t", index=False)

    # decision per prespec
    pc = p1[p1.dataset.isin(PRIMARY)]
    central = pc[pc.sigma == 0.5092].set_index("dataset").pi_sev
    sev = bool((central >= 0.5).all())
    below_all = bool((pc[pc.sigma > 0].pi_sev < 0.5).all())
    other_axes = [ax for ax, rr in zip(axes, rep_restr) if ax["replicated"]
                  and max(ax["comp_r2_adj_162694"], ax["comp_r2_adj_213621"]) < 0.5
                  and min(rr["restricted_congruence_162694"], rr["restricted_congruence_213621"]) >= 0.85]
    if sev:
        dec = "within-stage continuum variation is mainly unrecorded severity"
    elif below_all and other_axes:
        dec = "separate within-stage biology"
    else:
        dec = "scenario-dependent; no claim"
    res["decision"] = dec
    res["decision_inputs"] = dict(pi_sev_central=central.to_dict(), pi_below_0_5_all_sigma=below_all,
                                  qualifying_orthogonal_axes=other_axes)
    res["approx_severity_share_of_system_within_var"] = {
        ds: {str(s): p23[ds]["Fc"] * float(p1[(p1.dataset == ds) & (p1.sigma == s)].pi_sev.iloc[0]) for s in SIGMAS}
        for ds in PRIMARY}
    json.dump(res, open(f"{a.out}/results.json", "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
