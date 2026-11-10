#!/usr/bin/env python3
"""A4 step 5d: tau*, D = tau*_PC - tau*_PP with paired block-jackknife SE, BH families,
null calibration. Adapted copy of H22 06_summarize.py (bh, jk_cov, tau* from delete
values); H22 files are not edited.

tau*_c = coef_c * sd_c * M_tot / h2g (Gazal 2017); sd_c over common SNPs; M_tot = baselineLD
M_5_50 base; h2g = M_5_50 . coef from the same regression; delete values give the
jackknife variance and the within-regression covariance of PC and PP.
F_D: D across S3-eligible liver traits, BH. F_tau: tau* of PC and PP across eligible traits, BH.
Usage: 05d_ldsc_summarize.py <A4> <BASELINE_PREFIX>
"""
import json
import os
import re
import sys
from math import erfc, sqrt

import numpy as np
import pandas as pd

TRAITS = ["alt", "ast", "ggt", "pdff", "cirr", "nafld", "height"]
LIVER = TRAITS[:-1]


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def norm_p(z):
    return erfc(abs(z) / sqrt(2))


def jk_cov(a, b):
    B = len(a)
    return (B - 1) / B * np.sum((a - a.mean()) * (b - b.mean()))


class Ctx:
    def __init__(self, a4, base_prefix):
        self.a4 = a4
        bm = []
        for c in range(1, 23):
            with open(f"{base_prefix}{c}.l2.M_5_50") as fh:
                bm.append([float(x) for x in fh.readline().split()])
        self.base_m550 = np.sum(bm, axis=0)
        hdr = pd.read_csv(f"{base_prefix}1.l2.ldscore.gz", sep="\t", nrows=0).columns[3:]
        self.base_names = [h[:-2] for h in hdr]
        self.M_tot = float(self.base_m550[0])
        mm = [pd.read_csv(os.path.join(a4, f"work/ldsc/packed/M.{c}.tsv"), sep="\t") for c in range(1, 23)]
        self.m550 = pd.concat(mm).groupby("annotation", sort=False)["M_5_50"].sum()
        st = pd.concat([pd.read_csv(os.path.join(a4, f"work/ldsc/annot_stats/{c}.tsv"), sep="\t") for c in range(1, 23)])
        g = st.groupby("annotation", sort=False)[["n_common_in", "n_common_total"]].sum()
        p = g["n_common_in"] / g["n_common_total"]
        self.sd, self.prop_common = np.sqrt(p * (1 - p)), p

    def run(self, trait, model, cols):
        pre = os.path.join(self.a4, "work/ldsc/h2", trait, model)
        names = self.base_names + cols
        M = np.concatenate([self.base_m550, [self.m550[c] for c in cols]])
        log = open(pre + ".log").read()
        enr = None
        if os.path.exists(pre + ".results"):
            res = pd.read_csv(pre + ".results", sep="\t")
            coef = res["Coefficient"].to_numpy(float)
            ldsc_se = res["Coefficient_std_error"].to_numpy(float)
            enr = res.set_index(res["Category"].str.replace(r"L2_\d+$", "", regex=True))
        else:
            m = re.search(r"Coefficients:(.*?)Coefficient SE:(.*?)Lambda GC", log, re.S)
            coef = np.array(m.group(1).split(), float)
            ldsc_se = np.array(m.group(2).split(), float)
        assert len(coef) == len(names), (trait, model, len(coef), len(names))
        dele = np.loadtxt(pre + ".part_delete", ndmin=2)
        assert dele.shape[1] == len(names)
        h2 = float(M @ coef)
        h2_del = dele @ M
        out = dict(coef=coef, ldsc_se=ldsc_se, h2=h2, h2_del=h2_del, names=names, enr=enr,
                   n_snp=int(re.findall(r"(\d+) SNPs remain", log)[-1]), B=dele.shape[0])
        out["tau_star"], out["tau_star_del"] = {}, {}
        for c in cols:
            i = names.index(c)
            k = self.sd[c] * self.M_tot
            out["tau_star"][c] = coef[i] * k / h2
            out["tau_star_del"][c] = dele[:, i] * k / h2_del
        return out


def contrast(r, a, b):
    da, db = r["tau_star_del"][a], r["tau_star_del"][b]
    va, vb, cab = jk_cov(da, da), jk_cov(db, db), jk_cov(da, db)
    d = r["tau_star"][a] - r["tau_star"][b]
    se = np.sqrt(va + vb - 2 * cab)
    return dict(tau_star_PC=r["tau_star"][a], se_PC=np.sqrt(va), tau_star_PP=r["tau_star"][b], se_PP=np.sqrt(vb),
                jk_corr=cab / np.sqrt(va * vb), D=d, se_D=se, z_D=d / se, p_D=norm_p(d / se))


def main(a4, base_prefix):
    ctx = Ctx(a4, base_prefix)
    R = os.path.join(a4, "results")
    meta = pd.read_csv(os.path.join(a4, "work/ldsc/annotations.tsv"), sep="\t").set_index("annotation")
    s3 = json.load(open(os.path.join(a4, "work/ldsc/s3_checks.json")))
    rows = []
    for t in TRAITS:
        r = ctx.run(t, "BASE", [])
        se = np.sqrt(jk_cov(r["h2_del"], r["h2_del"]))
        rows.append(dict(trait=t, h2=r["h2"], h2_se_jk=se, h2_z=r["h2"] / se, h2_z_ldsc_log=s3[t]["z"],
                         eligible=s3[t]["eligible"], n_snp=r["n_snp"]))
    base = pd.DataFrame(rows)
    base.to_csv(os.path.join(R, "ldsc_baseline_h2.tsv"), sep="\t", index=False)
    elig = [t for t in LIVER if s3[t]["eligible"]]

    out, tau = [], []
    models = {"MAIN": ("ALL_HEP", "PC", "PP")}
    if bool(meta.loc["PCq", "testable"]):
        models["SEC"] = ("ALL_HEP", "PCq", "PPq")
    for model, (allc, pc, pp) in models.items():
        for t in TRAITS:
            r = ctx.run(t, model, [allc, pc, pp])
            row = dict(model=model, trait=t, eligible=t in elig, n_snp=r["n_snp"], h2=r["h2"])
            row.update(contrast(r, pc, pp))
            for c, lab in ((allc, "ALL_HEP"), (pc, "PC"), (pp, "PP")):
                i = r["names"].index(c)
                se_ts = np.sqrt(jk_cov(r["tau_star_del"][c], r["tau_star_del"][c]))
                e = r["enr"].loc[c] if r["enr"] is not None and c in r["enr"].index else None
                tau.append(dict(model=model, trait=t, eligible=t in elig, annotation=lab, n_lifted_peaks=int(meta.loc[c, "n_lifted"]),
                                prop_common_snps=float(ctx.prop_common[c]), tau=r["coef"][i], tau_se=r["ldsc_se"][i],
                                tau_z=r["coef"][i] / r["ldsc_se"][i], tau_star=r["tau_star"][c], tau_star_se=se_ts,
                                tau_star_z=r["tau_star"][c] / se_ts, tau_star_p=norm_p(r["tau_star"][c] / se_ts),
                                enrichment=None if e is None else float(e["Enrichment"]),
                                enrichment_se=None if e is None else float(e["Enrichment_std_error"]),
                                enrichment_p=None if e is None else float(e["Enrichment_p"])))
            out.append(row)
    D = pd.DataFrame(out)
    D["q_F_D"] = np.nan
    m = (D["model"] == "MAIN") & D["eligible"]
    D.loc[m, "q_F_D"] = bh(D.loc[m, "p_D"].values)
    TT = pd.DataFrame(tau)
    TT["q_F_tau"] = np.nan
    m = (TT["model"] == "MAIN") & TT["eligible"] & TT["annotation"].isin(["PC", "PP"])
    TT.loc[m, "q_F_tau"] = bh(TT.loc[m, "tau_star_p"].values)

    # null calibration
    draws = sorted(meta["draw"].dropna().astype(int).unique())
    nul = []
    for t in elig + ["height"]:
        for d in draws:
            r = ctx.run(t, f"NULL_{d}", ["ALL_HEP", f"NULL_PC_{d}", f"NULL_PP_{d}"])
            c = contrast(r, f"NULL_PC_{d}", f"NULL_PP_{d}")
            nul.append(dict(trait=t, draw=d, D=c["D"], se_D=c["se_D"], z_D=c["z_D"],
                            tau_star_PC=c["tau_star_PC"], tau_star_PP=c["tau_star_PP"]))
    N = pd.DataFrame(nul)
    N.to_csv(os.path.join(R, "ldsc_null_draws.tsv"), sep="\t", index=False)
    cal = N.groupby("trait").agg(n_draws=("z_D", "size"), frac_abs_z_ge_1_96=("z_D", lambda z: float(np.mean(np.abs(z) >= 1.96))),
                                 null_D_sd=("D", "std"), null_mean_se_D=("se_D", "mean")).reset_index()
    obs = D[D["model"] == "MAIN"].set_index("trait")
    cal["observed_z_D"] = cal["trait"].map(obs["z_D"])
    cal["observed_abs_z_rank_among_nulls"] = [int(np.sum(np.abs(N.loc[N["trait"] == t, "z_D"]) >= abs(obs.loc[t, "z_D"])))
                                              for t in cal["trait"]]
    cal.to_csv(os.path.join(R, "ldsc_null_calibration.tsv"), sep="\t", index=False)
    D.to_csv(os.path.join(R, "ldsc_zone_D.tsv"), sep="\t", index=False)
    TT.to_csv(os.path.join(R, "ldsc_zone_tau.tsv"), sep="\t", index=False)
    print(D[D["model"] == "MAIN"][["trait", "eligible", "tau_star_PC", "tau_star_PP", "D", "se_D", "z_D", "p_D", "q_F_D"]].to_string())
    print(cal.to_string())


if __name__ == "__main__":
    main(*sys.argv[1:3])
