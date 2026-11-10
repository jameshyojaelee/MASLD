#!/usr/bin/env python3
"""Follow-up step 5: tau* with block-jackknife SE, families F5, F5-GGT, F6, F6-cond, F7 (BH).

tau*_c = coef_c * sd_c * M_tot / h2g (Gazal 2017); sd_c over common SNPs; M_tot = baselineLD
M_5_50 of the base column; h2g = M_5_50 . coef of the same regression. Delete values give the
jackknife variances and, for F7, the AST-ALT covariance (identical SNPs and blocks).
Usage: 05_summarize.py <F> <BASELINE_PREFIX> <H22>
"""
import json
import os
import sys
from math import erfc, sqrt

import numpy as np
import pandas as pd

TRAITS = ["ast_panukbb", "ast_mvp", "alt_panukbb", "alt_mvp", "ggt_panukbb", "pdff", "height_panukbb"]
MODELS = ["J7", "RM_E122", "RM_E100", "RM_E108", "RM_E066", "RMJ"]


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def p_two(z):
    return erfc(abs(z) / sqrt(2))


def p_one(z):
    return 0.5 * erfc(z / sqrt(2))


def jk_cov(a, b):
    B = len(a)
    return (B - 1) / B * np.sum((a - a.mean()) * (b - b.mean()))


def m550(prefix):
    return np.sum([[float(x) for x in open(f"{prefix}{c}.l2.M_5_50").readline().split()] for c in range(1, 23)], axis=0)


def names_of(prefix):
    return [h[:-2] for h in pd.read_csv(f"{prefix}1.l2.ldscore.gz", sep="\t", nrows=0).columns[3:]]


class Ctx:
    def __init__(self, f, base, h22):
        self.f = f
        self.base_names, self.base_m550 = names_of(base), m550(base)
        self.M_tot = float(self.base_m550[0])
        st = pd.concat([pd.read_csv(os.path.join(h22, f"ldsc/annot_stats/{c}.tsv"), sep="\t") for c in range(1, 23)] +
                       [pd.read_csv(os.path.join(f, f"ldsc/annot_stats/rm.{c}.tsv"), sep="\t") for c in range(1, 23)])
        g = st.groupby("annotation")[["n_common_in", "n_common_total"]].sum()
        self.prop = g["n_common_in"] / g["n_common_total"]
        self.sd = np.sqrt(self.prop * (1 - self.prop))
        self.h22 = h22

    def run(self, trait, model):
        if model == "CT_endothelial_cell":
            prefix = os.path.join(self.h22, "ldsc/models", model, f"{model}.")
        else:
            prefix = os.path.join(self.f, "ldsc/models", model, f"{model}.")
        cols = names_of(prefix)
        names = self.base_names + cols
        M = np.concatenate([self.base_m550, m550(prefix)])
        pre = os.path.join(self.f, "ldsc/h2", trait, model)
        res = pd.read_csv(pre + ".results", sep="\t")
        coef = res["Coefficient"].to_numpy(float)
        assert len(coef) == len(names), (trait, model)
        dele = np.loadtxt(pre + ".part_delete", ndmin=2)
        h2, h2_del = float(M @ coef), dele @ M
        out = {}
        for c in cols:
            i = names.index(c)
            k = self.sd[c] * self.M_tot
            ts, ts_del = coef[i] * k / h2, dele[:, i] * k / h2_del
            se = sqrt(jk_cov(ts_del, ts_del))
            jk_coef_se = sqrt(jk_cov(dele[:, i], dele[:, i]))
            out[c] = dict(trait=trait, model=model, annotation=c, prop_common_snps=float(self.prop[c]),
                          tau=coef[i], tau_se_ldsc=res["Coefficient_std_error"].iloc[i],
                          tau_se_jk=jk_coef_se, tau_z=coef[i] / jk_coef_se, tau_star=ts, tau_star_se=se,
                          tau_star_z=ts / se, enrichment=res["Enrichment"].iloc[i],
                          enrichment_p=res["Enrichment_p"].iloc[i], h2=h2, ts_del=ts_del)
        return out


def main(f, base, h22):
    ctx = Ctx(f, base, h22)
    rows, rr = [], {}
    for t in TRAITS:
        for m in MODELS:
            rr[(t, m)] = ctx.run(t, m)
            rows += list(rr[(t, m)].values())
    tab = pd.DataFrame([{k: v for k, v in r.items() if k != "ts_del"} for r in rows])
    tab["p_one_sided"] = [p_one(z) for z in tab["tau_z"]]
    R = os.path.join(f, "results")
    tab.to_csv(os.path.join(R, "all_coefficients.tsv"), sep="\t", index=False)
    chk = float(np.nanmax(np.abs(tab["tau_se_jk"] / tab["tau_se_ldsc"] - 1)))

    def family(name, sel):
        d = pd.concat([tab[(tab.trait == t) & (tab.model == m) & (tab.annotation == a)] for t, m, a in sel])
        d = d.assign(family=name, q_BH=bh(d["p_one_sided"]))
        d.to_csv(os.path.join(R, f"{name}.tsv"), sep="\t", index=False)
        return d

    ast = ["ast_panukbb", "ast_mvp"]
    fam = {
        "F5": family("F5", [(t, "J7", "SPEC_endothelial_cell") for t in ast]),
        "F5_GGT": family("F5_GGT", [("ggt_panukbb", "J7", a) for a in ["SPEC_hepatocyte", "SPEC_cholangiocyte"]]),
        "F6": family("F6", [(t, "RM_" + x, "RM_" + x) for t in ast for x in ["E122", "E100", "E108", "E066"]]),
        "F6_cond": family("F6_cond", [(t, "RMJ", a) for t in ast for a in ["SPEC_endothelial_cell", "RM_E122"]]),
    }
    # F7: paired AST - ALT on identical SNPs, H22 CT_endothelial_cell model
    f7 = []
    for coh in ["panukbb", "mvp"]:
        a = ctx.run(f"pair_{coh}_ast", "CT_endothelial_cell")["SPEC_endothelial_cell"]
        b = ctx.run(f"pair_{coh}_alt", "CT_endothelial_cell")["SPEC_endothelial_cell"]
        d_del = a["ts_del"] - b["ts_del"]
        D = a["tau_star"] - b["tau_star"]
        se = sqrt(jk_cov(d_del, d_del))
        f7.append(dict(cohort=coh, tau_star_ast=a["tau_star"], se_ast=a["tau_star_se"], tau_star_alt=b["tau_star"],
                       se_alt=b["tau_star_se"], jk_corr=jk_cov(a["ts_del"], b["ts_del"]) / (a["tau_star_se"] * b["tau_star_se"]),
                       D=D, se_D=se, z_D=D / se, p_two_sided=p_two(D / se)))
    f7 = pd.DataFrame(f7)
    f7["q_BH"] = bh(f7["p_two_sided"])
    f7.to_csv(os.path.join(R, "F7.tsv"), sep="\t", index=False)
    summ = {"max_rel_diff_jk_vs_ldsc_coef_se": chk,
            **{k: {"n": len(v), "n_q_lt_0.05": int((v.q_BH < 0.05).sum())} for k, v in fam.items()},
            "F7": {"n": len(f7), "n_q_lt_0.05": int((f7.q_BH < 0.05).sum())}}
    json.dump(summ, open(os.path.join(R, "summary.json"), "w"), indent=2)
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:4])
