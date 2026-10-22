#!/usr/bin/env python
"""CGB8: gene properties of colocalized (G) vs disease-state (D) genes in the H57 universe.
Follows PRESPEC_CGB8.md (sha256 in PRESPEC_CGB8.md.sha256)."""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

SEED = 20261022
NMATCH = 1000
FEATS = {"F1_loeuf": "loeuf", "F2_abc_enh": "abc_liver_n_enh", "F3_length": "gene_length",
         "F4_density": "gene_density_1mb", "F5_eqtl_poscontrol": "gtex_liver_n_p1e5"}


def prepare(g, eqtl=None):
    if eqtl is not None:
        e = pd.read_csv(eqtl, sep="\t")[["g", "n_p1e5"]].rename(columns={"n_p1e5": "eqtl_n_p1e5"})
        g = g.merge(e, on="g", how="left")
    else:
        g["eqtl_n_p1e5"] = g.gtex_liver_n_p1e5
    deg = u.deg.fillna(False).astype(bool)
    gen = u.coloc_gen.astype(bool)
    u["set"] = np.select([gen & ~deg, deg & ~gen, gen & deg], ["G", "D", "GD"], "rest")
    u["F1_loeuf"] = u.loeuf
    u["F2_abc_enh"] = np.log1p(u.abc_liver_n_enh)
    u["F3_length"] = np.log(u.gene_length)
    u["F4_density"] = u.gene_density_1mb
    u["F5_eqtl_poscontrol"] = np.log1p(u.eqtl_n_p1e5)
    for f in FEATS:
        u[f] = (u[f] - u[f].mean()) / u[f].std(ddof=1)
    return u


def logit(df, a, b, feats):
    d = df[df.set.isin([a, b])].dropna(subset=feats + ["AveExpr"]).copy()
    d["y"] = (d.set == a).astype(int)
    fit = smf.logit("y ~ " + " + ".join(feats) + " + bs(AveExpr, df=4)", data=d).fit(disp=0)
    rows = []
    for f in feats:
        b_, se = fit.params[f], fit.bse[f]
        rows.append(dict(contrast=f"{a}_vs_{b}", feature=f, n_a=int(d.y.sum()), n_b=int((1 - d.y).sum()),
                         log_or=b_, or_=np.exp(b_), lo=np.exp(b_ - 1.96 * se), hi=np.exp(b_ + 1.96 * se),
                         p=fit.pvalues[f]))
    return pd.DataFrame(rows)


def matched(df, label, rng):
    dec = pd.qcut(df.AveExpr.rank(method="first"), 10, labels=False)
    df = df.assign(dec=dec)
    tgt = df[df.set == label]
    need = tgt.dec.value_counts()
    pools = {k: df.index[df.dec == k].values for k in need.index}
    rows = []
    for f in FEATS:
        obs = tgt[f].mean()
        null = np.empty(NMATCH)
        for i in range(NMATCH):
            idx = np.concatenate([rng.choice(pools[k], n, replace=False) for k, n in need.items()])
            null[i] = df.loc[idx, f].mean()
        diff = obs - null.mean()
        p = (1 + np.sum(np.abs(null - null.mean()) >= abs(diff))) / (NMATCH + 1)
        rows.append(dict(set=label, feature=f, n=int(len(tgt)), observed_mean_z=obs, matched_mean_z=null.mean(),
                         diff=diff, null_lo=np.quantile(null, 0.025) - null.mean(),
                         null_hi=np.quantile(null, 0.975) - null.mean(), p_emp=p))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--eqtl", help="per-gene eQTL counts for F5 (Broadaway, the COLOC r2 eQTL source)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=False)
    rng = np.random.default_rng(SEED)
    u = prepare(pd.read_csv(a.genes, sep="\t"), a.eqtl)
    res = {"set_sizes": u.set.value_counts().to_dict(),
           "missing_any_feature": {s: int(u[u.set == s][list(FEATS)].isna().any(axis=1).sum())
                                   for s in ["G", "D", "GD", "rest"]}}
    if (u.set == "G").sum() != 400 or (u.set == "D").sum() != 194:
        res["warning"] = "set sizes differ from 400/194"

    # positive control first
    pc = logit(u, "G", "D", ["F5_eqtl_poscontrol"])
    res["positive_control_F5"] = pc.iloc[0].to_dict()
    if not (pc.iloc[0].log_or > 0 and pc.iloc[0].p < 0.05):
        res["stop"] = "F5 positive control failed"
        json.dump(res, open(f"{a.out}/results.json", "w"), indent=1, default=float)
        print(json.dumps(res, indent=1, default=float))
        sys.exit(3)

    primary = logit(u, "G", "D", ["F1_loeuf", "F2_abc_enh", "F3_length", "F4_density"])
    primary["holm_p"] = multipletests(primary.p, method="holm")[1]
    sec = pd.concat([logit(u, "G", "rest", ["F1_loeuf", "F2_abc_enh", "F3_length", "F4_density"]),
                     logit(u, "D", "rest", ["F1_loeuf", "F2_abc_enh", "F3_length", "F4_density"])])
    mt = pd.concat([matched(u.dropna(subset=list(FEATS) + ["AveExpr"]).reset_index(drop=True), s, rng)
                    for s in ["G", "D"]])
    primary.to_csv(f"{a.out}/primary_G_vs_D.tsv", sep="\t", index=False)
    sec.to_csv(f"{a.out}/secondary_vs_rest.tsv", sep="\t", index=False)
    mt.to_csv(f"{a.out}/matched_check.tsv", sep="\t", index=False)
    u[["g", "gene_name", "set", "AveExpr", "loeuf", "abc_liver_n_enh", "gene_length", "gene_density_1mb",
       "gtex_liver_n_p1e5", "eqtl_n_p1e5"]].to_csv(f"{a.out}/gene_features.tsv", sep="\t", index=False)

    dec = []
    for _, r in primary.iterrows():
        m = mt[mt.feature == r.feature].set_index("set")
        same = np.sign(m.loc["G", "diff"] - m.loc["D", "diff"]) == np.sign(r.log_or)
        matched_p = min(m.loc["G", "p_emp"], m.loc["D", "p_emp"])
        dec.append(dict(feature=r.feature, holm_p=r.holm_p, matched_same_sign=bool(same),
                        matched_min_p=matched_p, passes=bool(r.holm_p < 0.05 and same and matched_p < 0.05)))
    res["decision_rows"] = dec
    res["decision"] = ("pursue: layers differ in gene type" if any(d["passes"] for d in dec)
                       else "no detected difference in constraint or regulatory complexity")
    res["primary"] = primary.to_dict(orient="records")
    json.dump(res, open(f"{a.out}/results.json", "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
