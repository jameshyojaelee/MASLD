#!/usr/bin/env python
"""CGB4 collation and prespecified decision."""
import argparse
import glob
import json
import os
import re

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
a = ap.parse_args()


def read_gsa(path):
    lines = [l for l in open(path) if not l.startswith("#")]
    rows = [re.split(r"\s+", l.strip()) for l in lines]
    return pd.DataFrame(rows[1:], columns=rows[0])


rows = []
for f in sorted(glob.glob(f"{a.out}/*.gsa.out")):
    trait, model = os.path.basename(f).split(".")[:2]
    d = read_gsa(f)
    d["trait"], d["model"] = trait, model
    rows.append(d)
t = pd.concat(rows)
for c in ["NGENES", "BETA", "BETA_STD", "SE", "P"]:
    t[c] = pd.to_numeric(t[c], errors="coerce")
t.to_csv(f"{a.out}/all_gsa.tsv", sep="\t", index=False)

res = {}
m1 = t[t.model.isin(["M1_cont", "M1_dis"])]
for prop in ["M1_cont", "M1_dis"]:
    p = m1[(m1.model == prop) & m1.trait.isin(["ukbb_alt", "ukbb_ast"])].sort_values("trait")
    holm = multipletests(p.P, method="holm")[1]
    for (_, r), h in zip(p.iterrows(), holm):
        rep = m1[(m1.model == prop) & (m1.trait == r.trait.replace("ukbb", "mvp") + "_eur")]
        rp = rep.iloc[0] if len(rep) else None
        res[f"{prop}_{r.trait}"] = dict(beta_std=r.BETA_STD, p=r.P, holm_p=h, n=r.NGENES,
                                        mvp_beta_std=None if rp is None else rp.BETA_STD,
                                        mvp_p=None if rp is None else rp.P,
                                        pursue=bool(h < 0.05 and rp is not None and rp.P < 0.05 and rp.BETA > 0))
res["descriptive_other_traits"] = m1[~m1.trait.isin(["ukbb_alt", "ukbb_ast"])][
    ["trait", "model", "BETA_STD", "P", "NGENES"]].to_dict(orient="records")

m2 = t[t.model == "M2_sys"].copy()
named = []
for tr in ["ukbb_alt", "ukbb_ast"]:
    s = m2[m2.trait == tr].copy()
    s["bh"] = multipletests(s.P, method="fdr_bh")[1]
    rep = m2[m2.trait == tr.replace("ukbb", "mvp") + "_eur"].set_index("VARIABLE")
    for _, r in s.iterrows():
        rp = rep.loc[r.VARIABLE] if r.VARIABLE in rep.index else None
        ok = bool(r.bh < 0.05 and rp is not None and rp.P < 0.05 and rp.BETA > 0)
        if r.bh < 0.1 or ok:
            named.append(dict(trait=tr, system=r.VARIABLE, ngenes=r.NGENES, beta_std=r.BETA_STD, p=r.P, bh=r.bh,
                              mvp_p=None if rp is None else rp.P, named=ok))
res["M2_systems_bh_lt_0.1_or_named"] = named
res["M2_n_named"] = int(sum(n["named"] for n in named))
json.dump(res, open(f"{a.out}/results.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
