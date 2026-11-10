#!/usr/bin/env python
"""BF3 depth-matched doublet test (DEVIATIONS.md 2026-10-23T23:33Z). Descriptive."""
import json
import sys

import numpy as np
import pandas as pd

DOUB, BUILD, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
SEED, NPAIR = 20261023, 500
rng = np.random.default_rng(SEED)
df = pd.read_csv(f"{DOUB}/nuclei_summary.tsv.gz", sep="\t")
dt = pd.read_csv(f"{BUILD}/donors.tsv", sep="\t").set_index("donor")
amb = pd.read_csv(f"{DOUB}/donor_ambient.tsv", sep="\t").set_index("donor")
SUM = ["depth", "H50", "C50", "F50", "bil", "str"]

hep = df.cell_type == "Hepatocytes"
clean = df.group == "hep_clean"
ref = df[clean].assign(hf=lambda d: d.H50 / d.depth).groupby("donor").agg(ref_depth=("depth", "median"),
                                                                           ref_hf=("hf", "median"))

syn = []
for (ds, r_), d in df.groupby(["dataset", "run"]):
    h = d[hep.loc[d.index]]
    for other, lab in [("Cholangiocytes", "hep+chol"), ("Fibroblasts", "hep+fib")]:
        o = d[d.cell_type == other]
        if len(h) and len(o):
            i, j = rng.integers(0, len(h), NPAIR), rng.integers(0, len(o), NPAIR)
            s = pd.DataFrame(h[SUM].values[i] + o[SUM].values[j], columns=SUM)
            s["donor"], s["dataset"], s["group"] = h.donor.values[i], ds, lab
            syn.append(s)
syn = pd.concat(syn, ignore_index=True)
tC = syn[syn.group == "hep+chol"].groupby("dataset").C50.quantile(0.05)
tF = syn[syn.group == "hep+fib"].groupby("dataset").F50.quantile(0.05)
RES = {"t_C50": tC.to_dict(), "t_F50": tF.to_dict(),
       "synthetic_share_above_t": {
           "hep+chol": float((syn[syn.group == "hep+chol"].C50 >= syn[syn.group == "hep+chol"].dataset.map(tC)).mean()),
           "hep+fib": float((syn[syn.group == "hep+fib"].F50 >= syn[syn.group == "hep+fib"].dataset.map(tF)).mean())}}

h = df[hep].join(ref, on="donor")
h["bin"] = np.floor(np.log2(h.depth / h.ref_depth) / 0.25).astype(int)
h["idr"] = (h.H50 / h.depth) / h.ref_hf
h["C_hi"] = h.C50 >= h.dataset.map(tC)
h["F_hi"] = h.F50 >= h.dataset.map(tF)
h["sgroup"] = h.sgroup.fillna("")
cl = h[h.group == "hep_clean"].groupby(["donor", "bin"]).agg(bg_C=("C_hi", "mean"), bg_F=("F_hi", "mean"),
                                                            bg_idr=("idr", "mean"), n_bg=("idr", "size"))
h = h.join(cl, on=["donor", "bin"])


def doublet_like(x, flag, bg):
    x = x.dropna(subset=[bg])
    return dict(n=int(len(x)), n_unmatched_dropped=None, share_above_t=float(x[flag].mean()),
                matched_background=float(x[bg].mean()),
                doublet_like_fraction=float((x[flag].mean() - x[bg].mean()) / 0.95))


out = {}
for ds, x in h.groupby("dataset"):
    o = {}
    for gname, col, flag, bg in [("hep_B1", "group", "C_hi", "bg_C"), ("hep_B2+", "group", "C_hi", "bg_C"),
                                 ("hep_S1", "sgroup", "F_hi", "bg_F"), ("hep_S2+", "sgroup", "F_hi", "bg_F")]:
        g = x[x[col] == gname]
        o[gname] = doublet_like(g, flag, bg)
        o[gname]["n_unmatched_dropped"] = int(g[bg].isna().sum())
        per = []
        for d_, gd in g.dropna(subset=[bg]).groupby("donor"):
            if len(gd) >= 5:
                per.append((gd[flag].mean() - gd[bg].mean()) / 0.95)
        o[gname]["donor_median_doublet_like_fraction"] = float(np.median(per)) if per else None
        o[gname]["n_donors_ge5"] = len(per)
    b2 = x[(x.group == "hep_B2+")]
    nd = b2[~b2.C_hi].dropna(subset=["bg_idr"])
    o["B2_below_tC"] = dict(n=int((~b2.C_hi).sum()), share_of_B2=float((~b2.C_hi).mean()),
                            identity_ratio_minus_matched=float((nd.idr - nd.bg_idr).mean()),
                            identity_ratio_minus_matched_median=float((nd.idr - nd.bg_idr).median()),
                            share_2plus_distinct_bil=float((b2[~b2.C_hi].bil_genes >= 2).mean()))
    b2d = b2.dropna(subset=["bg_idr"])
    o["B2_all_identity_ratio_minus_matched"] = float((b2d.idr - b2d.bg_idr).mean())
    sc = syn[(syn.group == "hep+chol") & (syn.dataset == ds)].join(ref, on="donor")
    sc["bin"] = np.floor(np.log2(sc.depth / sc.ref_depth) / 0.25).astype(int)
    sc = sc.join(cl[["bg_idr"]], on=["donor", "bin"]).dropna(subset=["bg_idr"])
    o["synthetic_hep+chol_identity_ratio_minus_matched"] = float(((sc.H50 / sc.depth) / sc.ref_hf - sc.bg_idr).mean())
    out[ds] = o
RES["by_dataset"] = out

# ambient expectation scaled: B2+ nuclei below t_C vs Poisson expectation of all B2+ nuclei
pd_ = []
for d_, x in h.groupby("donor"):
    b2 = x.group == "hep_B2+"
    pd_.append(dict(donor=d_, dataset=x.dataset.iloc[0], diagnosis=dt.loc[d_, "diagnosis"], n_hep=len(x),
                    B2=int(b2.sum()), B2_below_tC=int((b2 & ~x.C_hi).sum()), exp_B2=amb.loc[d_, "exp_B2"]))
pd_ = pd.DataFrame(pd_)
pd_["share_B2_below_tC"] = pd_.B2_below_tC / pd_.n_hep
pd_.to_csv(f"{OUT}/donor_B2.tsv", sep="\t", index=False)
pool = pd_.groupby("dataset")[["B2", "B2_below_tC", "exp_B2"]].sum()
pool["below_tC_obs_over_exp"] = pool.B2_below_tC / pool.exp_B2
RES["B2_below_tC_vs_ambient"] = pool.to_dict(orient="index")
RES["n_donors_B2_below_tC_gt_exp"] = int((pd_.B2_below_tC > pd_.exp_B2).sum())
RES["share_B2_below_tC_by_diagnosis_median"] = pd_.groupby(["dataset", "diagnosis"]).share_B2_below_tC.median().unstack().to_dict(orient="index")
RES["share_B2_below_tC_by_diagnosis_n"] = pd_.groupby(["dataset", "diagnosis"]).size().unstack().to_dict(orient="index")
json.dump(RES, open(f"{OUT}/doublet_matched.json", "w"), indent=1, default=float)
print(json.dumps(RES, indent=1, default=float))
