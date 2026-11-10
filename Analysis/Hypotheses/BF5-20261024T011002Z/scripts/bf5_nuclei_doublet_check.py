#!/usr/bin/env python
"""BF5 Part B, post hoc doublet check (DEVIATIONS.md 01:23Z). Same as bf5_nuclei.py after removing doublet-like nuclei.
Original question: do hepatocytes express the hepatocyte-made continuum genes broadly or in a subset,
and do gene-positive hepatocyte nuclei carry bile-duct transcripts more often than for comparable genes?"""
import argparse
import json
import os

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem
from statsmodels.stats.contingency_tables import StratifiedTable

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
BF3 = f"{ROOT}/Analysis/Hypotheses/BF3-20261023T154841Z"
BUILD = f"{BF3}/build-20261023T155101Z"
SUMMARY = f"{BF3}/doublet-20261023T232600Z/nuclei_summary.tsv.gz"
THRESH = f"{BF3}/doublet-matched-20261023T232809Z/doublet_matched.json"
ATLAS = f"{ROOT}/Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
CG0 = f"{ROOT}/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/genes.tsv"
TOPTYPE = f"{ROOT}/Analysis/Hypotheses/BF4-20261024T010648Z/out-20261024T011020Z/gene_top_celltype.tsv"
TARGETS = ["AKR1B10", "CA12", "ALDH3A1", "CIDEC", "EEF1A2", "PRAMEF10", "SFRP4"]
PRIMARY = ["AKR1B10", "CA12"]
DIAG = ["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"]
NCTRL, MINUMI, CHUNK = 20, 20, 100_000

# ---------------------------------------------------------------- controls: hepatocyte-made genes of similar rate
z = np.load(f"{BUILD}/donor_type_counts.npz", allow_pickle=True)  # BF3 build output; genes saved as object array
counts, donors, types, genes = z["counts"], z["donors"].astype(str), z["types"].astype(str), z["genes"].astype(str)
dt = pd.read_csv(f"{BUILD}/donors.tsv", sep="\t").set_index("donor").loc[donors]
ok = (dt.n_nuclei >= 500) & dt.diagnosis.isin(DIAG) & dt.sex.isin(["Male", "Female"])
counts, dt = counts[ok.values], dt[ok]
HEP = int(np.where(types == "Hepatocytes")[0][0])
hrate = pd.Series(counts[:, HEP, :].sum(0) / counts[:, HEP, :].sum(), index=genes)
gt = pd.read_csv(CG0, sep="\t")
sig = set(gt.gene_name[gt.signature_gene.fillna(False).astype(bool)].dropna())
tt = pd.read_csv(TOPTYPE, sep="\t")
pool = sorted(set(tt.gene[tt.hep_ratio >= 0.5]) - sig - set(TARGETS))
pool = [x for x in pool if hrate.get(x, 0) > 0]
CTRL = {}
for t in TARGETS:
    d = (np.log(hrate[pool]) - np.log(hrate[t])).abs().sort_values()
    CTRL[t] = list(d.index[:NCTRL])
ALL = TARGETS + sorted(set(sum(CTRL.values(), [])))
gidx = {x: i for i, x in enumerate(genes)}
cols = np.array([gidx[x] for x in ALL])

# ---------------------------------------------------------------- per-nucleus counts for hepatocyte nuclei
nu = pd.read_csv(SUMMARY, sep="\t")
nu = nu[(nu.cell_type == "Hepatocytes") & nu.donor.isin(dt.index)].sort_values("row").reset_index(drop=True)
th = json.load(open(THRESH))
doublet_like = (nu.C50 >= nu.dataset.map(th["t_C50"])) | (nu.F50 >= nu.dataset.map(th["t_F50"]))
print("doublet-like removed by cohort", nu[doublet_like].groupby("dataset").size().to_dict(),
      "of", nu.groupby("dataset").size().to_dict(), flush=True)
f = h5py.File(ATLAS, "r")
R = f["raw"]["X"]
indptr = R["indptr"][:]
G = len(genes)
assert f["raw"]["var"]["_index"].shape[0] == G
want = nu.row.values
Y = np.zeros((len(want), len(ALL)), dtype=np.float32)
for a0 in range(0, len(indptr) - 1, CHUNK):
    b0 = min(a0 + CHUNK, len(indptr) - 1)
    sel = np.where((want >= a0) & (want < b0))[0]
    if len(sel) == 0:
        continue
    lo, hi = indptr[a0], indptr[b0]
    X = sp.csr_matrix((R["data"][lo:hi], R["indices"][lo:hi], indptr[a0:b0 + 1] - lo), shape=(b0 - a0, G))
    Xs = X[want[sel] - a0]
    assert np.allclose(np.asarray(Xs.sum(1)).ravel(), nu.depth.values[sel])  # same nuclei as the BF3 summary
    Y[sel] = Xs[:, cols].toarray()
f.close()
C = pd.DataFrame(Y, columns=ALL)
nu = pd.concat([nu, C], axis=1)
nu = nu[~doublet_like.values].reset_index(drop=True)
nu["diagnosis"] = nu.donor.map(dt.diagnosis)
nu["bil_pos"] = nu.bil >= 1
nu["dbin"] = nu.groupby("donor").depth.transform(lambda s: pd.qcut(s.rank(method="first"), 10, labels=False))
print("hepatocyte nuclei", len(nu), "donors", nu.donor.nunique(), flush=True)

# ---------------------------------------------------------------- B1: subset or broad
b1 = []
for d, x in nu.groupby("donor"):
    idx = {}
    for gname in ALL:
        n = x[gname].values
        if n.sum() < MINUMI:
            continue
        lam = n.sum() / x.depth.sum()
        idx[gname] = (n >= 1).mean() / np.mean(1 - np.exp(-lam * x.depth.values))
    for t in TARGETS:
        if t not in idx:
            continue
        ci = [idx[c] for c in CTRL[t] if c in idx]
        if len(ci) < 5:
            continue
        b1.append(dict(donor=d, dataset=x.dataset.iloc[0], diagnosis=x.diagnosis.iloc[0], gene=t, index=idx[t],
                       control_median_index=float(np.median(ci)), n_controls=len(ci),
                       relative_index=idx[t] / float(np.median(ci)), umis=float(x[t].sum()),
                       share_positive=float((x[t] >= 1).mean())))
b1 = pd.DataFrame(b1)
b1.to_csv(f"{a.out}/b1_donor_index.tsv", sep="\t", index=False)
RES = {"controls": CTRL, "n_hepatocyte_nuclei": int(len(nu)), "n_donors": int(nu.donor.nunique()),
       "n_doublet_like_removed": int(doublet_like.sum())}
RES["B1"] = {t: dict(n_donors=int((b1.gene == t).sum()), median_relative_index=float(b1[b1.gene == t].relative_index.median()),
                     by_diagnosis=b1[b1.gene == t].groupby("diagnosis").relative_index.median().round(3).to_dict())
             for t in TARGETS if (b1.gene == t).any()}


# ---------------------------------------------------------------- B2 and B3: overlap with bile-duct transcripts, identity
def overlap(gname):
    pos = nu[gname] >= 1
    tabs, num, den = [], 0.0, 0.0
    for (_, _), s in nu.assign(pos=pos).groupby(["donor", "dbin"]):
        t = pd.crosstab(s.pos, s.bil_pos).reindex(index=[True, False], columns=[True, False], fill_value=0).values
        if t[0].sum() > 0 and t[1].sum() > 0:
            if t[:, 0].sum() > 0 and t[:, 1].sum() > 0:
                tabs.append(t.astype(float))
            num += len(s) * (s.hep_identity_ratio[s.pos].mean() - s.hep_identity_ratio[~s.pos].mean())
            den += len(s)
    if not tabs:  # a rare control gene in one cohort can have no informative stratum
        return dict(or_mh=np.nan, or_lo=np.nan, or_hi=np.nan, n_strata=0, n_positive=int(pos.sum()),
                    identity_diff=num / den if den else np.nan, donor_median_or=np.nan, n_donors_or=0)
    st = StratifiedTable(tabs)
    lo, hi = st.oddsratio_pooled_confint()
    per = []
    for d, s in nu.assign(pos=pos).groupby("donor"):
        if s.pos.sum() >= 5:
            t = pd.crosstab(s.pos, s.bil_pos).reindex(index=[True, False], columns=[True, False], fill_value=0).values + 0.5
            per.append(t[0, 0] * t[1, 1] / (t[0, 1] * t[1, 0]))
    return dict(or_mh=float(st.oddsratio_pooled), or_lo=float(lo), or_hi=float(hi), n_strata=len(tabs),
                n_positive=int(pos.sum()), identity_diff=num / den if den else np.nan,
                donor_median_or=float(np.median(per)) if per else np.nan, n_donors_or=len(per))


ALL_NU = nu
for scope in ["pooled", "GSE202379", "GSE244832"]:
    nu = ALL_NU if scope == "pooled" else ALL_NU[ALL_NU.dataset == scope]
    ov = {x: overlap(x) for x in ALL}
    pd.DataFrame(ov).T.to_csv(f"{a.out}/b2_b3_overlap_{scope}.tsv", sep="\t")
    RES[f"B2_B3_{scope}"] = {"n_nuclei": int(len(nu)), "n_donors": int(nu.donor.nunique())}
    for t in PRIMARY:
        c = pd.DataFrame({x: ov[x] for x in CTRL[t]}).T
        RES[f"B2_B3_{scope}"][t] = dict(ov[t], control_or_q95=float(c.or_mh.quantile(0.95)),
                                        control_or_median=float(c.or_mh.median()),
                                        or_above_control_q95=bool(ov[t]["or_mh"] > c.or_mh.quantile(0.95)),
                                        control_identity_median=float(c.identity_diff.median()),
                                        identity_pct_among_controls=float((c.identity_diff < ov[t]["identity_diff"]).mean()))
nu = ALL_NU

# per donor share of hepatocyte nuclei positive for the primary genes, for figures
sh = nu.groupby(["donor", "dataset", "diagnosis"])[PRIMARY].agg(lambda s: (s >= 1).mean()).reset_index()
sh.to_csv(f"{a.out}/donor_share_positive.tsv", sep="\t", index=False)
json.dump(RES, open(f"{a.out}/partB_doublet_check.json", "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in RES.items() if k != "controls"}, indent=1, default=float), flush=True)
