#!/usr/bin/env python
"""BF3 single-nucleus doublet check (DEVIATIONS.md 2026-10-23T23:26Z). Descriptive."""
import json
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

BUILD, OUT = sys.argv[1], sys.argv[2]
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
ATLAS = f"{ROOT}/Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
PAIR = {ds: f"{ROOT}/data/{ds}/metadata/donor_pairing.csv" for ds in ["GSE202379", "GSE244832"]}
SEED, NPAIR, CHUNK = 20261023, 500, 100_000
BIL, STR = ["KRT7", "KRT19", "EPCAM", "SOX9"], ["COL1A1", "COL1A2", "COL3A1", "DCN"]
IDENT = BIL + STR + ["LUM"]
KEEP_TYPES = ["Hepatocytes", "Cholangiocytes", "Fibroblasts"]
rng = np.random.default_rng(SEED)

# ---------------------------------------------------------------- programs from the build
z = np.load(f"{BUILD}/donor_type_counts.npz", allow_pickle=True)  # own build output; genes saved as object array
counts, donors, types, genes = z["counts"], z["donors"].astype(str), z["types"].astype(str), z["genes"].astype(str)
dt = pd.read_csv(f"{BUILD}/donors.tsv", sep="\t").set_index("donor").loc[donors]
ok = (dt.n_nuclei >= 500).values
counts, dt = counts[ok], dt[ok]
big = [c for c, t in enumerate(types) if dt[f"nuclei_{t}"].sum() >= 1000]
btypes = types[big].tolist()
rate = counts[:, big].sum(0) / counts[:, big].sum((0, 2))[:, None]
gidx = {g: i for i, g in enumerate(genes)}


def program(t):
    c = btypes.index(t)
    spec = rate[c] / (np.delete(rate, c, 0).max(0) + 1e-7)
    cand = (rate[c] >= 1e-5) & ~np.isin(genes, IDENT)
    return np.argsort(-np.where(cand, spec, -1))[:50]


PROG = {"H50": program("Hepatocytes"), "C50": program("Cholangiocytes"), "F50": program("Fibroblasts")}
RES = {"programs": {k: genes[v].tolist() for k, v in PROG.items()}}
COLS = BIL + STR + list(PROG)
ind = np.zeros((len(genes), len(COLS)))
for j, g in enumerate(BIL + STR):
    ind[gidx[g], j] = 1
for j, k in enumerate(PROG):
    ind[PROG[k], len(BIL) + len(STR) + j] = 1
ind = sp.csr_matrix(ind)

# donor-level marker fractions of all nuclear UMIs (all cell types)
Nd = counts.sum((1, 2))
fB = counts[:, :, [gidx[g] for g in BIL]].sum((1, 2)) / Nd
fS = counts[:, :, [gidx[g] for g in STR]].sum((1, 2)) / Nd
frac = pd.DataFrame({"fB": fB, "fS": fS}, index=dt.index)

# ---------------------------------------------------------------- per-nucleus summaries
f = h5py.File(ATLAS, "r")
obs = read_elem(f["obs"])
rmap = {}
for ds, p in PAIR.items():
    dp = pd.read_csv(p, dtype=str)
    for d, srrs in zip(dp.donor_id, dp.rna_srrs):
        for s in srrs.split(";"):
            rmap[s] = f"{ds}_{d}"
run = obs["sample"].astype(str).values
donor = pd.Series(run).map(rmap).values
keep = (obs.dataset.isin(list(PAIR)).values & obs.cell_type.astype(str).isin(KEEP_TYPES).values
        & pd.Series(donor).isin(dt.index).values)
R = f["raw"]["X"]
indptr = R["indptr"][:]
G = len(genes)
rows, depth, ngen, Y = [], [], [], []
for a in range(0, len(obs), CHUNK):
    b = min(a + CHUNK, len(obs))
    k = keep[a:b]
    if not k.any():
        continue
    lo, hi = indptr[a], indptr[b]
    X = sp.csr_matrix((R["data"][lo:hi], R["indices"][lo:hi], indptr[a:b + 1] - lo), shape=(b - a, G))[k]
    rows.append(np.where(k)[0] + a)
    depth.append(np.asarray(X.sum(1)).ravel())
    ngen.append(np.diff(X.indptr))
    Y.append((X @ ind).toarray())
f.close()
rows = np.concatenate(rows)
df = pd.DataFrame(np.vstack(Y), columns=COLS)
df.insert(0, "row", rows)
df["depth"], df["n_genes"] = np.concatenate(depth), np.concatenate(ngen)
df["run"], df["donor"] = run[rows], donor[rows]
df["dataset"] = obs.dataset.astype(str).values[rows]
df["cell_type"] = obs.cell_type.astype(str).values[rows]
df["conf"] = pd.to_numeric(obs.cell_type_conf.values[rows], errors="coerce")
assert np.allclose(df.depth, np.round(df.depth))
df["bil"], df["str"] = df[BIL].sum(1), df[STR].sum(1)
df["bil_genes"] = (df[BIL] > 0).sum(1)
df["str_genes"] = (df[STR] > 0).sum(1)

# ---------------------------------------------------------------- donor references and nucleus metrics
hep = df.cell_type == "Hepatocytes"
clean = hep & (df.bil == 0) & (df["str"] == 0)
ref = df[clean].assign(hf=lambda d: d.H50 / d.depth).groupby("donor").agg(ref_depth=("depth", "median"),
                                                                           ref_hf=("hf", "median"))


def score(d):
    d = d.join(ref, on="donor")
    d["depth_ratio"] = d.depth / d.ref_depth
    d["hep_identity_ratio"] = (d.H50 / d.depth) / d.ref_hf
    d["C50_frac"], d["F50_frac"] = d.C50 / d.depth, d.F50 / d.depth
    return d


df = score(df)
q_depth = df[clean].groupby("dataset").depth_ratio.quantile(0.95)
q_ident = df[clean].groupby("dataset").hep_identity_ratio.quantile(0.05)

# synthetic doublets within run
SUM = COLS + ["depth"]
syn = []
for (ds, r_), d in df.groupby(["dataset", "run"]):
    h = d[d.cell_type == "Hepatocytes"]
    for other, lab in [("Cholangiocytes", "synthetic_hep+chol"), ("Fibroblasts", "synthetic_hep+fib")]:
        o = d[d.cell_type == other]
        if len(h) == 0 or len(o) == 0:
            continue
        i, j = rng.integers(0, len(h), NPAIR), rng.integers(0, len(o), NPAIR)
        s = pd.DataFrame(h[SUM].values[i] + o[SUM].values[j], columns=SUM)
        s["donor"], s["dataset"], s["group"] = h.donor.values[i], ds, lab
        s["conf"] = np.nan
        syn.append(s)
syn = score(pd.concat(syn, ignore_index=True))

df["group"] = ""
g = df.group.copy()
g[hep & clean] = "hep_clean"
g[hep & (df.bil == 1)] = "hep_B1"
g[hep & (df.bil >= 2)] = "hep_B2+"
df["group"] = g
df["sgroup"] = np.select([hep & (df["str"] == 1), hep & (df["str"] >= 2)], ["hep_S1", "hep_S2+"], "")
df.loc[df.cell_type == "Cholangiocytes", "group"] = "cholangiocyte"
df.loc[df.cell_type == "Fibroblasts", "group"] = "fibroblast"


def summarize(d, ds):
    return dict(n=int(len(d)), median_depth_ratio=float(d.depth_ratio.median()),
                median_hep_identity_ratio=float(d.hep_identity_ratio.median()),
                median_C50_frac=float(d.C50_frac.median()), median_F50_frac=float(d.F50_frac.median()),
                median_conf=float(d.conf.median()) if d.conf.notna().any() else None,
                frac_depth_high=float((d.depth_ratio > q_depth[ds]).mean()),
                frac_identity_low=float((d.hep_identity_ratio < q_ident[ds]).mean()),
                frac_both=float(((d.depth_ratio > q_depth[ds]) & (d.hep_identity_ratio < q_ident[ds])).mean()),
                frac_2plus_distinct_bil=float((d.bil_genes >= 2).mean()) if "bil_genes" in d else None)


summ = {}
for ds in sorted(df.dataset.unique()):
    x = df[df.dataset == ds]
    s = {gname: summarize(x[x.group == gname], ds) for gname in ["hep_clean", "hep_B1", "hep_B2+", "cholangiocyte", "fibroblast"]}
    s.update({gname: summarize(x[x.sgroup == gname], ds) for gname in ["hep_S1", "hep_S2+"]})
    y = syn[syn.dataset == ds]
    s.update({gname: summarize(y[y.group == gname], ds) for gname in ["synthetic_hep+chol", "synthetic_hep+fib"]})
    s["n_hep"] = int((x.cell_type == "Hepatocytes").sum())
    summ[ds] = s
RES["groups"] = summ
RES["thresholds"] = {"depth_ratio_q95_clean": q_depth.to_dict(), "hep_identity_ratio_q05_clean": q_ident.to_dict()}

# ---------------------------------------------------------------- ambient calibration per donor
h = df[hep].join(frac, on="donor")
amb = []
for d_, x in h.groupby("donor"):
    D = x.depth.sum()
    rho = x["str"].sum() / (D * x.fS.iloc[0])
    lamB = rho * x.depth * x.fB.iloc[0]
    lamS = rho * x.depth * x.fS.iloc[0]
    p2 = lambda l: 1 - np.exp(-l) * (1 + l)
    b2 = x.bil >= 2
    amb.append(dict(donor=d_, dataset=x.dataset.iloc[0], diagnosis=dt.loc[d_, "diagnosis"], n_hep=len(x), rho=rho,
                    obs_bil=x.bil.sum(), exp_bil=lamB.sum(), obs_B2=int(b2.sum()), exp_B2=p2(lamB).sum(),
                    obs_S2=int((x["str"] >= 2).sum()), exp_S2=p2(lamS).sum(),
                    B2_median_depth_ratio=x.depth_ratio[b2].median() if b2.any() else np.nan,
                    B2_median_identity_ratio=x.hep_identity_ratio[b2].median() if b2.any() else np.nan,
                    B2_share=b2.mean()))
amb = pd.DataFrame(amb)
amb.to_csv(f"{OUT}/donor_ambient.tsv", sep="\t", index=False)
pool = amb.groupby("dataset")[["obs_bil", "exp_bil", "obs_B2", "exp_B2", "obs_S2", "exp_S2"]].sum()
pool["bil_obs_over_exp"] = pool.obs_bil / pool.exp_bil
pool["B2_obs_over_exp"] = pool.obs_B2 / pool.exp_B2
pool["S2_obs_over_exp"] = pool.obs_S2 / pool.exp_S2
RES["ambient_pooled"] = pool.to_dict(orient="index")
RES["ambient_donor_summary"] = dict(
    n_donors=int(len(amb)),
    median_rho=float(amb.rho.median()),
    median_bil_obs_over_exp=float((amb.obs_bil / amb.exp_bil).median()),
    n_donors_bil_obs_over_exp_gt1=int((amb.obs_bil > amb.exp_bil).sum()),
    median_B2_obs_over_exp=float((amb.obs_B2 / amb.exp_B2).median()),
    n_donors_B2_obs_over_exp_gt1=int((amb.obs_B2 > amb.exp_B2).sum()),
    median_S2_obs_over_exp=float((amb.obs_S2 / amb.exp_S2).median()),
    n_donors_B2_median_depth_ratio_gt1=int((amb.B2_median_depth_ratio > 1).sum()),
    median_B2_depth_ratio_over_donors=float(amb.B2_median_depth_ratio.median()),
    median_B2_identity_ratio_over_donors=float(amb.B2_median_identity_ratio.median()))
RES["B2_share_by_diagnosis"] = amb.groupby(["dataset", "diagnosis"]).B2_share.median().unstack().to_dict(orient="index")

df[["row", "run", "donor", "dataset", "cell_type", "conf", "depth", "n_genes", "bil", "str", "bil_genes",
        "str_genes", "H50", "C50", "F50", "depth_ratio", "hep_identity_ratio", "group", "sgroup"]].to_csv(
        f"{OUT}/nuclei_summary.tsv.gz", sep="\t", index=False)
json.dump(RES, open(f"{OUT}/doublet.json", "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in RES.items() if k != "programs"}, indent=1, default=float))
