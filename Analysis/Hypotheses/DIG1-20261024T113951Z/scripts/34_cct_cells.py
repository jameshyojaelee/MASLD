"""DIG1 job 34 (exploratory; follows job 33). GSE298719 Flex, 13 donors (6 healthy, 7 MASLD).

CCT chaperonin protein is higher in MASLD liver and rises with ballooning (jobs 30-32). Which cells
carry CCT RNA, and do hepatocytes themselves raise it with injury or disease?
Per-cell values: counts per 10k, log1p. Sets (fixed): CCT (TCP1, CCT2-8, CCT6A), prefoldin, tubulin,
actin, VIM/FLNA/MYH9, proliferation (MKI67, TOP2A, PCNA, MCM2, MCM3, MCM6).
Injury score in hepatocytes: MDB = KRT8, SQSTM1 (KRT18 not on the panel), z within donor.
Primary family (Holm over 2):
  S1  hepatocytes, within-donor rho(MDB, CCT) minus the mean over 200 matched random sets (hepatocyte
      mean and detection deciles, set genes excluded; seed 20261024), ranks residualised on log UMI and
      log genes; exact sign-flip over 13 donors.
  S2  hepatocytes, disease effect on donor-mean CCT (log1p CP10k), OLS ~ disease + sex, exact
      relabeling (1,716 splits).
Descriptive: donor-mean of every set per cell type; share of total CCT counts by cell type per
donor; disease effect per cell type and set (exact relabeling; BH over the table); same S1 statistic
for every set.
Usage: python -I 34_cct_cells.py <new_out_dir>
"""
import sys
from itertools import combinations, product
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
rng = np.random.default_rng(20261024)
N_NULL = 200

A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str})
dq = dq[dq.included].sort_values("gsm").reset_index(drop=True)
A = A[A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
genes = pd.Index(A.var_names)
SETS = {
    "CCT": ["TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8"],
    "prefoldin": ["PFDN1", "PFDN2", "PFDN4", "PFDN5", "PFDN6", "VBP1"],
    "tubulin": ["TUBA1A", "TUBA1B", "TUBA1C", "TUBA4A", "TUBB", "TUBB2A", "TUBB4B", "TUBB6"],
    "actin": ["ACTB", "ACTG1"],
    "VIM_FLNA_MYH9": ["VIM", "FLNA", "MYH9"],
    "proliferation": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM6"],
}
SETS = {k: [x for x in v if x in genes] for k, v in SETS.items()}
SETS = {k: v for k, v in SETS.items() if v}
INJ = [x for x in ["KRT8", "KRT18", "SQSTM1"] if x in genes]
pd.DataFrame([dict(set=k, n=len(v), genes=",".join(v)) for k, v in {**SETS, "MDB": INJ}.items()]).to_csv(OUT / "sets_used.tsv", sep="\t", index=False)
setg = sorted(set(x for v in SETS.values() for x in v) | set(INJ))

X = sparse.csr_matrix(A.X).astype(np.float32)
umi = np.asarray(X.sum(axis=1)).ravel()
lab = A.obs.label.astype(str).to_numpy()
gsm = A.obs.gsm.astype(str).to_numpy()
L = X.multiply(1e4 / umi[:, None]).tocsr(); L.data = np.log1p(L.data)
cidx = [genes.get_loc(x) for x in setg]
Ls = L[:, cidx].toarray(); Xs = X[:, cidx].toarray()
gi = {x: i for i, x in enumerate(setg)}
sc_all = {s: Ls[:, [gi[x] for x in v]].mean(1) for s, v in SETS.items()}

# cell-type means per donor and share of CCT counts
rows, share = [], []
cct_counts = Xs[:, [gi[x] for x in SETS["CCT"]]].sum(1)
for d in dq.gsm:
    m = gsm == d
    tot = cct_counts[m].sum()
    for t in np.unique(lab[m]):
        mt = m & (lab == t)
        rows.append(dict(gsm=d, label=t, n_cells=int(mt.sum()), **{s: float(sc_all[s][mt].mean()) for s in SETS}))
        share.append(dict(gsm=d, label=t, n_cells=int(mt.sum()), frac_cells=mt.sum() / m.sum(), frac_cct_counts=cct_counts[mt].sum() / tot))
ct = pd.DataFrame(rows).merge(dq[["gsm", "id", "disease", "sex_f"]], on="gsm")
ct.to_csv(OUT / "celltype_donor_means.tsv", sep="\t", index=False)
sh = pd.DataFrame(share).merge(dq[["gsm", "disease"]], on="gsm")
sh.to_csv(OUT / "cct_count_share.tsv", sep="\t", index=False)

n = len(dq)
SIGNS = np.array(list(product([1, -1], repeat=n)))
RELAB = []
for ix in combinations(range(n), int(dq.disease.sum())):
    v = np.zeros(n, int); v[list(ix)] = 1; RELAB.append(v)
RELAB = np.array(RELAB)


def signflip(v):
    v = np.asarray(v, float)
    return np.mean(np.abs((SIGNS * v).mean(1)) >= abs(v.mean()) - 1e-12)


def dis(y, dis_v, sex_v, relab):
    y = np.asarray(y, float)
    nn = len(y)
    beta = lambda v: np.linalg.lstsq(np.column_stack([np.ones(nn), v, sex_v]).astype(float), y, rcond=None)[0][1]
    b = beta(dis_v)
    return b, np.mean(np.abs([beta(v) for v in relab]) >= abs(b) - 1e-12)


# disease effect per cell type (only cell types present with >= 20 cells in every donor)
de = []
for t, o in ct.groupby("label"):
    o = o.set_index("gsm").reindex(dq.gsm)
    if o.n_cells.isna().any() or (o.n_cells < 20).any():
        continue
    for s in SETS:
        b, p = dis(o[s], dq.disease.to_numpy(), dq.sex_f.to_numpy(), RELAB)
        de.append(dict(label=t, set=s, beta=b, p_exact=p, mean_healthy=o[s][dq.disease.to_numpy() == 0].mean(), mean_masld=o[s][dq.disease.to_numpy() == 1].mean()))
de = pd.DataFrame(de)
de["bh"] = stats.false_discovery_control(de.p_exact)
de.to_csv(OUT / "disease_effect_by_celltype.tsv", sep="\t", index=False)

# hepatocytes: within-donor rho with MDB, against matched nulls
H = lab == "Hepatocytes"
XH = X[H]; LH = L[H]; gH = gsm[H]; umiH = umi[H]; ngH = np.asarray((XH > 0).sum(axis=1)).ravel()
mean_all = np.asarray(LH.mean(axis=0)).ravel(); det_all = np.asarray((XH > 0).mean(axis=0)).ravel()
keep = det_all >= 0.01
binm = pd.qcut(pd.Series(mean_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bind = pd.qcut(pd.Series(det_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
kg = np.array(genes)[keep]
bins = pd.Series(list(zip(binm, bind)), index=kg)
pool = {b: [x for x in idx if x not in set(setg)] for b, idx in bins.groupby(bins).groups.items()}
NULL = {s: [[rng.choice(pool[bins[x]]) for x in v if x in bins.index] for _ in range(N_NULL)] for s, v in SETS.items()}
need = sorted(set(setg) | set(x for v in NULL.values() for r in v for x in r))
LdH = LH[:, [genes.get_loc(x) for x in need]].toarray()
gj = {x: i for i, x in enumerate(need)}


def resid_rank(M, Z):
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


wd = []
for d in dq.gsm:
    m = gH == d
    D = LdH[m]; sd = D.std(0); sd[sd == 0] = np.nan
    Zg = (D - D.mean(0)) / sd
    sc = lambda gs: np.nanmean(Zg[:, [gj[x] for x in gs]], 1)
    Zc = np.column_stack([np.ones(m.sum()), stats.rankdata(np.log(umiH[m])), stats.rankdata(np.log(ngH[m]))])
    ri = resid_rank(sc(INJ)[:, None], Zc)[:, 0]
    for s, v in SETS.items():
        ro = resid_rank(sc(v)[:, None], Zc)[:, 0]
        Rn = resid_rank(np.column_stack([sc(r) for r in NULL[s]]), Zc)
        nm = np.mean([np.corrcoef(ri, Rn[:, j])[0, 1] for j in range(Rn.shape[1])])
        wd.append(dict(gsm=d, set=s, n_cells=int(m.sum()), rho=np.corrcoef(ri, ro)[0, 1], null_mean=nm))
    print(d, int(m.sum()), flush=True)
wd = pd.DataFrame(wd).merge(dq[["gsm", "id", "disease"]], on="gsm")
wd["excess"] = wd.rho - wd.null_mean
wd.to_csv(OUT / "hep_within_donor_rho.tsv", sep="\t", index=False)
summ = []
for s, o in wd.groupby("set"):
    o = o.set_index("gsm").reindex(dq.gsm)
    summ.append(dict(set=s, mean_rho=o.rho.mean(), mean_null=o.null_mean.mean(), mean_excess=o.excess.mean(),
                     n_excess_positive=int((o.excess > 0).sum()), p_excess=signflip(o.excess)))
summ = pd.DataFrame(summ); summ["bh"] = stats.false_discovery_control(summ.p_excess)
summ.to_csv(OUT / "hep_within_donor_summary.tsv", sep="\t", index=False)

s1 = summ.set_index("set").loc["CCT"]
hep = ct[ct.label == "Hepatocytes"].set_index("gsm").reindex(dq.gsm)
b2, p2 = dis(hep.CCT, dq.disease.to_numpy(), dq.sex_f.to_numpy(), RELAB)
prim = pd.DataFrame([dict(id="S1", test="hepatocytes: within-donor rho(MDB, CCT) minus matched null", est=s1.mean_excess,
                          n_positive=int(s1.n_excess_positive), p=s1.p_excess),
                     dict(id="S2", test="hepatocytes: disease effect on donor-mean CCT (log1p CP10k)", est=b2, n_positive=np.nan, p=p2)])
o = np.argsort(prim.p.to_numpy()); ph = np.maximum.accumulate(prim.p.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 220)
print(prim.round(4).to_string())
print(ct.groupby("label")[list(SETS)].mean().round(3).to_string())
print(sh.groupby(["label", "disease"]).frac_cct_counts.mean().unstack().round(3).to_string())
print(de.sort_values("p_exact").round(3).to_string())
print(summ.round(3).to_string())
