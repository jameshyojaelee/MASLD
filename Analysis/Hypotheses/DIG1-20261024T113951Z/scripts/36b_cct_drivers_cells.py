"""DIG1 job 36b (exploratory; follows job 34). GSE298719 Flex hepatocytes, 13 donors.

Hepatocytes with more KRT8/SQSTM1 RNA carry more CCT RNA (job 34). Which transcriptional program
travels with hepatocyte CCT? Candidates (fixed):
  YAP   CCN1 CCN2 ANKRD1 AMOTL2 AXL NUAK2 LATS2 [Unverified recall of common YAP/TAZ targets]
  HSF1  HSPA1A HSPA1B DNAJB1 HSPH1 BAG3 HSPA4L HSP90AA1 (heat-shock response)
  ductular KRT7 KRT19 SOX9 EPCAM SPP1 (hepatocyte-to-ductular reprogramming)
  MDB   KRT8 KRT18 SQSTM1; CI as job 24; ATF4 targets as job 24.
Statistic: within donor, rho(set, CCT) with ranks residualised on log UMI and log genes, minus the
mean rho of the set with 200 random gene sets matched to CCT (hepatocyte mean and detection
deciles; seed 20261024). Exact sign-flip over 13 donors. Also the same after adding the MDB score
to the residualisation (does the link hold among cells with equal KRT8/SQSTM1?).
Primary family (Holm over 2): Y1 YAP, Y2 HSF1 (excess, no MDB adjustment).
Descriptive: all sets, both adjustments, BH; MASLD vs healthy donor means (exact relabeling).
Usage: python -I 36b_cct_drivers_cells.py <new_out_dir>
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
hep = A[(A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
del A
genes = pd.Index(hep.var_names); g = list(genes)
X = sparse.csr_matrix(hep.X).astype(np.float32)
gsm = hep.obs.gsm.astype(str).to_numpy()
umi = np.asarray(X.sum(axis=1)).ravel(); ngene = np.asarray((X > 0).sum(axis=1)).ravel()
L = X.multiply(1e4 / umi[:, None]).tocsr(); L.data = np.log1p(L.data)
SETS = {
    "CCT": ["TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8"],
    "YAP": ["CCN1", "CCN2", "ANKRD1", "AMOTL2", "AXL", "NUAK2", "LATS2"],
    "HSF1": ["HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "BAG3", "HSPA4L", "HSP90AA1"],
    "ductular": ["KRT7", "KRT19", "SOX9", "EPCAM", "SPP1"],
    "MDB": ["KRT8", "KRT18", "SQSTM1"],
    "CI": [x for x in g if x.startswith("NDUF") and x[4:5] in "ABSCV" and x[5:6].isdigit() and x != "NDUFA4"],
    "ATF4_targets": ["PHGDH", "PSAT1", "PSPH", "SHMT2", "MTHFD2", "MTHFD1L", "ALDH18A1", "PYCR1", "ASNS", "DDIT3", "TRIB3", "CHAC1"],
    "tubulin": ["TUBA1A", "TUBA1B", "TUBA1C", "TUBA4A", "TUBB", "TUBB2A", "TUBB4B", "TUBB6"],
}
SETS = {k: [x for x in v if x in genes] for k, v in SETS.items()}
SETS = {k: v for k, v in SETS.items() if v}
pd.DataFrame([dict(set=k, n=len(v), genes=",".join(v)) for k, v in SETS.items()]).to_csv(OUT / "sets_used.tsv", sep="\t", index=False)
setg = set(x for v in SETS.values() for x in v)
mean_all = np.asarray(L.mean(axis=0)).ravel(); det_all = np.asarray((X > 0).mean(axis=0)).ravel()
keep = det_all >= 0.01
binm = pd.qcut(pd.Series(mean_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bind = pd.qcut(pd.Series(det_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bins = pd.Series(list(zip(binm, bind)), index=np.array(g)[keep])
pool = {b: [x for x in idx if x not in setg] for b, idx in bins.groupby(bins).groups.items()}
NULL = [[rng.choice(pool[bins[x]]) for x in SETS["CCT"] if x in bins.index] for _ in range(N_NULL)]
need = sorted(setg | set(x for r in NULL for x in r))
Ld = L[:, [genes.get_loc(x) for x in need]].toarray(); gi = {x: i for i, x in enumerate(need)}


def resid_rank(M, Z):
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


rows = []
for d in dq.gsm:
    m = gsm == d
    D = Ld[m]; sd = D.std(0); sd[sd == 0] = np.nan
    Zg = (D - D.mean(0)) / sd
    sc = lambda gs: np.nanmean(Zg[:, [gi[x] for x in gs]], 1)
    base = np.column_stack([np.ones(m.sum()), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m]))])
    for adj in ["depth", "depth+MDB"]:
        Zc = base if adj == "depth" else np.column_stack([base, stats.rankdata(sc(SETS["MDB"]))])
        rc = resid_rank(sc(SETS["CCT"])[:, None], Zc)[:, 0]
        Rn = resid_rank(np.column_stack([sc(r) for r in NULL]), Zc)
        for s, v in SETS.items():
            if s == "CCT" or (adj == "depth+MDB" and s == "MDB"):
                continue
            rs = resid_rank(sc(v)[:, None], Zc)[:, 0]
            nm = np.mean([np.corrcoef(rs, Rn[:, j])[0, 1] for j in range(Rn.shape[1])])
            rows.append(dict(gsm=d, adjust=adj, set=s, n_cells=int(m.sum()), rho=np.corrcoef(rs, rc)[0, 1], null_mean=nm))
    print(d, int(m.sum()), flush=True)
wd = pd.DataFrame(rows).merge(dq[["gsm", "id", "disease"]], on="gsm")
wd["excess"] = wd.rho - wd.null_mean
wd.to_csv(OUT / "within_donor_rho.tsv", sep="\t", index=False)
n = len(dq)
SIGNS = np.array(list(product([1, -1], repeat=n)))
signflip = lambda v: np.mean(np.abs((SIGNS * np.asarray(v, float)).mean(1)) >= abs(np.mean(v)) - 1e-12)
summ = []
for (adj, s), o in wd.groupby(["adjust", "set"]):
    o = o.set_index("gsm").reindex(dq.gsm)
    summ.append(dict(adjust=adj, set=s, mean_rho=o.rho.mean(), mean_null=o.null_mean.mean(), mean_excess=o.excess.mean(),
                     n_excess_positive=int((o.excess > 0).sum()), p_excess=signflip(o.excess)))
summ = pd.DataFrame(summ); summ["bh"] = stats.false_discovery_control(summ.p_excess)
summ.to_csv(OUT / "within_donor_summary.tsv", sep="\t", index=False)
prim = summ[(summ.adjust == "depth") & summ.set.isin(["YAP", "HSF1"])].copy()
prim["id"] = prim.set.map({"YAP": "Y1", "HSF1": "Y2"}); prim = prim.sort_values("id").reset_index(drop=True)
o = np.argsort(prim.p_excess.to_numpy()); ph = np.maximum.accumulate(prim.p_excess.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
# donor means, disease effect
RELAB = []
for ix in combinations(range(n), int(dq.disease.sum())):
    v = np.zeros(n, int); v[list(ix)] = 1; RELAB.append(v)
RELAB = np.array(RELAB)
dm = pd.DataFrame([{**dict(gsm=d), **{s: Ld[gsm == d][:, [gi[x] for x in v]].mean() for s, v in SETS.items()}} for d in dq.gsm])
de = []
for s in SETS:
    y = dm[s].to_numpy()
    beta = lambda v: np.linalg.lstsq(np.column_stack([np.ones(n), v, dq.sex_f]).astype(float), y, rcond=None)[0][1]
    b = beta(dq.disease.to_numpy())
    de.append(dict(set=s, beta=b, p_exact=np.mean(np.abs([beta(v) for v in RELAB]) >= abs(b) - 1e-12)))
de = pd.DataFrame(de); de.to_csv(OUT / "disease_effect_donor_means.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 200)
print(prim.round(4).to_string()); print(summ.round(3).to_string()); print(de.round(4).to_string())
