"""DIG1 job 39 (exploratory; robustness check for jobs 34, 36b, 37b). GSE298719 Flex hepatocytes, 13 donors.

Hepatocytes with more KRT8/SQSTM1 RNA (MDB score) carry more CCT, heat-shock and chaperone RNA
and less PPARGC1A. If the MDB score and these genes both vary with lobular position, the
couplings could be zonation rather than injury. Test: add zonation scores to the residualisation.
Zonation scores (fixed; mean within-donor z): central = CYP2E1, GLUL, CYP1A2, OAT, CYP3A4;
portal = CPS1, HAL, ASS1, SDS, PCK1, ALDOB [standard markers, recalled].
Statistic as job 34: within donor, rho(MDB, test) with ranks residualised on log UMI and log genes
(depth) or depth + ranks of the central and portal scores (depth+zonation), minus the mean over 200
random sets matched gene by gene on hepatocyte mean x detection deciles (set, MDB and zonation
genes excluded; seed 20261024). Exact sign-flip over 13 donors.
Primary family (Holm over 3), depth+zonation: Z1 CCT, Z2 HSF1 set (HSPA1A HSPA1B DNAJB1 HSPH1 BAG3
HSPA4L HSP90AA1), Z3 HSPB1.
Descriptive (BH): all tests, both adjustments; the same statistic within central-score tertiles
(depth only); rho of MDB and each test with the zonation scores.
Usage: python -I 39_zonation_check_cells.py <new_out_dir>
"""
import sys
from itertools import product
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
ZON = {"central": ["CYP2E1", "GLUL", "CYP1A2", "OAT", "CYP3A4"], "portal": ["CPS1", "HAL", "ASS1", "SDS", "PCK1", "ALDOB"]}
MDB = ["KRT8", "KRT18", "SQSTM1"]
TESTS = {
    "CCT": ["TCP1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8"],
    "HSF1": ["HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "BAG3", "HSPA4L", "HSP90AA1"],
    "chaperone_network": ["HSPA8", "HSP90AA1", "HSP90AB1", "STIP1", "AHSA1", "CHORDC1", "DNAJA1", "HSPA4", "HSPA4L", "HSPH1"],
    "tubulin": ["TUBA1A", "TUBA1B", "TUBA1C", "TUBA4A", "TUBB", "TUBB2A", "TUBB4B", "TUBB6"],
    "actin": ["ACTB", "ACTG1"],
    "CI": [x for x in g if x.startswith("NDUF") and x[4:5] in "ABSCV" and x[5:6].isdigit() and x != "NDUFA4"],
    "mito_chaperones": ["HSPD1", "HSPE1"],
    "HSPB1": ["HSPB1"], "HSPA8": ["HSPA8"], "HSP90AA1": ["HSP90AA1"], "HSPA1A": ["HSPA1A"], "PPARGC1A": ["PPARGC1A"], "CPT1A": ["CPT1A"],
}
TESTS = {k: [x for x in v if x in genes] for k, v in TESTS.items()}
TESTS = {k: v for k, v in TESTS.items() if v}
ZON = {k: [x for x in v if x in genes] for k, v in ZON.items()}
MDB = [x for x in MDB if x in genes]
pd.DataFrame([dict(set=k, n=len(v), genes=",".join(v)) for k, v in {**TESTS, **ZON, "MDB": MDB}.items()]).to_csv(OUT / "sets_used.tsv", sep="\t", index=False)
excl = set(x for v in TESTS.values() for x in v) | set(x for v in ZON.values() for x in v) | set(MDB)
mean_all = np.asarray(L.mean(axis=0)).ravel(); det_all = np.asarray((X > 0).mean(axis=0)).ravel()
keep = det_all >= 0.01
binm = pd.qcut(pd.Series(mean_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bind = pd.qcut(pd.Series(det_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bins = pd.Series(list(zip(binm, bind)), index=np.array(g)[keep])
pool = {b: [x for x in idx if x not in excl] for b, idx in bins.groupby(bins).groups.items()}
NULL = {t: [[rng.choice(pool[bins[x]]) for x in v if x in bins.index] for _ in range(N_NULL)] for t, v in TESTS.items()}
need = sorted(excl | set(x for r in NULL.values() for s in r for x in s))
Ld = L[:, [genes.get_loc(x) for x in need]].toarray(); gi = {x: i for i, x in enumerate(need)}


def resid_rank(M, Z):
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


def cor_cols(a, M):
    a = (a - a.mean()) / a.std(); sd = M.std(0); sd[sd == 0] = np.nan
    return ((M - M.mean(0)) / sd * a[:, None]).mean(0)


rows, zc = [], []
for d in dq.gsm:
    m = gsm == d
    D = Ld[m]; sd = D.std(0); sd[sd == 0] = np.nan
    Zg = (D - D.mean(0)) / sd
    sc = lambda gs: np.nanmean(Zg[:, [gi[x] for x in gs]], 1)
    mdb, cen, por = sc(MDB), sc(ZON["central"]), sc(ZON["portal"])
    base = np.column_stack([np.ones(m.sum()), stats.rankdata(np.log(umi[m])), stats.rankdata(np.log(ngene[m]))])
    zadj = np.column_stack([base, stats.rankdata(cen), stats.rankdata(por)])
    zc.append(dict(gsm=d, test="MDB", rho_central=stats.spearmanr(mdb, cen)[0], rho_portal=stats.spearmanr(mdb, por)[0]))
    for t, v in TESTS.items():
        zc.append(dict(gsm=d, test=t, rho_central=stats.spearmanr(sc(v), cen)[0], rho_portal=stats.spearmanr(sc(v), por)[0]))
    tert = pd.qcut(stats.rankdata(cen), 3, labels=False)
    strata = [("all", "depth", np.ones(m.sum(), bool), base), ("all", "depth+zonation", np.ones(m.sum(), bool), zadj)]
    strata += [(f"central_tertile_{k + 1}", "depth", tert == k, base) for k in range(3)]
    for stratum, adj, sel, Zfull in strata:
        Zr = Zfull[sel]
        ri = resid_rank(mdb[sel][:, None], Zr)[:, 0]
        for t, v in TESTS.items():
            ro = resid_rank(sc(v)[sel][:, None], Zr)[:, 0]
            Rn = resid_rank(np.column_stack([sc(r)[sel] for r in NULL[t]]), Zr)
            rows.append(dict(gsm=d, stratum=stratum, adjust=adj, test=t, n_cells=int(sel.sum()), rho=np.corrcoef(ri, ro)[0, 1],
                             null_mean=np.nanmean(cor_cols(ri, Rn))))
    print(d, int(m.sum()), flush=True)
wd = pd.DataFrame(rows).merge(dq[["gsm", "id", "disease"]], on="gsm")
wd["excess"] = wd.rho - wd.null_mean
wd.to_csv(OUT / "within_donor_rho.tsv", sep="\t", index=False)
zc = pd.DataFrame(zc); zc.to_csv(OUT / "zonation_correlations.tsv", sep="\t", index=False)
n = len(dq)
SIGNS = np.array(list(product([1, -1], repeat=n)))
signflip = lambda v: np.mean(np.abs((SIGNS * np.asarray(v, float)).mean(1)) >= abs(np.mean(v)) - 1e-12)
summ = []
for (st, adj, t), o in wd.groupby(["stratum", "adjust", "test"]):
    o = o.set_index("gsm").reindex(dq.gsm)
    summ.append(dict(stratum=st, adjust=adj, test=t, mean_rho=o.rho.mean(), mean_null=o.null_mean.mean(), mean_excess=o.excess.mean(),
                     n_excess_positive=int((o.excess > 0).sum()), p_excess=signflip(o.excess)))
summ = pd.DataFrame(summ)
summ["bh"] = summ.groupby(["stratum", "adjust"]).p_excess.transform(stats.false_discovery_control)
summ.to_csv(OUT / "within_donor_summary.tsv", sep="\t", index=False)
prim = summ[(summ.stratum == "all") & (summ.adjust == "depth+zonation") & summ.test.isin(["CCT", "HSF1", "HSPB1"])].copy()
prim["id"] = prim.test.map({"CCT": "Z1", "HSF1": "Z2", "HSPB1": "Z3"}); prim = prim.sort_values("id").reset_index(drop=True)
o = np.argsort(prim.p_excess.to_numpy()); ph = np.maximum.accumulate(prim.p_excess.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 220); pd.set_option("display.max_rows", 300)
print(prim.round(4).to_string())
print(summ.pivot_table(index="test", columns=["stratum", "adjust"], values="mean_excess").round(3).to_string())
print(summ.pivot_table(index="test", columns=["stratum", "adjust"], values="n_excess_positive").to_string())
print(zc.groupby("test")[["rho_central", "rho_portal"]].mean().round(3).to_string())
