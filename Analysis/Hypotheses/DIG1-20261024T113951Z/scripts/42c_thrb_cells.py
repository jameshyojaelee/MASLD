"""DIG1 job 42c (exploratory; follows job 42b). GSE298719 Flex, 13 donors (6 healthy, 7 MASLD).
Which cells express THRB, does hepatocyte THRB change with disease, and inside one liver do hepatocytes
with more THRB RNA carry more thyroid-hormone target RNA (TR4 = CPT1A, ANGPTL4, DIO1, ME1;
doi 10.1371/journal.pone.0240338)? Do injured (KRT8/SQSTM1-high) hepatocytes carry less THRB?
Per-cell values: counts per 10k, log1p. Scores: mean within-donor z.
Within-donor statistic: rho(reference, test) with ranks residualised on log UMI and log genes, minus
the mean rho over 200 random tests matched gene by gene on hepatocyte mean x detection decile bins
(TH and MDB genes excluded; seed 20261024). Exact sign-flip over 13 donors.
Couplings (descriptive, BH over the table): C1 ref THRB, test TR4; C2 ref MDB (KRT8, SQSTM1), test
THRB; C3 ref MDB, test TR4; C4 ref central zonation score, test THRB; C5 ref THRB, each TR4 gene.
Disease effect: OLS donor mean ~ disease + sex, exact relabeling (1,716 splits), per cell type.
No primary family: detection of THRB in Flex hepatocytes is unknown before this run.
Usage: python -I 42c_thrb_cells.py <new_out_dir>
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
genes = pd.Index(A.var_names); g = list(genes)
TR4 = ["CPT1A", "ANGPTL4", "DIO1", "ME1"]
TH = ["THRB", "THRA", "DIO1", "DIO2", "DIO3", "SLC10A1", "SLC16A2", "SLC16A10", "CPT1A", "ANGPTL4", "ME1", "THRSP", "PPARGC1A", "LDLR", "CYP7A1", "SHBG"]
MDB = ["KRT8", "SQSTM1"]
CEN = ["CYP2E1", "GLUL", "CYP1A2", "OAT", "CYP3A4"]
TH = [x for x in TH if x in genes]; TR4 = [x for x in TR4 if x in genes]
MDB = [x for x in MDB if x in genes]; CEN = [x for x in CEN if x in genes]
X = sparse.csr_matrix(A.X).astype(np.float32)
umi = np.asarray(X.sum(axis=1)).ravel()
lab = A.obs.label.astype(str).to_numpy(); gsm = A.obs.gsm.astype(str).to_numpy()
L = X.multiply(1e4 / umi[:, None]).tocsr(); L.data = np.log1p(L.data)
setg = sorted(set(TH) | set(MDB) | set(CEN))
cidx = [genes.get_loc(x) for x in setg]
Ls = L[:, cidx].toarray(); Xs = X[:, cidx].toarray(); gi = {x: i for i, x in enumerate(setg)}

rows, share = [], []
for d in dq.gsm:
    m = gsm == d
    tot = {x: Xs[m, gi[x]].sum() for x in TH}
    for t in np.unique(lab[m]):
        mt = m & (lab == t)
        rows.append(dict(gsm=d, label=t, n_cells=int(mt.sum()), **{x: float(Ls[mt, gi[x]].mean()) for x in TH},
                         **{f"det_{x}": float((Xs[mt, gi[x]] > 0).mean()) for x in TH}))
        share.append(dict(gsm=d, label=t, frac_cells=mt.sum() / m.sum(), **{f"frac_{x}_counts": (Xs[mt, gi[x]].sum() / v if v > 0 else np.nan) for x, v in tot.items()}))
ct = pd.DataFrame(rows).merge(dq[["gsm", "id", "disease", "sex_f"]], on="gsm"); ct.to_csv(OUT / "celltype_donor_means.tsv", sep="\t", index=False)
sh = pd.DataFrame(share).merge(dq[["gsm", "disease"]], on="gsm"); sh.to_csv(OUT / "count_share.tsv", sep="\t", index=False)

n = len(dq)
SIGNS = np.array(list(product([1, -1], repeat=n)))
RELAB = []
for ix in combinations(range(n), int(dq.disease.sum())):
    v = np.zeros(n, int); v[list(ix)] = 1; RELAB.append(v)
RELAB = np.array(RELAB)
signflip = lambda v: np.mean(np.abs((SIGNS * np.asarray(v, float)).mean(1)) >= abs(np.mean(v)) - 1e-12)


def dis(y, dis_v, sex_v):
    y = np.asarray(y, float); nn = len(y)
    beta = lambda v: np.linalg.lstsq(np.column_stack([np.ones(nn), v, sex_v]).astype(float), y, rcond=None)[0][1]
    b = beta(dis_v)
    return b, np.mean(np.abs([beta(v) for v in RELAB]) >= abs(b) - 1e-12)


de = []
for t, o in ct.groupby("label"):
    o = o.set_index("gsm").reindex(dq.gsm)
    if o.n_cells.isna().any() or (o.n_cells < 20).any():
        continue
    for x in TH:
        b, p = dis(o[x], dq.disease.to_numpy(), dq.sex_f.to_numpy())
        de.append(dict(label=t, gene=x, beta=b, p_exact=p, mean_healthy=o[x][dq.disease.to_numpy() == 0].mean(), mean_masld=o[x][dq.disease.to_numpy() == 1].mean()))
de = pd.DataFrame(de); de["bh"] = de.groupby("label").p_exact.transform(stats.false_discovery_control)
de.to_csv(OUT / "disease_effect_by_celltype.tsv", sep="\t", index=False)

# hepatocyte within-donor couplings
H = lab == "Hepatocytes"
XH = X[H]; LH = L[H]; gH = gsm[H]; umiH = umi[H]; ngH = np.asarray((XH > 0).sum(axis=1)).ravel()
mean_all = np.asarray(LH.mean(axis=0)).ravel(); det_all = np.asarray((XH > 0).mean(axis=0)).ravel()
keep = det_all >= 0.01
binm = pd.qcut(pd.Series(mean_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bind = pd.qcut(pd.Series(det_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bins = pd.Series(list(zip(binm, bind)), index=np.array(g)[keep])
pool = {b: [x for x in idx if x not in set(setg)] for b, idx in bins.groupby(bins).groups.items()}
pd.DataFrame(dict(gene=setg, hep_mean=[mean_all[genes.get_loc(x)] for x in setg], hep_detection=[det_all[genes.get_loc(x)] for x in setg],
                  binned=[x in bins.index for x in setg])).to_csv(OUT / "hep_gene_detection.tsv", sep="\t", index=False)
TESTS = {"THRB": ["THRB"], "TR4": TR4, **{x: [x] for x in TR4}}
TESTS = {k: [x for x in v if x in bins.index] for k, v in TESTS.items()}
TESTS = {k: v for k, v in TESTS.items() if v}
NULL = {t: [[rng.choice(pool[bins[x]]) for x in v] for _ in range(N_NULL)] for t, v in TESTS.items()}
COUP = [("C1", "THRB", "TR4"), ("C2", "MDB", "THRB"), ("C3", "MDB", "TR4"), ("C4", "central", "THRB")] + [("C5", "THRB", x) for x in TR4]
COUP = [c for c in COUP if c[2] in TESTS and (c[1] in ("MDB", "central") or c[1] in TESTS)]
need = sorted(set(setg) | set(x for r in NULL.values() for s in r for x in s))
Ld = LH[:, [genes.get_loc(x) for x in need]].toarray(); gj = {x: i for i, x in enumerate(need)}


def resid_rank(M, Zc):
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return R - Zc @ np.linalg.lstsq(Zc, R, rcond=None)[0]


def cor_cols(a, M):
    a = (a - a.mean()) / a.std(); sd = M.std(0); sd[sd == 0] = np.nan
    return ((M - M.mean(0)) / sd * a[:, None]).mean(0)


wd = []
for d in dq.gsm:
    m = gH == d
    D = Ld[m]; sd = D.std(0); sd[sd == 0] = np.nan
    Zg = (D - D.mean(0)) / sd
    sc = lambda gs: np.nanmean(Zg[:, [gj[x] for x in gs]], 1)
    Zc = np.column_stack([np.ones(m.sum()), stats.rankdata(np.log(umiH[m])), stats.rankdata(np.log(ngH[m]))])
    refs = {"MDB": sc(MDB), "central": sc(CEN), **{k: sc(v) for k, v in TESTS.items()}}
    for cid, rname, tname in COUP:
        ref = refs[rname]
        if not np.isfinite(ref).all() or np.nanstd(ref) == 0:
            continue
        ri = resid_rank(ref[:, None], Zc)[:, 0]
        ro = resid_rank(sc(TESTS[tname])[:, None], Zc)[:, 0]
        Rn = resid_rank(np.column_stack([np.nan_to_num(sc(r)) for r in NULL[tname]]), Zc)
        wd.append(dict(gsm=d, coupling=cid, ref=rname, test=tname, n_cells=int(m.sum()), rho=np.corrcoef(ri, ro)[0, 1], null_mean=np.nanmean(cor_cols(ri, Rn))))
    print(d, int(m.sum()), flush=True)
wd = pd.DataFrame(wd).merge(dq[["gsm", "id", "disease"]], on="gsm")
wd["excess"] = wd.rho - wd.null_mean
wd.to_csv(OUT / "hep_within_donor_rho.tsv", sep="\t", index=False)
summ = []
for (cid, r, t), o in wd.groupby(["coupling", "ref", "test"]):
    o = o.set_index("gsm").reindex(dq.gsm)
    if o.excess.isna().any():
        continue
    summ.append(dict(coupling=cid, ref=r, test=t, mean_rho=o.rho.mean(), mean_null=o.null_mean.mean(), mean_excess=o.excess.mean(),
                     n_excess_positive=int((o.excess > 0).sum()), p_excess=signflip(o.excess)))
summ = pd.DataFrame(summ); summ["bh"] = stats.false_discovery_control(summ.p_excess)
summ.to_csv(OUT / "hep_within_donor_summary.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 220); pd.set_option("display.max_rows", 500)
print(pd.read_csv(OUT / "hep_gene_detection.tsv", sep="\t").round(4).to_string())
print(summ.round(4).to_string())
print(de[de.label == "Hepatocytes"].round(3).to_string())
print(ct.groupby("label")[TH].mean().round(3).to_string())
print(sh.groupby("label")[[c for c in sh.columns if c.startswith("frac_")]].mean().round(3).to_string())
