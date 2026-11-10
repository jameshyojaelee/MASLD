"""DIG1 job 37b (exploratory; follows job 36b). GSE298719 Flex, 13 donors (6 healthy, 7 MASLD).

A published NASH biopsy study reports low HSP27 (HSPB1) in ballooned hepatocytes and lower HSPB1
mRNA with ballooning (doi 10.1038/srep22528). Job 36b found the HSF1-target set rising with
KRT8/SQSTM1 in hepatocytes. Gene by gene, in hepatocytes: which heat-shock genes rise with the
KRT8/SQSTM1 injury score and with disease, and which fall? Where is HSPB1 expressed?
Per-cell values: counts per 10k, log1p.
Within-donor statistic per gene: rho(gene, MDB score) with ranks residualised on log UMI and log
genes, minus the mean rho of 200 random genes from the same hepatocyte mean x detection decile bin
(set genes excluded; one draw per bin, seed 20261024). MDB = mean within-donor z of KRT8, SQSTM1.
Exact sign-flip over 13 donors.
Disease effect: OLS donor mean ~ disease + sex, exact relabeling (1,716 splits).
Primary family (Holm over 2), direct check of the published claim, two-sided:
  H1  hepatocytes, HSPB1 within-donor coupling to MDB (excess over matched genes).
  H2  hepatocytes, HSPB1 disease effect on donor mean.
Descriptive (BH within each table): every gene, both statistics; donor-mean log1p CP10k of every
gene by cell type; share of HSPB1 and CRYAB counts by cell type.
Usage: python -I 37b_hsp_genes_cells.py <new_out_dir>
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
HSP = ["HSPB1", "CRYAB", "HSPB8", "HSPA1A", "HSPA1B", "HSPA6", "HSPA8", "HSPA4", "HSPA4L", "HSPH1", "DNAJB1", "DNAJA1",
       "BAG3", "HSP90AA1", "HSP90AB1", "STIP1", "CHORDC1", "AHSA1", "FKBP4", "HSPD1", "HSPE1", "HSF1", "UCHL1", "CPT1A",
       "PPARGC1A", "HSPA5"]
HSP = [x for x in HSP if x in genes]
INJ = [x for x in ["KRT8", "KRT18", "SQSTM1"] if x in genes]
pd.DataFrame(dict(gene=HSP + INJ, role=["test"] * len(HSP) + ["MDB"] * len(INJ))).to_csv(OUT / "genes_used.tsv", sep="\t", index=False)
setg = sorted(set(HSP) | set(INJ))

X = sparse.csr_matrix(A.X).astype(np.float32)
umi = np.asarray(X.sum(axis=1)).ravel()
lab = A.obs.label.astype(str).to_numpy()
gsm = A.obs.gsm.astype(str).to_numpy()
L = X.multiply(1e4 / umi[:, None]).tocsr(); L.data = np.log1p(L.data)
cidx = [genes.get_loc(x) for x in setg]
Ls = L[:, cidx].toarray(); Xs = X[:, cidx].toarray()
gi = {x: i for i, x in enumerate(setg)}

# cell-type donor means and count shares
rows, share = [], []
for d in dq.gsm:
    m = gsm == d
    tot = {x: Xs[m, gi[x]].sum() for x in ["HSPB1", "CRYAB"] if x in gi}
    for t in np.unique(lab[m]):
        mt = m & (lab == t)
        rows.append(dict(gsm=d, label=t, n_cells=int(mt.sum()), **{x: float(Ls[mt, gi[x]].mean()) for x in HSP}))
        share.append(dict(gsm=d, label=t, frac_cells=mt.sum() / m.sum(), **{f"frac_{x}_counts": Xs[mt, gi[x]].sum() / v for x, v in tot.items()}))
ct = pd.DataFrame(rows).merge(dq[["gsm", "id", "disease", "sex_f"]], on="gsm")
ct.to_csv(OUT / "celltype_donor_means.tsv", sep="\t", index=False)
sh = pd.DataFrame(share).merge(dq[["gsm", "disease"]], on="gsm")
sh.to_csv(OUT / "count_share.tsv", sep="\t", index=False)

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


# disease effect per gene, hepatocytes and other cell types with >= 20 cells in every donor
de = []
for t, o in ct.groupby("label"):
    o = o.set_index("gsm").reindex(dq.gsm)
    if o.n_cells.isna().any() or (o.n_cells < 20).any():
        continue
    for x in HSP:
        b, p = dis(o[x], dq.disease.to_numpy(), dq.sex_f.to_numpy(), RELAB)
        de.append(dict(label=t, gene=x, beta=b, p_exact=p, mean_healthy=o[x][dq.disease.to_numpy() == 0].mean(),
                       mean_masld=o[x][dq.disease.to_numpy() == 1].mean()))
de = pd.DataFrame(de)
de["bh"] = de.groupby("label").p_exact.transform(stats.false_discovery_control)
de.to_csv(OUT / "disease_effect_by_celltype.tsv", sep="\t", index=False)

# hepatocytes: within-donor rho with MDB, against matched single genes
H = lab == "Hepatocytes"
XH = X[H]; LH = L[H]; gH = gsm[H]; umiH = umi[H]; ngH = np.asarray((XH > 0).sum(axis=1)).ravel()
mean_all = np.asarray(LH.mean(axis=0)).ravel(); det_all = np.asarray((XH > 0).mean(axis=0)).ravel()
keep = det_all >= 0.01
binm = pd.qcut(pd.Series(mean_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
bind = pd.qcut(pd.Series(det_all[keep]).rank(method="first"), 10, labels=False).to_numpy()
kg = np.array(genes)[keep]
bins = pd.Series(list(zip(binm, bind)), index=kg)
pool = {b: [x for x in idx if x not in set(setg)] for b, idx in bins.groupby(bins).groups.items()}
test = [x for x in HSP if x in bins.index]
pd.DataFrame(dict(gene=HSP, hep_mean=[mean_all[genes.get_loc(x)] for x in HSP], hep_detection=[det_all[genes.get_loc(x)] for x in HSP],
                  tested=[x in bins.index for x in HSP])).to_csv(OUT / "hep_gene_detection.tsv", sep="\t", index=False)
NB = {b: list(rng.choice(pool[b], size=min(N_NULL, len(pool[b])), replace=False)) for b in sorted(set(bins[test]))}
need = sorted(set(setg) | set(x for v in NB.values() for x in v))
LdH = LH[:, [genes.get_loc(x) for x in need]].toarray()
gj = {x: i for i, x in enumerate(need)}


def resid_rank(M, Z):
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]


def cor_cols(a, M):
    a = (a - a.mean()) / a.std(); M = (M - M.mean(0)) / M.std(0)
    return (M * a[:, None]).mean(0)


wd = []
for d in dq.gsm:
    m = gH == d
    D = LdH[m]; sd = D.std(0); sd[sd == 0] = np.nan
    Zc = np.column_stack([np.ones(m.sum()), stats.rankdata(np.log(umiH[m])), stats.rankdata(np.log(ngH[m]))])
    zi = np.nanmean(np.column_stack([(D[:, gj[x]] - D[:, gj[x]].mean()) / sd[gj[x]] for x in INJ]), 1)
    ri = resid_rank(zi[:, None], Zc)[:, 0]
    nullrho = {}
    for b, v in NB.items():
        cols = [gj[x] for x in v if np.isfinite(sd[gj[x]])]
        nullrho[b] = np.nanmean(cor_cols(ri, resid_rank(D[:, cols], Zc)))
    for x in test:
        if not np.isfinite(sd[gj[x]]):
            continue
        ro = resid_rank(D[:, [gj[x]]], Zc)[:, 0]
        wd.append(dict(gsm=d, gene=x, n_cells=int(m.sum()), rho=np.corrcoef(ri, ro)[0, 1], null_mean=nullrho[bins[x]]))
    print(d, int(m.sum()), flush=True)
wd = pd.DataFrame(wd).merge(dq[["gsm", "id", "disease"]], on="gsm")
wd["excess"] = wd.rho - wd.null_mean
wd.to_csv(OUT / "hep_within_donor_rho.tsv", sep="\t", index=False)
summ = []
for x, o in wd.groupby("gene"):
    o = o.set_index("gsm").reindex(dq.gsm)
    if o.excess.isna().any():
        continue
    summ.append(dict(gene=x, mean_rho=o.rho.mean(), mean_null=o.null_mean.mean(), mean_excess=o.excess.mean(),
                     n_excess_positive=int((o.excess > 0).sum()), p_excess=signflip(o.excess)))
summ = pd.DataFrame(summ); summ["bh"] = stats.false_discovery_control(summ.p_excess)
summ.to_csv(OUT / "hep_within_donor_summary.tsv", sep="\t", index=False)

s1 = summ.set_index("gene").loc["HSPB1"]
h2 = de[(de.label == "Hepatocytes") & (de.gene == "HSPB1")].iloc[0]
prim = pd.DataFrame([dict(id="H1", test="hepatocytes: within-donor rho(MDB, HSPB1) minus matched genes", est=s1.mean_excess,
                          n_positive=int(s1.n_excess_positive), p=s1.p_excess),
                     dict(id="H2", test="hepatocytes: disease effect on donor-mean HSPB1 (log1p CP10k)", est=h2.beta, n_positive=np.nan, p=h2.p_exact)])
o = np.argsort(prim.p.to_numpy()); ph = np.maximum.accumulate(prim.p.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 220); pd.set_option("display.max_rows", 500)
print(prim.round(4).to_string())
print(summ.sort_values("mean_excess").round(3).to_string())
print(de[de.label == "Hepatocytes"].sort_values("beta").round(3).to_string())
print(ct.groupby("label")[["HSPB1", "CRYAB", "HSPA1A", "BAG3", "HSP90AA1"]].mean().round(3).to_string())
print(sh.groupby(["label", "disease"])[[c for c in sh.columns if c.startswith("frac_") and c.endswith("_counts")]].mean().round(3).to_string())
