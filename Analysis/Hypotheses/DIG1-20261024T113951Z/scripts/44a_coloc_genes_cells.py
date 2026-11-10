"""DIG1 job 44a (exploratory; follows jobs 41c and 43). GSE298719 Flex, 13 donors (6 healthy, 7 MASLD).
For the colocalized genes with a genetic direction (41c table: SuSiE top pair, tier 1-2, one gene per
1 Mb locus), two hepatocyte statistics per gene:
  coupling  within-donor rho(gene, injury score) with ranks residualised on log UMI and log genes,
            minus the mean rho of 200 random genes from the same hepatocyte mean x detection decile
            bin (one draw per bin, seed 20261024). Injury score = mean within-donor z of KRT8, KRT18,
            SQSTM1 (as job 37b code). Exact sign-flip over 13 donors; BH over genes.
  disease   OLS donor mean ~ disease + sex, exact relabeling (1,716 splits); BH over genes.
Input for 44c. Usage: python -I 44a_coloc_genes_cells.py <new_out_dir>
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
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
GT = pd.read_csv(RUN / "out_direction_posthoc_controls/gene_directions_susie_top_dedup.tsv", sep="\t")
HSP = [x for x in GT.gene.astype(str).unique() if x not in ("KRT8", "KRT18", "SQSTM1")]
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

# hepatocyte donor means
rows = []
for d in dq.gsm:
    mt = (gsm == d) & (lab == "Hepatocytes")
    rows.append(dict(gsm=d, label="Hepatocytes", n_cells=int(mt.sum()), **dict(zip(HSP, Ls[mt][:, [gi[x] for x in HSP]].mean(0).astype(float)))))
ct = pd.DataFrame(rows).merge(dq[["gsm", "id", "disease", "sex_f"]], on="gsm")
ct.to_csv(OUT / "hep_donor_means.tsv", sep="\t", index=False)

n = len(dq)
SIGNS = np.array(list(product([1, -1], repeat=n)))
RELAB = []
for ix in combinations(range(n), int(dq.disease.sum())):
    v = np.zeros(n, int); v[list(ix)] = 1; RELAB.append(v)
RELAB = np.array(RELAB)


def signflip(v):
    v = np.asarray(v, float)
    return np.mean(np.abs((SIGNS * v).mean(1)) >= abs(v.mean()) - 1e-12)


# disease effect per gene in hepatocytes, exact relabeling (vectorised)
o = ct.set_index("gsm").reindex(dq.gsm)
Y = o[HSP].to_numpy(float); sexv = dq.sex_f.to_numpy(float)
row1 = lambda v: np.linalg.pinv(np.column_stack([np.ones(n), v, sexv]))[1]
b_obs = row1(dq.disease.to_numpy(float)) @ Y
B_rel = np.array([row1(v.astype(float)) for v in RELAB]) @ Y
p_ex = np.mean(np.abs(B_rel) >= np.abs(b_obs)[None, :] - 1e-12, axis=0)
dmask = dq.disease.to_numpy() == 1
de = pd.DataFrame(dict(gene=HSP, beta=b_obs, p_exact=p_ex, mean_healthy=Y[~dmask].mean(0), mean_masld=Y[dmask].mean(0)))
de["bh"] = stats.false_discovery_control(de.p_exact)
de.to_csv(OUT / "hep_disease_effect.tsv", sep="\t", index=False)

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

pd.set_option("display.width", 220); pd.set_option("display.max_rows", 500)
print(summ.sort_values("p_excess").head(40).round(4).to_string()); print(de.sort_values("p_exact").head(30).round(4).to_string())
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
