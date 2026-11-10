"""DIG1 job 9 (exploratory). GSE298719 Flex hepatocytes. Two questions.

Q1. Does MASLD move hepatocytes along an axis that already exists among healthy hepatocytes, or
    along a new one? PCA on healthy-donor hepatocytes only (analytic Pearson residuals, theta 100,
    2,000 most variable genes among healthy cells, 3,000 cells sampled per donor, seed 20261024).
    The disease vector is the edgeR hepatocyte logFC (out_pseudobulk_v2). Report the Spearman rho
    between the disease vector and each healthy PC loading (PC1-10), and the fraction of the
    disease vector's variance (over the 2,000 genes, centered) captured by the healthy PC subspace.
    Null: the same quantity for 1,000 random unit vectors (sanity) and for 200 donor relabelings
    of the pseudobulk-free cell mean difference (the observed donor split is one of 1,716).
Q2. Within donor, which genes co-vary with GNMT across hepatocytes? Partial Spearman of GNMT with
    every gene detected in >= 10% of hepatocytes, depth (log UMI) removed from ranks, per donor;
    Fisher z mean over 13 donors; one-sample t; BH over genes. Then: do the top within-donor GNMT
    partners overlap the genes that fall in MASLD (edgeR FDR < 0.10, down)?
Donor is the unit for every inference. Usage: python -I 09_hep_axes.py
"""
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
OUT = RUN / "out_hep_axes"
OUT.mkdir(exist_ok=False)
rng = np.random.default_rng(20261024)

A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str})
dq = dq[dq.included].sort_values("gsm").reset_index(drop=True)
hep = A[(A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
del A
genes = pd.Index(hep.var_names)
X = sparse.csr_matrix(hep.X).astype(np.float64)
gsm = hep.obs.gsm.astype(str).to_numpy()
umi = np.asarray(X.sum(axis=1)).ravel()
det = np.asarray((X > 0).mean(axis=0)).ravel()

# sample cells per donor
take = np.concatenate([rng.choice(np.where(gsm == g)[0], min(3000, (gsm == g).sum()), replace=False) for g in dq.gsm])
Xs, gs, us = X[take], gsm[take], umi[take]
healthy = set(dq.gsm[dq.disease == 0])
hmask = np.isin(gs, list(healthy))


def pearson_resid(Xc, u, cols, gene_frac):
    mu = np.outer(u, gene_frac[cols])
    Y = Xc[:, cols].toarray()
    r = (Y - mu) / np.sqrt(mu + mu ** 2 / 100.0)
    lim = np.sqrt(len(u))
    return np.clip(r, -lim, lim)


keep = np.where(det >= 0.01)[0]
gene_frac_h = np.asarray(Xs[hmask].sum(axis=0)).ravel() / us[hmask].sum()
Rh = pearson_resid(Xs[hmask], us[hmask], keep, gene_frac_h)
var = Rh.var(0)
top = keep[np.argsort(var)[::-1][:2000]]
Rh = pearson_resid(Xs[hmask], us[hmask], top, gene_frac_h)
Rh -= Rh.mean(0)
U, S, Vt = np.linalg.svd(Rh, full_matrices=False)
L = Vt[:10]                                   # loadings, PCs x genes
evr = (S[:10] ** 2) / (S ** 2).sum()
sym = genes[top]

de = pd.read_csv(RUN / "out_pseudobulk_v2/gene_DE_edgeR_all_lineages.tsv.gz", sep="\t")
de = de[de.lkey == "GSE298719:hepatocytes"].drop_duplicates("symbol").set_index("symbol")
ok = np.array([s in de.index for s in sym])
dv = de.logFC.reindex(sym[ok]).to_numpy()
dvc = dv - dv.mean()
Lk = L[:, ok]
rows = [dict(pc=i + 1, var_explained_healthy=evr[i], spearman_with_disease_logFC=stats.spearmanr(Lk[i], dv)[0]) for i in range(10)]
# fraction of disease-vector variance in healthy PC1-10 subspace (Lk rows ~orthonormal on subset)
Q, _ = np.linalg.qr(Lk.T)
frac_obs = np.sum((Q.T @ dvc) ** 2) / np.sum(dvc ** 2)
rand = []
for _ in range(1000):
    v = rng.standard_normal(len(dvc)); v -= v.mean()
    rand.append(np.sum((Q.T @ v) ** 2) / np.sum(v ** 2))
pd.DataFrame(rows).to_csv(OUT / "healthy_pc_vs_disease_vector.tsv", sep="\t", index=False)
pd.DataFrame([dict(n_genes=int(ok.sum()), frac_disease_var_in_healthy_pc1_10=frac_obs, random_vector_median=np.median(rand),
                   random_vector_p95=np.quantile(rand, 0.95))]).to_csv(OUT / "subspace_fraction.tsv", sep="\t", index=False)
pd.DataFrame(L.T, index=sym, columns=[f"PC{i+1}" for i in range(10)]).to_csv(OUT / "healthy_pc_loadings.tsv.gz", sep="\t")

# project all sampled cells (healthy-gene-frac scaling) onto healthy PCs; donor means and spread
Ra = pearson_resid(Xs, us, top, gene_frac_h) - pearson_resid(Xs[hmask], us[hmask], top, gene_frac_h).mean(0)
P = Ra @ L.T
proj = pd.DataFrame(P, columns=[f"PC{i+1}" for i in range(10)]); proj["gsm"] = gs
dm = proj.groupby("gsm").agg(["mean", "std"])
dm.columns = [f"{a}_{b}" for a, b in dm.columns]
dm = dq[["gsm", "id", "disease", "fibrosis"]].set_index("gsm").join(dm)
# residual: disease-axis score of each cell = projection on the unit disease vector, restricted to the
# part orthogonal to healthy PCs
dfull = np.zeros(len(top)); dfull[ok] = dvc
d_orth = dfull - Q @ (Q.T @ dvc) if False else None
dunit = np.zeros(len(top)); dunit[ok] = dvc / np.linalg.norm(dvc)
Qf = np.zeros((len(top), Q.shape[1])); Qf[ok] = Q
dperp = dunit - Qf @ (Qf.T @ dunit)
proj_perp = Ra @ (dperp / np.linalg.norm(dperp))
dm["disease_axis_orthogonal_mean"] = pd.Series(proj_perp).groupby(gs).mean()
dm["disease_axis_orthogonal_sd"] = pd.Series(proj_perp).groupby(gs).std()
dm.to_csv(OUT / "donor_projection_summary.tsv", sep="\t")

n = len(dq)
RELAB = []
for ix in combinations(range(n), int(dq.disease.sum())):
    v = np.zeros(n, int); v[list(ix)] = 1; RELAB.append(v)
RELAB = np.array(RELAB)


def eff(y):
    def beta(v):
        Xd = np.column_stack([np.ones(n), v, dq.sex_f]).astype(float)
        return np.nan if np.linalg.matrix_rank(Xd) < 3 else np.linalg.lstsq(Xd, y, rcond=None)[0][1]
    b = beta(dq.disease.to_numpy()); nb = np.array([beta(v) for v in RELAB])
    return dict(beta=b, p_exact=np.nanmean(np.abs(nb) >= abs(b) - 1e-12))


er = []
for c in [f"PC{i+1}_mean" for i in range(10)] + [f"PC{i+1}_std" for i in range(10)] + ["disease_axis_orthogonal_mean", "disease_axis_orthogonal_sd"]:
    er.append(dict(statistic=c, **eff(dm.loc[dq.gsm, c].to_numpy())))
pd.DataFrame(er).to_csv(OUT / "disease_effect_on_projections.tsv", sep="\t", index=False)

# Q2: within-donor GNMT partners
g_idx = genes.get_loc("GNMT")
cand = np.where(det >= 0.10)[0]
cand = cand[cand != g_idx]
Z = []
for g in dq.gsm:
    ix = np.where(gsm == g)[0]
    lu = np.log(umi[ix])
    C = np.column_stack([np.ones(len(ix)), lu])
    def rres(M):
        Rk = stats.rankdata(M, axis=0)
        return Rk - C @ np.linalg.lstsq(C, Rk, rcond=None)[0]
    ry = rres(np.asarray(X[ix, g_idx].todense()).ravel() / umi[ix])
    rho = np.empty(len(cand))
    for s in range(0, len(cand), 500):
        cols = cand[s:s + 500]
        M = X[ix][:, cols].toarray() / umi[ix][:, None]
        rx = rres(M)
        rho[s:s + 500] = (rx * ry[:, None]).sum(0) / np.sqrt((rx ** 2).sum(0) * (ry ** 2).sum())
    Z.append(np.arctanh(np.clip(rho, -0.999, 0.999)))
Z = np.array(Z)
t = stats.ttest_1samp(Z, 0)
res = pd.DataFrame(dict(gene=genes[cand], detection=det[cand], mean_z=Z.mean(0), n_pos=(Z > 0).sum(0),
                        mean_z_healthy=Z[dq.disease.to_numpy() == 0].mean(0), mean_z_masld=Z[dq.disease.to_numpy() == 1].mean(0),
                        t=t.statistic, p=t.pvalue))
res["bh_q"] = stats.false_discovery_control(res.p.fillna(1))
res = res.merge(de[["logFC", "FDR"]].rename(columns={"logFC": "disease_logFC", "FDR": "disease_FDR"}), left_on="gene", right_index=True, how="left")
res.sort_values("mean_z", ascending=False).to_csv(OUT / "gnmt_within_donor_partners.tsv.gz", sep="\t", index=False)
# overlap: are positive GNMT partners enriched among MASLD-down genes, negative among MASLD-up?
r2 = res.dropna(subset=["disease_logFC"])
r2["partner"] = np.where((r2.bh_q < 0.05) & (r2.mean_z > 0), "pos", np.where((r2.bh_q < 0.05) & (r2.mean_z < 0), "neg", "none"))
r2["dz"] = np.where((r2.disease_FDR < 0.10) & (r2.disease_logFC < 0), "down", np.where((r2.disease_FDR < 0.10) & (r2.disease_logFC > 0), "up", "ns"))
pd.crosstab(r2.partner, r2.dz).to_csv(OUT / "partner_vs_disease_crosstab.tsv", sep="\t")
pd.DataFrame([dict(spearman_meanz_vs_disease_logFC=stats.spearmanr(r2.mean_z, r2.disease_logFC)[0], n=len(r2))]).to_csv(
    OUT / "partner_vs_disease_rho.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).round(3).to_string(), frac_obs, np.quantile(rand, 0.95), flush=True)
