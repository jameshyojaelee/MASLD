"""DIG1 job 18b (exploratory). Is the GNMT within-cell axis special, after removing generic
co-detection structure? GSE298719 Flex hepatocytes; anchors from job 17 (GNMT, 40 null, 40
changed). Within-cell partner z per anchor, two covariate sets for the per-donor rank residuals:
  base   : log UMI (job 9)
  quality: log UMI, log genes detected, fraction of UMI in the 50 most abundant hepatocyte genes,
           mitochondrial fraction
from all 13 donors and from the 6 healthy donors. Agreement of each anchor's partner z (anchor
excluded) with: GSE298719 logFC, GSE256398 logFC, bulk between-unit partner rho and protein
partner rho (job 18a). Spearman and partial Spearman, the latter with gene ranks residualised on
Flex hepatocyte detection and log mean CP10k (plus bulk mean logCPM for bulk). GNMT is reported as
a percentile among null and changed anchors. Descriptive. Usage: python -I 18b_anchor_agreement.py
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
BETW = RUN / "out_anchor_between"
OUT = RUN / "out_anchor_agreement"
OUT.mkdir(exist_ok=False)
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str})
dq = dq[dq.included].sort_values("gsm").reset_index(drop=True)
anc = pd.read_csv(BETW / "anchors.tsv", sep="\t")
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad", backed="r")
genes = pd.Index(A.var_names)
m = (A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()
S = A[m].to_memory()
X = S.X.tocsr().astype(np.float64)
gsm = S.obs.gsm.astype(str).to_numpy()
umi = np.asarray(X.sum(axis=1)).ravel()
ngen = np.asarray((X > 0).sum(axis=1)).ravel()
gsum = np.asarray(X.sum(axis=0)).ravel()
top50 = np.argsort(gsum)[::-1][:50]
dom = np.asarray(X[:, top50].sum(axis=1)).ravel() / umi
mt = np.where(genes.str.startswith("MT-"))[0]
mtf = np.asarray(X[:, mt].sum(axis=1)).ravel() / umi if len(mt) else np.zeros(len(umi))
det = np.asarray((X > 0).mean(axis=0)).ravel()
cand = np.where(det >= 0.10)[0]
a_idx = np.array([genes.get_loc(g) for g in anc.anchor])
cols = np.union1d(cand, a_idx)
pos = {c: i for i, c in enumerate(cols)}
Z = {"base": [], "quality": []}
dis = []
for g in dq.gsm:
    ix = np.where(gsm == g)[0]
    M = X[ix][:, cols].toarray() / umi[ix][:, None]
    R0 = stats.rankdata(M, axis=0)
    del M
    for cov, C in [("base", np.column_stack([np.ones(len(ix)), np.log(umi[ix])])),
                   ("quality", np.column_stack([np.ones(len(ix)), np.log(umi[ix]), np.log(ngen[ix]), dom[ix], mtf[ix]]))]:
        Rk = R0 - C @ np.linalg.lstsq(C, R0, rcond=None)[0]
        Rk /= np.sqrt((Rk ** 2).sum(0)) + 1e-12
        Z[cov].append(np.arctanh(np.clip(Rk.T @ Rk[:, [pos[c] for c in a_idx]], -0.999, 0.999)))
        del Rk
    dis.append(int(dq.set_index("gsm").disease[g]))
    del R0
dis = np.array(dis)
zg = pd.Index(genes[cols])
de = pd.read_csv(RUN / "out_pseudobulk_v2/gene_DE_edgeR_all_lineages.tsv.gz", sep="\t").drop_duplicates(["lkey", "symbol"])
Fl = de[de.lkey == "GSE298719:hepatocytes"].set_index("symbol").logFC
Vl = de[de.lkey == "GSE256398:hepatocytes"].set_index("symbol").logFC
Bp = pd.read_csv(BETW / "bulk_partner_rho.tsv.gz", sep="\t").set_index("gene")
Pp = pd.read_csv(BETW / "protein_partner_rho.tsv.gz", sep="\t").set_index("gene")
bmean = pd.read_csv(BETW / "bulk_gene_means.tsv.gz", sep="\t").set_index("gene").bulk_mean_logcpm
gcov = pd.DataFrame({"det": det, "lmean": np.log(gsum / umi.sum() * 1e4 + 1e-6)}, index=genes)


def pspear(a, b, covs):
    ok = a.notna() & b.notna() & covs.notna().all(axis=1)
    a, b, c = a[ok], b[ok], covs[ok]
    if len(a) < 100:
        return np.nan, np.nan, int(len(a))
    s = stats.spearmanr(a, b)[0]
    C = np.column_stack([np.ones(len(a))] + [stats.rankdata(c[k]) for k in c.columns])
    ra, rb = stats.rankdata(a), stats.rankdata(b)
    ra = ra - C @ np.linalg.lstsq(C, ra, rcond=None)[0]
    rb = rb - C @ np.linalg.lstsq(C, rb, rcond=None)[0]
    return s, float(np.corrcoef(ra, rb)[0, 1]), int(len(a))


rows = []
cand_names = genes[cand]
for k, (ga, kind) in enumerate(zip(anc.anchor, anc.kind)):
    for cov in Z:
        Zc = np.array(Z[cov])
        for sub, msk in [("all", np.ones(len(dis), bool)), ("healthy", dis == 0)]:
            z = pd.Series(Zc[msk][:, :, k].mean(0), index=zg).drop(ga)
            z = z[z.index.isin(cand_names)]
            gc = gcov.reindex(z.index)
            for tgt, vec, extra in [("flex_logFC", Fl, None), ("v3_logFC", Vl, None),
                                    ("bulk_partner", Bp[ga] if ga in Bp.columns else None, bmean),
                                    ("protein_partner", Pp[ga] if ga in Pp.columns else None, None)]:
                if vec is None:
                    continue
                v = vec.reindex(z.index)
                cv = gc if extra is None else gc.assign(bm=extra.reindex(z.index))
                s, ps, n = pspear(z, v, cv)
                rows.append(dict(anchor=ga, kind=kind if isinstance(kind, str) else "null", covariates=cov, donors=sub,
                                 target=tgt, n_genes=n, spearman=s, partial_spearman=ps))
res = pd.DataFrame(rows)
res.to_csv(OUT / "anchor_agreement.tsv", sep="\t", index=False)
summ = []
for (cov, sub, tgt), g in res.groupby(["covariates", "donors", "target"]):
    real = g[g.anchor == "GNMT"]
    if real.empty:
        continue
    for stat in ["spearman", "partial_spearman"]:
        x = real[stat].iloc[0]
        for kd in ["null", "changed"]:
            ref = g[g.kind == kd][stat].dropna()
            summ.append(dict(covariates=cov, donors=sub, target=tgt, stat=stat, ref=kd, n_ref=len(ref), gnmt=x,
                             ref_median=ref.median(), ref_p05=ref.quantile(0.05), ref_p95=ref.quantile(0.95),
                             frac_ref_more_extreme=float((ref.abs() >= abs(x)).mean())))
summ = pd.DataFrame(summ)
summ.to_csv(OUT / "gnmt_vs_control_anchors.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(summ.round(3).to_string(index=False))
