"""DIG1 job 2 (exploratory; EXPLORATORY_PLAN.md question A4). GSE298719 hepatocytes, cell level.

Is the Ductular injury rise a hepatocyte subpopulation or a shift of all hepatocytes, and what
else marks Ductular-high hepatocytes (zonation, cholangiocyte identity, doublet signs)?
Cells are summarised per donor; donors are the unit for every disease comparison.
Usage: python -I 02_hepatocyte_cells.py
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
HOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
OUT = RUN / "out_cells"
OUT.mkdir(exist_ok=False)
DUCT = "hotspot_hepatocytes_48f39dd4d817a10e"
TOP_FRAC = 0.05

PANELS = {
    "cholangiocyte_identity": ["KRT19", "KRT7", "EPCAM", "CFTR", "SLC4A2", "ANXA4", "SPP1"],
    "biliary_hybrid": ["SOX9", "CD24", "KRT23", "PROM1", "ONECUT1", "HNF1B"],
    "periportal": ["CPS1", "ASS1", "HAL", "SDS", "PCK1", "HAMP"],
    "pericentral": ["GLUL", "CYP2E1", "CYP1A2", "OAT", "CYP3A4"],
    "hepatocyte_identity": ["ALB", "APOA1", "TF", "HNF4A", "SERPINA1"],
    "proliferation": ["MKI67", "TOP2A"],
    "stromal_identity": ["COL1A1", "COL1A2", "COL3A1", "DCN"],
    "immune_endothelial": ["PTPRC", "CD74", "PECAM1"],
}


def log_normalize(counts, library):
    scale = np.divide(1e4, library, out=np.zeros_like(library, dtype=float), where=library > 0)
    norm = (sparse.diags(scale.astype(np.float32)) @ counts).tocsr()
    norm.data = np.log1p(norm.data)
    return norm


A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str})
dq = dq[dq.included].sort_values("gsm").reset_index(drop=True)
hep = A[(A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
del A
X = sparse.csr_matrix(hep.X)
genes = pd.Index(hep.var_names)
umi = np.asarray(X.sum(axis=1)).ravel().astype(float)
norm = log_normalize(X, umi)
obs = hep.obs[["gsm", "umi", "n_genes"]].reset_index(drop=True)

realised = pd.read_csv(V5 / "outputs/job3/realised_gene_set_check.tsv", sep="\t").realised_genes[0].split(";")
mem = pd.read_csv(HOT / "program_membership_v2.tsv", sep="\t")
mem = mem[(mem.program_uid == DUCT) & mem.mapped_symbol.isin(realised)].groupby("mapped_symbol").original_l1_weight.sum()
assert len(mem) == 18 and all(g in genes for g in mem.index)
members = set(mem.index)
w = mem / mem.sum()

# Per-cell scores. z over all included hepatocytes, so donors share one scale.
M = norm[:, [genes.get_loc(g) for g in w.index]].toarray()
Z = (M - M.mean(0)) / M.std(0)
obs["duct_score"] = Z @ w.to_numpy()
raw_member = np.asarray(X[:, [genes.get_loc(g) for g in w.index]].sum(axis=1)).ravel()
obs["member_umi_frac"] = raw_member / umi


def panel_score(name):
    gs = [g for g in PANELS[name] if g in genes and g not in members]
    m = norm[:, [genes.get_loc(g) for g in gs]].toarray()
    return ((m - m.mean(0)) / np.where(m.std(0) > 0, m.std(0), 1)).mean(1), gs


panel_genes = {}
for p in PANELS:
    obs[p], panel_genes[p] = panel_score(p)
obs["zonation_pp_minus_pc"] = obs.periportal - obs.pericentral
pd.Series({p: ";".join(g) for p, g in panel_genes.items()}).to_csv(OUT / "panel_genes_used.tsv", sep="\t", header=["genes"])

design = dq[["gsm", "disease", "sex_f"]].reset_index(drop=True)
n = len(design)
RELAB = []
for ix in combinations(range(n), int(design.disease.sum())):
    v = np.zeros(n, int); v[list(ix)] = 1; RELAB.append(v)
RELAB = np.array(RELAB)
OBS_COL = int(np.where((RELAB == design.disease.to_numpy()).all(1))[0][0])


def disease_effect(y):
    """OLS beta for disease (~ disease + sex), HC3 CI, exact relabeling p (two-sided)."""
    y = np.asarray(y, float)
    def beta(v):
        Xd = np.column_stack([np.ones(n), v, design.sex_f]).astype(float)
        return np.linalg.lstsq(Xd, y, rcond=None)[0][1]
    Xd = np.column_stack([np.ones(n), design.disease, design.sex_f]).astype(float)
    XtXi = np.linalg.inv(Xd.T @ Xd); b = XtXi @ Xd.T @ y; e = y - Xd @ b
    h = np.einsum("ij,jk,ik->i", Xd, XtXi, Xd); C = XtXi @ Xd.T
    se = np.sqrt(np.sum(C[1] ** 2 * e ** 2 / (1 - h) ** 2)); df = n - 3
    nb = np.array([beta(v) for v in RELAB])
    return dict(beta=b[1], ci_low=b[1] - stats.t.ppf(0.975, df) * se, ci_high=b[1] + stats.t.ppf(0.975, df) * se,
                p_exact=np.mean(np.abs(nb) >= abs(b[1]) - 1e-12), n_relabelings=len(RELAB))


# ---- 1. distribution shape per donor ------------------------------------------------------
healthy = design.gsm[design.disease == 0]
thr = np.median([np.quantile(obs.duct_score[obs.gsm == g], 0.95) for g in healthy])
qs = np.arange(0.05, 1.0, 0.05)
rows = []
for g in design.gsm:
    s = obs.duct_score[obs.gsm == g].to_numpy()
    r = dict(gsm=g, n_hep=len(s), mean=s.mean(), frac_above_healthy_p95=np.mean(s > thr),
             member_umi_frac_mean=obs.member_umi_frac[obs.gsm == g].mean())
    r.update({f"q{int(round(q * 100)):02d}": np.quantile(s, q) for q in qs})
    rows.append(r)
donor = design.merge(pd.DataFrame(rows), on="gsm")
donor.to_csv(OUT / "donor_score_distribution.tsv", sep="\t", index=False)
eff = []
for col in ["mean", "frac_above_healthy_p95", "member_umi_frac_mean"] + [f"q{int(round(q * 100)):02d}" for q in qs]:
    y = donor[col].to_numpy(float)
    if col == "frac_above_healthy_p95":
        y = np.log(np.clip(y, 1e-3, None) / np.clip(1 - y, 1e-3, None))
    eff.append(dict(statistic=col, scale="logit" if col == "frac_above_healthy_p95" else "score",
                    sd_healthy=donor.loc[donor.disease == 0, col].std(), **disease_effect(y)))
pd.DataFrame(eff).to_csv(OUT / "disease_effect_by_quantile.tsv", sep="\t", index=False)
pd.Series(dict(healthy_p95_threshold=thr)).to_csv(OUT / "threshold.tsv", sep="\t", header=["value"])

# ---- 2. what marks Ductular-high hepatocytes (top 5% within donor vs rest) -----------------
panel_cols = list(PANELS) + ["zonation_pp_minus_pc", "umi", "n_genes"]
prow, corr = [], []
for g in design.gsm:
    o = obs[obs.gsm == g]
    top = o.duct_score >= o.duct_score.quantile(1 - TOP_FRAC)
    r = dict(gsm=g, n_top=int(top.sum()))
    for c in panel_cols:
        v = np.log10(o[c]) if c in ("umi", "n_genes") else o[c]
        r[c] = v[top].mean() - v[~top].mean()
    prow.append(r)
    corr.append(dict(gsm=g, **{f"rho_{c}": stats.spearmanr(o.duct_score, o[c])[0] for c in list(PANELS) + ["zonation_pp_minus_pc"]},
                     rho_pp_vs_pc=stats.spearmanr(o.periportal, o.pericentral)[0]))
pt = design.merge(pd.DataFrame(prow), on="gsm")
pt.to_csv(OUT / "top5_vs_rest_by_donor.tsv", sep="\t", index=False)
cr = design.merge(pd.DataFrame(corr), on="gsm")
cr.to_csv(OUT / "cell_correlations_by_donor.tsv", sep="\t", index=False)
summ = []
for c in panel_cols:
    v = pt[c].to_numpy()
    summ.append(dict(feature=c, unit="log10 ratio" if c in ("umi", "n_genes") else "z difference",
                     median_all=np.median(v), n_positive=int((v > 0).sum()), n=len(v),
                     sign_test_p=stats.binomtest(int((v > 0).sum()), len(v)).pvalue,
                     median_healthy=np.median(v[pt.disease == 0]), median_disease=np.median(v[pt.disease == 1])))
for c in [k for k in cr.columns if k.startswith("rho_")]:
    v = cr[c].to_numpy()
    summ.append(dict(feature=c, unit="within-donor Spearman", median_all=np.median(v), n_positive=int((v > 0).sum()),
                     n=len(v), sign_test_p=stats.binomtest(int((v > 0).sum()), len(v)).pvalue,
                     median_healthy=np.median(v[cr.disease == 0]), median_disease=np.median(v[cr.disease == 1])))
pd.DataFrame(summ).to_csv(OUT / "top5_and_correlation_summary.tsv", sep="\t", index=False)

# ---- 3. unbiased: genes co-enriched in Ductular-high hepatocytes (members excluded) ---------
det = np.asarray((X > 0).mean(axis=0)).ravel()
keep = np.where((det >= 0.05) & ~genes.isin(list(members)))[0]
cp = (sparse.diags((1e4 / umi).astype(np.float32)) @ X[:, keep]).tocsr()
lr = []
for g in design.gsm:
    idx = np.where(obs.gsm == g)[0]
    s = obs.duct_score.to_numpy()[idx]
    top = s >= np.quantile(s, 1 - TOP_FRAC)
    a = np.asarray(cp[idx[top]].mean(axis=0)).ravel(); b = np.asarray(cp[idx[~top]].mean(axis=0)).ravel()
    lr.append(np.log2((a + 0.1) / (b + 0.1)))
L = np.array(lr)                                   # donors x genes
t_all = stats.ttest_1samp(L, 0).statistic
dis = design.disease.to_numpy() == 1
res = pd.DataFrame(dict(gene=genes[keep], detection=det[keep], mean_log2ratio=L.mean(0), t_one_sample=t_all,
                        n_donors_up=(L > 0).sum(0), mean_log2ratio_healthy=L[~dis].mean(0),
                        mean_log2ratio_disease=L[dis].mean(0)))
res["p"] = 2 * stats.t.sf(np.abs(res.t_one_sample), n - 1)
res["bh_q"] = stats.false_discovery_control(res.p.fillna(1))
res.sort_values("t_one_sample", ascending=False).to_csv(OUT / "coenriched_genes_top5_vs_rest.tsv.gz", sep="\t", index=False)
obs.to_csv(OUT / "hepatocyte_cell_scores.tsv.gz", sep="\t", index=False)
print("done", flush=True)
