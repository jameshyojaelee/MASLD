"""DIG1 job 7 (exploratory; follows C and job 6). GSE298719 Flex hepatocytes, cell level.

Q1. Is lower GNMT/MAT1A in MASLD hepatocytes a shift of all hepatocytes (every quantile lower) or a
    subset? Q2. Do SOX9+/EPCAM+ ('plastic') hepatocytes differ from other hepatocytes in the same
    donor and UMI decile in methionine-cycle, zonation, identity and sterol-synthesis genes?
Q3. Does the GNMT/MAT1A disease effect remain when plastic cells are removed? Q4. Is it present in
    both periportal-high and pericentral-high hepatocytes?
Donor is the unit throughout. Disease effects: OLS ~ disease + sex, exact relabeling over 1,716
splits. Per-cell expression: counts per 10k, log1p. Usage: python -I 07_gnmt_plastic_cells.py
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
OUT = RUN / "out_gnmt_cells"
OUT.mkdir(exist_ok=False)
GENES = {
    "methionine": ["GNMT", "MAT1A", "BHMT", "AHCY", "CBS"],
    "periportal": ["CPS1", "HAL", "SDS", "ASS1", "PCK1"],
    "pericentral": ["GLUL", "CYP2E1", "OAT", "CYP1A2"],
    "identity": ["ALB", "APOA1", "SERPINA1", "TF", "HP"],
    "sterol": ["HMGCR", "SQLE", "LDLR", "INSIG1", "SCD", "FASN"],
    "biliary": ["SOX9", "EPCAM", "CFTR", "KRT7", "SPP1", "KRT23", "CD24", "ITGB6"],
    "stress": ["CDKN1A", "GDF15", "KRT8", "SOX4", "LCN2"],
}
MIN_PLASTIC = 20

A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str})
dq = dq[dq.included].sort_values("gsm").reset_index(drop=True)
hep = A[(A.obs.label == "Hepatocytes").to_numpy() & A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
del A
genes = pd.Index(hep.var_names)
allg = [g for gs in GENES.values() for g in gs if g in genes]
pd.Series({k: ";".join(g for g in v if g in genes) for k, v in GENES.items()}).to_csv(OUT / "genes_used.tsv", sep="\t", header=["genes"])
X = sparse.csr_matrix(hep.X)
umi = np.asarray(X.sum(axis=1)).ravel().astype(float)
raw = X[:, [genes.get_loc(g) for g in allg]].toarray()
cell = pd.DataFrame(np.log1p(raw / umi[:, None] * 1e4), columns=allg)
cell["gsm"] = hep.obs.gsm.astype(str).to_numpy()
cell["umi"] = umi
cell["plastic"] = (raw[:, allg.index("SOX9")] + raw[:, allg.index("EPCAM")]) > 0
zp = lambda gs: cell[[g for g in gs if g in allg]].apply(lambda v: (v - v.mean()) / v.std()).mean(1)
cell["pp_minus_pc"] = zp(GENES["periportal"]) - zp(GENES["pericentral"])
cell["umi_decile"] = cell.groupby("gsm").umi.transform(lambda u: pd.qcut(u.rank(method="first"), 10, labels=False))
cell["zone_tertile"] = cell.groupby("gsm").pp_minus_pc.transform(lambda z: pd.qcut(z.rank(method="first"), 3, labels=False))

n = len(dq)
RELAB = []
for ix in combinations(range(n), int(dq.disease.sum())):
    v = np.zeros(n, int); v[list(ix)] = 1; RELAB.append(v)
RELAB = np.array(RELAB)


def eff(y):
    y = np.asarray(y, float)
    def beta(v):
        Xd = np.column_stack([np.ones(n), v, dq.sex_f]).astype(float)
        return np.nan if np.linalg.matrix_rank(Xd) < 3 else np.linalg.lstsq(Xd, y, rcond=None)[0][1]
    b = beta(dq.disease.to_numpy())
    nb = np.array([beta(v) for v in RELAB])
    return dict(beta=b, p_exact=np.nanmean(np.abs(nb) >= abs(b) - 1e-12),
                mean_healthy=y[dq.disease == 0].mean(), mean_disease=y[dq.disease == 1].mean())


# Q1 and Q3/Q4: donor summaries under several cell subsets, then disease effect
QS = [0.1, 0.25, 0.5, 0.75, 0.9]
rows = []
for gname in ["GNMT", "MAT1A", "BHMT", "ALB", "SERPINA1", "HMGCR", "SQLE"]:
    subsets = {"all": cell, "non_plastic": cell[~cell.plastic]}
    subsets.update({f"zone_tertile_{k}": cell[cell.zone_tertile == k] for k in range(3)})
    for sname, s in subsets.items():
        means = s.groupby("gsm")[gname].mean().reindex(dq.gsm)
        rows.append(dict(gene=gname, subset=sname, statistic="mean", **eff(means.to_numpy())))
        if sname == "all":
            for q in QS:
                qq = s.groupby("gsm")[gname].quantile(q).reindex(dq.gsm)
                rows.append(dict(gene=gname, subset=sname, statistic=f"q{int(q * 100):02d}", **eff(qq.to_numpy())))
            det = s.groupby("gsm")[gname].apply(lambda v: np.mean(v > 0)).reindex(dq.gsm)
            rows.append(dict(gene=gname, subset=sname, statistic="frac_detected", **eff(det.to_numpy())))
pd.DataFrame(rows).to_csv(OUT / "disease_effect_by_subset.tsv", sep="\t", index=False)

# Q2: plastic vs other hepatocytes within donor and UMI decile
prow = []
for g, o in cell.groupby("gsm"):
    p = o.plastic.to_numpy()
    r = dict(gsm=g, n_hep=len(o), n_plastic=int(p.sum()), frac_plastic=p.mean())
    if p.sum() >= MIN_PLASTIC:
        for gname in allg + ["pp_minus_pc"]:
            d = [od[gname][od.plastic].mean() - od[gname][~od.plastic].mean()
                 for _, od in o.groupby("umi_decile") if od.plastic.sum() >= 3 and (~od.plastic).sum() >= 3]
            r[gname] = np.mean(d) if d else np.nan
    prow.append(r)
pt = dq[["gsm", "id", "disease", "fibrosis"]].merge(pd.DataFrame(prow), on="gsm")
pt.to_csv(OUT / "plastic_vs_other_by_donor.tsv", sep="\t", index=False)
summ = []
for gname in allg + ["pp_minus_pc"]:
    if gname not in pt:
        continue
    v = pt[gname].dropna()
    vd = pt.loc[pt.disease == 1, gname].dropna()
    summ.append(dict(gene=gname, n_donors=len(v), median_diff=v.median(), n_positive=int((v > 0).sum()),
                     sign_test_p=stats.binomtest(int((v > 0).sum()), len(v)).pvalue if len(v) else np.nan,
                     n_masld=len(vd), median_diff_masld=vd.median(), n_positive_masld=int((vd > 0).sum())))
pd.DataFrame(summ).to_csv(OUT / "plastic_vs_other_summary.tsv", sep="\t", index=False)
fr = cell.groupby("gsm").plastic.mean().reindex(dq.gsm)
pd.DataFrame([dict(statistic="logit_frac_plastic", **eff(np.log(fr / (1 - fr)).to_numpy()))]).to_csv(
    OUT / "plastic_fraction_effect.tsv", sep="\t", index=False)
print(pd.DataFrame(summ).round(3).to_string(), flush=True)
