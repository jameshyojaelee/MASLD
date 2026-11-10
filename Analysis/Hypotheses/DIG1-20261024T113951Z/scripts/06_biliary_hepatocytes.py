"""DIG1 job 6 (exploratory; follows A2). Are biliary transcripts in MASLD hepatocytes carried by
hepatocytes or by mixed droplets? GSE298719 Flex, CellTypist-labelled hepatocytes.

Per cell: raw counts of biliary markers (KRT7, KRT19, EPCAM, SOX9, CFTR, SPP1, ANXA4, CD24) and of
non-parenchymal markers (PTPRC, PECAM1, CD74, COL1A1, LYVE1, STAB2). A 'biliary+' hepatocyte has
>= 2 biliary UMIs from >= 2 distinct biliary genes. For each donor: fraction biliary+, and for
biliary+ vs biliary- cells matched on UMI decile, the hepatocyte identity per UMI (ALB, APOA1, TF,
SERPINA1, HP, APOC3) and the non-parenchymal fraction. A doublet predicts biliary+ cells with
higher UMI and diluted hepatocyte identity; a hepatocyte state predicts normal UMI and intact
hepatocyte identity. Donor is the unit; exact relabeling over 1,716 splits for disease contrasts.
Usage: python -I 06_biliary_hepatocytes.py
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
OUT = RUN / "out_biliary"
OUT.mkdir(exist_ok=False)
BIL = ["KRT7", "KRT19", "EPCAM", "SOX9", "CFTR", "SPP1", "ANXA4", "CD24"]
NPC = ["PTPRC", "PECAM1", "CD74", "COL1A1", "LYVE1", "STAB2", "DCN", "CD163"]
HEP = ["ALB", "APOA1", "TF", "SERPINA1", "HP", "APOC3"]

A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str})
dq = dq[dq.included].sort_values("gsm").reset_index(drop=True)
genes = pd.Index(A.var_names)
for gs in (BIL, NPC, HEP):
    missing = [g for g in gs if g not in genes]
    print("missing:", missing, flush=True)
BIL, NPC, HEP = ([g for g in gs if g in genes] for gs in (BIL, NPC, HEP))
keep = A.obs.gsm.isin(dq.gsm).to_numpy()
A = A[keep].copy()
X = sparse.csr_matrix(A.X)
umi = np.asarray(X.sum(axis=1)).ravel().astype(float)
lab = A.obs.label.astype(str).to_numpy()
gsm = A.obs.gsm.astype(str).to_numpy()


def sub(gs):
    return X[:, [genes.get_loc(g) for g in gs]].toarray()


B, N, H = sub(BIL), sub(NPC), sub(HEP)
bil_umi, bil_ngenes = B.sum(1), (B > 0).sum(1)
cell = pd.DataFrame(dict(gsm=gsm, label=lab, umi=umi, bil_umi=bil_umi, bil_ngenes=bil_ngenes,
                         npc_frac=N.sum(1) / umi, hep_frac=H.sum(1) / umi))
cell["biliary_pos"] = (cell.bil_umi >= 2) & (cell.bil_ngenes >= 2)
for g in BIL:
    cell[g] = B[:, BIL.index(g)]

# Reference: per-UMI biliary and hepatocyte content in true cholangiocytes and hepatocytes
ref = cell.groupby("label").agg(n=("umi", "size"), umi_median=("umi", "median"), bil_per_1k=("bil_umi", lambda v: 0),
                                frac_biliary_pos=("biliary_pos", "mean"), hep_frac_median=("hep_frac", "median"),
                                npc_frac_median=("npc_frac", "median"))
ref["bil_per_1k"] = cell.assign(r=cell.bil_umi / cell.umi * 1e3).groupby("label").r.mean()
ref.sort_values("n", ascending=False).to_csv(OUT / "label_reference_profiles.tsv", sep="\t")

hep = cell[cell.label == "Hepatocytes"].copy()
hep["umi_decile"] = hep.groupby("gsm").umi.transform(lambda u: pd.qcut(u.rank(method="first"), 10, labels=False))
rows = []
for g, o in hep.groupby("gsm"):
    pos = o.biliary_pos.to_numpy()
    r = dict(gsm=g, n_hep=len(o), n_biliary_pos=int(pos.sum()), frac_biliary_pos=pos.mean())
    # decile-matched contrasts: mean over deciles that have both kinds of cell
    d_umi, d_hep, d_npc = [], [], []
    for _, od in o.groupby("umi_decile"):
        p = od.biliary_pos.to_numpy()
        if p.sum() >= 3 and (~p).sum() >= 3:
            d_hep.append(np.log2(od.hep_frac[p].median() / od.hep_frac[~p].median()))
            d_npc.append(od.npc_frac[p].mean() - od.npc_frac[~p].mean())
    r.update(log2_umi_pos_vs_neg=np.log2(o.umi[pos].median() / o.umi[~pos].median()) if pos.sum() >= 3 else np.nan,
             log2_hepfrac_pos_vs_neg_decile=np.mean(d_hep) if d_hep else np.nan,
             npcfrac_diff_pos_vs_neg_decile=np.mean(d_npc) if d_npc else np.nan, n_deciles_used=len(d_hep),
             npc_pos_among_biliary_pos=np.mean(o.npc_frac[pos] > 0) if pos.sum() else np.nan,
             npc_pos_among_biliary_neg=np.mean(o.npc_frac[~pos] > 0))
    for gname in BIL:
        r[f"frac_{gname}_pos"] = np.mean(o[gname] > 0)
    rows.append(r)
donor = dq[["gsm", "id", "group", "disease", "sex_f", "fibrosis"]].merge(pd.DataFrame(rows), on="gsm")
chol_frac = cell.assign(c=cell.label == "Cholangiocytes").groupby("gsm").c.mean()
donor["cholangiocyte_frac_all_cells"] = donor.gsm.map(chol_frac)
donor.to_csv(OUT / "biliary_hepatocytes_by_donor.tsv", sep="\t", index=False)

n = len(donor)
RELAB = []
for ix in combinations(range(n), int(donor.disease.sum())):
    v = np.zeros(n, int); v[list(ix)] = 1; RELAB.append(v)
RELAB = np.array(RELAB)


def eff(y):
    def beta(v):
        Xd = np.column_stack([np.ones(n), v, donor.sex_f]).astype(float)
        if np.linalg.matrix_rank(Xd) < 3:
            return np.nan
        return np.linalg.lstsq(Xd, y, rcond=None)[0][1]
    b = beta(donor.disease.to_numpy())
    nb = np.array([beta(v) for v in RELAB])
    return dict(beta=b, p_exact=np.nanmean(np.abs(nb) >= abs(b) - 1e-12))


res = []
for col in ["frac_biliary_pos", "log2_umi_pos_vs_neg", "log2_hepfrac_pos_vs_neg_decile", "npcfrac_diff_pos_vs_neg_decile"]:
    y = donor[col].to_numpy(float)
    if col == "frac_biliary_pos":
        y = np.log(np.clip(y, 1e-4, None) / np.clip(1 - y, 1e-4, None))
    ok = ~np.isnan(y)
    summ = dict(statistic=col, median_healthy=np.nanmedian(donor[col][donor.disease == 0]),
                median_disease=np.nanmedian(donor[col][donor.disease == 1]), n_nonmissing=int(ok.sum()))
    if col == "frac_biliary_pos":
        summ.update(eff(y))
    else:
        v = y[ok]
        summ.update(n_positive=int((v > 0).sum()), sign_test_p=stats.binomtest(int((v > 0).sum()), len(v)).pvalue)
    res.append(summ)
pd.DataFrame(res).to_csv(OUT / "biliary_summary.tsv", sep="\t", index=False)
cell[cell.label.isin(["Hepatocytes", "Cholangiocytes"])].to_csv(OUT / "cells_hep_chol_markers.tsv.gz", sep="\t", index=False)
print(pd.DataFrame(res).to_string(), flush=True)
