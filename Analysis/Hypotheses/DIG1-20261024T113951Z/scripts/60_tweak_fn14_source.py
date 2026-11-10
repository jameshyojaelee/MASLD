"""DIG1 job 60 (exploratory, post hoc to the receptor screen after job 59). Which cells make the ligands that
could drive the inflamed hepatocyte state, and do donors with more state cells have more of them?
Background: of 15 candidate receptors, only TNFRSF12A (Fn14, the TWEAK receptor) rises in CCL20+ hepatocytes in
all four cohorts, with the non-canonical NF-kB genes MAP3K14, RELB and NFKB2 and the TWEAK-inducible KRT23.
Cohorts: Flex GSE298719 (13 donors, V5 labels) and atlas snRNA GSE202379 (46) and GSE244832 (18). GSE256398 has
hepatocytes only and is not used.
A. Lineage source of each ligand: share of its UMIs per cell label, median over MASLD donors (job 57 method).
B. Donor level: y = fraction of hepatocytes with CCL20 > 0 (job 56 definition) and hepatocyte CPM of TNFRSF12A;
   x = ligand CPM in each lineage pseudobulk (donors with >= 20 cells of that lineage) and in all non-hepatocyte
   cells pooled. Spearman rho, and partial rho after regressing ranks on disease stage and log median hepatocyte
   UMI. BH within cohort and y over all x. Descriptive; donor is the unit; cross-sectional.
Usage: python -I 60_tweak_fn14_source.py <new_out_dir>
"""
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
V5 = ROOT / "Analysis/Hypotheses/V5-20261021T122723Z"
ATLAS = ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
DONOR_META = ROOT / ("Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/"
                     "metadata/donor_metadata_stage_corrected_donor.tsv")
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
LIGANDS = ["TNFSF12", "TNF", "LTA", "LTB", "TNFSF14", "CD40LG", "TNFSF13B", "TNFSF10", "TNFSF15", "IL1A", "IL1B", "OSM",
           "IL17A", "IL22", "IFNG", "CXCL8"]
HEP_GENES = ["CCL20", "TNFRSF12A", "KRT23"]
LINEAGES = ["Macrophages", "Endothelial cells", "Fibroblasts", "T cells", "Cholangiocytes", "Mono+mono derived cells",
            "Resident NK", "Circulating NK/NKT", "Plasma cells", "B cells"]
STAGE = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2, "Cirrhosis": 3}


def cat(o, k):
    g = o[k]; c = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]); return c[g["codes"][:]]


def lineage_source(tag, X, gidx, label, donor, disease):
    rows = []
    for d in sorted(set(donor)):
        md = np.where(donor == d)[0]; Xd = X[md]; Ld = label[md]
        for lab in sorted(set(Ld)):
            m = np.where(Ld == lab)[0]
            s = np.asarray(Xd[m][:, list(gidx.values())].sum(axis=0)).ravel()
            for g, v in zip(gidx, s):
                rows.append(dict(cohort=tag, donor=d, disease=disease[d], label=lab, gene=g, umis=float(v), n_cells=len(m)))
    R = pd.DataFrame(rows)
    R["share"] = R.umis / R.groupby(["donor", "gene"]).umis.transform("sum")
    return R


def donor_table(tag, X, gnames, label, donor, stage):
    gi = {}
    for i, g in enumerate(gnames):
        gi.setdefault(g, i)
    rows = []
    for d in sorted(set(donor)):
        md = np.where(donor == d)[0]; Xd = X[md]; Ld = label[md]
        h = np.where(Ld == "Hepatocytes")[0]
        if len(h) < 30:
            continue
        Xh = Xd[h]; umi = np.asarray(Xh.sum(axis=1)).ravel(); tot = umi.sum()
        rec = dict(cohort=tag, donor=d, stage=stage[d], n_hep=len(h), median_umi=float(np.median(umi)),
                   frac_CCL20=float((Xh[:, gi["CCL20"]].toarray().ravel() > 0).mean()))
        for g in HEP_GENES:
            if g in gi:
                rec[f"hep_cpm_{g}"] = float(Xh[:, gi[g]].sum() / tot * 1e6)
        nh = np.where(Ld != "Hepatocytes")[0]
        groups = {lab: np.where(Ld == lab)[0] for lab in LINEAGES} | {"nonhep": nh}
        for lab, idx in groups.items():
            if len(idx) < 20:
                continue
            Xl = Xd[idx]; t = Xl.sum()
            for g in LIGANDS:
                if g in gi:
                    rec[f"{lab}|{g}"] = float(Xl[:, gi[g]].sum() / t * 1e6)
        rows.append(rec)
    return pd.DataFrame(rows)


def partial_rho(x, y, Z):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    B = np.column_stack([np.ones(len(x))] + [stats.rankdata(z) for z in Z.T])
    ex = rx - B @ np.linalg.lstsq(B, rx, rcond=None)[0]; ey = ry - B @ np.linalg.lstsq(B, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def bh(p):
    p = np.asarray(p, float); q = np.full(len(p), np.nan); ok = np.isfinite(p)
    if ok.sum():
        pv = p[ok]; o = np.argsort(pv); r = pv[o] * ok.sum() / np.arange(1, ok.sum() + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]; tmp = np.empty(ok.sum()); tmp[o] = np.minimum(r, 1); q[ok] = tmp
    return q


S, T = [], []
# Flex
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str}); dq = dq[dq.included]
A = A[A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
dis_f = dict(zip(dq.gsm, dq.disease.astype(int)))
gn = list(A.var_names); gidx = {g: gn.index(g) for g in LIGANDS + HEP_GENES if g in gn}
Xf = sparse.csr_matrix(A.X); lab = A.obs.label.astype(str).to_numpy(); don = A.obs.gsm.astype(str).to_numpy()
S.append(lineage_source("flex", Xf, gidx, lab, don, dis_f))
T.append(donor_table("flex", Xf, gn, lab, don, dis_f))
del A, Xf
# atlas
mapping = {}
for ds in ["GSE202379", "GSE244832"]:
    t = pd.read_csv(ROOT / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
    for r in t.itertuples(index=False):
        for s in str(r.rna_srrs).split(";"):
            if s.strip():
                mapping[s.strip()] = f"{ds}_{r.donor_id}"
meta = pd.read_csv(DONOR_META, sep="\t", dtype=str)
meta = meta[meta.dataset.isin(["GSE202379", "GSE244832"])]
dis_a = {r.sample: int(r.disease_stage_coarse != "Healthy") for r in meta.itertuples()}
stage_a = {r.sample: STAGE[r.disease_stage_coarse] for r in meta.itertuples()}
with h5py.File(ATLAS, "r") as f:
    o = f["obs"]
    ds_all = cat(o, "dataset"); sm_all = cat(o, "sample"); ct_all = cat(o, "cell_type")
    gnames = [x.decode() for x in f["raw"]["var"]["_index"][:]]
    rows = np.where(np.isin(ds_all, ["GSE202379", "GSE244832"]))[0]
    Xr = f["raw"]["X"]; indptr = Xr["indptr"][:]
    blocks = []
    for b0 in range(0, len(indptr) - 1, 100_000):
        b1 = min(b0 + 100_000, len(indptr) - 1)
        sel = rows[(rows >= b0) & (rows < b1)]
        if len(sel) == 0:
            continue
        lo, hi = indptr[b0], indptr[b1]
        B = sparse.csr_matrix((Xr["data"][lo:hi].astype(np.float32), Xr["indices"][lo:hi], indptr[b0:b1 + 1] - lo), shape=(b1 - b0, len(gnames)))
        blocks.append(B[sel - b0])
    X = sparse.vstack(blocks).tocsr(); del blocks
cds = ds_all[rows]; cct = ct_all[rows]; cdon = np.array([mapping[s] for s in sm_all[rows]])
gidx = {g: gnames.index(g) for g in LIGANDS + HEP_GENES if g in gnames}
for tag in ["GSE244832", "GSE202379"]:
    m = np.where(cds == tag)[0]
    S.append(lineage_source(tag, X[m], gidx, cct[m], cdon[m], dis_a))
    T.append(donor_table(tag, X[m], gnames, cct[m], cdon[m], stage_a))
S = pd.concat(S, ignore_index=True); T = pd.concat(T, ignore_index=True)
S.to_csv(OUT / "ligand_lineage_source_by_donor.tsv.gz", sep="\t", index=False)
T.to_csv(OUT / "donor_ligand_table.tsv", sep="\t", index=False)
SS = (S[S.disease == 1].groupby(["cohort", "label", "gene"]).agg(median_share=("share", "median"), n_donors=("donor", "nunique"),
                                                                 total_umis=("umis", "sum")).reset_index())
SS.to_csv(OUT / "summary_ligand_source.tsv", sep="\t", index=False)
L = []
for co, D in T.groupby("cohort"):
    for y in ["frac_CCL20", "hep_cpm_TNFRSF12A"]:
        res = []
        for x in [c for c in D.columns if "|" in c]:
            dd = D[[x, y, "stage", "median_umi"]].dropna()
            if len(dd) < 8 or dd[x].nunique() < 3:
                continue
            rho, p = stats.spearmanr(dd[x], dd[y])
            res.append(dict(cohort=co, y=y, x=x, lineage=x.split("|")[0], ligand=x.split("|")[1], n=len(dd), rho=rho, p=p,
                            partial_rho=partial_rho(dd[x].to_numpy(), dd[y].to_numpy(), dd[["stage", "median_umi"]].to_numpy())))
        if res:
            R = pd.DataFrame(res); R["q_bh"] = bh(R.p); L.append(R)
L = pd.concat(L, ignore_index=True)
L.to_csv(OUT / "donor_ligand_correlations.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
P = SS[SS.gene.isin(["TNFSF12", "TNF", "LTB", "TNFSF14", "IL1B", "OSM"])].pivot_table(index=["gene", "label"], columns="cohort", values="median_share")
print(P[P.max(axis=1) > 0.03].round(3).to_string())
print(L.sort_values("p").head(40).round(4).to_string())
print(L[L.ligand == "TNFSF12"].round(4).to_string())
