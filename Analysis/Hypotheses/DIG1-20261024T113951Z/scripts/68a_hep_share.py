"""DIG1 job 68a (post hoc to job 68). Hepatocyte share of each gene's UMIs, for every gene, in MASLD donors.
Per donor: hepatocyte-labelled UMIs / all UMIs of the gene; median over MASLD donors per cohort.
Cohorts and labels as job 64: Flex GSE298719 (V5 labels) and atlas snRNA GSE202379, GSE244832 (atlas cell_type).
Used by job 68b to keep proteins whose RNA comes mainly from hepatocytes. Descriptive.
Usage: python -I 68a_hep_share.py <new_out_dir>
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




def hep_share(tag, X, genes, label, donor, disease):
    ishep = np.array(["epatocyte" in str(l) for l in label])
    H, T, ds = [], [], []
    for d in sorted(set(donor)):
        if not disease.get(d, 0):
            continue
        m = np.where(donor == d)[0]; Xd = X[m]
        T.append(np.asarray(Xd.sum(axis=0)).ravel()); H.append(np.asarray(Xd[ishep[m]].sum(axis=0)).ravel()); ds.append(d)
    H = np.vstack(H); T = np.vstack(T)
    with np.errstate(invalid="ignore", divide="ignore"):
        sh = np.where(T > 0, H / T, np.nan)
    return pd.DataFrame(dict(cohort=tag, gene=genes, median_hep_share=np.nanmedian(sh, axis=0), n_donors_detected=(T > 0).sum(axis=0),
                             n_donors=len(ds), total_umis=T.sum(axis=0), hep_labels=";".join(sorted(set(label[ishep])))))


S = []
# Flex
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str}); dq = dq[dq.included]
A = A[A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
dis_f = dict(zip(dq.gsm, dq.disease.astype(int)))
gn = np.array(A.var_names)
Xf = sparse.csr_matrix(A.X); lab = A.obs.label.astype(str).to_numpy(); don = A.obs.gsm.astype(str).to_numpy()
S.append(hep_share("flex", Xf, gn, lab, don, dis_f))
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
for tag in ["GSE244832", "GSE202379"]:
    m = np.where(cds == tag)[0]
    S.append(hep_share(tag, X[m], np.array(gnames), cct[m], cdon[m], dis_a))
S = pd.concat(S, ignore_index=True)
S.to_csv(OUT / "gene_hep_share.tsv.gz", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
print(S.groupby("cohort").agg(n=("gene", "size"), n_donors=("n_donors", "max"), labels=("hep_labels", "first")).to_string())
for g in ["ALB", "APOA1", "CYP2E1", "PALLD", "ANXA2", "PLOD2", "PLIN2", "AKR1B10", "KRT7", "SOD2", "NFKB2", "ICAM1", "PTPRC", "COL1A1"]:
    print(g, S[S.gene == g][["cohort", "median_hep_share", "n_donors_detected"]].round(3).values.tolist())
