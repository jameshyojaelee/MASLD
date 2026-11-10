"""DIG1 job 64 (exploratory, post hoc to job 61). Which cells make GDF15 and the p53-target genes that rise in
CCL20+ hepatocytes? Koo 2017 (Liver Int, doi 10.1111/liv.13587) found palmitate raised GDF15 mRNA in Kupffer
cells but not hepatocytes. Method of job 60 part A: share of each gene's UMIs per cell label, per donor; median
over MASLD donors and over healthy donors; cell share of each label for comparison. Descriptive only.
Cohorts: Flex GSE298719 (V5 labels) and atlas snRNA GSE202379, GSE244832 (atlas cell_type).
Usage: python -I 64_p53_gene_source.py <new_out_dir>
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
LIGANDS = ["GDF15", "CDKN1A", "MDM2", "TNFRSF10B", "TNFRSF10D", "TP53", "PLK3", "BCL2L1", "SERPINE1", "FGF21", "MYC", "ATF3"]
HEP_GENES = []
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


S = []
# Flex
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str}); dq = dq[dq.included]
A = A[A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
dis_f = dict(zip(dq.gsm, dq.disease.astype(int)))
gn = list(A.var_names); gidx = {g: gn.index(g) for g in LIGANDS + HEP_GENES if g in gn}
Xf = sparse.csr_matrix(A.X); lab = A.obs.label.astype(str).to_numpy(); don = A.obs.gsm.astype(str).to_numpy()
S.append(lineage_source("flex", Xf, gidx, lab, don, dis_f))
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
S = pd.concat(S, ignore_index=True)
S.to_csv(OUT / "gene_lineage_source_by_donor.tsv.gz", sep="\t", index=False)
cs = S.drop_duplicates(["cohort", "donor", "label"])[["cohort", "donor", "disease", "label", "n_cells"]].copy()
cs["cell_share"] = cs.n_cells / cs.groupby("donor").n_cells.transform("sum")
S = S.merge(cs[["donor", "label", "cell_share"]], on=["donor", "label"])
SS = (S.groupby(["cohort", "disease", "label", "gene"]).agg(median_share=("share", "median"), median_cell_share=("cell_share", "median"),
                                                            n_donors=("donor", "nunique"), total_umis=("umis", "sum")).reset_index())
SS.to_csv(OUT / "summary_gene_source.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
P = SS[SS.disease == 1].pivot_table(index=["gene", "label"], columns="cohort", values="median_share")
print(P[P.max(axis=1) > 0.05].round(3).to_string())
C = SS[(SS.disease == 1) & (SS.gene == "GDF15")].pivot_table(index="label", columns="cohort", values="median_cell_share")
print(C.round(3).to_string())
