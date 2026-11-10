"""DIG1 job 57 (exploratory, post hoc to jobs 53-56). Does the inflamed hepatocyte state carry bulk MASLD rises?
Several top bulk rises (CCL20 12th, TRIM31 30th, FCAMR 48th, IL32 76th of 6,390 up genes in the 774-unit
release) are CCL20-state genes. Two questions, per donor:
A. Lineage source: share of each gene's UMIs by cell label (Flex GSE298719 V5 labels, 13 donors; atlas
   cell_type for GSE244832, 18 donors, and GSE202379, 46 donors). Median over MASLD and over healthy donors.
B. Concentration inside hepatocytes: hepatocytes ranked by the CORE NF-kB score (mean within-donor z of log1p
   CP10k of the 12 CORE genes, the target gene left out), residualised on ranks of log UMI and log genes.
   Enrichment = (target UMIs in the top 5% / all hepatocyte target UMIs) / (all UMIs in the top 5% / all
   hepatocyte UMIs). 1 = no concentration. Null: the same for 100 random genes with detection matched to the
   target (0.67-1.5x; seed 20261028); excess = log2 enrichment - median log2 null. Donors with >= 20 target UMIs.
Descriptive only; donor is the unit; no test family is claimed.
Usage: python -I 57_state_bulk_source.py <new_out_dir>
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
CORE = ["CXCL1", "CXCL8", "BIRC3", "TNFAIP3", "ICAM1", "ZC3H12A", "RELB", "NFKBIA", "NFKB2", "MAP3K14", "UBD", "FCAMR"]
TARGETS = ["CCL20", "TRIM31", "FCAMR", "IL32", "CXCL10", "CXCL8", "CXCL1", "LOXL4", "KRT23", "UBD", "FSTL3", "BIRC3", "PIGR", "SAA2",
           "SPP1", "SOX9", "ALB", "APOA1", "CYP2E1", "HP", "IGHA1", "COL1A1", "PTPRC"]
N_NULL, TOP = 100, 0.05
rng = np.random.default_rng(20261028)


def pct(x):
    return (stats.rankdata(x) - 0.5) / len(x)


def cat(o, k):
    g = o[k]; c = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]); return c[g["codes"][:]]


def lineage_source(tag, X, gidx, label, donor, disease):
    rows = []
    for d in sorted(set(donor)):
        md = np.where(donor == d)[0]
        Xd = X[md]; Ld = label[md]
        for lab in sorted(set(Ld)):
            m = np.where(Ld == lab)[0]
            s = np.asarray(Xd[m][:, list(gidx.values())].sum(axis=0)).ravel()
            for g, v in zip(gidx, s):
                rows.append(dict(cohort=tag, donor=d, disease=disease[d], label=lab, gene=g, umis=float(v), n_cells=len(m),
                                 lineage_umis=float(Xd[m].sum())))
    R = pd.DataFrame(rows)
    R["share"] = R.umis / R.groupby(["donor", "gene"]).umis.transform("sum")
    R["cell_share"] = R.n_cells / R.groupby(["donor", "gene"]).n_cells.transform("sum")
    R["umi_share"] = R.lineage_umis / R.groupby(["donor", "gene"]).lineage_umis.transform("sum")
    return R


def concentration(tag, X, gnames, donor, disease):
    X = sparse.csr_matrix(X)
    genes = {}
    for i, g in enumerate(gnames):
        genes.setdefault(g, i)
    det = np.asarray((X > 0).mean(axis=0)).ravel()
    core = [g for g in CORE if g in genes]
    targets = [g for g in TARGETS if g in genes]
    nulls = {}
    for t in targets:
        cand = np.where((det >= 0.67 * det[genes[t]]) & (det <= 1.5 * det[genes[t]]))[0]
        cand = np.setdiff1d(cand, [genes[g] for g in set(core) | set(targets)])
        nulls[t] = rng.choice(cand, size=min(N_NULL, len(cand)), replace=False)
    rows = []
    for d in sorted(set(donor)):
        m = np.where(donor == d)[0]; Xd = X[m]
        umi = np.asarray(Xd.sum(axis=1)).ravel(); ng = np.asarray((Xd > 0).sum(axis=1)).ravel()
        L = np.log1p(Xd[:, [genes[g] for g in core]].toarray() * (1e4 / umi[:, None]))
        sd = L.std(0); sd[sd == 0] = np.nan; Z = (L - L.mean(0)) / sd
        B = np.column_stack([np.ones(len(m)), stats.rankdata(np.log(umi)), stats.rankdata(np.log(ng))])
        for t in targets:
            keep = [j for j, g in enumerate(core) if g != t]
            sc = np.nan_to_num(np.nanmean(Z[:, keep], 1))
            r = sc - B @ np.linalg.lstsq(B, sc, rcond=None)[0]
            top = pct(r) >= 1 - TOP
            umi_ratio = umi[top].sum() / umi.sum()
            def enr(cols):
                C = Xd[:, cols]
                tot = np.asarray(C.sum(axis=0)).ravel(); tt = np.asarray(C[top].sum(axis=0)).ravel()
                with np.errstate(divide="ignore", invalid="ignore"):
                    return np.where(tot >= 20, np.log2((tt / tot) / umi_ratio), np.nan), tot
            e_t, tot_t = enr([genes[t]])
            e_n, _ = enr(list(nulls[t]))
            rows.append(dict(cohort=tag, donor=d, disease=disease[d], gene=t, n_hep=len(m), target_umis=float(tot_t[0]),
                             top_umi_share=umi_ratio, log2_enrich=float(e_t[0]), null_median=float(np.nanmedian(e_n)),
                             excess=float(e_t[0] - np.nanmedian(e_n)) if np.isfinite(e_t[0]) else np.nan))
        print(tag, d, flush=True)
    return pd.DataFrame(rows)


# Flex
A = ad.read_h5ad(V5 / "derived/gse298719_qc_cells.h5ad")
dq = pd.read_csv(V5 / "outputs/donor_qc.tsv", sep="\t", dtype={"id": str}); dq = dq[dq.included]
A = A[A.obs.gsm.isin(dq.gsm).to_numpy()].copy()
dis_f = dict(zip(dq.gsm, dq.disease.astype(int)))
gn = list(A.var_names); gidx = {g: gn.index(g) for g in TARGETS if g in gn}
Xf = sparse.csr_matrix(A.X); lab = A.obs.label.astype(str).to_numpy(); don = A.obs.gsm.astype(str).to_numpy()
S = [lineage_source("flex", Xf, gidx, lab, don, dis_f)]
hm = lab == "Hepatocytes"
K = [concentration("flex", Xf[np.where(hm)[0]], gn, don[hm], dis_f)]
del A, Xf

# atlas cohorts
mapping = {}
for ds in ["GSE202379", "GSE244832"]:
    t = pd.read_csv(ROOT / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
    for r in t.itertuples(index=False):
        for s in str(r.rna_srrs).split(";"):
            if s.strip():
                mapping[s.strip()] = f"{ds}_{r.donor_id}"
meta = pd.read_csv(DONOR_META, sep="\t", dtype=str)
dis_a = {r.sample: int(r.disease_stage_coarse != "Healthy") for r in meta.itertuples() if r.dataset in ("GSE202379", "GSE244832")}
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
gidx = {}
for g in TARGETS:
    if g in gnames:
        gidx[g] = gnames.index(g)
for tag in ["GSE244832", "GSE202379"]:
    m = cds == tag
    S.append(lineage_source(tag, X[np.where(m)[0]], gidx, cct[m], cdon[m], dis_a))
    hm = m & (cct == "Hepatocytes")
    K.append(concentration(tag, X[np.where(hm)[0]], gnames, cdon[hm], dis_a))
S = pd.concat(S, ignore_index=True); K = pd.concat(K, ignore_index=True)
S.to_csv(OUT / "lineage_source_by_donor.tsv.gz", sep="\t", index=False)
K.to_csv(OUT / "hepatocyte_concentration_by_donor.tsv", sep="\t", index=False)
SS = (S.groupby(["cohort", "disease", "label", "gene"]).agg(median_share=("share", "median"), median_cell_share=("cell_share", "median"),
                                                          median_umi_share=("umi_share", "median"), n_donors=("donor", "nunique")).reset_index())
SS.to_csv(OUT / "summary_lineage_source.tsv", sep="\t", index=False)
KS = (K.groupby(["cohort", "disease", "gene"]).agg(n_donors=("excess", lambda v: int(np.isfinite(v).sum())), median_log2_enrich=("log2_enrich", "median"),
                                                 median_null=("null_median", "median"), median_excess=("excess", "median"),
                                                 n_pos=("excess", lambda v: int((v > 0).sum()))).reset_index())
KS.to_csv(OUT / "summary_concentration.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as fh:
    import importlib.metadata as im
    fh.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
h = SS[SS.label == "Hepatocytes"].pivot_table(index="gene", columns=["cohort", "disease"], values="median_share")
print("hepatocyte share of each gene's UMIs (median over donors)"); print(h.round(3).to_string())
top = SS[SS.disease == 1].sort_values("median_share", ascending=False).groupby(["cohort", "gene"]).head(2)
print(top.sort_values(["gene", "cohort"]).round(3).to_string())
print(KS.round(3).to_string())
