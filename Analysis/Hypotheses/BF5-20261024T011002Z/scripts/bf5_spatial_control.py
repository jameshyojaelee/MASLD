#!/usr/bin/env python
"""BF5 Part A (PRESPEC_BF5.md): are the hepatocyte-made continuum-up genes more central in the lobule than other
hepatocyte genes of similar detection? Spot scoring and lobule position as in BF1."""
import argparse
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--draws", type=int, default=1000)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CG0 = f"{ROOT}/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z/genes.tsv"
R9D = (f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
       "catalog_r9d-20261016T105059Z/build_out/catalog_r9d_gene_table.parquet")
SP = f"{ROOT}/Analysis/Spatial/results/preprocessed"
DMAP = f"{ROOT}/Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv"
GROUPS = f"{ROOT}/figures/candidates/continuum-biology-20261024T002544Z/data/visium_up_gene_groups.tsv"
TOPTYPE = f"{ROOT}/Analysis/Hypotheses/BF4-20261024T010648Z/out-20261024T011020Z/gene_top_celltype.tsv"
H88_MARKERS = {"CYP2E1", "GLUL", "OAT", "CYP1A2", "CYP3A4", "CPS1", "HAL", "SDS", "ASS1", "PCK1",
               "KRT7", "KRT19", "EPCAM", "SOX9", "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"}
SEED = 20261024
rng = np.random.default_rng(SEED)

g = pd.read_csv(CG0, sep="\t")
g["signature_gene"] = g.signature_gene.fillna(False).astype(bool)
r9 = pd.read_parquet(R9D, columns=["ensembl_id", "zon_protein_call_weiss2026", "zon_mrna_call_yakubovsky2026"]
                     ).drop_duplicates("ensembl_id")
g = g.merge(r9, left_on="g", right_on="ensembl_id", how="left")
sym = g.gene_name
called_any = set(sym[g.zon_protein_call_weiss2026.isin(["central", "portal"]) |
                     g.zon_mrna_call_yakubovsky2026.isin(["central", "portal"])].dropna())
CEN = set(sym[(~g.signature_gene) & (g.zon_protein_call_weiss2026 == "central")].dropna())
POR = set(sym[(~g.signature_gene) & (g.zon_protein_call_weiss2026 == "portal")].dropna())
grp = pd.read_csv(GROUPS, sep="\t")
TARGET = set(grp.gene[grp.group == "Hepatocytes"])
tt = pd.read_csv(TOPTYPE, sep="\t")
POOL = (set(tt.gene[tt.hep_ratio >= 0.5]) - set(sym[g.signature_gene].dropna()) - called_any - H88_MARKERS - TARGET)
NEED = CEN | POR | TARGET | POOL
print("targets", sorted(TARGET), "pool", len(POOL), flush=True)


def section_test(X, syms, tot_obs, ngen_obs):
    tot = np.asarray(X.sum(1)).ravel()
    keep = np.isin(syms, list(NEED))
    uu, inv = np.unique(syms[keep], return_inverse=True)
    M = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(uu)))
    C = (X[:, keep] @ M).toarray().astype(np.float64)
    det = (C > 0).mean(0)
    Lg = np.log1p(C / tot[:, None] * 1e4)
    del C
    sd = Lg.std(0)
    Zm = (Lg - Lg.mean(0)) / np.where(sd > 0, sd, 1)
    del Lg
    cov = np.column_stack([np.ones(len(tot)), np.log1p(tot_obs), np.log1p(ngen_obs)])
    B, *_ = np.linalg.lstsq(cov, Zm, rcond=None)
    R = Zm - cov @ B
    del Zm
    ok = det >= 0.01
    idx = {x: i for i, x in enumerate(uu) if ok[i]}
    ps = -(R[:, [idx[x] for x in CEN if x in idx]].mean(1) - R[:, [idx[x] for x in POR if x in idx]].mean(1))
    ps = (ps - ps.mean()) / ps.std()

    def corr(cols):
        s = R[:, cols].mean(1)
        return float(np.mean(ps * (s - s.mean()) / s.std()))

    it = np.array([idx[x] for x in TARGET if x in idx])
    if len(it) == 0:
        return None
    ip = np.array([idx[x] for x in POOL if x in idx])
    cand = np.concatenate([it, ip])
    dec = pd.qcut(pd.Series(det[cand]).rank(method="first"), 10, labels=False).values
    dmap = dict(zip(cand, dec))
    need = pd.Series([dmap[i] for i in it]).value_counts()
    pools = {d: ip[[dmap[i] == d for i in ip]] for d in need.index}
    null = np.array([corr(np.concatenate([rng.choice(pools[d], need[d], replace=len(pools[d]) < need[d])
                                          for d in need.index])) for _ in range(a.draws)])
    r_obs = corr(it)
    return dict(n_target=len(it), targets=",".join(sorted(uu[it])), n_pool=len(ip), r_obs=r_obs,
                null_mean_r=float(null.mean()), diff_fisher_z=float(np.arctanh(r_obs) - np.arctanh(null).mean()),
                pct_in_null=float((null < r_obs).mean()))


rows = []
dmap_ = pd.read_csv(DMAP, sep="\t")
for arr in dmap_.array_id[dmap_.group == "primary"]:
    A = ad.read_h5ad(f"{SP}/HRA007511_starsolo_hexscaled/{arr}.h5ad")
    X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
    o = section_test(X, A.var.gene_symbol.astype(str).values, A.obs.total_counts.values, A.obs.n_genes.values)
    if o:
        rows.append(dict(o, section=arr, dataset="HMSMA"))
        print(rows[-1], flush=True)
    del A, X
A = ad.read_h5ad(f"{SP}/merged_spatial_vu_raw.h5ad")
Xall = A.layers["counts"]
Xall = Xall.tocsr() if sp.issparse(Xall) else sp.csr_matrix(Xall)
for s in A.obs.sample_id.unique():
    m = (A.obs.sample_id == s).values
    ob = A.obs[m]
    o = section_test(Xall[m], A.var_names.astype(str).values, ob.total_counts.values, ob.n_genes_by_counts.values)
    if o:
        rows.append(dict(o, section=str(s), dataset="Vu 2025"))
        print(rows[-1], flush=True)
del A, Xall

v = pd.DataFrame(rows)
v.to_csv(f"{a.out}/partA_sections.tsv", sep="\t", index=False)
RES = {}
for ds, d in v.groupby("dataset"):
    x = d.diff_fisher_z.values
    flips = rng.choice([-1, 1], size=(100000, len(x))) @ x / len(x)
    RES[ds] = dict(n_sections=len(x), mean_diff_fisher_z=float(x.mean()), n_below_null=int((x < 0).sum()),
                   signflip_p=float((1 + np.sum(np.abs(flips) >= abs(x.mean()))) / 100001),
                   median_r_obs=float(d.r_obs.median()), median_null_r=float(d.null_mean_r.median()))
RES["reading_met"] = all(RES[ds]["mean_diff_fisher_z"] < 0 and RES[ds]["signflip_p"] < 0.05 for ds in ["HMSMA", "Vu 2025"])
json.dump(RES, open(f"{a.out}/partA.json", "w"), indent=1)
print(json.dumps(RES, indent=1), flush=True)
