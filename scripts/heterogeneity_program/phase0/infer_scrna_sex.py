#!/usr/bin/env python
# ─────────────────────────────────────────────────────────────────────────────
# Phase 0e — scRNA donor sex inference (XIST / Y-gene k-means)
#
# GATE-G is inoperative for scRNA because donor sex is 100% blank. We infer it
# the same way the bulk pipeline does (01_sample_qc.R): per-donor pseudobulk of
# XIST vs a Y-chromosome gene, log1p, 2-means; female = higher XIST. XIST/DDX3Y
# were stripped from the dialogue pseudobulk (confounder genes), so we read the
# full atlas h5ad. We use per-donor MEAN expression on X as-is (the F/M signal is
# strongly bimodal and robust to normalization — no CPM needed for a binary call).
#
# Output: RNA-seq/results/heterogeneity_program/phase0/scrna_donor_sex.tsv
#         (sample, xist, y_gene, y_name, lxist, ly, inferred_sex, n_cells)
# Env: spatial (scanpy). Compute node only (40GB h5ad).
# ─────────────────────────────────────────────────────────────────────────────
import sys, numpy as np, pandas as pd, anndata as ad
from sklearn.cluster import KMeans

BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
H5   = f"{BASE}/Analysis/SingleCell/results_gpu_v2/pseudotime/.atlas_patched.h5ad"
OUT  = f"{BASE}/RNA-seq/results/heterogeneity_program/phase0/scrna_donor_sex.tsv"

print(f"Opening {H5} (backed) ...", flush=True)
A = ad.read_h5ad(H5, backed="r")
print(f"  n_obs={A.n_obs:,}  n_vars={A.n_vars:,}", flush=True)
print(f"  obs columns: {list(A.obs.columns)}", flush=True)

# donor column
donor_col = next((c for c in ["sample","donor","sample_id","donor_id","patient",
                              "orig.ident","Sample","Donor"] if c in A.obs.columns), None)
assert donor_col, "no donor column found in obs"
print(f"  donor column = '{donor_col}'  ({A.obs[donor_col].nunique()} donors)", flush=True)

# gene index resolution (var_names may be symbols or carry a gene_name column)
def gene_series():
    for col in ["gene_name","gene_symbol","symbol","features","_index"]:
        if col in A.var.columns:
            return pd.Index(A.var[col].astype(str))
    return pd.Index(A.var_names.astype(str))
gs = gene_series()
def find_gene(names):
    for nm in names:
        hit = np.where(gs == nm)[0]
        if len(hit): return int(hit[0]), nm
    return None, None

xist_idx, xist_nm = find_gene(["XIST"])
y_idx,   y_nm     = find_gene(["DDX3Y","RPS4Y1","UTY","EIF1AY","KDM5D","USP9Y"])
assert xist_idx is not None, "XIST not found in atlas var"
assert y_idx is not None, "no Y-chromosome gene found in atlas var"
print(f"  XIST @ {xist_idx}; Y gene = {y_nm} @ {y_idx}", flush=True)

# pull the two gene columns into memory (695k × 2 — small)
sub = A[:, [xist_idx, y_idx]].to_memory()
X = sub.X
X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
df = pd.DataFrame({"sample": A.obs[donor_col].astype(str).values,
                   "xist": np.asarray(X[:, 0]).ravel(),
                   "y":    np.asarray(X[:, 1]).ravel()})
g = df.groupby("sample").agg(xist=("xist","mean"), y=("y","mean"), n_cells=("xist","size"))
g["lxist"] = np.log1p(g["xist"]); g["ly"] = np.log1p(g["y"])
print(f"  per-donor table: {g.shape[0]} donors", flush=True)

km = KMeans(n_clusters=2, n_init=25, random_state=42).fit(g[["lxist","ly"]].values)
g["cluster"] = km.labels_
female_cluster = g.groupby("cluster")["lxist"].mean().idxmax()
g["inferred_sex"] = np.where(g["cluster"] == female_cluster, "F", "M")
g["y_name"] = y_nm

out = g.reset_index()[["sample","xist","y","y_name","lxist","ly","n_cells","inferred_sex"]]
out.to_csv(OUT, sep="\t", index=False)
print("\nSex distribution:\n", out["inferred_sex"].value_counts().to_string(), flush=True)
print(f"\nWrote {OUT}  ({len(out)} donors)", flush=True)
