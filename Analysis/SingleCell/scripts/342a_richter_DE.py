#!/usr/bin/env python
"""
342a_richter_DE.py — Phase 1.1: derive polyploid signature from Richter 2021 snRNA-seq2.

Strategy:
- Pseudobulk per (individual × ploidy) — paired design within 6 mice
- DE: 4n vs 2n hepatocytes using DESeq2 (via pydeseq2)
- Use individual as covariate
- Extract top 100 up + top 100 down genes at FDR<0.05, |log2FC|>0.5
- Save per-gene table + ranked signature

Inputs:  data/external/richter2021/richter2021_hepatocytes.h5ad
Outputs: data/ploidy_signatures/richter2021_de_full.tsv
         data/ploidy_signatures/richter2021_polyploid_{up,down}.tsv
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from pathlib import Path

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT",
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
os.chdir(ROOT)

H5AD = ROOT / "data/external/richter2021/richter2021_hepatocytes.h5ad"
OUT_DIR = ROOT / "data/ploidy_signatures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"[{pd.Timestamp.now()}] Loading Richter h5ad", flush=True)
a = ad.read_h5ad(H5AD)
print(f"  shape: {a.shape}", flush=True)

# Filter: passing QC + hepatocytes + canonical 2n/4n
mask = (
    (a.obs["passing_QC"].astype(str) == "pass")
    & (a.obs["inferred_cell_label"].astype(str) == "Hepatocytes")
    & a.obs["Ploidy"].astype(str).isin(["2n", "4n"])
)
a = a[mask].copy()
print(f"[{pd.Timestamp.now()}] After filtering: {a.shape}", flush=True)
print(f"  Ploidy: {a.obs['Ploidy'].value_counts().to_dict()}", flush=True)
print(f"  Individuals: {a.obs['ID.Individual'].value_counts().to_dict()}", flush=True)

# Pseudobulk: sum counts per (individual × ploidy)
print(f"[{pd.Timestamp.now()}] Building pseudobulk", flush=True)
a.obs["pb_group"] = (a.obs["ID.Individual"].astype(str).str.replace("[^A-Za-z0-9]", "_", regex=True)
                    + "@@" + a.obs["Ploidy"].astype(str))
groups = a.obs["pb_group"].unique()
X = a.X if not sp.issparse(a.X) else a.X.toarray()
# scaled (Smart-seq2 floats); convert to integer counts by rounding for DESeq2
print(f"  X dtype: {X.dtype}  min/max: {X.min():.2f}/{X.max():.2f}", flush=True)
counts_int = np.rint(X).astype(np.int64)

pb_data = {}
pb_meta = []
for g in groups:
    cells = a.obs["pb_group"] == g
    n = int(cells.sum())
    if n < 5:
        print(f"  SKIP {g} (n={n} cells < 5)", flush=True)
        continue
    sub = counts_int[cells.values]
    pb_data[g] = sub.sum(axis=0)
    ind, ploidy = g.rsplit("@@", 1)
    pb_meta.append({"sample": g, "individual": ind, "ploidy": ploidy, "n_cells": n})

pb_counts = pd.DataFrame(pb_data, index=a.var_names).T  # samples × genes
pb_meta = pd.DataFrame(pb_meta).set_index("sample")
print(f"  pseudobulk samples: {pb_counts.shape}", flush=True)
print(pb_meta.to_string())

# Filter genes (min expression for DESeq2)
keep_genes = (pb_counts > 5).sum(axis=0) >= max(3, pb_counts.shape[0] // 3)
pb_counts = pb_counts.loc[:, keep_genes]
print(f"  genes after filter (>5 in >=1/3 samples): {pb_counts.shape[1]}", flush=True)

# Map Ensembl → gene symbol
gene_meta = a.var[["gene_name", "biotype"]].copy()
gene_meta = gene_meta.loc[pb_counts.columns]

# DESeq2
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    print(f"[{pd.Timestamp.now()}] Running pyDESeq2", flush=True)
    inference = None
    dds = DeseqDataSet(
        counts=pb_counts.astype(int),
        metadata=pb_meta,
        design_factors=["individual", "ploidy"],
        ref_level=["ploidy", "2n"],
        refit_cooks=True,
        quiet=True,
    )
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["ploidy", "4n", "2n"], quiet=True)
    stat.summary()
    res = stat.results_df.copy()
    res["gene_symbol"] = gene_meta.loc[res.index, "gene_name"].values
    res["biotype"] = gene_meta.loc[res.index, "biotype"].values
    res["gene_ensembl"] = res.index
    res = res.sort_values("padj")
    res.to_csv(OUT_DIR / "richter2021_de_full.tsv", sep="\t", index=False)
    print(f"  wrote richter2021_de_full.tsv ({len(res)} genes)", flush=True)
except Exception as e:
    print(f"  pyDESeq2 failed ({e}); falling back to Wilcoxon rank-sum on log1p-normalized pseudobulk", flush=True)
    # Fallback: normalize pseudobulk and do simple rank tests
    from scipy.stats import ranksums
    norm = np.log1p(pb_counts.div(pb_counts.sum(axis=1), axis=0) * 1e6)
    rows = []
    for g in norm.columns:
        v2 = norm.loc[pb_meta["ploidy"] == "2n", g].values
        v4 = norm.loc[pb_meta["ploidy"] == "4n", g].values
        if v2.size < 2 or v4.size < 2:
            continue
        s, p = ranksums(v4, v2)
        rows.append({"gene_ensembl": g, "log2FC": np.log2((v4.mean() + 1e-6) / (v2.mean() + 1e-6)),
                     "stat": s, "pvalue": p,
                     "gene_symbol": gene_meta.loc[g, "gene_name"], "biotype": gene_meta.loc[g, "biotype"]})
    res = pd.DataFrame(rows).sort_values("pvalue")
    from statsmodels.stats.multitest import multipletests
    res["padj"] = multipletests(res["pvalue"], method="fdr_bh")[1]
    res.to_csv(OUT_DIR / "richter2021_de_full.tsv", sep="\t", index=False)
    print(f"  fallback wrote {len(res)} rows", flush=True)

# Extract signatures
res = pd.read_csv(OUT_DIR / "richter2021_de_full.tsv", sep="\t")
LFC = "log2FoldChange" if "log2FoldChange" in res.columns else "log2FC"
res = res.dropna(subset=[LFC, "padj"])
sig_up = res[(res["padj"] < 0.05) & (res[LFC] > 0.5)].sort_values(LFC, ascending=False).head(150)
sig_dn = res[(res["padj"] < 0.05) & (res[LFC] < -0.5)].sort_values(LFC, ascending=True).head(150)
print(f"\n=== Richter signature ===", flush=True)
print(f"  up_in_4n: {len(sig_up)} genes", flush=True)
print(f"  down_in_4n: {len(sig_dn)} genes", flush=True)
print(f"  top 10 up (gene_symbol): {sig_up['gene_symbol'].head(10).tolist()}", flush=True)
print(f"  top 10 down (gene_symbol): {sig_dn['gene_symbol'].head(10).tolist()}", flush=True)

sig_up[["gene_ensembl", "gene_symbol", LFC, "padj"]].to_csv(
    OUT_DIR / "richter2021_polyploid_up.tsv", sep="\t", index=False)
sig_dn[["gene_ensembl", "gene_symbol", LFC, "padj"]].to_csv(
    OUT_DIR / "richter2021_polyploid_down.tsv", sep="\t", index=False)

# Sanity overlap with mechanism markers
mm = pd.read_csv(ROOT / "data/ploidy_signatures/mechanism_markers.tsv", sep="\t")
mm_human = set(mm["gene"].str.upper())
# Mouse equivalents: capitalize first letter only
mm_mouse = {g[0].upper() + g[1:].lower() for g in mm["gene"]}
print(f"\n=== Mechanism marker overlap (human/mouse symbols) ===", flush=True)
up_syms_upper = set(sig_up["gene_symbol"].dropna().astype(str).str.upper())
dn_syms_upper = set(sig_dn["gene_symbol"].dropna().astype(str).str.upper())
overlap_up = mm_human & up_syms_upper
overlap_dn = mm_human & dn_syms_upper
print(f"  mechanism markers in UP-in-4n: {sorted(overlap_up)}", flush=True)
print(f"  mechanism markers in DOWN-in-4n: {sorted(overlap_dn)}", flush=True)
print(f"  total mechanism overlap: {len(overlap_up | overlap_dn)} of {len(mm_human)} markers", flush=True)
print(f"\n[{pd.Timestamp.now()}] DONE", flush=True)
