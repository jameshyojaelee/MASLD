#!/usr/bin/env python3
"""Aggregate per-cell bulk-NMF (k=6) and cNMF (k=16) program activations by
(program x cell_type x disease_stage_coarse).

Disease stage mapping (matches Script 309):
  Healthy            -> 0
  NAFLD / MASL       -> 1 (Steatosis)
  NASH  / MASH       -> 2 (Steatohepatitis)
  Cirrhotic          -> 3 (Cirrhosis)
  MASLD (unlabeled)  -> recover from GSE244832 donor_pairing if dataset matches,
                        else NaN (cell dropped).

Outputs:
  results_gpu_v2/disease_signatures/
    bulk_nmf_celltype_stage_mean.csv  (rows=cell_type x stage, cols=P1..P6)
    cnmf_celltype_stage_mean.csv      (rows=cell_type x stage, cols=P1..P16)
"""

import os
import re
import h5py
import numpy as np
import pandas as pd
import scanpy as sc

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SC_DIR = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2")
ATLAS = os.path.join(SC_DIR, "integrated_atlas.h5ad")
SCORED = os.path.join(SC_DIR, "disease_signatures/scored_atlas.h5ad")
CNMF_USAGES = os.path.join(
    SC_DIR, "mcp/cnmf_runs/global/global.usages.k_16.dt_0_03.consensus.txt"
)
NMF_ATLAS_CSV = os.path.join(BASE, "RNA-seq/results/subtypes/nmf_ksweep_program_atlas.csv")
DONOR_PAIRING = os.path.join(BASE, "data/GSE244832/metadata/donor_pairing.csv")
OUT_DIR = os.path.join(SC_DIR, "disease_signatures")
OUT_BULK = os.path.join(OUT_DIR, "bulk_nmf_celltype_stage_mean.csv")
OUT_CNMF = os.path.join(OUT_DIR, "cnmf_celltype_stage_mean.csv")

STAGE_LABELS = {0: "Healthy", 1: "Steatosis", 2: "Steatohepatitis", 3: "Cirrhosis"}


def read_obs_categorical(f, key):
    obj = f[f"obs/{key}"]
    if "codes" in obj:
        cats = obj["categories"][:]
        if cats.dtype.kind == "O":
            cats = np.array([c.decode() if isinstance(c, bytes) else str(c)
                             for c in cats])
        elif cats.dtype.kind == "S":
            cats = cats.astype(str)
        return cats[obj["codes"][:]]
    vals = obj[:]
    if vals.dtype.kind in ("S", "O"):
        return np.array([v.decode() if isinstance(v, bytes) else str(v)
                         for v in vals])
    return vals


def build_stage_vector(samples, datasets, conditions):
    """Map each cell to disease stage 0..3 (NaN if unknown)."""
    cond_to_stage = {
        "Healthy": 0, "Mixed": 0,
        "NAFLD": 1, "MASL": 1,
        "NASH": 2, "MASH": 2,
        "Cirrhotic": 3,
    }
    stage = np.array([cond_to_stage.get(c, np.nan) for c in conditions], dtype=float)

    # Recover GSE244832 stage from donor_pairing: SRR -> condition (NORMAL/MASL/MASH)
    dp = pd.read_csv(DONOR_PAIRING)
    srr_to_stage = {}
    g244_map = {"NORMAL": 0, "MASL": 1, "MASH": 2}
    for _, row in dp.iterrows():
        cond = g244_map.get(row["condition"])
        if cond is None or pd.isna(row["rna_srrs"]):
            continue
        for srr in str(row["rna_srrs"]).split(";"):
            srr = srr.strip()
            if srr:
                srr_to_stage[srr] = cond

    g244_mask = (datasets == "GSE244832") & np.isnan(stage)
    n_recovered = 0
    for i in np.where(g244_mask)[0]:
        s = srr_to_stage.get(samples[i])
        if s is not None:
            stage[i] = s
            n_recovered += 1
    print(f"  Recovered stage for {n_recovered:,} GSE244832 MASLD cells via donor_pairing")
    return stage


def aggregate(scores, cell_types, stages, programs):
    """Group-mean scores by (cell_type, stage). scores: (n_cells, n_programs)."""
    keep = ~np.isnan(stages)
    scores = scores[keep]
    cell_types = np.asarray(cell_types)[keep]
    stages = stages[keep].astype(int)

    df = pd.DataFrame(scores, columns=programs)
    df["cell_type"] = cell_types
    df["stage"] = [STAGE_LABELS[s] for s in stages]
    grp = df.groupby(["cell_type", "stage"]).mean(numeric_only=True)
    n = df.groupby(["cell_type", "stage"]).size().rename("n_cells")
    grp = grp.join(n)
    return grp.reset_index()


# ---------------------------------------------------------------------------
# 1. Cell metadata + stage mapping (from scored_atlas.h5ad obs)
# ---------------------------------------------------------------------------
print(f"Reading metadata from {SCORED}")
with h5py.File(SCORED, "r") as f:
    barcodes = read_obs_categorical(f, "_index")
    cell_types = read_obs_categorical(f, "cell_type")
    samples = read_obs_categorical(f, "sample")
    datasets = read_obs_categorical(f, "dataset")
    conditions = read_obs_categorical(f, "condition")

print(f"  {len(barcodes):,} cells")
stage_vec = build_stage_vector(samples, datasets, conditions)
n_staged = (~np.isnan(stage_vec)).sum()
print(f"  {n_staged:,} cells have a disease stage label "
      f"({100 * n_staged / len(barcodes):.1f}%)")

meta_df = pd.DataFrame({
    "barcode": barcodes,
    "cell_type": cell_types,
    "stage": stage_vec,
}).set_index("barcode")

# ---------------------------------------------------------------------------
# 2. cNMF per-cell usages -> aggregate
# ---------------------------------------------------------------------------
print(f"\nReading cNMF usages: {CNMF_USAGES}")
cnmf = pd.read_csv(CNMF_USAGES, sep="\t", index_col=0)
cnmf.columns = [f"P{c}" for c in cnmf.columns]
print(f"  shape: {cnmf.shape[0]:,} cells x {cnmf.shape[1]} programs")

common = cnmf.index.intersection(meta_df.index)
print(f"  joined: {len(common):,} of {cnmf.shape[0]:,} cells matched")
cnmf = cnmf.loc[common]
m = meta_df.loc[common]

cnmf_agg = aggregate(cnmf.values, m["cell_type"].values, m["stage"].values,
                     programs=list(cnmf.columns))
cnmf_agg.to_csv(OUT_CNMF, index=False)
print(f"  wrote {OUT_CNMF}  ({len(cnmf_agg)} rows)")

# ---------------------------------------------------------------------------
# 3. Bulk NMF k=6 program signatures -> score on integrated atlas -> aggregate
# ---------------------------------------------------------------------------
print(f"\nReading bulk NMF k=6 program atlas: {NMF_ATLAS_CSV}")
nmf_atlas = pd.read_csv(NMF_ATLAS_CSV)
k6 = nmf_atlas[nmf_atlas["k"] == 6][["program", "label", "top10_genes"]]
print(f"  k=6 programs: {len(k6)}")

program_genes = {}
for _, row in k6.iterrows():
    genes = [g.strip() for g in str(row["top10_genes"]).split(",") if g.strip()]
    program_genes[row["program"]] = genes

print(f"\nOpening atlas in backed mode: {ATLAS}")
adata = sc.read_h5ad(ATLAS, backed="r")
print(f"  shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# Map gene IDs (versioned Ensembl) to atlas var_names
var_lookup = set(adata.var_names)
gene_id_col = None
for c in ("gene_ids", "gene_id", "ensembl_id"):
    if c in adata.var.columns:
        gene_id_col = c
        break
ens_to_var = {}
if gene_id_col is not None:
    for vn, gid in zip(adata.var_names, adata.var[gene_id_col]):
        gid = str(gid)
        ens_to_var[gid] = vn
        ens_to_var[gid.split(".")[0]] = vn

def map_gene(g):
    if g in var_lookup:
        return g
    return ens_to_var.get(g) or ens_to_var.get(g.split(".")[0])

mapped_program_genes = {}
all_genes = []
for prog, genes in program_genes.items():
    mapped = [map_gene(g) for g in genes]
    mapped = [g for g in mapped if g is not None]
    mapped_program_genes[prog] = mapped
    all_genes.extend(mapped)
    print(f"  {prog}: {len(mapped)}/{len(genes)} genes mapped")

unique_genes = sorted(set(all_genes))
print(f"  {len(unique_genes)} unique program genes to load from atlas")

# Slice the relevant gene columns. AnnData backed mode supports column slicing.
gene_idx = [adata.var_names.get_loc(g) for g in unique_genes]
print("  Reading gene-subset expression matrix into memory...")
X_sub = adata.X[:, gene_idx]
if hasattr(X_sub, "toarray"):
    X_sub = X_sub.toarray()
X_sub = np.asarray(X_sub, dtype=np.float32)
print(f"  dense shape: {X_sub.shape}, mem={X_sub.nbytes / 1e9:.2f} GB")

# Library-size normalize + log1p (manual, on the gene subset means we need
# total counts per cell — read once)
print("  Reading per-cell library size...")
# Use the existing per-cell total if available; else compute from full row sum
# (which would densify each row -- skip and use a constant target for log1p
# normalization across the gene subset).
# Simpler: compute z-score per gene across cells, then mean across program's
# genes for each cell. This is robust to library-size differences and matches
# how disease_signature_scoring.py builds its consensus z-score component.
mu = X_sub.mean(axis=0)
sigma = X_sub.std(axis=0)
sigma[sigma == 0] = 1.0
Z = (X_sub - mu) / sigma  # z-score per gene
gene_to_col = {g: i for i, g in enumerate(unique_genes)}

program_scores = {}
for prog, genes in mapped_program_genes.items():
    if not genes:
        program_scores[prog] = np.full(adata.n_obs, np.nan)
        continue
    cols = [gene_to_col[g] for g in genes if g in gene_to_col]
    program_scores[prog] = Z[:, cols].mean(axis=1)

obs_names_atlas = np.array(adata.obs_names)
score_df = pd.DataFrame(program_scores, index=obs_names_atlas)
score_df = score_df.loc[meta_df.index.intersection(score_df.index)]
m2 = meta_df.loc[score_df.index]

bulk_agg = aggregate(score_df.values, m2["cell_type"].values, m2["stage"].values,
                     programs=list(score_df.columns))
bulk_agg.to_csv(OUT_BULK, index=False)
print(f"  wrote {OUT_BULK}  ({len(bulk_agg)} rows)")
