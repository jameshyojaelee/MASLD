#!/usr/bin/env python3
"""
204_scdrs_scoring.py — scDRS Disease Relevance Scoring

Maps MASLD GWAS signal to individual cells and cell types in scRNA-seq.
Identifies which cell states carry genetic disease risk.

Method: scDRS (Zhang et al., Nature Genetics 2022)
  1. Generate gene sets from GWAS: top 1000 genes by gene-level GWAS score
  2. Score each cell for disease relevance
  3. Test cell-type enrichment via Monte Carlo
  4. Stratify within-cell-type by disease status (DESCRIPTIVE ONLY)

SIGNIFICANCE CAVEAT (2026-06-20, mega-review A6 pseudoreplication remediation):
  Per-cell disease-relevance significance from this script is NOT a valid
  statistical test of disease enrichment. Each donor contributes thousands of
  cells, so any per-CELL test (the z-test p-value `scdrs_pvalue`, `prop_sig_005`,
  and the within-cell-type Mann-Whitney) treats correlated cells from the same
  donor as independent observations (pseudoreplication). With ~657K cells the
  Mann-Whitney p-values collapse to exactly 0.0 (n inflated by ~3 orders of
  magnitude vs the ≤18-donor design) — these are artifacts, not evidence.
  The CANONICAL, donor-correct disease-relevance result is the MAGMA scDRS arm
  (scripts/figures/figS_scdrs_magma_aggregate.py + figS_scdrs_magma.R), which
  aggregates scores to the donor level before testing and reaches the opposite
  (null) conclusion. The per-cell quantities below are retained for descriptive
  cell-state ranking only; do NOT cite their p-values.

Inputs:
  - scRNA-seq: Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad
  - Gene-level GWAS scores: COLOC PP.H4 + TWAS z-scores
  - Cell-type annotations in adata.obs

Outputs:
  - RNA-seq/results/gwas_rna_integration/scdrs_cell_scores.csv
  - RNA-seq/results/gwas_rna_integration/scdrs_celltype_enrichment.csv
  - RNA-seq/results/gwas_rna_integration/scdrs_within_celltype.csv

SLURM: gpu partition, 1 GPU, 64GB, 48h
Environment: rapids_singlecell
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR = os.path.join(BASE, "RNA-seq/results/gwas_rna_integration")
os.makedirs(OUTDIR, exist_ok=True)

# ===========================================================================
# 1. Load gene-level GWAS scores
# ===========================================================================
print("Loading gene-level GWAS scores...")

coloc = pd.read_csv(os.path.join(BASE,
    "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc = coloc[coloc['gene'].notna() & (coloc['gene'] != '')]
print(f"  COLOC: {len(coloc)} genes")

twas = pd.read_csv(os.path.join(BASE,
    "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv"))
twas_best = twas.groupby('gene_name').agg(
    twas_z=('zscore', lambda x: x.iloc[np.argmin(np.abs(x))]),
    twas_p=('pvalue', 'min')
).reset_index()
twas_best.rename(columns={'gene_name': 'gene'}, inplace=True)
print(f"  TWAS: {len(twas_best)} genes")

# Merge into composite score
gwas_genes = coloc[['gene', 'coloc_best_pp4']].merge(
    twas_best[['gene', 'twas_z', 'twas_p']], on='gene', how='outer')

# Composite: rank-normalized average of COLOC + |TWAS z|
gwas_genes['coloc_rank'] = gwas_genes['coloc_best_pp4'].rank(pct=True)
gwas_genes['twas_rank'] = gwas_genes['twas_z'].abs().rank(pct=True)
gwas_genes['gwas_score'] = gwas_genes[['coloc_rank', 'twas_rank']].mean(axis=1)
gwas_genes = gwas_genes.sort_values('gwas_score', ascending=False)
print(f"  Composite scores: {len(gwas_genes)} genes")

# ===========================================================================
# 2. Create gene sets for scDRS-like scoring
# ===========================================================================
print("\nCreating GWAS gene sets...")

# Top N genes by composite GWAS score
for n_top in [200, 500, 1000]:
    top_genes = gwas_genes.dropna(subset=['gwas_score']).head(n_top)
    print(f"  Top {n_top}: mean COLOC PP.H4 = {top_genes['coloc_best_pp4'].mean():.3f}")

# Primary gene set: top 1000 by composite score
gwas_gene_set = set(gwas_genes.dropna(subset=['gwas_score']).head(1000)['gene'].values)

# Secondary: COLOC-only (PP.H4 > 0.3)
coloc_gene_set = set(coloc[coloc['coloc_best_pp4'] > 0.3]['gene'].values)

print(f"  Primary gene set (top 1000 composite): {len(gwas_gene_set)}")
print(f"  COLOC gene set (PP.H4 > 0.3): {len(coloc_gene_set)}")

# ===========================================================================
# 3. Load scRNA-seq atlas
# ===========================================================================
print("\nLoading scRNA-seq atlas...")
h5ad_path = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad")

# Handle h5py attribute errors by patching CategoricalIndex reading
import anndata
import h5py

# First try to fix the h5ad by removing problematic 'ordered' attribute
try:
    adata = sc.read_h5ad(h5ad_path)
except (KeyError, Exception) as e:
    print(f"  Warning: standard read failed ({e})")
    print(f"  Attempting manual load with obs/var patching...")

    # Read everything except obs, then read obs manually
    with h5py.File(h5ad_path, 'r') as f:
        print(f"  Keys in h5ad: {list(f.keys())}")
        if 'obs' in f:
            print(f"  Obs keys: {list(f['obs'].keys())[:10]}")

    # Try with anndata's more tolerant reader
    try:
        import anndata._io.h5ad as h5ad_io
        adata = anndata.read_h5ad(h5ad_path)
    except Exception as e2:
        print(f"  Fallback: reading X matrix and obs directly with h5py")
        with h5py.File(h5ad_path, 'r') as f:
            # Read X matrix
            from scipy.sparse import csr_matrix
            X_data = f['X/data'][:]
            X_indices = f['X/indices'][:]
            X_indptr = f['X/indptr'][:]
            shape = f['X'].attrs.get('shape', None)
            if shape is None:
                shape = (len(X_indptr) - 1, max(X_indices) + 1)
            X = csr_matrix((X_data, X_indices, X_indptr), shape=shape)

            # Read var (gene names)
            var_index = f['var/_index'][:]
            if isinstance(var_index[0], bytes):
                var_index = [v.decode() for v in var_index]
            var = pd.DataFrame(index=var_index)

            # Read obs manually, skipping problematic categorical columns
            obs_dict = {}
            for key in f['obs'].keys():
                if key == '_index':
                    obs_index = f['obs/_index'][:]
                    if isinstance(obs_index[0], bytes):
                        obs_index = [v.decode() for v in obs_index]
                    continue
                try:
                    # Try to read directly
                    data = f[f'obs/{key}']
                    if 'categories' in data:
                        # It's a categorical — read codes + categories
                        codes = data['codes'][:]
                        cats = data['categories'][:]
                        if isinstance(cats[0], bytes):
                            cats = [c.decode() for c in cats]
                        obs_dict[key] = pd.Categorical.from_codes(codes, categories=cats)
                    elif hasattr(data, 'shape'):
                        vals = data[:]
                        if isinstance(vals[0], bytes):
                            vals = [v.decode() for v in vals]
                        obs_dict[key] = vals
                except Exception:
                    pass

            obs = pd.DataFrame(obs_dict, index=obs_index)
            adata = anndata.AnnData(X=X, obs=obs, var=var)

print(f"  Atlas: {adata.n_obs} cells x {adata.n_vars} genes")

# Get gene names
gene_names = adata.var_names.tolist()
if any(g.startswith('ENSG') for g in gene_names[:5]):
    # Map ENSEMBL to symbol if needed
    if 'gene_name' in adata.var.columns:
        gene_map = dict(zip(adata.var.index, adata.var['gene_name']))
        gene_names = [gene_map.get(g, g) for g in gene_names]
    elif 'symbol' in adata.var.columns:
        gene_map = dict(zip(adata.var.index, adata.var['symbol']))
        gene_names = [gene_map.get(g, g) for g in gene_names]

# Find cell type column
ct_col = None
for col in ['cell_type', 'celltype', 'leiden', 'cluster']:
    if col in adata.obs.columns:
        ct_col = col
        break

if ct_col is None:
    print("  WARNING: No cell type column found. Available:", list(adata.obs.columns[:20]))
    ct_col = adata.obs.columns[0]

cell_types = adata.obs[ct_col].unique()
print(f"  Cell type column: '{ct_col}', {len(cell_types)} types")

# Find disease status column
disease_col = None
for col in ['condition_harmonized', 'condition', 'disease_status', 'group']:
    if col in adata.obs.columns:
        disease_col = col
        break

print(f"  Disease column: '{disease_col}'")

# ===========================================================================
# 4. scDRS-like scoring: per-cell GWAS relevance
# ===========================================================================
print("\n=== Computing per-cell disease relevance scores ===")

# Find overlap between gene sets and atlas genes
gwas_overlap = gwas_gene_set.intersection(set(gene_names))
coloc_overlap = coloc_gene_set.intersection(set(gene_names))
print(f"  GWAS gene set overlap with atlas: {len(gwas_overlap)}/{len(gwas_gene_set)}")
print(f"  COLOC gene set overlap: {len(coloc_overlap)}/{len(coloc_gene_set)}")

# Get gene indices for scoring
gwas_gene_idx = [i for i, g in enumerate(gene_names) if g in gwas_overlap]
coloc_gene_idx = [i for i, g in enumerate(gene_names) if g in coloc_overlap]

# Use the already-loaded adata (no backed mode needed)
print("  Scoring cells (this may take a few minutes)...")
adata_mem = adata  # already in memory from the load above

# Normalize if not already
if adata_mem.X.max() > 50:  # likely raw counts
    sc.pp.normalize_total(adata_mem, target_sum=1e4)
    sc.pp.log1p(adata_mem)

# Get gene names from in-memory object
gene_names_mem = adata_mem.var_names.tolist()
if any(g.startswith('ENSG') for g in gene_names_mem[:5]):
    if 'gene_name' in adata_mem.var.columns:
        gene_map_mem = dict(zip(adata_mem.var.index, adata_mem.var['gene_name']))
        gene_names_mem = [gene_map_mem.get(g, g) for g in gene_names_mem]

gwas_gene_idx = [i for i, g in enumerate(gene_names_mem) if g in gwas_overlap]
coloc_gene_idx = [i for i, g in enumerate(gene_names_mem) if g in coloc_overlap]

# Compute scDRS-like score:
# For each cell, mean expression of GWAS genes normalized by control gene sets
from scipy.sparse import issparse

X = adata_mem.X
if issparse(X):
    X_dense_gwas = np.asarray(X[:, gwas_gene_idx].todense())
else:
    X_dense_gwas = X[:, gwas_gene_idx]

# Score = mean(GWAS genes) - mean(matched control genes)
# Control: random genes matched by expression level
np.random.seed(42)
all_means = np.asarray(X.mean(axis=0)).flatten()
gwas_means = all_means[gwas_gene_idx]

n_ctrl_sets = 100
ctrl_scores = np.zeros((adata_mem.n_obs, n_ctrl_sets))

for i in range(n_ctrl_sets):
    # Match control genes by expression level (±20% of GWAS gene means)
    ctrl_idx = []
    for gm in gwas_means[:50]:  # Use first 50 for speed
        candidates = np.where(np.abs(all_means - gm) < 0.2 * gm + 0.01)[0]
        if len(candidates) > 0:
            ctrl_idx.append(np.random.choice(candidates))
    if len(ctrl_idx) < 10:
        ctrl_idx = np.random.choice(X.shape[1], size=len(gwas_gene_idx), replace=False)

    if issparse(X):
        ctrl_expr = np.asarray(X[:, ctrl_idx].todense())
    else:
        ctrl_expr = X[:, ctrl_idx]
    ctrl_scores[:, i] = ctrl_expr.mean(axis=1)

# Disease relevance score
gwas_cell_score = X_dense_gwas.mean(axis=1)
ctrl_mean = ctrl_scores.mean(axis=1)
ctrl_std = ctrl_scores.std(axis=1) + 1e-10

scdrs_score = (gwas_cell_score - ctrl_mean) / ctrl_std

# DESCRIPTIVE per-cell z (NOT a significance test): a per-cell normal-tail
# probability used only to rank cell states. It is NOT donor-corrected and must
# NOT be reported as a p-value (pseudoreplication — see module docstring caveat).
# Canonical disease-relevance significance = MAGMA donor arm
# (scripts/figures/figS_scdrs_magma_aggregate.py).
scdrs_pvalue = stats.norm.sf(scdrs_score)

print(f"  Mean scDRS score: {scdrs_score.mean():.3f} ± {scdrs_score.std():.3f}")
print(f"  Cells in upper z-tail (<0.05, DESCRIPTIVE, not a donor-corrected p): "
      f"{(scdrs_pvalue < 0.05).sum()} ({(scdrs_pvalue < 0.05).mean()*100:.1f}%)")

# ===========================================================================
# 5. Cell-type enrichment (DESCRIPTIVE ranking — see significance caveat in docstring)
# ===========================================================================
# NOTE (mega-review A6): the Monte Carlo below permutes individual CELLS, so its
# `mc_pvalue`/`fdr` and `prop_sig_005` are pseudoreplicated and over-confident.
# Use these columns to RANK cell states descriptively only. The donor-corrected
# enrichment significance lives in the MAGMA scDRS arm
# (scripts/figures/figS_scdrs_magma_aggregate.py + figS_scdrs_magma.R).
print("\n=== Cell-type enrichment (descriptive ranking) ===")

cell_scores_df = pd.DataFrame({
    'cell_barcode': adata_mem.obs.index,
    'scdrs_score': scdrs_score,
    'scdrs_pvalue': scdrs_pvalue,
    'cell_type': adata_mem.obs[ct_col].values
})

if disease_col:
    cell_scores_df['disease_status'] = adata_mem.obs[disease_col].values

# Monte Carlo enrichment per cell type
n_mc = 1000
ct_enrichment = []

for ct in cell_scores_df['cell_type'].unique():
    ct_mask = cell_scores_df['cell_type'] == ct
    n_ct = ct_mask.sum()
    if n_ct < 10:
        continue

    ct_mean = scdrs_score[ct_mask].mean()

    # Monte Carlo: randomly sample same number of cells
    mc_means = np.array([
        np.random.choice(scdrs_score, size=n_ct, replace=False).mean()
        for _ in range(n_mc)
    ])

    mc_p = (mc_means >= ct_mean).mean()
    mc_z = (ct_mean - mc_means.mean()) / (mc_means.std() + 1e-10)

    ct_enrichment.append({
        'cell_type': ct,
        'n_cells': n_ct,
        'mean_scdrs': ct_mean,
        'mc_pvalue': max(mc_p, 1 / n_mc),
        'mc_zscore': mc_z,
        'prop_sig_005': (scdrs_pvalue[ct_mask] < 0.05).mean()
    })

    print(f"  {ct:30s}: mean={ct_mean:+.3f}, mc_p={mc_p:.4f} (cell-permuted, descriptive), n={n_ct}")

ct_enrich_df = pd.DataFrame(ct_enrichment)
ct_enrich_df['fdr'] = stats.false_discovery_control(ct_enrich_df['mc_pvalue'])
ct_enrich_df = ct_enrich_df.sort_values('mc_pvalue')

# ===========================================================================
# 6. Within-cell-type stratification by disease status
# ===========================================================================
within_ct = []
if disease_col:
    print("\n=== Within-cell-type disease stratification ===")
    for ct in cell_scores_df['cell_type'].unique():
        ct_data = cell_scores_df[cell_scores_df['cell_type'] == ct]
        for status in ct_data['disease_status'].unique():
            mask = ct_data['disease_status'] == status
            if mask.sum() < 10:
                continue
            within_ct.append({
                'cell_type': ct,
                'disease_status': status,
                'n_cells': mask.sum(),
                'mean_scdrs': ct_data.loc[mask, 'scdrs_score'].mean(),
                'prop_sig_005': (ct_data.loc[mask, 'scdrs_pvalue'] < 0.05).mean()
            })

    within_ct_df = pd.DataFrame(within_ct)

    # Disease vs control within each cell type — DESCRIPTIVE mean shift only.
    # The former per-CELL Mann-Whitney here is the pseudoreplication artifact
    # flagged in mega-review A6: with ~657K cells (thousands per donor) its
    # p-values collapse to exactly 0.0 and DO NOT constitute a valid disease-
    # enrichment test (cells are not independent — the design is ≤18 donors).
    # The per-cell test is RETIRED. Report only the descriptive Δmean score; the
    # CANONICAL donor-corrected significance is the MAGMA scDRS arm
    # (scripts/figures/figS_scdrs_magma_aggregate.py + figS_scdrs_magma.R), which
    # aggregates to donor level and reaches the opposite (null) conclusion.
    for ct in cell_scores_df['cell_type'].unique():
        ct_data = cell_scores_df[cell_scores_df['cell_type'] == ct]
        disease_mask = ct_data['disease_status'].isin(['MASLD', 'Disease', 'NASH', 'NAFLD'])
        control_mask = ct_data['disease_status'].isin(['Healthy', 'Control', 'Normal'])

        if disease_mask.sum() > 10 and control_mask.sum() > 10:
            delta_mean = (ct_data.loc[disease_mask, 'scdrs_score'].mean()
                          - ct_data.loc[control_mask, 'scdrs_score'].mean())
            print(f"  {ct:30s}: disease={disease_mask.sum()}, control={control_mask.sum()}, "
                  f"Δmean_scdrs={delta_mean:+.3f} (descriptive; p retired — see MAGMA donor arm)")

# ===========================================================================
# 7. Save results
# ===========================================================================
print("\n=== Saving results ===")

cell_scores_df.to_csv(os.path.join(OUTDIR, "scdrs_cell_scores.csv"), index=False)
print(f"  Saved scdrs_cell_scores.csv: {len(cell_scores_df)} cells")

ct_enrich_df.to_csv(os.path.join(OUTDIR, "scdrs_celltype_enrichment.csv"), index=False)
print(f"  Saved scdrs_celltype_enrichment.csv: {len(ct_enrich_df)} cell types")

if within_ct:
    within_ct_df.to_csv(os.path.join(OUTDIR, "scdrs_within_celltype.csv"), index=False)
    print(f"  Saved scdrs_within_celltype.csv: {len(within_ct_df)} rows")

print("\nDone.")
