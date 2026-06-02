#!/usr/bin/env python3
"""
204_scdrs_scoring_cpu.py — scDRS-like Disease Relevance Scoring (CPU version)

Maps MASLD GWAS signal to individual cells and cell types via h5py direct read
(bypasses anndata 'ordered' attribute bug). Runs on CPU partition.

SLURM: cpu partition, 16 CPUs, 128GB, 48h
Environment: spatial
"""

import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import h5py
from scipy import stats, sparse

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

twas = pd.read_csv(os.path.join(BASE,
    "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv"))
twas_best = twas.groupby('gene_name').agg(
    twas_z=('zscore', lambda x: x.iloc[np.argmin(np.abs(x))]),
    twas_p=('pvalue', 'min')
).reset_index().rename(columns={'gene_name': 'gene'})

gwas_genes = coloc[['gene', 'coloc_best_pp4']].merge(
    twas_best[['gene', 'twas_z', 'twas_p']], on='gene', how='outer')
gwas_genes['coloc_rank'] = gwas_genes['coloc_best_pp4'].rank(pct=True)
gwas_genes['twas_rank'] = gwas_genes['twas_z'].abs().rank(pct=True)
gwas_genes['gwas_score'] = gwas_genes[['coloc_rank', 'twas_rank']].mean(axis=1)
gwas_genes = gwas_genes.sort_values('gwas_score', ascending=False).dropna(subset=['gwas_score'])

gwas_gene_set = set(gwas_genes.head(1000)['gene'].values)
print(f"  GWAS gene set: {len(gwas_gene_set)} genes")

# ===========================================================================
# 2. Load scRNA-seq atlas via h5py (bypass anndata bug)
# ===========================================================================
print("\nLoading scRNA-seq atlas via h5py...")
h5path = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad")

f = h5py.File(h5path, 'r')

# Gene names
var_idx = f['var/_index'][:]
gene_names = [v.decode() if isinstance(v, bytes) else v for v in var_idx]
n_genes = len(gene_names)

# Cell type
ct_codes = f['obs/cell_type/codes'][:]
ct_cats = f['obs/cell_type/categories'][:]
ct_cats = [c.decode() if isinstance(c, bytes) else c for c in ct_cats]
cell_types = pd.Categorical.from_codes(ct_codes, categories=ct_cats)

# Condition (for disease status)
try:
    cond_codes = f['obs/condition_harmonized/codes'][:]
    cond_cats = f['obs/condition_harmonized/categories'][:]
    cond_cats = [c.decode() if isinstance(c, bytes) else c for c in cond_cats]
    conditions = pd.Categorical.from_codes(cond_codes, categories=cond_cats)
except Exception:
    conditions = None
    print("  WARNING: condition_harmonized not readable")

# Cell barcodes
obs_idx = f['obs/_index'][:]
barcodes = [v.decode() if isinstance(v, bytes) else v for v in obs_idx]

n_cells = len(barcodes)
print(f"  Atlas: {n_cells} cells x {n_genes} genes")
print(f"  Cell types: {ct_cats}")

# Find GWAS gene indices
gwas_overlap = gwas_gene_set.intersection(set(gene_names))
gwas_gene_idx = np.array([i for i, g in enumerate(gene_names) if g in gwas_overlap])
print(f"  GWAS genes in atlas: {len(gwas_gene_idx)}/{len(gwas_gene_set)}")

# ===========================================================================
# 3. Score cells in chunks (memory-efficient)
# ===========================================================================
print("\nScoring cells...")

# Read X matrix structure
X_grp = f['X']
is_sparse = 'data' in X_grp  # CSR format

CHUNK_SIZE = 50000  # Process 50K cells at a time
n_chunks = (n_cells + CHUNK_SIZE - 1) // CHUNK_SIZE

gwas_scores = np.zeros(n_cells)
all_gene_means = np.zeros(n_genes)

# First pass: compute gene means for control matching
print("  Computing gene means...")
if is_sparse:
    data = X_grp['data'][:]
    indices = X_grp['indices'][:]
    indptr = X_grp['indptr'][:]
    X_sparse = sparse.csr_matrix((data, indices, indptr), shape=(n_cells, n_genes))

    all_gene_means = np.asarray(X_sparse.mean(axis=0)).flatten()
else:
    # Dense matrix — read in chunks
    for chunk_i in range(n_chunks):
        start = chunk_i * CHUNK_SIZE
        end = min(start + CHUNK_SIZE, n_cells)
        chunk = X_grp[start:end, :]
        all_gene_means += chunk.sum(axis=0)
    all_gene_means /= n_cells

print(f"  Gene mean range: {all_gene_means.min():.3f} - {all_gene_means.max():.3f}")

# Second pass: compute per-cell GWAS gene scores + control scores
print("  Computing per-cell scores...")
np.random.seed(42)

# Select control genes (matched by expression level)
gwas_means = all_gene_means[gwas_gene_idx]
n_ctrl = min(len(gwas_gene_idx), 200)
ctrl_idx_sets = []
for _ in range(50):  # 50 control sets
    ctrl_idx = []
    for gm in gwas_means[:n_ctrl]:
        candidates = np.where(np.abs(all_gene_means - gm) < 0.2 * abs(gm) + 0.01)[0]
        candidates = candidates[~np.isin(candidates, gwas_gene_idx)]
        if len(candidates) > 0:
            ctrl_idx.append(np.random.choice(candidates))
    ctrl_idx_sets.append(np.array(ctrl_idx))

if is_sparse:
    # CSR: efficient row slicing
    gwas_expr = np.asarray(X_sparse[:, gwas_gene_idx].mean(axis=1)).flatten()

    ctrl_scores = np.zeros((n_cells, len(ctrl_idx_sets)))
    for i, ctrl_idx in enumerate(ctrl_idx_sets):
        if len(ctrl_idx) > 0:
            ctrl_scores[:, i] = np.asarray(X_sparse[:, ctrl_idx].mean(axis=1)).flatten()

    ctrl_mean = ctrl_scores.mean(axis=1)
    ctrl_std = ctrl_scores.std(axis=1) + 1e-10

    gwas_scores = (gwas_expr - ctrl_mean) / ctrl_std
else:
    for chunk_i in range(n_chunks):
        start = chunk_i * CHUNK_SIZE
        end = min(start + CHUNK_SIZE, n_cells)
        chunk = X_grp[start:end, :]

        gwas_expr = chunk[:, gwas_gene_idx].mean(axis=1)
        ctrl_vals = np.array([chunk[:, ci].mean(axis=1) for ci in ctrl_idx_sets if len(ci) > 0])
        ctrl_mean = ctrl_vals.mean(axis=0)
        ctrl_std = ctrl_vals.std(axis=0) + 1e-10

        gwas_scores[start:end] = (gwas_expr - ctrl_mean) / ctrl_std

        if chunk_i % 5 == 0:
            print(f"    Chunk {chunk_i+1}/{n_chunks}")

f.close()

# P-values
scdrs_pvalue = stats.norm.sf(gwas_scores)
print(f"  Mean score: {gwas_scores.mean():.3f} ± {gwas_scores.std():.3f}")
print(f"  Cells with p < 0.05: {(scdrs_pvalue < 0.05).sum()} ({(scdrs_pvalue < 0.05).mean()*100:.1f}%)")

# ===========================================================================
# 4. Cell-type enrichment (Monte Carlo)
# ===========================================================================
print("\n=== Cell-type enrichment ===")

cell_types_arr = np.array(cell_types)
n_mc = 1000
ct_enrichment = []

for ct in ct_cats:
    ct_mask = cell_types_arr == ct
    n_ct = ct_mask.sum()
    if n_ct < 10:
        continue

    ct_mean = gwas_scores[ct_mask].mean()

    mc_means = np.array([
        np.random.choice(gwas_scores, size=n_ct, replace=False).mean()
        for _ in range(n_mc)
    ])

    mc_p = max((mc_means >= ct_mean).mean(), 1 / n_mc)
    mc_z = (ct_mean - mc_means.mean()) / (mc_means.std() + 1e-10)

    ct_enrichment.append({
        'cell_type': ct,
        'n_cells': int(n_ct),
        'mean_scdrs': ct_mean,
        'mc_pvalue': mc_p,
        'mc_zscore': mc_z,
        'prop_sig_005': float((scdrs_pvalue[ct_mask] < 0.05).mean())
    })

    print(f"  {ct:30s}: mean={ct_mean:+.3f}, p={mc_p:.4f}, n={n_ct}")

ct_enrich_df = pd.DataFrame(ct_enrichment)
ct_enrich_df['fdr'] = stats.false_discovery_control(ct_enrich_df['mc_pvalue'])
ct_enrich_df = ct_enrich_df.sort_values('mc_pvalue')

# ===========================================================================
# 5. Within-cell-type disease stratification
# ===========================================================================
within_ct = []
if conditions is not None:
    print("\n=== Within-cell-type disease stratification ===")
    cond_arr = np.array(conditions)

    for ct in ct_cats:
        ct_mask = cell_types_arr == ct
        disease_mask = ct_mask & np.isin(cond_arr, ['MASLD', 'Disease', 'NASH', 'NAFLD', 'Cirrhosis'])
        control_mask = ct_mask & np.isin(cond_arr, ['Healthy', 'Control', 'Normal'])

        if disease_mask.sum() > 10 and control_mask.sum() > 10:
            stat, p = stats.mannwhitneyu(
                gwas_scores[disease_mask],
                gwas_scores[control_mask],
                alternative='two-sided'
            )
            within_ct.append({
                'cell_type': ct,
                'n_disease': int(disease_mask.sum()),
                'n_control': int(control_mask.sum()),
                'mean_disease': float(gwas_scores[disease_mask].mean()),
                'mean_control': float(gwas_scores[control_mask].mean()),
                'mannwhitney_p': float(p)
            })
            print(f"  {ct:30s}: D={disease_mask.sum()}, C={control_mask.sum()}, p={p:.4f}")

# ===========================================================================
# 6. Save results
# ===========================================================================
print("\n=== Saving results ===")

# Cell-level scores (subsample to keep file manageable — 1.2M rows is large)
cell_df = pd.DataFrame({
    'cell_barcode': barcodes,
    'scdrs_score': gwas_scores,
    'scdrs_pvalue': scdrs_pvalue,
    'cell_type': cell_types_arr
})
if conditions is not None:
    cell_df['disease_status'] = np.array(conditions)

cell_df.to_csv(os.path.join(OUTDIR, "scdrs_cell_scores.csv.gz"),
               index=False, compression='gzip')
print(f"  Saved scdrs_cell_scores.csv.gz: {len(cell_df)} cells")

ct_enrich_df.to_csv(os.path.join(OUTDIR, "scdrs_celltype_enrichment.csv"), index=False)
print(f"  Saved scdrs_celltype_enrichment.csv: {len(ct_enrich_df)} cell types")

if within_ct:
    within_ct_df = pd.DataFrame(within_ct)
    within_ct_df.to_csv(os.path.join(OUTDIR, "scdrs_within_celltype.csv"), index=False)
    print(f"  Saved scdrs_within_celltype.csv: {len(within_ct_df)} rows")

print("\nDone.")
