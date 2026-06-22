#!/usr/bin/env python3
"""
204_scdrs_scoring_lean.py — scDRS-like scoring, memory-lean version

Reads only GWAS gene columns from the 38GB h5ad using CSR row slicing,
avoiding full matrix load. Processes cells in chunks.

SIGNIFICANCE CAVEAT (2026-06-20, mega-review A6 pseudoreplication remediation):
  Per-cell disease-relevance significance from this script is NOT a valid
  statistical test of disease enrichment. Each donor contributes thousands of
  cells, so any per-CELL test (the z-tail `scdrs_pvalue`, `prop_sig`, the cell-
  permuted `mc_pvalue`/`fdr`, and the within-cell-type Mann-Whitney) treats
  correlated cells from the same donor as independent observations
  (pseudoreplication). With ~657K cells the Mann-Whitney p-values collapse to
  exactly 0.0 (n inflated by ~3 orders of magnitude vs the ≤18-donor design) —
  these are artifacts, not evidence. The CANONICAL, donor-correct disease-
  relevance result is the MAGMA scDRS arm
  (scripts/figures/figS_scdrs_magma_aggregate.py + figS_scdrs_magma.R), which
  aggregates scores to the donor level before testing and reaches the opposite
  (null) conclusion. The per-cell quantities below are retained for descriptive
  cell-state ranking only; do NOT cite their p-values.

SLURM: cpu partition, 4 CPUs, 48G, 48h
Environment: spatial
"""

import os, sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import h5py
from scipy import stats, sparse

# Force unbuffered output for SLURM logs
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)

BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR = os.path.join(BASE, "RNA-seq/results/gwas_rna_integration")
os.makedirs(OUTDIR, exist_ok=True)

# ===========================================================================
# 1. GWAS gene set
# ===========================================================================
print("Loading GWAS gene scores...")

coloc = pd.read_csv(os.path.join(BASE,
    "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc = coloc[coloc['gene'].notna() & (coloc['gene'] != '')]

twas = pd.read_csv(os.path.join(BASE,
    "RNA-seq/results/causal_inference/twas_multi_gwas_combined.csv"))
twas_best = twas.groupby('gene_name').agg(
    twas_p=('pvalue', 'min')).reset_index().rename(columns={'gene_name': 'gene'})

gwas_genes = coloc[['gene', 'coloc_best_pp4']].merge(
    twas_best[['gene', 'twas_p']], on='gene', how='outer')
gwas_genes['score'] = gwas_genes['coloc_best_pp4'].rank(pct=True).fillna(0)
gwas_genes = gwas_genes.sort_values('score', ascending=False)
gwas_gene_set = set(gwas_genes.head(1000)['gene'].values)
print(f"  GWAS gene set: {len(gwas_gene_set)}")

# ===========================================================================
# 2. Read atlas metadata + gene names only (no X matrix)
# ===========================================================================
print("Reading atlas metadata...")
h5path = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad")
f = h5py.File(h5path, 'r')

# Gene names
var_idx = f['var/_index'][:]
gene_names = [v.decode() if isinstance(v, bytes) else v for v in var_idx]
n_genes = len(gene_names)

# Cell metadata
ct_codes = f['obs/cell_type/codes'][:]
ct_cats = [c.decode() if isinstance(c, bytes) else c for c in f['obs/cell_type/categories'][:]]
cell_types = pd.Categorical.from_codes(ct_codes, categories=ct_cats)

try:
    cond_codes = f['obs/condition_harmonized/codes'][:]
    cond_cats = [c.decode() if isinstance(c, bytes) else c
                 for c in f['obs/condition_harmonized/categories'][:]]
    conditions = pd.Categorical.from_codes(cond_codes, categories=cond_cats)
except Exception:
    conditions = None

obs_idx = f['obs/_index'][:]
barcodes = [v.decode() if isinstance(v, bytes) else v for v in obs_idx]
n_cells = len(barcodes)

print(f"  {n_cells} cells x {n_genes} genes")

# GWAS gene indices
gwas_idx = sorted([i for i, g in enumerate(gene_names) if g in gwas_gene_set])
print(f"  GWAS genes found: {len(gwas_idx)}/{len(gwas_gene_set)}")

# Expression-matched control genes using row-density from indptr
# For CSR, indptr gives per-row nnz. We need per-column nnz.
# Instead of scanning all indices (too slow for 38GB), use a sampling approach:
# read a few chunks of data, compute column sums on a subset, use as proxy.
np.random.seed(42)
gwas_set = set(gwas_idx)
non_gwas_idx_arr = np.array([i for i in range(n_genes) if i not in gwas_set])

print("  Computing expression-matched controls via sampled column means...")
X_grp = f['X']
indptr = X_grp['indptr'][:]  # (n_cells+1,) — small, ~5MB

# Sample 100K cells for expression proxy (read their data from CSR)
sample_n = min(100000, n_cells)
sample_rows = np.sort(np.random.choice(n_cells, sample_n, replace=False))
col_sum = np.zeros(n_genes, dtype=np.float64)

for si in range(0, len(sample_rows), 10000):
    batch = sample_rows[si:si+10000]
    for ri in batch:
        ptr_s = int(indptr[ri])
        ptr_e = int(indptr[ri + 1])
        if ptr_e > ptr_s:
            cols = X_grp['indices'][ptr_s:ptr_e]
            vals = X_grp['data'][ptr_s:ptr_e]
            col_sum[cols] += vals

col_mean = col_sum / sample_n

# Match each GWAS gene to a non-GWAS gene with similar mean expression
gwas_means_est = col_mean[gwas_idx]
non_gwas_means = col_mean[non_gwas_idx_arr]

ctrl_idx = []
used = set()
for g_mean in gwas_means_est:
    diffs = np.abs(non_gwas_means - g_mean)
    for u in used:
        pos = np.searchsorted(non_gwas_idx_arr, u)
        if pos < len(diffs) and non_gwas_idx_arr[pos] == u:
            diffs[pos] = 1e12
    best = np.argmin(diffs)
    ctrl_idx.append(int(non_gwas_idx_arr[best]))
    used.add(int(non_gwas_idx_arr[best]))

ctrl_idx = sorted(ctrl_idx)
print(f"  Matched {len(ctrl_idx)} expression-matched control genes")

# Combined column indices we need
needed_cols = sorted(set(gwas_idx + list(ctrl_idx)))
col_to_local = {c: i for i, c in enumerate(needed_cols)}
gwas_local = [col_to_local[c] for c in gwas_idx]
ctrl_local = [col_to_local[c] for c in ctrl_idx]

print(f"  Will read {len(needed_cols)} gene columns total")

# ===========================================================================
# 3. Stream CSR row-by-row, extracting only needed columns
# ===========================================================================
print("Scoring cells via CSR streaming...")

X_grp = f['X']
# indptr already loaded above; data and indices read in chunks

CHUNK = 50000
n_chunks = (n_cells + CHUNK - 1) // CHUNK

gwas_scores_raw = np.zeros(n_cells, dtype=np.float32)
ctrl_scores_raw = np.zeros(n_cells, dtype=np.float32)

needed_set = set(needed_cols)

for ci in range(n_chunks):
    s = ci * CHUNK
    e = min(s + CHUNK, n_cells)

    # Read only the data/indices slice for these rows
    ptr_s = int(indptr[s])
    ptr_e = int(indptr[e])

    data_chunk = X_grp['data'][ptr_s:ptr_e]
    indices_chunk = X_grp['indices'][ptr_s:ptr_e]
    local_indptr = indptr[s:e+1] - ptr_s

    # Build small CSR for this chunk
    chunk_csr = sparse.csr_matrix(
        (data_chunk, indices_chunk, local_indptr),
        shape=(e - s, n_genes)
    )

    # Extract needed columns only
    sub = chunk_csr[:, needed_cols]  # (chunk_size x len(needed_cols))
    sub_dense = np.asarray(sub.todense())  # This is small: chunk x ~1500

    gwas_scores_raw[s:e] = sub_dense[:, gwas_local].mean(axis=1)
    ctrl_scores_raw[s:e] = sub_dense[:, ctrl_local].mean(axis=1)

    if ci % 5 == 0:
        print(f"  Chunk {ci+1}/{n_chunks} ({e}/{n_cells} cells)")

f.close()

# Z-normalize: (GWAS mean - control mean) / control std across cells
# Use a global normalization approach
ctrl_mean = ctrl_scores_raw.mean()
ctrl_std = ctrl_scores_raw.std() + 1e-10
gwas_mean_global = gwas_scores_raw.mean()

# Per-cell z-score: how much more does this cell express GWAS genes vs control?
scdrs_score = (gwas_scores_raw - ctrl_scores_raw)
score_mean = scdrs_score.mean()
score_std = scdrs_score.std() + 1e-10
scdrs_z = (scdrs_score - score_mean) / score_std

# DESCRIPTIVE per-cell z-tail (NOT a significance test): used only to rank cell
# states. It is NOT donor-corrected and must NOT be reported as a p-value
# (pseudoreplication — see module docstring caveat). Canonical disease-relevance
# significance = MAGMA donor arm (scripts/figures/figS_scdrs_magma_aggregate.py).
scdrs_pvalue = stats.norm.sf(scdrs_z)
print(f"\nMean z: {scdrs_z.mean():.3f} ± {scdrs_z.std():.3f}")
print(f"Cells in upper z-tail (<0.05, DESCRIPTIVE, not a donor-corrected p): "
      f"{(scdrs_pvalue < 0.05).sum()} ({(scdrs_pvalue < 0.05).mean()*100:.1f}%)")

# ===========================================================================
# 4. Cell-type enrichment (DESCRIPTIVE ranking — see significance caveat in docstring)
# ===========================================================================
# NOTE (mega-review A6): the Monte Carlo below permutes individual CELLS, so its
# `mc_pvalue`/`fdr` and `prop_sig` are pseudoreplicated and over-confident. Use
# these columns to RANK cell states descriptively only. The donor-corrected
# enrichment significance lives in the MAGMA scDRS arm
# (scripts/figures/figS_scdrs_magma_aggregate.py + figS_scdrs_magma.R).
print("\n=== Cell-type enrichment (descriptive ranking) ===")
ct_arr = np.array(cell_types)
n_mc = 1000
ct_results = []

for ct in ct_cats:
    mask = ct_arr == ct
    n_ct = mask.sum()
    if n_ct < 10:
        continue
    ct_mean = scdrs_z[mask].mean()
    mc = np.array([np.random.choice(scdrs_z, n_ct, replace=False).mean() for _ in range(n_mc)])
    mc_p = max((mc >= ct_mean).mean(), 1/n_mc)
    ct_results.append({
        'cell_type': ct, 'n_cells': int(n_ct),
        'mean_scdrs': float(ct_mean), 'mc_pvalue': float(mc_p),
        'prop_sig': float((scdrs_pvalue[mask] < 0.05).mean())
    })
    print(f"  {ct:30s}: z={ct_mean:+.3f}  p={mc_p:.4f} (cell-permuted, descriptive)  n={n_ct}")

ct_df = pd.DataFrame(ct_results)
ct_df['fdr'] = stats.false_discovery_control(ct_df['mc_pvalue'])
ct_df = ct_df.sort_values('mc_pvalue')

# ===========================================================================
# 5. Within-cell-type disease vs control
# ===========================================================================
within_results = []
if conditions is not None:
    print("\n=== Disease vs control within cell types (descriptive Δmean only) ===")
    cond_arr = np.array(conditions)
    disease_labels = {'MASLD', 'Disease', 'NASH', 'NAFLD', 'Cirrhosis'}
    control_labels = {'Healthy', 'Control', 'Normal'}

    # Disease vs control within each cell type — DESCRIPTIVE mean shift only.
    # The former per-CELL Mann-Whitney here is the pseudoreplication artifact
    # flagged in mega-review A6: with ~657K cells (thousands per donor) its
    # p-values collapse to exactly 0.0 and DO NOT constitute a valid disease-
    # enrichment test (cells are not independent — the design is ≤18 donors).
    # The per-cell test is RETIRED. Report only the descriptive Δmean score; the
    # CANONICAL donor-corrected significance is the MAGMA scDRS arm
    # (scripts/figures/figS_scdrs_magma_aggregate.py + figS_scdrs_magma.R), which
    # aggregates to donor level and reaches the opposite (null) conclusion.
    for ct in ct_cats:
        ct_mask = ct_arr == ct
        d_mask = ct_mask & np.isin(cond_arr, list(disease_labels))
        c_mask = ct_mask & np.isin(cond_arr, list(control_labels))
        if d_mask.sum() > 10 and c_mask.sum() > 10:
            mean_d = float(scdrs_z[d_mask].mean())
            mean_c = float(scdrs_z[c_mask].mean())
            within_results.append({
                'cell_type': ct, 'n_disease': int(d_mask.sum()), 'n_control': int(c_mask.sum()),
                'mean_disease': mean_d,
                'mean_control': mean_c,
                'delta_mean_scdrs': mean_d - mean_c
            })
            print(f"  {ct:30s}: D_mean={mean_d:+.3f} C_mean={mean_c:+.3f} "
                  f"Δmean={mean_d - mean_c:+.3f} (descriptive; p retired — see MAGMA donor arm)")

# ===========================================================================
# 6. Save
# ===========================================================================
print("\n=== Saving ===")

# Cell scores (compressed — 1.2M rows)
cell_df = pd.DataFrame({
    'cell_barcode': barcodes, 'scdrs_score': scdrs_z,
    'scdrs_pvalue': scdrs_pvalue, 'cell_type': ct_arr
})
if conditions is not None:
    cell_df['condition'] = np.array(conditions)
cell_df.to_csv(os.path.join(OUTDIR, "scdrs_cell_scores.csv.gz"), index=False, compression='gzip')
print(f"  scdrs_cell_scores.csv.gz: {len(cell_df)} cells")

ct_df.to_csv(os.path.join(OUTDIR, "scdrs_celltype_enrichment.csv"), index=False)
print(f"  scdrs_celltype_enrichment.csv: {len(ct_df)} types")

if within_results:
    pd.DataFrame(within_results).to_csv(os.path.join(OUTDIR, "scdrs_within_celltype.csv"), index=False)
    print(f"  scdrs_within_celltype.csv: {len(within_results)} types")

print("\nDone.")
