#!/usr/bin/env python
"""
Build 3 reference signatures for D2 Reversal arm sensitivity analysis.

ref_a: hepatocyte subtype-based (Disease-Progressor + Disease-Associated vs
       Healthy + Disease-Neutral), pseudobulk from scRNA atlas.
ref_b: bulk 5-cohort Dream Disease-vs-Control (reformat existing CSV).
ref_c: F0 vs F4 hepatocyte stratified (disease_stage_coarse Healthy vs Cirrhosis;
       falls back to Steatohepatitis if Cirrhosis too sparse).

Output:  Analysis/Perturbation/data/reference_signatures/
Schema:  gene, healthy_mean, diseased_mean, LFC, pval, padj
"""

import os
import sys
import numpy as np
import pandas as pd
import h5py
from scipy import sparse
from scipy.stats import ttest_ind_from_stats

PROJECT = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design'
H5AD = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad'
META_MAP_CSV = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/meta_subtype_mapping.csv'
DREAM_CSV = f'{PROJECT}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv'
GENCODE_TSV = f'{PROJECT}/data/gencode_v49_gene_metadata.tsv.gz'
OUT_DIR = f'{PROJECT}/Analysis/Perturbation/data/reference_signatures'

PSEUDOCOUNT = 1e-6  # for log2 ratio of group means in log-normalised space


def log(msg):
    print(f'[ref_curator] {msg}', flush=True)


def load_h5ad_pieces():
    """Return (X csr, gene_symbols, obs_df) without triggering uns parsing."""
    f = h5py.File(H5AD, 'r')
    log(f'opened {H5AD}')

    # X (csr_matrix)
    shape = tuple(f['X'].attrs['shape'])
    log(f'loading X csr {shape} ...')
    X = sparse.csr_matrix(
        (f['X/data'][:], f['X/indices'][:], f['X/indptr'][:]),
        shape=shape,
    )
    log(f'X loaded; nnz={X.nnz:,}')

    # genes (symbols)
    gene_bytes = f['var/gene'][:]
    gene_symbols = np.array([g.decode() if isinstance(g, bytes) else g for g in gene_bytes])
    log(f'n_genes={len(gene_symbols)}')

    # obs columns we need
    obs = {}
    for k in ['hepatocyte_subtype', 'disease_stage_coarse', 'sample']:
        g = f[f'obs/{k}']
        if isinstance(g, h5py.Group) and 'categories' in g.keys():
            cats = g['categories'][:]
            cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
            codes = g['codes'][:]
            arr = np.where(codes >= 0, cats[codes.clip(min=0)], None)
            obs[k] = arr
        else:
            arr = g[:]
            obs[k] = np.array([x.decode() if isinstance(x, bytes) else x for x in arr])
    obs_df = pd.DataFrame(obs)
    f.close()
    return X, gene_symbols, obs_df


def welch_test_sparse(X, mask_a, mask_b):
    """Per-gene Welch t-test of X[mask_a] vs X[mask_b] (sparse-aware).
    Returns mean_a, mean_b, t, p, df.
    """
    Xa = X[mask_a]
    Xb = X[mask_b]
    n_a, n_b = Xa.shape[0], Xb.shape[0]
    # Means
    m_a = np.asarray(Xa.mean(axis=0)).ravel()
    m_b = np.asarray(Xb.mean(axis=0)).ravel()
    # Variances via E[X^2] - E[X]^2 (population), then convert to sample var
    sqa = Xa.multiply(Xa).sum(axis=0)
    sqb = Xb.multiply(Xb).sum(axis=0)
    e2_a = np.asarray(sqa).ravel() / n_a
    e2_b = np.asarray(sqb).ravel() / n_b
    var_a = e2_a - m_a ** 2
    var_b = e2_b - m_b ** 2
    # sample variance (n/(n-1)) correction
    var_a = var_a * n_a / max(n_a - 1, 1)
    var_b = var_b * n_b / max(n_b - 1, 1)
    var_a = np.clip(var_a, 0, None)
    var_b = np.clip(var_b, 0, None)
    sd_a = np.sqrt(var_a)
    sd_b = np.sqrt(var_b)
    t, p = ttest_ind_from_stats(m_a, sd_a, n_a, m_b, sd_b, n_b, equal_var=False)
    return m_a, m_b, t, p, n_a, n_b


def bh_padj(p):
    p = np.asarray(p, dtype=float)
    mask = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(n) + 1)
    # enforce monotone
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    inv = np.empty_like(order)
    inv[order] = np.arange(n)
    out[mask] = adj[inv]
    return out


def build_ref_a(X, genes, obs_df, meta_map):
    """Disease-Progressor + Disease-Associated vs Healthy + Disease-Neutral."""
    # Map hep subtype int → meta_subtype
    obs_df = obs_df.copy()
    obs_df['hepatocyte_subtype_int'] = pd.to_numeric(obs_df['hepatocyte_subtype'], errors='coerce')
    sub2meta = dict(zip(meta_map['subtype'].astype(int), meta_map['meta_subtype']))
    obs_df['meta_subtype'] = obs_df['hepatocyte_subtype_int'].map(sub2meta)
    counts = obs_df['meta_subtype'].value_counts(dropna=False)
    log(f'ref_a meta_subtype counts:\n{counts.to_string()}')

    diseased_labels = ['Disease-Progressor', 'Disease-Associated']
    healthy_labels = ['Healthy', 'Disease-Neutral']
    mask_d = obs_df['meta_subtype'].isin(diseased_labels).to_numpy()
    mask_h = obs_df['meta_subtype'].isin(healthy_labels).to_numpy()
    log(f'ref_a diseased n={mask_d.sum():,}  healthy n={mask_h.sum():,}')

    m_d, m_h, t, p, n_d, n_h = welch_test_sparse(X, mask_d, mask_h)
    # LFC: log2((m_d + eps) / (m_h + eps)) — X is already log1p
    lfc = np.log2((m_d + PSEUDOCOUNT) / (m_h + PSEUDOCOUNT))
    padj = bh_padj(p)

    df = pd.DataFrame({
        'gene': genes,
        'healthy_mean': m_h,
        'diseased_mean': m_d,
        'LFC': lfc,
        'pval': p,
        'padj': padj,
    })
    df.attrs['n_diseased'] = int(n_d)
    df.attrs['n_healthy'] = int(n_h)
    df.attrs['diseased_labels'] = ','.join(diseased_labels)
    df.attrs['healthy_labels'] = ','.join(healthy_labels)
    return df


def build_ref_c(X, genes, obs_df):
    """F0 vs F4 stratified: Healthy vs Cirrhosis (fallback Steatohepatitis)."""
    obs_df = obs_df.copy()
    counts = obs_df['disease_stage_coarse'].value_counts(dropna=False)
    log(f'ref_c disease_stage_coarse counts:\n{counts.to_string()}')

    healthy_label = 'Healthy'
    cirr_n = (obs_df['disease_stage_coarse'] == 'Cirrhosis').sum()
    if cirr_n >= 100:
        diseased_label = 'Cirrhosis'
    else:
        diseased_label = 'Steatohepatitis'
        log(f'Cirrhosis sparse (n={cirr_n}); falling back to Steatohepatitis')

    mask_d = (obs_df['disease_stage_coarse'] == diseased_label).to_numpy()
    mask_h = (obs_df['disease_stage_coarse'] == healthy_label).to_numpy()
    log(f'ref_c {diseased_label} n={mask_d.sum():,}  {healthy_label} n={mask_h.sum():,}')

    m_d, m_h, t, p, n_d, n_h = welch_test_sparse(X, mask_d, mask_h)
    lfc = np.log2((m_d + PSEUDOCOUNT) / (m_h + PSEUDOCOUNT))
    padj = bh_padj(p)

    df = pd.DataFrame({
        'gene': genes,
        'healthy_mean': m_h,
        'diseased_mean': m_d,
        'LFC': lfc,
        'pval': p,
        'padj': padj,
    })
    df.attrs['n_diseased'] = int(n_d)
    df.attrs['n_healthy'] = int(n_h)
    df.attrs['diseased_label'] = diseased_label
    df.attrs['healthy_label'] = healthy_label
    return df


def build_ref_b():
    """Reformat dream Disease-vs-Control (5-cohort canonical)."""
    log(f'reading {DREAM_CSV}')
    d = pd.read_csv(DREAM_CSV)
    log(f'dream rows={len(d):,}')

    # gene column in dream is ENSG.version; add symbol via gencode
    gen = pd.read_csv(GENCODE_TSV, sep='\t')
    gen_map = dict(zip(gen['gene_id'], gen['gene_name']))
    d['gene_id'] = d['gene']
    d['gene'] = d['gene_id'].map(gen_map).fillna(d['gene_id'])

    out = pd.DataFrame({
        'gene': d['gene'],
        'gene_id': d['gene_id'],
        'healthy_mean': np.nan,
        'diseased_mean': np.nan,
        'LFC': d['logFC'],
        'pval': d['P.Value'],
        'padj': d['padj'],
        'AveExpr': d['AveExpr'],
    })
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    meta_map = pd.read_csv(META_MAP_CSV)

    log('=== loading h5ad ===')
    X, genes, obs_df = load_h5ad_pieces()

    log('=== ref_a (subtype-based) ===')
    ref_a = build_ref_a(X, genes, obs_df, meta_map)
    ref_a_path = os.path.join(OUT_DIR, 'ref_a_subtype.csv')
    ref_a.to_csv(ref_a_path, index=False)
    log(f'wrote {ref_a_path}  n_genes={len(ref_a)}  '
        f'n_diseased={ref_a.attrs["n_diseased"]}  n_healthy={ref_a.attrs["n_healthy"]}')

    log('=== ref_c (F0 vs F4) ===')
    ref_c = build_ref_c(X, genes, obs_df)
    ref_c_path = os.path.join(OUT_DIR, 'ref_c_F0vF4.csv')
    ref_c.to_csv(ref_c_path, index=False)
    log(f'wrote {ref_c_path}  n_genes={len(ref_c)}  '
        f'n_diseased={ref_c.attrs["n_diseased"]}  n_healthy={ref_c.attrs["n_healthy"]} '
        f'(labels d={ref_c.attrs["diseased_label"]}, h={ref_c.attrs["healthy_label"]})')

    log('=== ref_b (bulk 5-cohort) ===')
    ref_b = build_ref_b()
    ref_b_path = os.path.join(OUT_DIR, 'ref_b_bulk5cohort.csv')
    ref_b.to_csv(ref_b_path, index=False)
    log(f'wrote {ref_b_path}  n_genes={len(ref_b)}')

    # Provenance / summary for README
    summary_path = os.path.join(OUT_DIR, '_build_summary.tsv')
    rows = [
        ('ref_a_subtype', len(ref_a), ref_a.attrs['n_diseased'], ref_a.attrs['n_healthy'],
         f'Diseased={{{ref_a.attrs["diseased_labels"]}}} Healthy={{{ref_a.attrs["healthy_labels"]}}}'),
        ('ref_b_bulk5cohort', len(ref_b), np.nan, np.nan, 'Dream 5-cohort Disease-vs-Control (847 samples)'),
        ('ref_c_F0vF4', len(ref_c), ref_c.attrs['n_diseased'], ref_c.attrs['n_healthy'],
         f'Diseased={ref_c.attrs["diseased_label"]} Healthy={ref_c.attrs["healthy_label"]}'),
    ]
    pd.DataFrame(rows, columns=['ref', 'n_genes', 'n_diseased', 'n_healthy', 'definition']).to_csv(
        summary_path, sep='\t', index=False)
    log(f'wrote {summary_path}')


if __name__ == '__main__':
    main()
