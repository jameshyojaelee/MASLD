#!/usr/bin/env python3
"""
312d_crossmodal_deepening.py
Workstream 4: Fix empty GWAS-ATAC TF activity + spatial niche + plasma secretome

Part 1: Score 24 disease regulons on 100K hepatocyte subset (all 37K genes)
        and aggregate per meta-subtype. Wilcoxon tests for 12 GWAS-disrupted TFs.
Part 2: Correlate meta-subtype spatial scores with fibrosis markers.
Part 3: Cross-reference Progressor markers with Olink plasma panel.

Environment: spatial (scanpy, anndata, scipy, numpy, pandas)

# CHANGED 2026-04-22 (T0.11): mannwhitneyu now on pseudobulk (sample-level means)
# not per-cell, to avoid pseudo-replication (was n~23k cells, now n~246 samples).
"""

import os
import sys
import logging
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats, sparse

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design'

# Donor-aware spatial statistics helper (lives in the Spatial scripts dir).
# Part 2 below uses spearman_by_donor to replace the spot-level pooled Spearman
# (pseudoreplication fix F221/F222/F225).
sys.path.insert(0, f'{PROJECT}/Analysis/Spatial/scripts')
from spatial_stats import spearman_by_donor

# Inputs
PSEUDOTIME_H5AD = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/pseudotime/Hepatocytes_subset.h5ad'
METADATA_CSV = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv'
META_MAPPING_CSV = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/meta_subtype_mapping.csv'
REGULON_CSV = f'{PROJECT}/Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv'
SPATIAL_SCORES_CSV = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/crossmodal/spatial/subtype_spatial_scores.csv'
SPATIAL_H5AD = f'{PROJECT}/Analysis/Spatial/results/preprocessed/merged_spatial.h5ad'
MARKERS_CSV = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/subtype_markers.csv'
OLINK_TXT = f'{PROJECT}/Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt'

# Outputs
OUT_GWAS = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/crossmodal/gwas_atac'
OUT_SPATIAL = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/crossmodal/spatial'
OUT_PLASMA = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/crossmodal/plasma'

# GWAS-disrupted TFs (12 with motif disruption evidence)
GWAS_TFS = [
    'HNF4A', 'RORA', 'THRB', 'NR1H4', 'CEBPB', 'FOXA1',
    'FOXA2', 'RXRA', 'AR', 'ESR1', 'PPARA', 'KLF15'
]

# Fibrosis markers for spatial niche analysis
FIBROSIS_MARKERS = ['ACTA2', 'COL1A1', 'COL3A1', 'COL1A2', 'TGFB1', 'TIMP1']

# Meta-subtype Progressor subtypes
PROGRESSOR_SUBTYPES = [1, 6, 7, 14, 17, 31]

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger(__name__)

for d in [OUT_GWAS, OUT_SPATIAL, OUT_PLASMA]:
    os.makedirs(d, exist_ok=True)


# ===================================================================
# PART 1: TF Regulon Activity Scoring
# ===================================================================
def part1_tf_regulon_activity():
    log.info('='*60)
    log.info('PART 1: TF Regulon Activity Scoring')
    log.info('='*60)

    # Load 100K hepatocyte subset (all genes)
    log.info(f'Loading hepatocyte subset from {PSEUDOTIME_H5AD}')
    adata = sc.read_h5ad(PSEUDOTIME_H5AD)
    log.info(f'  Shape: {adata.shape[0]} cells x {adata.shape[1]} genes')
    log.info(f'  obs columns: {list(adata.obs.columns[:10])}...')

    # Ensure var_names are gene symbols
    log.info(f'  First 5 var_names: {list(adata.var_names[:5])}')

    # Load subtype metadata (657K cells)
    log.info('Loading subtype metadata...')
    meta = pd.read_csv(METADATA_CSV, index_col=0)
    log.info(f'  Metadata: {meta.shape[0]} cells, columns: {list(meta.columns[:8])}')

    # Load meta-subtype mapping
    meta_map = pd.read_csv(META_MAPPING_CSV)
    log.info(f'  Meta-subtype mapping: {meta_map.shape[0]} subtypes')

    # Match cell barcodes between 100K subset and full metadata
    common = adata.obs.index.intersection(meta.index)
    log.info(f'  Barcode overlap: {len(common)} / {adata.shape[0]} cells in subset')

    if len(common) < 1000:
        log.warning('Very low barcode overlap -- checking index formats')
        log.info(f'  adata index samples: {list(adata.obs.index[:3])}')
        log.info(f'  meta index samples: {list(meta.index[:3])}')

    # Transfer subtype labels
    adata.obs['hepatocyte_subtype'] = np.nan
    adata.obs.loc[common, 'hepatocyte_subtype'] = meta.loc[common, 'hepatocyte_subtype'].values

    # Build meta-subtype lookup
    subtype_to_meta = dict(zip(meta_map['subtype'].astype(int), meta_map['meta_subtype']))
    adata.obs['meta_subtype'] = adata.obs['hepatocyte_subtype'].map(
        lambda x: subtype_to_meta.get(int(x), 'Unknown') if pd.notna(x) else np.nan
    )

    n_labelled = adata.obs['meta_subtype'].notna().sum()
    log.info(f'  Cells with meta-subtype labels: {n_labelled}')

    if n_labelled < 1000:
        log.error('Too few labelled cells for meaningful analysis. Aborting Part 1.')
        return

    # Filter to labelled cells
    adata_labelled = adata[adata.obs['meta_subtype'].notna()].copy()
    log.info(f'  Working with {adata_labelled.shape[0]} labelled cells')

    # Load disease regulons
    log.info('Loading disease regulons...')
    regulons = pd.read_csv(REGULON_CSV)
    log.info(f'  {len(regulons)} regulons loaded')

    # Parse target genes and score each regulon
    all_genes = set(adata_labelled.var_names)
    tf_scores = {}

    for _, row in regulons.iterrows():
        tf = row['tf_name']
        targets_raw = str(row['target_genes']).split(';')
        targets_in_data = [g for g in targets_raw if g in all_genes]

        log.info(f'  {tf}: {len(targets_in_data)}/{len(targets_raw)} targets found in data')

        if len(targets_in_data) < 2:
            log.warning(f'    Skipping {tf} -- fewer than 2 valid target genes')
            continue

        score_name = f'{tf}_activity'
        try:
            sc.tl.score_genes(adata_labelled, targets_in_data, score_name=score_name)
            tf_scores[tf] = score_name
        except Exception as e:
            log.warning(f'    score_genes failed for {tf}: {e}')

    log.info(f'  Successfully scored {len(tf_scores)} regulons')

    # --- Aggregate per fine-grained subtype ---
    scored_tfs = list(tf_scores.keys())
    score_cols = [tf_scores[tf] for tf in scored_tfs]

    subtype_agg = adata_labelled.obs.groupby('hepatocyte_subtype')[score_cols].mean()
    subtype_agg.columns = scored_tfs
    subtype_agg.index.name = 'subtype'
    subtype_agg.to_csv(f'{OUT_GWAS}/subtype_tf_activity.csv')
    log.info(f'  Saved subtype_tf_activity.csv ({subtype_agg.shape})')

    # --- Aggregate per meta-subtype ---
    meta_agg = adata_labelled.obs.groupby('meta_subtype')[score_cols].mean()
    meta_agg.columns = scored_tfs
    meta_agg.index.name = 'meta_subtype'
    meta_agg.to_csv(f'{OUT_GWAS}/metasubtype_tf_activity.csv')
    log.info(f'  Saved metasubtype_tf_activity.csv ({meta_agg.shape})')
    log.info(f'\n  Meta-subtype TF activity:\n{meta_agg.round(4).to_string()}')

    # --- Wilcoxon tests: Disease-Progressor vs Healthy for GWAS-disrupted TFs ---
    # T0.11 FIX (2026-04-22): aggregate to sample-level pseudobulk BEFORE mannwhitneyu
    # to avoid pseudo-replication. Previously treated ~23k cells as independent obs
    # (p-values massively inflated); now uses one value per (sample, meta_subtype).
    log.info('Running Wilcoxon tests (Progressor vs Healthy) for GWAS-disrupted TFs...')
    if 'sample' not in adata_labelled.obs.columns:
        raise RuntimeError("obs column 'sample' missing -- required for pseudobulk aggregation")

    prog_mask = adata_labelled.obs['meta_subtype'] == 'Disease-Progressor'
    healthy_mask = adata_labelled.obs['meta_subtype'] == 'Healthy'
    n_prog_cells = int(prog_mask.sum())
    n_healthy_cells = int(healthy_mask.sum())
    log.info(f'  Cell-level counts -- Progressor: {n_prog_cells}, Healthy: {n_healthy_cells}')

    # Build a single long-form DataFrame with all scored TF columns, sample, meta_subtype.
    pb_cols = [tf_scores[tf] for tf in scored_tfs]
    pb_df = adata_labelled.obs[['sample', 'meta_subtype'] + pb_cols].copy()
    # Keep only Progressor / Healthy rows for the contrast.
    pb_df = pb_df[pb_df['meta_subtype'].isin(['Disease-Progressor', 'Healthy'])].copy()
    # Sample-level mean of each TF activity score per (sample, meta_subtype).
    pseudobulk = (
        pb_df.groupby(['sample', 'meta_subtype'], observed=True)[pb_cols]
        .mean()
        .reset_index()
        .dropna(subset=pb_cols, how='all')
    )
    n_prog_samples = int((pseudobulk['meta_subtype'] == 'Disease-Progressor').sum())
    n_healthy_samples = int((pseudobulk['meta_subtype'] == 'Healthy').sum())
    log.info(f'  Pseudobulk sample counts -- Progressor: {n_prog_samples}, Healthy: {n_healthy_samples}')

    def _pb_vals(tf_col, subtype):
        return (
            pseudobulk.loc[pseudobulk['meta_subtype'] == subtype, tf_col]
            .dropna()
            .values
            .astype(float)
        )

    wilcox_results = []
    for tf in GWAS_TFS:
        if tf not in tf_scores:
            wilcox_results.append({
                'tf': tf, 'gwas_disrupted': True,
                'n_progressor_samples': n_prog_samples,
                'n_healthy_samples': n_healthy_samples,
                'n_progressor_cells': n_prog_cells,
                'n_healthy_cells': n_healthy_cells,
                'mean_progressor': np.nan, 'mean_healthy': np.nan,
                'delta': np.nan,
                'wilcoxon_stat': np.nan, 'pvalue': np.nan,
                'note': 'regulon not scored (too few targets)'
            })
            continue

        col = tf_scores[tf]
        prog_vals = _pb_vals(col, 'Disease-Progressor')
        healthy_vals = _pb_vals(col, 'Healthy')

        if len(prog_vals) < 3 or len(healthy_vals) < 3:
            stat, pval = np.nan, np.nan
            note = f'insufficient samples (prog={len(prog_vals)}, healthy={len(healthy_vals)})'
        else:
            try:
                stat, pval = stats.mannwhitneyu(prog_vals, healthy_vals, alternative='two-sided')
                note = ''
            except Exception:
                stat, pval = np.nan, np.nan
                note = 'mannwhitneyu failed'

        wilcox_results.append({
            'tf': tf, 'gwas_disrupted': True,
            'n_progressor_samples': len(prog_vals),
            'n_healthy_samples': len(healthy_vals),
            'n_progressor_cells': n_prog_cells,
            'n_healthy_cells': n_healthy_cells,
            'mean_progressor': float(np.nanmean(prog_vals)) if len(prog_vals) else np.nan,
            'mean_healthy': float(np.nanmean(healthy_vals)) if len(healthy_vals) else np.nan,
            'delta': (float(np.nanmean(prog_vals)) - float(np.nanmean(healthy_vals)))
                     if (len(prog_vals) and len(healthy_vals)) else np.nan,
            'wilcoxon_stat': stat, 'pvalue': pval,
            'note': note
        })

    # Also test non-GWAS TFs for comparison
    for tf in scored_tfs:
        if tf in GWAS_TFS:
            continue
        col = tf_scores[tf]
        prog_vals = _pb_vals(col, 'Disease-Progressor')
        healthy_vals = _pb_vals(col, 'Healthy')

        if len(prog_vals) < 3 or len(healthy_vals) < 3:
            stat, pval = np.nan, np.nan
            note = f'insufficient samples (prog={len(prog_vals)}, healthy={len(healthy_vals)})'
        else:
            try:
                stat, pval = stats.mannwhitneyu(prog_vals, healthy_vals, alternative='two-sided')
                note = ''
            except Exception:
                stat, pval = np.nan, np.nan
                note = 'mannwhitneyu failed'

        wilcox_results.append({
            'tf': tf, 'gwas_disrupted': False,
            'n_progressor_samples': len(prog_vals),
            'n_healthy_samples': len(healthy_vals),
            'n_progressor_cells': n_prog_cells,
            'n_healthy_cells': n_healthy_cells,
            'mean_progressor': float(np.nanmean(prog_vals)) if len(prog_vals) else np.nan,
            'mean_healthy': float(np.nanmean(healthy_vals)) if len(healthy_vals) else np.nan,
            'delta': (float(np.nanmean(prog_vals)) - float(np.nanmean(healthy_vals)))
                     if (len(prog_vals) and len(healthy_vals)) else np.nan,
            'wilcoxon_stat': stat, 'pvalue': pval,
            'note': note
        })

    wilcox_df = pd.DataFrame(wilcox_results)
    # BH correction
    valid = wilcox_df['pvalue'].notna()
    wilcox_df['padj'] = np.nan
    if valid.sum() > 0:
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(wilcox_df.loc[valid, 'pvalue'], method='fdr_bh')
        wilcox_df.loc[valid, 'padj'] = padj

    wilcox_df = wilcox_df.sort_values('pvalue')
    wilcox_df.to_csv(f'{OUT_GWAS}/gwas_tf_wilcoxon.csv', index=False)
    log.info(f'  Saved gwas_tf_wilcoxon.csv ({len(wilcox_df)} TFs)')

    # Summary
    sig = wilcox_df[(wilcox_df['padj'] < 0.05) & (wilcox_df['gwas_disrupted'])]
    log.info(f'  GWAS-disrupted TFs significant (padj<0.05): {len(sig)}/{len(GWAS_TFS)}')
    if len(sig) > 0:
        log.info(f'    {sig[["tf","mean_progressor","mean_healthy","delta","padj"]].to_string(index=False)}')

    del adata, adata_labelled
    import gc; gc.collect()
    log.info('Part 1 complete.\n')


# ===================================================================
# PART 2: Spatial Niche Co-localization
# ===================================================================
def part2_spatial_niche():
    log.info('='*60)
    log.info('PART 2: Spatial Niche Co-localization')
    log.info('='*60)

    # Load spatial ssGSEA scores (subtype per spot)
    log.info('Loading spatial subtype scores...')
    spatial_scores = pd.read_csv(SPATIAL_SCORES_CSV)
    log.info(f'  Shape: {spatial_scores.shape}')

    # Load meta-subtype mapping to build meta-subtype spatial scores
    meta_map = pd.read_csv(META_MAPPING_CSV)
    subtype_to_meta = dict(zip(meta_map['subtype'].astype(int), meta_map['meta_subtype']))

    # Build meta-subtype scores by averaging constituent subtypes
    meta_subtypes = sorted(set(subtype_to_meta.values()))
    meta_spatial = pd.DataFrame(index=spatial_scores.index)

    for ms in meta_subtypes:
        member_ids = [k for k, v in subtype_to_meta.items() if v == ms]
        member_cols = [f'Hep_{i}' for i in member_ids if f'Hep_{i}' in spatial_scores.columns]
        if member_cols:
            meta_spatial[ms] = spatial_scores[member_cols].mean(axis=1)
            log.info(f'  {ms}: averaged {len(member_cols)} subtypes')
        else:
            log.warning(f'  {ms}: no matching columns found')

    # Load spatial h5ad for fibrosis marker expression
    log.info(f'Loading spatial h5ad from {SPATIAL_H5AD}...')
    sp_adata = sc.read_h5ad(SPATIAL_H5AD)
    log.info(f'  Shape: {sp_adata.shape}')

    # Extract fibrosis marker expression
    available_markers = [m for m in FIBROSIS_MARKERS if m in sp_adata.var_names]
    log.info(f'  Fibrosis markers found: {available_markers} ({len(available_markers)}/{len(FIBROSIS_MARKERS)})')

    if len(available_markers) == 0:
        log.error('No fibrosis markers found in spatial data. Aborting Part 2.')
        return

    # Get expression matrix for these markers.
    # NOTE: merged_spatial.h5ad .X is library-size-normalized log expression
    # (preprocessed by the spatial pipeline), so it is safe for Spearman here.
    marker_idx = [list(sp_adata.var_names).index(m) for m in available_markers]
    X = sp_adata.X
    if sparse.issparse(X):
        marker_expr = pd.DataFrame(
            X[:, marker_idx].toarray(),
            columns=available_markers,
            index=sp_adata.obs.index
        )
    else:
        marker_expr = pd.DataFrame(
            X[:, marker_idx],
            columns=available_markers,
            index=sp_adata.obs.index
        )

    # Per-spot donor labels (sample_id) for donor-blocked correlation.
    # Each donor belongs to exactly one dataset, so blocking by donor also
    # stratifies by dataset (Vu vs GSE192741) — addresses the cross-cohort
    # Simpson's-paradox confound (F222) as well as pseudoreplication (F221).
    # Defensive: fall back to individual/dataset if sample_id is absent.
    donor_col = next((c for c in ('sample_id', 'individual', 'dataset')
                      if c in sp_adata.obs.columns), None)
    if donor_col is None:
        log.error('No donor column (sample_id/individual/dataset) in spatial '
                  'h5ad; cannot run donor-blocked correlation. Aborting Part 2.')
        return
    log.info(f'  Donor-blocking correlation by obs column: {donor_col}')
    donor_labels = pd.Series(
        sp_adata.obs[donor_col].astype(str).values,
        index=sp_adata.obs.index,
        name=donor_col,
    )

    log.info(f'  Marker expression shape: {marker_expr.shape}')

    # Align spatial scores with marker expression (same spot order)
    # spatial_scores has no index column -- use positional alignment if row counts match
    n_scores = len(spatial_scores)
    n_spots = sp_adata.shape[0]
    log.info(f'  Spatial scores: {n_scores} rows, Spatial h5ad: {n_spots} spots')

    if n_scores != n_spots:
        log.warning(f'  Row count mismatch ({n_scores} vs {n_spots}). Using min rows.')
        n_use = min(n_scores, n_spots)
        meta_spatial = meta_spatial.iloc[:n_use]
        marker_expr = marker_expr.iloc[:n_use]
        donor_labels = donor_labels.iloc[:n_use]
    else:
        log.info('  Row counts match -- proceeding with positional alignment.')

    # Compute DONOR-AWARE Spearman correlations: each meta-subtype vs each
    # fibrosis marker. The previous implementation pooled all 24,058 spots as
    # independent observations (n_spots=24058 in every row), which is severe
    # pseudoreplication — the true biological n is ~15 donors, all from 2
    # datasets (F221). spearman_by_donor computes rho WITHIN each donor
    # (sample_id), then tests the per-donor rho distribution against 0 with a
    # Wilcoxon signed-rank test, so the reported p-value is donor-level, not
    # spot-level. Donor-blocking also removes the Vu-vs-GSE192741 cross-cohort
    # confound (F222), since each donor lies entirely within one dataset.
    # Column meaning changes (names kept for the in-script summary consumer):
    #   spearman_rho -> MEAN of per-donor Spearman rho (was pooled spot rho)
    #   pvalue       -> donor-level Wilcoxon p (was pooled spot Spearman p)
    #   n_donors     -> number of donors contributing a rho (replaces n_spots)
    results = []
    for ms in meta_subtypes:
        if ms not in meta_spatial.columns:
            continue
        ms_vals = meta_spatial[ms].values.astype(float)
        for marker in available_markers:
            m_vals = marker_expr[marker].values.astype(float)
            res = spearman_by_donor(ms_vals, m_vals, donor_labels.values,
                                    min_per_donor=20)
            results.append({
                'meta_subtype': ms,
                'fibrosis_marker': marker,
                'spearman_rho': res['mean_rho'],   # mean of per-donor rho
                'pvalue': res['pval'],             # donor-level Wilcoxon p
                'n_donors': res['n_donors'],
                'method': res['method'],
            })

    niche_df = pd.DataFrame(results)

    # BH correction over the donor-level p-values (F225: FDR now applied to
    # valid donor-level tests, not pseudoreplicated spot-level p-values).
    niche_df['padj'] = np.nan
    if len(niche_df) > 0 and niche_df['pvalue'].notna().any():
        from statsmodels.stats.multitest import multipletests
        mask = niche_df['pvalue'].notna()
        niche_df.loc[mask, 'padj'] = multipletests(
            niche_df.loc[mask, 'pvalue'], method='fdr_bh')[1]

    niche_df = niche_df.sort_values('pvalue')
    niche_df.to_csv(f'{OUT_SPATIAL}/niche_colocalization.csv', index=False)
    log.info(f'  Saved niche_colocalization.csv ({len(niche_df)} rows)')

    # Summary: Progressor correlations
    prog = niche_df[niche_df['meta_subtype'] == 'Disease-Progressor']
    if len(prog) > 0:
        log.info(f'\n  Progressor spatial-fibrosis correlations:')
        log.info(f'{prog[["fibrosis_marker","spearman_rho","padj"]].to_string(index=False)}')

    del sp_adata
    import gc; gc.collect()
    log.info('Part 2 complete.\n')


# ===================================================================
# PART 3: Plasma Secretome Analysis
# ===================================================================
def part3_plasma_secretome():
    log.info('='*60)
    log.info('PART 3: Plasma Secretome Analysis')
    log.info('='*60)

    # Load subtype markers
    log.info('Loading subtype markers...')
    markers = pd.read_csv(MARKERS_CSV)
    log.info(f'  Shape: {markers.shape}, columns: {list(markers.columns)}')

    # Load meta-subtype mapping
    meta_map = pd.read_csv(META_MAPPING_CSV)
    subtype_to_meta = dict(zip(meta_map['subtype'].astype(int), meta_map['meta_subtype']))

    # Get Progressor marker genes (union of top 50 from each Progressor subtype)
    progressor_markers = set()
    for st in PROGRESSOR_SUBTYPES:
        st_markers = markers[markers['subtype'] == st].head(50)
        progressor_markers.update(st_markers['names'].tolist())
        log.info(f'  Subtype {st}: {len(st_markers)} markers extracted')

    log.info(f'  Total unique Progressor markers: {len(progressor_markers)}')

    # Load Olink protein list
    log.info('Loading Olink protein panel...')
    olink = pd.read_csv(OLINK_TXT, sep='\t', nrows=0)  # just headers
    # Actually need first column values
    olink_full = pd.read_csv(OLINK_TXT, sep='\t')
    olink_proteins = set(olink_full.iloc[:, 0].dropna().astype(str).tolist())
    log.info(f'  Olink proteins: {len(olink_proteins)}')

    # Cross-reference
    secreted = progressor_markers.intersection(olink_proteins)
    log.info(f'  Progressor markers in Olink panel: {len(secreted)}/{len(progressor_markers)}')

    # Build output table with marker info
    secreted_rows = []
    for gene in sorted(secreted):
        # Get which Progressor subtypes this gene is a marker for
        gene_subtypes = markers[(markers['names'] == gene) &
                                (markers['subtype'].isin(PROGRESSOR_SUBTYPES))]
        subtypes_str = ';'.join(gene_subtypes['subtype'].astype(str).tolist())
        best_score = gene_subtypes['scores'].max() if len(gene_subtypes) > 0 else np.nan
        best_lfc = gene_subtypes['logfoldchanges'].max() if len(gene_subtypes) > 0 else np.nan

        secreted_rows.append({
            'gene': gene,
            'progressor_subtypes': subtypes_str,
            'n_progressor_subtypes': len(gene_subtypes),
            'best_marker_score': best_score,
            'best_logfc': best_lfc,
            'in_olink_panel': True
        })

    secreted_df = pd.DataFrame(secreted_rows)
    if len(secreted_df) > 0:
        secreted_df = secreted_df.sort_values('best_marker_score', ascending=False)

    secreted_df.to_csv(f'{OUT_PLASMA}/secreted_progressor_markers.csv', index=False)
    log.info(f'  Saved secreted_progressor_markers.csv ({len(secreted_df)} genes)')

    if len(secreted_df) > 0:
        log.info(f'\n  Top secreted Progressor markers:')
        log.info(f'{secreted_df.head(20).to_string(index=False)}')

    log.info('Part 3 complete.\n')


# ===================================================================
# MAIN
# ===================================================================
if __name__ == '__main__':
    log.info('312d_crossmodal_deepening.py')
    log.info(f'Project root: {PROJECT}')

    # Run all three parts
    part1_tf_regulon_activity()
    part2_spatial_niche()
    part3_plasma_secretome()

    log.info('='*60)
    log.info('All parts complete.')
    log.info(f'Outputs:')
    log.info(f'  GWAS-ATAC: {OUT_GWAS}/')
    log.info(f'  Spatial:   {OUT_SPATIAL}/')
    log.info(f'  Plasma:    {OUT_PLASMA}/')
    log.info('='*60)
