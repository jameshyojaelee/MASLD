#!/usr/bin/env python3
"""
312e_hnf4a_zonation_resolved.py
C1 HNF4A zonation re-test (ATAC improvement plan).

HYPOTHESIS: HNF4A is a periportal (Zone 1) marker. The bulk hepatocyte
Wilcoxon test (Script 312d) averages across zones and may mask a
periportal-specific Progressor-vs-Healthy difference. Re-run regulon
activity Wilcoxon WITHIN zone strata (Periportal vs Pericentral) for
all 12 GWAS-disrupted TFs.

Inputs:
  - Hepatocyte 100K subset h5ad (Hepatocytes_subset.h5ad, full gene set)
  - Subtype metadata + meta-subtype mapping (Script 308/312 outputs)
  - Disease regulons CSV (Script 312d source: disease_regulons.csv)

Method:
  1. Score periportal (PP) and pericentral (PC) signatures with
     scanpy.tl.score_genes. HNF4A is dropped from PP markers to avoid
     circularity (it is BOTH the TF we test and a PP zonation marker).
  2. Compute zone_score = PP_score - PC_score per cell.
  3. Approach A (primary): tercile split on zone_score -> PP / Mid / PC.
  4. Approach B (sensitivity): sharp boundary using per-cell quantiles
     of PP_score and PC_score (q75 / q25).
  5. For each of the 12 GWAS-disrupted TFs (same panel as 312d) and each
     zone, run sample-level pseudobulk Mann-Whitney U (Progressor vs
     Healthy) and BH-adjust over the 12 TFs PER ZONE.

Outputs (results_gpu_v2/hepatocyte_subtypes/zonation_resolved/):
  - zone_assignments.csv        (cell-level zone labels + scores)
  - periportal_tf_wilcoxon.csv  (Wilcoxon results, approach A PP cells)
  - pericentral_tf_wilcoxon.csv (Wilcoxon results, approach A PC cells)
  - mid_tf_wilcoxon.csv         (Wilcoxon results, approach A Mid cells)
  - periportal_tf_wilcoxon_sharp.csv (Approach B PP)
  - pericentral_tf_wilcoxon_sharp.csv (Approach B PC)
  - zone_summary.csv            (per-zone cell + sample counts)

Environment: spatial (scanpy, anndata, scipy, numpy, pandas)
"""

import os
import sys
import logging
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design'

# Inputs (same as 312d -- the 100K subset has full gene set and `sample` col)
PSEUDOTIME_H5AD = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/pseudotime/Hepatocytes_subset.h5ad'
METADATA_CSV = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv'
META_MAPPING_CSV = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/meta_subtype_mapping.csv'
REGULON_CSV = f'{PROJECT}/Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv'

# Output
OUT_DIR = f'{PROJECT}/Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/zonation_resolved'
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Marker panels (Halpern 2017 + Aizarani 2019, human symbols)
# HNF4A intentionally dropped from PP -- it is BOTH a zonation marker and
# one of the TFs we test for regulon-activity differences (circularity).
# ---------------------------------------------------------------------------
PP_MARKERS = [
    # Urea cycle / nitrogen metabolism
    'CPS1', 'OTC', 'ASS1', 'ASL', 'ARG1',
    # Gluconeogenesis
    'PCK1', 'G6PC', 'FBP1',
    # Secretory / acute phase
    'ALB', 'TTR', 'C3', 'HAL',
    # Sulfation / IGF axis
    'SULT1A1', 'IGFBP1',
    # Bile acid uptake
    'SLCO1B1', 'SLC10A1',
    # NOTE: HNF4A omitted (circular -- it is a TF we test)
]

PC_MARKERS = [
    # Drug metabolism (CYP450)
    'CYP2E1', 'CYP1A2', 'CYP3A4',
    # Glutamine synthesis / Wnt targets
    'GLUL', 'OAT', 'AXIN2', 'LGR5',
    # Bile acid synthesis
    'CYP7A1', 'CYP8B1',
    # Lipogenesis
    'ACLY', 'SCD',
    # Other PC
    'AKR1D1', 'ALDH1A1',
]

# 12 GWAS-disrupted TFs (same panel as 312d)
GWAS_TFS = [
    'HNF4A', 'RORA', 'THRB', 'NR1H4', 'CEBPB', 'FOXA1',
    'FOXA2', 'RXRA', 'AR', 'ESR1', 'PPARA', 'KLF15'
]

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_and_label_atlas():
    """Load 100K hepatocyte subset and transfer meta-subtype labels."""
    log.info(f'Loading hepatocyte atlas from {PSEUDOTIME_H5AD}')
    adata = sc.read_h5ad(PSEUDOTIME_H5AD)
    log.info(f'  Shape: {adata.shape[0]} cells x {adata.shape[1]} genes')
    log.info(f'  obs columns: {list(adata.obs.columns)[:12]}')

    log.info('Loading subtype metadata + meta-subtype mapping...')
    meta = pd.read_csv(METADATA_CSV, index_col=0)
    meta_map = pd.read_csv(META_MAPPING_CSV)
    log.info(f'  Subtype metadata: {meta.shape}')
    log.info(f'  Meta-subtype mapping: {meta_map.shape}')

    common = adata.obs.index.intersection(meta.index)
    log.info(f'  Barcode overlap: {len(common)} / {adata.shape[0]} subset cells')

    if len(common) < 1000:
        raise RuntimeError(
            f'Too few barcode overlap ({len(common)}); aborting.')

    adata.obs['hepatocyte_subtype'] = np.nan
    adata.obs.loc[common, 'hepatocyte_subtype'] = (
        meta.loc[common, 'hepatocyte_subtype'].values
    )
    subtype_to_meta = dict(zip(
        meta_map['subtype'].astype(int), meta_map['meta_subtype']))
    adata.obs['meta_subtype'] = adata.obs['hepatocyte_subtype'].map(
        lambda x: subtype_to_meta.get(int(x), 'Unknown')
        if pd.notna(x) else np.nan
    )
    n_lab = adata.obs['meta_subtype'].notna().sum()
    log.info(f'  Cells labelled with meta_subtype: {n_lab}')

    if 'sample' not in adata.obs.columns:
        raise RuntimeError(
            "obs column 'sample' missing -- required for pseudobulk.")

    adata = adata[adata.obs['meta_subtype'].notna()].copy()
    log.info(f'  After filtering to labelled cells: {adata.shape}')
    return adata


def score_zonation(adata):
    """Score PP + PC signatures and derive zone_score = PP - PC."""
    log.info('Scoring periportal + pericentral signatures...')
    available_pp = [g for g in PP_MARKERS if g in adata.var_names]
    available_pc = [g for g in PC_MARKERS if g in adata.var_names]
    log.info(f'  PP markers found: {len(available_pp)}/{len(PP_MARKERS)} '
             f'-- missing: {set(PP_MARKERS) - set(available_pp)}')
    log.info(f'  PC markers found: {len(available_pc)}/{len(PC_MARKERS)} '
             f'-- missing: {set(PC_MARKERS) - set(available_pc)}')

    if len(available_pp) < 5 or len(available_pc) < 5:
        raise RuntimeError(
            'Too few zonation markers in data; aborting.')

    sc.tl.score_genes(adata, available_pp, score_name='PP_score')
    sc.tl.score_genes(adata, available_pc, score_name='PC_score')
    adata.obs['zone_score'] = (
        adata.obs['PP_score'].astype(float)
        - adata.obs['PC_score'].astype(float)
    )
    log.info(f'  PP_score: mean={adata.obs["PP_score"].mean():.3f}, '
             f'std={adata.obs["PP_score"].std():.3f}')
    log.info(f'  PC_score: mean={adata.obs["PC_score"].mean():.3f}, '
             f'std={adata.obs["PC_score"].std():.3f}')
    log.info(f'  zone_score: mean={adata.obs["zone_score"].mean():.3f}, '
             f'std={adata.obs["zone_score"].std():.3f}')


def assign_zones(adata):
    """Approach A (tercile) + Approach B (sharp boundary)."""
    zs = adata.obs['zone_score'].astype(float).values
    pp = adata.obs['PP_score'].astype(float).values
    pc = adata.obs['PC_score'].astype(float).values

    # Approach A: tercile on zone_score
    q_lo, q_hi = np.quantile(zs, [1.0 / 3, 2.0 / 3])
    zone_A = np.where(zs <= q_lo, 'Pericentral',
                      np.where(zs >= q_hi, 'Periportal', 'Mid'))
    adata.obs['zone_A'] = pd.Categorical(
        zone_A, categories=['Periportal', 'Mid', 'Pericentral']
    )
    log.info(f'  Approach A tercile cutoffs: q33={q_lo:.3f}, q67={q_hi:.3f}')
    log.info(f'  Approach A counts: '
             f'{adata.obs["zone_A"].value_counts().to_dict()}')

    # Approach B: sharp boundary (PP high & PC low; PC high & PP low)
    pp_q75 = np.quantile(pp, 0.75)
    pp_q25 = np.quantile(pp, 0.25)
    pc_q75 = np.quantile(pc, 0.75)
    pc_q25 = np.quantile(pc, 0.25)
    zone_B = np.full(len(zs), 'Mid', dtype=object)
    zone_B[(pp > pp_q75) & (pc < pc_q25)] = 'Periportal'
    zone_B[(pc > pc_q75) & (pp < pp_q25)] = 'Pericentral'
    adata.obs['zone_B'] = pd.Categorical(
        zone_B, categories=['Periportal', 'Mid', 'Pericentral']
    )
    log.info(f'  Approach B PP cutoff: PP>{pp_q75:.3f} & PC<{pc_q25:.3f}')
    log.info(f'  Approach B PC cutoff: PC>{pc_q75:.3f} & PP<{pp_q25:.3f}')
    log.info(f'  Approach B counts: '
             f'{adata.obs["zone_B"].value_counts().to_dict()}')


def score_tf_regulons(adata):
    """Score each GWAS-TF regulon (target gene set) via score_genes."""
    log.info('Loading + scoring TF regulons...')
    regulons = pd.read_csv(REGULON_CSV)
    log.info(f'  Disease regulons: {len(regulons)}')

    all_genes = set(adata.var_names)
    tf_to_score_col = {}
    for _, row in regulons.iterrows():
        tf = row['tf_name']
        targets = [g for g in str(row['target_genes']).split(';')
                   if g in all_genes]
        if len(targets) < 2:
            log.warning(f'  {tf}: only {len(targets)} valid targets; skipped.')
            continue
        col = f'{tf}_activity'
        try:
            sc.tl.score_genes(adata, targets, score_name=col)
            tf_to_score_col[tf] = col
            log.info(f'  {tf}: scored ({len(targets)} targets)')
        except Exception as e:
            log.warning(f'  {tf}: score_genes failed: {e}')
    return tf_to_score_col


def wilcoxon_by_zone(adata, tf_to_score_col, zone_col, zone_label, tag):
    """Sample-level pseudobulk Mann-Whitney U for each TF, within one zone."""
    sub = adata.obs[
        (adata.obs[zone_col] == zone_label)
        & (adata.obs['meta_subtype'].isin(
            ['Disease-Progressor', 'Healthy']))
    ].copy()
    n_cells = len(sub)
    n_prog_cells = int((sub['meta_subtype'] == 'Disease-Progressor').sum())
    n_healthy_cells = int((sub['meta_subtype'] == 'Healthy').sum())

    if n_cells < 50:
        log.warning(f'  [{tag}/{zone_label}] only {n_cells} cells; '
                    f'skipping Wilcoxon.')
        return pd.DataFrame(columns=[
            'tf', 'zone', 'approach', 'n_progressor_samples',
            'n_healthy_samples', 'n_progressor_cells', 'n_healthy_cells',
            'mean_progressor', 'mean_healthy', 'delta', 'wilcoxon_stat',
            'pvalue', 'padj', 'note'
        ])

    cols = list(tf_to_score_col.values())
    pseudobulk = (
        sub[['sample', 'meta_subtype'] + cols]
        .groupby(['sample', 'meta_subtype'], observed=True)[cols]
        .mean()
        .reset_index()
        .dropna(subset=cols, how='all')
    )
    n_prog_samples = int(
        (pseudobulk['meta_subtype'] == 'Disease-Progressor').sum()
    )
    n_healthy_samples = int(
        (pseudobulk['meta_subtype'] == 'Healthy').sum()
    )
    log.info(f'  [{tag}/{zone_label}] cells: Prog={n_prog_cells}, '
             f'Healthy={n_healthy_cells} | samples: Prog={n_prog_samples}, '
             f'Healthy={n_healthy_samples}')

    rows = []
    for tf in GWAS_TFS:
        if tf not in tf_to_score_col:
            rows.append({
                'tf': tf, 'zone': zone_label, 'approach': tag,
                'n_progressor_samples': n_prog_samples,
                'n_healthy_samples': n_healthy_samples,
                'n_progressor_cells': n_prog_cells,
                'n_healthy_cells': n_healthy_cells,
                'mean_progressor': np.nan, 'mean_healthy': np.nan,
                'delta': np.nan,
                'wilcoxon_stat': np.nan, 'pvalue': np.nan, 'padj': np.nan,
                'note': 'regulon not scored (too few targets)'
            })
            continue

        col = tf_to_score_col[tf]
        prog_v = (pseudobulk.loc[
            pseudobulk['meta_subtype'] == 'Disease-Progressor', col]
            .dropna().values.astype(float))
        h_v = (pseudobulk.loc[
            pseudobulk['meta_subtype'] == 'Healthy', col]
            .dropna().values.astype(float))

        if len(prog_v) < 3 or len(h_v) < 3:
            stat, pval = np.nan, np.nan
            note = f'insufficient samples (prog={len(prog_v)}, ' \
                   f'healthy={len(h_v)})'
        else:
            try:
                stat, pval = stats.mannwhitneyu(
                    prog_v, h_v, alternative='two-sided')
                note = ''
            except Exception as e:
                stat, pval = np.nan, np.nan
                note = f'mannwhitneyu failed: {e}'

        mp = float(np.nanmean(prog_v)) if len(prog_v) else np.nan
        mh = float(np.nanmean(h_v)) if len(h_v) else np.nan
        rows.append({
            'tf': tf, 'zone': zone_label, 'approach': tag,
            'n_progressor_samples': int(len(prog_v)),
            'n_healthy_samples': int(len(h_v)),
            'n_progressor_cells': n_prog_cells,
            'n_healthy_cells': n_healthy_cells,
            'mean_progressor': mp, 'mean_healthy': mh,
            'delta': (mp - mh) if (np.isfinite(mp) and np.isfinite(mh))
                     else np.nan,
            'wilcoxon_stat': stat, 'pvalue': pval, 'padj': np.nan,
            'note': note
        })

    df = pd.DataFrame(rows)
    valid = df['pvalue'].notna()
    if valid.sum() > 0:
        _, padj, _, _ = multipletests(df.loc[valid, 'pvalue'],
                                      method='fdr_bh')
        df.loc[valid, 'padj'] = padj
    df = df.sort_values('pvalue', na_position='last')
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info('=' * 60)
    log.info('312e_hnf4a_zonation_resolved.py')
    log.info('=' * 60)

    adata = load_and_label_atlas()
    score_zonation(adata)
    assign_zones(adata)
    tf_to_score_col = score_tf_regulons(adata)
    log.info(f'  Successfully scored {len(tf_to_score_col)} TF regulons.')

    # Save cell-level zone assignments
    cols_to_save = [
        'sample', 'meta_subtype', 'hepatocyte_subtype',
        'PP_score', 'PC_score', 'zone_score', 'zone_A', 'zone_B',
    ]
    cols_to_save = [c for c in cols_to_save if c in adata.obs.columns]
    adata.obs[cols_to_save].to_csv(f'{OUT_DIR}/zone_assignments.csv')
    log.info(f'Wrote zone_assignments.csv ({adata.shape[0]} cells)')

    # Approach A (tercile)
    log.info('\n--- Approach A: tercile on zone_score ---')
    pp_A = wilcoxon_by_zone(
        adata, tf_to_score_col, 'zone_A', 'Periportal', 'tercile')
    mid_A = wilcoxon_by_zone(
        adata, tf_to_score_col, 'zone_A', 'Mid', 'tercile')
    pc_A = wilcoxon_by_zone(
        adata, tf_to_score_col, 'zone_A', 'Pericentral', 'tercile')

    pp_A.to_csv(f'{OUT_DIR}/periportal_tf_wilcoxon.csv', index=False)
    mid_A.to_csv(f'{OUT_DIR}/mid_tf_wilcoxon.csv', index=False)
    pc_A.to_csv(f'{OUT_DIR}/pericentral_tf_wilcoxon.csv', index=False)
    log.info(f'\n  Periportal (Approach A) top 5:\n'
             f'{pp_A.head(5)[["tf","mean_progressor","mean_healthy","delta","pvalue","padj"]].to_string(index=False)}')
    log.info(f'\n  Pericentral (Approach A) top 5:\n'
             f'{pc_A.head(5)[["tf","mean_progressor","mean_healthy","delta","pvalue","padj"]].to_string(index=False)}')

    # Approach B (sharp boundary)
    log.info('\n--- Approach B: sharp boundary on PP_score and PC_score ---')
    pp_B = wilcoxon_by_zone(
        adata, tf_to_score_col, 'zone_B', 'Periportal', 'sharp')
    pc_B = wilcoxon_by_zone(
        adata, tf_to_score_col, 'zone_B', 'Pericentral', 'sharp')
    pp_B.to_csv(f'{OUT_DIR}/periportal_tf_wilcoxon_sharp.csv', index=False)
    pc_B.to_csv(f'{OUT_DIR}/pericentral_tf_wilcoxon_sharp.csv', index=False)
    log.info(f'\n  Periportal (Approach B) top 5:\n'
             f'{pp_B.head(5)[["tf","mean_progressor","mean_healthy","delta","pvalue","padj"]].to_string(index=False)}')

    # Per-zone summary
    summary_rows = []
    for zone_col, tag in [('zone_A', 'tercile'), ('zone_B', 'sharp')]:
        for z in ['Periportal', 'Mid', 'Pericentral']:
            sub = adata.obs[
                (adata.obs[zone_col] == z)
                & (adata.obs['meta_subtype'].isin(
                    ['Disease-Progressor', 'Healthy']))
            ]
            ns = (sub.groupby(['sample', 'meta_subtype'],
                              observed=True).size().reset_index())
            summary_rows.append({
                'approach': tag, 'zone': z,
                'n_cells': len(sub),
                'n_prog_cells': int(
                    (sub['meta_subtype'] == 'Disease-Progressor').sum()),
                'n_healthy_cells': int(
                    (sub['meta_subtype'] == 'Healthy').sum()),
                'n_prog_samples': int(
                    (ns['meta_subtype'] == 'Disease-Progressor').sum()),
                'n_healthy_samples': int(
                    (ns['meta_subtype'] == 'Healthy').sum()),
            })
    pd.DataFrame(summary_rows).to_csv(
        f'{OUT_DIR}/zone_summary.csv', index=False)
    log.info('Wrote zone_summary.csv')

    # HNF4A spotlight summary
    log.info('\n' + '=' * 60)
    log.info('HNF4A SPOTLIGHT')
    log.info('=' * 60)
    for zone_df, tag in [
        (pp_A, 'Periportal/tercile'),
        (pc_A, 'Pericentral/tercile'),
        (mid_A, 'Mid/tercile'),
        (pp_B, 'Periportal/sharp'),
        (pc_B, 'Pericentral/sharp'),
    ]:
        row = zone_df[zone_df['tf'] == 'HNF4A']
        if len(row) == 0:
            log.info(f'  {tag}: HNF4A not present.')
            continue
        r = row.iloc[0]
        log.info(f'  {tag}: HNF4A '
                 f'delta={r["delta"]:.4f}, '
                 f'pvalue={r["pvalue"]}, padj={r["padj"]} '
                 f'(prog mean={r["mean_progressor"]:.4f}, '
                 f'healthy mean={r["mean_healthy"]:.4f})')

    log.info('\nDone.')


if __name__ == '__main__':
    main()
