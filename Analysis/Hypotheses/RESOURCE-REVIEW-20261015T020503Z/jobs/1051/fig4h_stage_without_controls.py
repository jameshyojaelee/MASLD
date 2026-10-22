"""Recompute the Figure 4H x-axis (system fibrosis-stage coefficient) with and without
the 68 GSE213621 healthy controls that the stage code places at level 0.
Read-only inputs; writes only to this job directory. Model per feature and cohort:
score ~ fibrosis_stage_z + inferred_sex (OLS), fixed-effect inverse-variance meta over
GSE162694 and GSE213621, then the published collection-balanced system aggregation
(median within system x collection, mean over collections)."""
import sys, numpy as np, pandas as pd
from scipy.stats import spearmanr
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT = ROOT + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1051/'
ML = ROOT + 'RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z/'
SRC = ROOT + 'figures/main/fig4_singlecell_programs/source_tables/'
collections = sys.argv[1].split(',') if len(sys.argv) > 1 else ['hallmark', 'kegg', 'reactome', 'go_bp', 'go_mf', 'go_cc']
man = pd.read_csv(ROOT + 'figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv', sep='\t', keep_default_na=False)
man = man[man.dataset.isin(['GSE162694', 'GSE213621']) & (man.fibrosis_stage != 'NA')].copy()
man['fibrosis_stage'] = man.fibrosis_stage.astype(float)
man['male'] = (man.inferred_sex == 'M').astype(float)
print(man.groupby(['dataset', 'group_binary']).size(), flush=True)

def cohort_fit(Y, meta):
    """Y: donors x features (DataFrame aligned to meta rows). Returns beta, se for fz."""
    fz = (meta.fibrosis_stage - meta.fibrosis_stage.mean()) / meta.fibrosis_stage.std()
    X = np.column_stack([np.ones(len(meta)), fz.values, meta.male.values])
    Yv = Y.values
    ok = np.isfinite(Yv).all(axis=0)
    XtX_inv = np.linalg.inv(X.T @ X)
    B = XtX_inv @ X.T @ np.nan_to_num(Yv)
    resid = np.nan_to_num(Yv) - X @ B
    dfree = len(meta) - X.shape[1]
    s2 = (resid ** 2).sum(axis=0) / dfree
    se = np.sqrt(s2 * XtX_inv[1, 1])
    b = B[1].copy(); b[~ok] = np.nan; se[~ok] = np.nan
    return pd.Series(b, index=Y.columns), pd.Series(se, index=Y.columns)

def feature_meta(wide, roster_filter):
    res = {}
    for ds in ['GSE162694', 'GSE213621']:
        m = man[man.dataset == ds]
        if ds == 'GSE213621' and roster_filter == 'disease_only':
            m = m[m.group_binary == 'Disease']
        m = m[m.sample_id.isin(wide.index)]
        res[ds] = cohort_fit(wide.loc[m.sample_id], m)
    (b1, s1), (b2, s2) = res['GSE162694'], res['GSE213621']
    w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
    return (b1 * w1 + b2 * w2) / (w1 + w2)

rows = []
for col in collections:
    d = pd.read_csv(ML + f'pathway_tf/pathways/{col}/donor_scores.tsv.gz', sep='\t', usecols=['set_id', 'sample_id', 'pathway_score', 'dataset'])
    d = d[d.dataset.isin(['GSE162694', 'GSE213621'])]
    wide = d.pivot(index='sample_id', columns='set_id', values='pathway_score')
    del d
    for roster in ['as_published', 'disease_only']:
        b = feature_meta(wide, roster)
        rows.append(pd.DataFrame({'feature_id': b.index, 'collection': col, 'roster': roster, 'stage_beta': b.values}))
    print(col, wide.shape, flush=True)
h = pd.read_csv(ML + 'programs/hotspot/hotspot_participant_scores.tsv.gz', sep='\t', usecols=['program_uid', 'sample_id', 'dataset', 'outcome_z'])
h = h[h.dataset.isin(['GSE162694', 'GSE213621'])]
wide = h.pivot(index='sample_id', columns='program_uid', values='outcome_z')
for roster in ['as_published', 'disease_only']:
    b = feature_meta(wide, roster)
    rows.append(pd.DataFrame({'feature_id': b.index, 'collection': 'hotspot', 'roster': roster, 'stage_beta': b.values}))
fe = pd.concat(rows)
fe.to_csv(OUT + 'feature_stage_beta_with_without_controls.tsv.gz', sep='\t', index=False)

# reproduction against published feature-level stage betas
pub_p = pd.read_csv(SRC + 'corrected_774/stage_continuum_comparison/pathway_stage_vs_continuum.tsv.gz', sep='\t').rename(columns={'set_id': 'feature_id'})
pub_h = pd.read_csv(SRC + 'corrected_774/stage_continuum_comparison/hotspot_stage_vs_continuum.tsv', sep='\t').rename(columns={'program_uid': 'feature_id'})
pub = pd.concat([pub_p[['feature_id', 'beta']], pub_h[['feature_id', 'beta']]])
chk = fe[fe.roster == 'as_published'].merge(pub, on='feature_id')
diff = (chk.stage_beta - chk.beta).abs()
print('feature reproduction: n', len(chk), 'max abs diff', diff.max(), 'n diff>1e-6', int((diff > 1e-6).sum()), flush=True)

coords = pd.read_csv(SRC + 'current_candidate/molecular_systems_atlas_20260824/geometry/latent_coordinates.tsv.gz', sep='\t')
coords = coords[coords.component_role == 'main_system_component'][['feature_id', 'collection', 'community_id']]
src4h = pd.read_csv(SRC + 'corrected_774/fig4h_single_outline/fig4h_source.tsv', sep='\t').set_index('community_id')
summary = []
for roster in ['as_published', 'disease_only']:
    x = coords.merge(fe[fe.roster == roster], on=['feature_id', 'collection'], how='left')
    fam = x.dropna(subset=['stage_beta']).groupby(['community_id', 'collection']).stage_beta.median()
    sysb = fam.groupby('community_id').mean()
    sysb = sysb.loc[src4h.index]
    cb = src4h.continuum_beta
    rho = spearmanr(sysb, cb)[0]
    same = int((np.sign(sysb) == np.sign(cb)).sum())
    if roster == 'as_published':
        print('system stage_beta reproduction vs fig4h_source: max abs diff', float((sysb - src4h.stage_beta).abs().max()), flush=True)
    summary.append({'roster': roster, 'n_GSE213621': int((man.dataset == 'GSE213621').sum() if roster == 'as_published' else ((man.dataset == 'GSE213621') & (man.group_binary == 'Disease')).sum()),
                    'spearman_rho_stage_vs_continuum': rho, 'same_direction_of_43': same,
                    'median_system_stage_beta': float(sysb.median())})
    sysb.rename('stage_beta_' + roster).to_csv(OUT + f'system_stage_beta_{roster}.tsv', sep='\t')
s = pd.DataFrame(summary)
print(s.to_string(), flush=True)
s.to_csv(OUT + 'fig4h_stage_without_controls_summary.tsv', sep='\t', index=False)
