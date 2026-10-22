"""Blind re-derivation of Figure 4H (system stage vs continuum coefficients) with and
without the 68 GSE213621 controls, from release files only.

Inputs (read-only): five-cohort manifest, corrected774 continuum axes, corrected774
molecular-layer pathway donor scores and Hotspot participant scores, frozen system
coordinates. Model per feature and cohort (OLS, classical SE):
  stage:     score ~ stage_z + sex            (stage_z standardized in the fit sample)
  continuum: score ~ axis_z + factor(stage) + sex   (axis_z standardized in the fit sample)
Fixed-effect inverse-variance meta over GSE162694 and GSE213621. Systems: median within
system x collection, mean over collections (equal-collection).
"""
import sys, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr

ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT = sys.argv[1]
ML = ROOT + 'RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z/'
HAC = ROOT + 'RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261009T171017Z/'
SRC = ROOT + 'figures/main/fig4_singlecell_programs/source_tables/'
COH = ['GSE162694', 'GSE213621']

man = pd.read_csv(ROOT + 'figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv', sep='\t', keep_default_na=False)
man = man[man.dataset.isin(COH) & (man.fibrosis_stage != 'NA') & (man.inferred_sex.isin(['M', 'F']))].copy()
man['stage'] = man.fibrosis_stage.astype(float)
man['male'] = (man.inferred_sex == 'M').astype(float)
census = man.groupby(['dataset', 'group_binary', 'stage']).size().rename('n').reset_index()
census.to_csv(OUT + 'roster_census_recomputed.tsv', sep='\t', index=False)
print(census.to_string(), flush=True)

proj = pd.read_csv(HAC + 'projection/participant_scores.tsv', sep='\t')
uns = pd.read_csv(HAC + 'unsupervised/participant_scores.tsv', sep='\t')
axes = pd.concat([proj[['sample_id', 'dataset', 'axis_id', 'axis_raw']], uns[['sample_id', 'dataset', 'axis_id', 'axis_raw']]])
axes = axes[axes.axis_id.isin(['fixed_projection', 'signature_pc1'])].drop_duplicates(['sample_id', 'axis_id'])
axw = axes.pivot(index='sample_id', columns='axis_id', values='axis_raw')
man = man.merge(axw, left_on='sample_id', right_index=True, how='left')
print('axis missing:', man[['fixed_projection', 'signature_pc1']].isna().sum().to_dict(), flush=True)

ROSTERS = {
    'as_published': man,
    'disease_only': man[man.group_binary == 'Disease'],
}
for r, m in ROSTERS.items():
    print(r, m.groupby('dataset').size().to_dict(), 'total', len(m), flush=True)


def zs(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std(ddof=1)


def ols(Y, X, col):
    """Y n x p (no NaN), X n x k. returns beta, se for column col."""
    XtX_inv = np.linalg.inv(X.T @ X)
    B = XtX_inv @ X.T @ Y
    R = Y - X @ B
    df = X.shape[0] - X.shape[1]
    s2 = (R ** 2).sum(axis=0) / df
    return B[col], np.sqrt(s2 * XtX_inv[col, col])


def fit_cohort(wide, m, model, axis='fixed_projection'):
    m = m[m.sample_id.isin(wide.index)]
    Y = wide.loc[m.sample_id.values]
    if model == 'stage':
        X = np.column_stack([np.ones(len(m)), zs(m.stage), m.male.values])
    else:
        st = pd.get_dummies(m.stage.astype(int).astype(str), drop_first=True).astype(float).values
        X = np.column_stack([np.ones(len(m)), zs(m[axis]), st, m.male.values])
    Yv = Y.values.astype(float)
    good = np.isfinite(Yv).all(axis=0)
    b = np.full(Yv.shape[1], np.nan); s = np.full(Yv.shape[1], np.nan)
    if good.any():
        bb, ss = ols(Yv[:, good], X, 1)
        b[good] = bb; s[good] = ss
    # features with some missing values: fit on complete rows
    for j in np.where(~good)[0]:
        ok = np.isfinite(Yv[:, j])
        if ok.sum() >= 20:
            Xo = X[ok]
            if model == 'stage':
                Xo = np.column_stack([np.ones(ok.sum()), zs(m.stage.values[ok]), m.male.values[ok]])
            else:
                st = pd.get_dummies(m.stage.values[ok].astype(int).astype(str), drop_first=True).astype(float).values
                Xo = np.column_stack([np.ones(ok.sum()), zs(m[axis].values[ok]), st, m.male.values[ok]])
            bb, ss = ols(Yv[ok, j][:, None], Xo, 1)
            b[j] = bb[0]; s[j] = ss[0]
    return pd.Series(b, index=Y.columns), pd.Series(s, index=Y.columns), int(len(m))


def meta(wide, m, model, axis='fixed_projection'):
    bs, ss, ns = [], [], {}
    for ds in COH:
        b, s, n = fit_cohort(wide, m[m.dataset == ds], model, axis)
        bs.append(b); ss.append(s); ns[ds] = n
    w = [1 / s ** 2 for s in ss]
    return (bs[0] * w[0] + bs[1] * w[1]) / (w[0] + w[1]), ns


rows = []
nlog = []
collections = ['hallmark', 'kegg', 'reactome', 'go_bp', 'go_mf', 'go_cc']
tables = {}
for col in collections:
    d = pd.read_csv(ML + f'pathway_tf/pathways/{col}/donor_scores.tsv.gz', sep='\t', usecols=['set_id', 'sample_id', 'pathway_score', 'dataset'])
    d = d[d.dataset.isin(COH)]
    tables[col] = d.pivot(index='sample_id', columns='set_id', values='pathway_score')
    del d
    print('loaded', col, tables[col].shape, flush=True)
h = pd.read_csv(ML + 'programs/hotspot/hotspot_participant_scores.tsv.gz', sep='\t', usecols=['program_uid', 'sample_id', 'dataset', 'outcome_z'])
h = h[h.dataset.isin(COH)]
tables['hotspot'] = h.pivot(index='sample_id', columns='program_uid', values='outcome_z')
print('loaded hotspot', tables['hotspot'].shape, flush=True)

for col, wide in tables.items():
    for r, m in ROSTERS.items():
        out = {}
        for model, axis in [('stage', None), ('cont_fixed_projection', 'fixed_projection'), ('cont_signature_pc1', 'signature_pc1')]:
            b, ns = meta(wide, m, 'stage' if model == 'stage' else 'cont', axis or 'fixed_projection')
            out[model] = b
            nlog.append({'collection': col, 'roster': r, 'model': model, **ns})
        df = pd.DataFrame(out)
        df['feature_id'] = df.index; df['collection'] = col; df['roster'] = r
        rows.append(df.reset_index(drop=True))
fe = pd.concat(rows, ignore_index=True)
fe.to_csv(OUT + 'feature_stage_continuum_betas_recomputed.tsv.gz', sep='\t', index=False)
pd.DataFrame(nlog).drop_duplicates().to_csv(OUT + 'fit_sample_sizes.tsv', sep='\t', index=False)

# reproduction against published feature tables
pub_p = pd.read_csv(SRC + 'corrected_774/stage_continuum_comparison/pathway_stage_vs_continuum.tsv.gz', sep='\t').rename(columns={'set_id': 'feature_id'})
pub_h = pd.read_csv(SRC + 'corrected_774/stage_continuum_comparison/hotspot_stage_vs_continuum.tsv', sep='\t').rename(columns={'program_uid': 'feature_id'})
pub_h['collection'] = 'hotspot'
pub = pd.concat([pub_p[['feature_id', 'collection', 'beta', 'continuum_beta']], pub_h[['feature_id', 'collection', 'beta', 'continuum_beta']]])
chk = fe[fe.roster == 'as_published'].merge(pub, on=['feature_id', 'collection'])
repro = {
    'n_features_matched': len(chk),
    'stage_max_abs_diff': float((chk.stage - chk.beta).abs().max()),
    'cont_fp_max_abs_diff': float((chk.cont_fixed_projection - chk.continuum_beta).abs().max()),
}
print('feature reproduction', repro, flush=True)

# feature-layer rho per collection (features with finite published stage+continuum)
lay = []
for col in collections + ['hotspot']:
    keep = pub[(pub.collection == col) & pub.beta.notna() & pub.continuum_beta.notna()][['feature_id']]
    for r in ROSTERS:
        x = fe[(fe.roster == r) & (fe.collection == col)].merge(keep, on='feature_id')
        x = x.dropna(subset=['stage', 'cont_fixed_projection'])
        lay.append({'collection': col, 'roster': r, 'n': len(x),
                    'rho': spearmanr(x.stage, x.cont_fixed_projection)[0],
                    'same': int((np.sign(x.stage) == np.sign(x.cont_fixed_projection)).sum())})
lay = pd.DataFrame(lay)
lay.to_csv(OUT + 'feature_layer_rho_recomputed.tsv', sep='\t', index=False)
print(lay.to_string(), flush=True)

# systems
coords = pd.read_csv(SRC + 'current_candidate/molecular_systems_atlas_20260824/geometry/latent_coordinates.tsv.gz', sep='\t')
coords = coords[coords.component_role == 'main_system_component'][['feature_id', 'collection', 'community_id']]
src4h = pd.read_csv(SRC + 'corrected_774/fig4h_single_outline/fig4h_source.tsv', sep='\t').set_index('community_id')
sysres = {}
summ = []
for r in ROSTERS:
    x = coords.merge(fe[fe.roster == r], on=['feature_id', 'collection'], how='left')
    out = {}
    for v in ['stage', 'cont_fixed_projection', 'cont_signature_pc1']:
        fam = x.dropna(subset=[v]).groupby(['community_id', 'collection'])[v].median()
        out[v] = fam.groupby('community_id').mean()
    s = pd.DataFrame(out).loc[src4h.index]
    sysres[r] = s
    for yv in ['cont_fixed_projection', 'cont_signature_pc1']:
        summ.append({'roster': r, 'y': yv, 'n_systems': len(s),
                     'spearman_rho': spearmanr(s.stage, s[yv])[0],
                     'same_direction': int((np.sign(s.stage) == np.sign(s[yv])).sum()),
                     'median_stage': float(s.stage.median()), 'median_cont': float(s[yv].median())})
summ = pd.DataFrame(summ)
print(summ.to_string(), flush=True)
summ.to_csv(OUT + 'system_summary_recomputed.tsv', sep='\t', index=False)
a = sysres['as_published']
sysrep = {'stage_vs_fig4h_source_max_abs_diff': float((a.stage - src4h.stage_beta).abs().max()),
          'cont_vs_fig4h_source_max_abs_diff': float((a.cont_fixed_projection - src4h.continuum_beta).abs().max())}
print('system reproduction', sysrep, flush=True)
both = pd.concat({k: v for k, v in sysres.items()}, axis=1)
both.columns = ['_'.join(c) for c in both.columns]
both.to_csv(OUT + 'system_betas_recomputed.tsv', sep='\t')
# hotspot programs 116 (S4J)
json.dump({'feature_reproduction': repro, 'system_reproduction': sysrep}, open(OUT + 'reproduction.json', 'w'), indent=1)
