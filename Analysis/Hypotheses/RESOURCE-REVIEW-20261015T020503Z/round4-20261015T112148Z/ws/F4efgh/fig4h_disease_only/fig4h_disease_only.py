"""Figure 4H and S4J on the stated estimand: GSE162694 disease units (109) and GSE213621
disease units (258); the 68 GSE213621 controls, which the five-cohort manifest codes as
fibrosis_stage 0, are removed from BOTH axes.

Extends RESOURCE-REVIEW jobs/1051 (stage axis only) to the continuum axis.
Stage model per feature and cohort (as 1051, reproduces published betas):
    score ~ stage_z + male                         (OLS, classical SE)
Continuum model per feature and cohort (as 20_pathway_continuum.R pathway_fit_matrix):
    outcome_z ~ factor(stage) + sex + axis_z        (outcome and axis standardized over the kept units)
Fixed-effect inverse-variance meta over the two cohorts; systems = collection-balanced
aggregation (median within system x collection, mean over collections).
Rosters: as_published (326 GSE213621 units incl. 68 controls) to reproduce the panel, and
disease_only (258). Stage axis also given with outcomes re-standardized over the kept units.
No randomness. Read-only inputs."""
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT = sys.argv[1].rstrip('/') + '/'
ML = ROOT + 'RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z/'
SRC = ROOT + 'figures/main/fig4_singlecell_programs/source_tables/'
REG = ROOT + 'RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261011T152106Z/projection/score_registry.tsv'
COLL = ['hallmark', 'kegg', 'reactome', 'go_bp', 'go_mf', 'go_cc']
COH = ['GSE162694', 'GSE213621']

man = pd.read_csv(ROOT + 'figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv', sep='\t', keep_default_na=False)
man = man[man.dataset.isin(COH) & (man.fibrosis_stage != 'NA')].copy()
man['fibrosis_stage'] = man.fibrosis_stage.astype(float)
man['male'] = (man.inferred_sex == 'M').astype(float)
reg = pd.read_csv(REG, sep='\t')[['sample_id', 'fixed_projection_raw', 'signature_pc1_raw']]
man = man.merge(reg, on='sample_id', how='left')
assert man.fixed_projection_raw.notna().all()
census = man.groupby(['dataset', 'group_binary', 'fibrosis_stage']).size()
print(census, flush=True)
census.rename('n').to_csv(OUT + 'roster_census.tsv', sep='\t')


def zs(x):
    return (x - x.mean()) / x.std(ddof=1)


def ols(X, Y):
    XtX_inv = np.linalg.inv(X.T @ X)
    B = XtX_inv @ X.T @ Y
    resid = Y - X @ B
    s2 = (resid ** 2).sum(axis=0) / (X.shape[0] - X.shape[1])
    return B, s2, XtX_inv


def cohort_fits(Y, meta, standardize_outcome=True):
    """Y: units x features (no missing). Returns dict of beta/se Series per model."""
    out = {}
    n = len(meta)
    fz = zs(meta.fibrosis_stage).values
    Xs = np.column_stack([np.ones(n), fz, meta.male.values])
    for tag, Yv in [('stage', Y.values), ('stage_restd', zs(Y).values)]:
        B, s2, XI = ols(Xs, Yv)
        out[tag] = (B[1], np.sqrt(s2 * XI[1, 1]))
    # 20_pathway_continuum.R standardizes pathway outcomes over the kept units;
    # 31_hotspot_continuum.R uses outcome_z = score_raw unchanged.
    Yz = zs(Y).values if standardize_outcome else Y.values
    stage_dum = pd.get_dummies(meta.fibrosis_stage.astype(int), prefix='s', drop_first=True).astype(float).values
    for axis in ['fixed_projection', 'signature_pc1']:
        az = zs(meta[axis + '_raw']).values
        Xc = np.column_stack([np.ones(n), stage_dum, meta.male.values, az])
        assert np.linalg.matrix_rank(Xc) == Xc.shape[1]
        B, s2, XI = ols(Xc, Yz)
        out['cont_' + axis] = (B[-1], np.sqrt(s2 * XI[-1, -1]))
    return {k: (pd.Series(b, index=Y.columns), pd.Series(s, index=Y.columns)) for k, (b, s) in out.items()}


def feature_meta(wide, roster, standardize_outcome=True):
    res = {}
    for ds in COH:
        m = man[man.dataset == ds]
        if roster == 'disease_only':
            m = m[m.group_binary == 'Disease']
        m = m[m.sample_id.isin(wide.index)].reset_index(drop=True)
        Y = wide.loc[m.sample_id]
        ok = np.isfinite(Y.values).all(axis=0) & (Y.std(ddof=1).values > 0)
        Y = Y.loc[:, ok]
        res[ds] = cohort_fits(Y, m, standardize_outcome)
    rows = {}
    for k in res[COH[0]]:
        (b1, s1), (b2, s2) = res[COH[0]][k], res[COH[1]][k]
        b1, s1 = b1.align(b2, join='inner')[0], s1.reindex(b1.index)
        b2, s2 = b2.reindex(b1.index), s2.reindex(b1.index)
        w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
        rows[k] = (b1 * w1 + b2 * w2) / (w1 + w2)
    return pd.DataFrame(rows)


blocks = []
for col in COLL:
    d = pd.read_csv(ML + f'pathway_tf/pathways/{col}/donor_scores.tsv.gz', sep='\t', usecols=['set_id', 'sample_id', 'pathway_score', 'dataset'])
    d = d[d.dataset.isin(COH)]
    wide = d.pivot(index='sample_id', columns='set_id', values='pathway_score')
    del d
    for roster in ['as_published', 'disease_only']:
        fm = feature_meta(wide, roster)
        fm['feature_id'], fm['collection'], fm['roster'] = fm.index, col, roster
        blocks.append(fm.reset_index(drop=True))
    print(col, wide.shape, flush=True)
h = pd.read_csv(ML + 'programs/hotspot/hotspot_participant_scores.tsv.gz', sep='\t', usecols=['program_uid', 'sample_id', 'dataset', 'outcome_z'])
h = h[h.dataset.isin(COH)]
wide = h.pivot(index='sample_id', columns='program_uid', values='outcome_z')
for roster in ['as_published', 'disease_only']:
    fm = feature_meta(wide, roster, standardize_outcome=False)
    fm['feature_id'], fm['collection'], fm['roster'] = fm.index, 'hotspot', roster
    blocks.append(fm.reset_index(drop=True))
fe = pd.concat(blocks, ignore_index=True)
fe.to_csv(OUT + 'feature_stage_continuum_betas_by_roster.tsv.gz', sep='\t', index=False)

# Reproduction of published feature-level values (as_published roster)
pub_p = pd.read_csv(SRC + 'corrected_774/stage_continuum_comparison/pathway_stage_vs_continuum.tsv.gz', sep='\t').rename(columns={'set_id': 'feature_id'})
pub_h = pd.read_csv(SRC + 'corrected_774/stage_continuum_comparison/hotspot_stage_vs_continuum.tsv', sep='\t').rename(columns={'program_uid': 'feature_id'})
pub = pd.concat([pub_p[['feature_id', 'beta', 'continuum_beta']], pub_h[['feature_id', 'beta', 'continuum_beta']]])
chk = fe[fe.roster == 'as_published'].merge(pub, on='feature_id')
rep = []
for k, ref in [('stage', 'beta'), ('cont_fixed_projection', 'continuum_beta'), ('cont_signature_pc1', 'continuum_beta')]:
    dif = (chk[k] - chk[ref]).abs()
    rep.append(dict(level='feature', ours=k, published=ref, n=len(chk), max_abs_diff=float(dif.max()), n_diff_gt_1e6=int((dif > 1e-6).sum())))

coords = pd.read_csv(SRC + 'current_candidate/molecular_systems_atlas_20260824/geometry/latent_coordinates.tsv.gz', sep='\t')
coords = coords[coords.component_role == 'main_system_component'][['feature_id', 'collection', 'community_id']]
src4h = pd.read_csv(SRC + 'corrected_774/fig4h_single_outline/fig4h_source.tsv', sep='\t').set_index('community_id')

sysrows, summ = [], []
for roster in ['as_published', 'disease_only']:
    x = coords.merge(fe[fe.roster == roster], on=['feature_id', 'collection'], how='left')
    S = {}
    for k in ['stage', 'stage_restd', 'cont_fixed_projection', 'cont_signature_pc1']:
        fam = x.dropna(subset=[k]).groupby(['community_id', 'collection'])[k].median()
        S[k] = fam.groupby('community_id').mean()
    S = pd.DataFrame(S).loc[src4h.index]
    S['roster'] = roster
    sysrows.append(S.reset_index())
    if roster == 'as_published':
        for k, ref in [('stage', 'stage_beta'), ('cont_fixed_projection', 'continuum_beta'), ('cont_signature_pc1', 'continuum_beta')]:
            rep.append(dict(level='system', ours=k, published=ref, n=43, max_abs_diff=float((S[k] - src4h[ref]).abs().max()), n_diff_gt_1e6=int(((S[k] - src4h[ref]).abs() > 1e-6).sum())))
    for xk in ['stage', 'stage_restd']:
        for yk in ['cont_fixed_projection', 'cont_signature_pc1']:
            summ.append(dict(level='systems_43', roster=roster, x=xk, y=yk,
                             spearman_rho=spearmanr(S[xk], S[yk])[0],
                             same_direction=int((np.sign(S[xk]) == np.sign(S[yk])).sum()), n=43,
                             median_x=float(S[xk].median()), median_y=float(S[yk].median())))
    hp = fe[(fe.roster == roster) & (fe.collection == 'hotspot')].dropna(subset=['stage', 'cont_fixed_projection'])
    for xk in ['stage', 'stage_restd']:
        for yk in ['cont_fixed_projection', 'cont_signature_pc1']:
            summ.append(dict(level='hotspot_programs', roster=roster, x=xk, y=yk,
                             spearman_rho=spearmanr(hp[xk], hp[yk])[0],
                             same_direction=int((np.sign(hp[xk]) == np.sign(hp[yk])).sum()), n=len(hp),
                             median_x=float(hp[xk].median()), median_y=float(hp[yk].median())))
pd.DataFrame(rep).to_csv(OUT + 'reproduction_check.tsv', sep='\t', index=False)
pd.concat(sysrows).to_csv(OUT + 'system_stage_continuum_by_roster.tsv', sep='\t', index=False)
s = pd.DataFrame(summ)
s.to_csv(OUT + 'summary_stage_vs_continuum_by_roster.tsv', sep='\t', index=False)
print(pd.DataFrame(rep).to_string(), flush=True)
print(s.to_string(), flush=True)
