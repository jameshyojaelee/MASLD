"""Blind re-derivation of composition-adjusted continuum coefficients (S4 focal forest,
R05-04/05/06, R13-05). Written from the prespecification text, not from the producing script.

Unit: genetic-identity cluster; cohorts GSE162694 and GSE213621; disease-only; recorded
stage and inferred sex. Per cohort and feature: OLS of outcome_z on axis_z (cohort z of
fixed_projection_raw) + categorical stage + sex [+ ILR composition]; HC3 SE; t interval.
Meta: fixed-effect inverse variance with HC3 SE; normal p. BH within 43 systems and 117
programs (family size fixed). Arms A, B, C, D1, D2, D3 x fraction sets earlier/rebuilt, plus U.
"""
import sys
import numpy as np, pandas as pd
from scipy import stats
from scipy.linalg import helmert

np.random.seed(20260923)
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT = sys.argv[1]
COH = ['GSE162694', 'GSE213621']
man = pd.read_csv(ROOT + 'figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv', sep='\t', keep_default_na=False)
man = man[man.dataset.isin(COH) & (man.group_binary == 'Disease') & (man.fibrosis_stage != 'NA') & man.inferred_sex.isin(['M', 'F'])].copy()
man['stage'] = man.fibrosis_stage.astype(int)
ax = pd.read_csv(ROOT + 'RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261011T152106Z/projection/score_registry.tsv', sep='\t')
man = man.merge(ax[['sample_id', 'fixed_projection_raw']], on='sample_id', how='inner')
print(man.groupby('dataset').size().to_dict(), flush=True)

hs = pd.read_csv(ROOT + 'RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261011T152106Z/programs/hotspot/hotspot_participant_scores.tsv.gz', sep='\t', usecols=['program_uid', 'sample_id', 'dataset', 'outcome_z'])
hs = hs[hs.dataset.isin(COH)].pivot(index='sample_id', columns='program_uid', values='outcome_z')
sy = pd.read_csv(ROOT + 'RNA-seq/results/histology_anchored_continuum/molecular_systems_inference/hac-system-inference-corrected774-20261009T171017Z/tables/participant_system_scores.tsv.gz', sep='\t')
sy = sy[sy.dataset.isin(COH) & (sy.aggregation_method == 'collection_balanced')].pivot(index='sample_id', columns='community_id', values='system_score_z')
FAM = {'programs': (hs, 117), 'systems': (sy, 43)}
print('programs', hs.shape, 'systems', sy.shape, flush=True)

FR = {}
for fs, base in [('earlier', 'Analysis/Deconvolution/results/'), ('rebuilt', 'Analysis/Deconvolution/corrected_774_20261009T141331Z/results/')]:
    parts = []
    for ds in COH:
        f = pd.read_csv(ROOT + base + f'{ds}/{ds}_bayesprism_proportions.tsv', sep='\t', index_col=0)
        parts.append(f)
    FR[fs] = pd.concat(parts)
TYPES = list(FR['rebuilt'].columns)
assert len(TYPES) == 16 and list(FR['earlier'].columns) == TYPES
man = man[man.sample_id.isin(FR['earlier'].index) & man.sample_id.isin(FR['rebuilt'].index)]
print('with fractions', man.groupby('dataset').size().to_dict(), flush=True)


def fail_rate(frac):
    out = {}
    for fs in ['earlier', 'rebuilt']:
        for ds in COH:
            ids = man[man.dataset == ds].sample_id
            out[(fs, ds)] = (frac[fs].loc[ids] < 1e-6).mean()
    return pd.DataFrame(out)


def closure(df):
    return df.div(df.sum(axis=1), axis=0)


# arm B pooled list
rates16 = fail_rate(FR)
failing16 = [t for t in TYPES if (rates16.loc[t] > 0.20).any()]
print('failing 16-type parts:', failing16, flush=True)

LINEAGE = {
    'hepatocytes': ['Hepatocytes'], 'cholangiocytes': ['Cholangiocytes'], 'endothelial': ['Endothelial cells'],
    'fibroblasts': ['Fibroblasts'], 'macrophage/monocyte': ['Macrophages', 'Mono+mono derived cells'],
    'other myeloid': ['cDC1s', 'cDC2s', 'pDCs', 'Neutrophils', 'Basophils'],
    'T/NK': ['T cells', 'Circulating NK/NKT', 'Resident NK'], 'B/plasma': ['B cells', 'Plasma cells'],
}
NEAREST = {'cholangiocytes': 'hepatocytes', 'other myeloid': 'macrophage/monocyte', 'macrophage/monocyte': 'other myeloid',
           'B/plasma': 'T/NK', 'T/NK': 'B/plasma', 'fibroblasts': 'endothelial', 'endothelial': 'fibroblasts',
           'hepatocytes': 'cholangiocytes'}


def group_sum(frac_df, groups):
    return pd.DataFrame({g: frac_df[m].sum(axis=1) for g, m in groups.items()})


groups_C = {k: list(v) for k, v in LINEAGE.items()}
while True:
    gf = {fs: group_sum(FR[fs], groups_C) for fs in FR}
    r = fail_rate(gf)
    mx = r.max(axis=1)
    bad = mx[mx > 0.20]
    if bad.empty:
        break
    worst = bad.idxmax()
    tgt = NEAREST[worst]
    # nearest may itself have been merged: find the group containing tgt's members
    tgt_key = next(k for k in groups_C if set(LINEAGE.get(tgt, [])) <= set(groups_C[k])) if tgt not in groups_C else tgt
    new_key = f'{tgt_key} + {worst}'
    groups_C[new_key] = groups_C.pop(tgt_key) + groups_C.pop(worst)
    print('arm C merge', worst, '->', tgt_key, flush=True)
print('arm C parts', groups_C, flush=True)


def comp_matrix(arm, fs, ids):
    f = FR[fs].loc[ids, TYPES].astype(float)
    if arm == 'A':
        x = np.maximum(f.values, 1e-6)
    elif arm == 'D1':
        x = np.maximum(f.values, 1e-4)
    elif arm == 'D2':
        x = np.maximum(f.values, 1e-3)
    elif arm == 'D3':
        v = f.values.copy(); d = 1e-4
        for i in range(v.shape[0]):
            low = v[i] < d; k = low.sum()
            if k:
                hi = ~low
                v[i, hi] = v[i, hi] * (1 - k * d) / v[i, hi].sum()
                v[i, low] = d
        x = v
    elif arm == 'B':
        keep = [t for t in TYPES if t not in failing16]
        g = f[keep].copy(); g['other'] = f[failing16].sum(axis=1)
        x = np.maximum(closure(g).values, 1e-6)
    elif arm == 'C':
        g = group_sum(f, groups_C)
        x = np.maximum(closure(g).values, 1e-6)
    D = x.shape[1]
    ilr = np.log(x) @ (-helmert(D).T)
    c = ilr - ilr.mean(axis=0)
    U, S, Vt = np.linalg.svd(c, full_matrices=False)
    keep = S > 1e-10 * S[0]
    return U[:, keep] * S[keep], D, int(keep.sum())


def design(m, comp=None):
    z = (m.fixed_projection_raw - m.fixed_projection_raw.mean()) / m.fixed_projection_raw.std(ddof=1)
    st = pd.get_dummies(m.stage.astype(str), drop_first=True).astype(float).values
    sex = (m.inferred_sex == 'M').astype(float).values[:, None]
    cols = [np.ones((len(m), 1)), z.values[:, None], st, sex]
    if comp is not None:
        cols.append(comp)
    return np.hstack(cols)


def fit_hc3(Y, X):
    inv = np.linalg.inv(X.T @ X)
    B = inv @ X.T @ Y
    R = Y - X @ B
    h = np.einsum('ij,jk,ik->i', X, inv, X)
    df = X.shape[0] - X.shape[1]
    if h.max() >= 1 - 1e-10 or np.linalg.matrix_rank(X) < X.shape[1]:
        return None
    se = np.empty(Y.shape[1]); ols_se = np.sqrt((R ** 2).sum(0) / df * inv[1, 1])
    for j in range(Y.shape[1]):
        u = R[:, j] / (1 - h)
        meat = (X * u[:, None]).T @ (X * u[:, None])
        V = inv @ meat @ inv
        se[j] = np.sqrt(V[1, 1])
    return B[1], se, ols_se, df, h.max()


def bh(p, n):
    p = np.asarray(p, float); ok = np.isfinite(p); q = np.full(len(p), np.nan)
    pv = p[ok]; o = np.argsort(pv); r = np.empty(len(pv)); r[o] = np.arange(1, len(pv) + 1)
    qq = pv * n / r
    qs = np.minimum.accumulate(qq[o][::-1])[::-1]
    out = np.empty(len(pv)); out[o] = np.minimum(qs, 1)
    q[ok] = out
    return q


fits = [(('U', 'none'))] + [(a, fs) for a in ['A', 'B', 'C', 'D1', 'D2', 'D3'] for fs in ['earlier', 'rebuilt']]
coh_rows, meta_rows, dims = [], [], []
for arm, fs in fits:
    for fam, (wide, nfam) in FAM.items():
        per = {}
        for ds in COH:
            m = man[man.dataset == ds].reset_index(drop=True)
            m = m[m.sample_id.isin(wide.index)].reset_index(drop=True)
            comp = None
            if arm != 'U':
                comp, D, kdim = comp_matrix(arm, fs, m.sample_id.values)
                dims.append({'arm': arm, 'fraction_set': fs, 'dataset': ds, 'family': fam, 'n_parts': D, 'dims_kept': kdim})
            X = design(m, comp)
            Y = wide.loc[m.sample_id.values].values.astype(float)
            good = np.isfinite(Y).all(axis=0)  # features with any missing score are left non-estimable
            b = np.full(Y.shape[1], np.nan); se = b.copy(); ose = b.copy(); df = X.shape[0] - X.shape[1]; hmax = 1.0
            res = fit_hc3(Y[:, good], X)
            if res is not None:
                b[good], se[good], ose[good], df, hmax = res
            tcrit = stats.t.ppf(0.975, df)
            p = 2 * stats.t.sf(np.abs(b / se), df)
            per[ds] = pd.DataFrame({'feature_id': wide.columns, 'beta': b, 'hc3_se': se, 'ols_se': ose,
                                    'hc3_low': b - tcrit * se, 'hc3_high': b + tcrit * se, 'hc3_p': p,
                                    'n': len(m), 'df': df, 'max_leverage': hmax})
            per[ds]['hc3_q'] = bh(p, nfam)
            per[ds]['dataset'] = ds; per[ds]['arm'] = arm; per[ds]['fraction_set'] = fs; per[ds]['family'] = fam
            coh_rows.append(per[ds])
        a, b2 = per[COH[0]].set_index('feature_id'), per[COH[1]].set_index('feature_id')
        w1, w2 = 1 / a.hc3_se ** 2, 1 / b2.hc3_se ** 2
        mb = (a.beta * w1 + b2.beta * w2) / (w1 + w2)
        mse = np.sqrt(1 / (w1 + w2))
        mp = 2 * stats.norm.sf(np.abs(mb / mse))
        md = pd.DataFrame({'beta': mb, 'se': mse, 'low': mb - 1.96 * mse, 'high': mb + 1.96 * mse, 'p': mp,
                           'beta_GSE162694': a.beta, 'beta_GSE213621': b2.beta})
        md['q'] = bh(md.p.values, nfam)
        md['call'] = md.q < 0.05
        md['arm'] = arm; md['fraction_set'] = fs; md['family'] = fam
        meta_rows.append(md.reset_index())
coh = pd.concat(coh_rows, ignore_index=True)
meta = pd.concat(meta_rows, ignore_index=True)
coh.to_csv(OUT + 'cohort_effects_recomputed.tsv.gz', sep='\t', index=False)
meta.to_csv(OUT + 'meta_effects_recomputed.tsv.gz', sep='\t', index=False)
pd.DataFrame(dims).to_csv(OUT + 'ilr_dims.tsv', sep='\t', index=False)
summ = meta.groupby(['arm', 'fraction_set', 'family']).agg(calls=('call', 'sum'), estimable=('beta', lambda s: int(np.isfinite(s).sum()))).reset_index()
print(summ.to_string(), flush=True)
summ.to_csv(OUT + 'calls_by_arm.tsv', sep='\t', index=False)
FOC = {'Stromal ECM': 'hotspot_hepatocytes_f05c535ae5bbc0b9', 'Ductular injury': 'hotspot_hepatocytes_48f39dd4d817a10e'}
foc = meta[meta.feature_id.isin(FOC.values())].copy()
foc['program'] = foc.feature_id.map({v: k for k, v in FOC.items()})
print(foc[['program', 'arm', 'fraction_set', 'beta', 'low', 'high', 'q', 'beta_GSE162694', 'beta_GSE213621']].to_string(), flush=True)
foc.to_csv(OUT + 'focal_meta_recomputed.tsv', sep='\t', index=False)
fc = coh[coh.feature_id.isin(FOC.values())].copy()
fc['program'] = fc.feature_id.map({v: k for k, v in FOC.items()})
print(fc[['program', 'dataset', 'arm', 'fraction_set', 'n', 'beta', 'hc3_low', 'hc3_high', 'hc3_q']].to_string(), flush=True)
fc.to_csv(OUT + 'focal_cohort_recomputed.tsv', sep='\t', index=False)
