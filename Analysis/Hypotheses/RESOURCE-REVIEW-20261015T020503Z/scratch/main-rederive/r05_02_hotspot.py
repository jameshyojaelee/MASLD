# Independent re-derivation of R05-02 on the 116/117 Hotspot programs (main loop).
import numpy as np, pandas as pd
from scipy.stats import spearmanr
SRC = 'figures/main/fig4_singlecell_programs/source_tables/corrected_774/stage_continuum_comparison/hotspot_stage_vs_continuum.tsv'
ML = 'RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z/programs/hotspot/hotspot_participant_scores.tsv.gz'
MAN = 'figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv'
pub = pd.read_csv(SRC, sep='\t')
s = pd.read_csv(ML, sep='\t', usecols=['program_uid', 'sample_id', 'dataset', 'inferred_sex', 'fibrosis_stage', 'outcome_z'])
man = pd.read_csv(MAN, sep='\t', keep_default_na=False)[['sample_id', 'group_binary']]
s = s.merge(man, on='sample_id', how='left')
s = s[s.fibrosis_stage.notna() & s.dataset.isin(['GSE162694', 'GSE213621'])]  # published n_cohorts = 2
print(s.drop_duplicates('sample_id').groupby(['dataset', 'group_binary']).size())

def fit(d):
    z = (d.fibrosis_stage - d.fibrosis_stage.mean()) / d.fibrosis_stage.std()
    X = np.column_stack([np.ones(len(d)), z, (d.inferred_sex == 'M').astype(float)])
    y = d.outcome_z.values
    B, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ B
    s2 = r @ r / (len(d) - 3)
    return B[1], np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])

def meta(df):
    out = {}
    for pid, g in df.groupby('program_uid'):
        bs = [fit(c) for _, c in g.groupby('dataset')]
        w = np.array([1 / se**2 for _, se in bs]); b = np.array([b for b, _ in bs])
        out[pid] = (b * w).sum() / w.sum()
    return pd.Series(out)

full = meta(s)
dis = meta(s[~((s.dataset == 'GSE213621') & (s.group_binary == 'Control'))])
p = pub[pub.estimable.astype(str) == 'True'].set_index('program_uid')
p = p.join(full.rename('mine_full')).join(dis.rename('mine_dis'))
print('programs', len(p), 'max |mine_full - published beta|', float((p.mine_full - p.beta).abs().max()))
for col in ['beta', 'mine_full', 'mine_dis']:
    rho = spearmanr(p[col], p.continuum_beta).correlation
    same = int((np.sign(p[col]) == np.sign(p.continuum_beta)).sum())
    print(f'{col}: rho vs continuum {rho:.4f}; same direction {same}/{len(p)}')
