"""Ancestry-stratified Fig 2D coding/splice contrast using the adopted producer's
own functions (snapshot of fig2_multiplicity_correction.py, job 22218816).
Input: the producer's gated, coordinate-restored table. Read-only."""
import sys
from pathlib import Path
ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
RUN = ROOT / 'figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816'
sys.path.insert(0, str(RUN / 'snapshot/scripts/figures'))
import numpy as np, pandas as pd
import fig2_multiplicity_correction as m

a = pd.read_csv(RUN / 'multiplicity/corrected_locus_coordinates.tsv', sep='\t', dtype={'locus': str})
a['coding_mass'] = a[m.CODING_COLS].sum(axis=1)
a['arch_total'] = a[m.ARCH_COLS].sum(axis=1)

def run(sub, label, fit_id='coding_mass'):
    sub = m.assign_clusters(sub)
    tab = m.cluster_table(sub, 'coding_mass', 'arch_total')
    r = m.analyse(tab, m.fit_rng(fit_id))
    n_inst = sub.groupby('trait_scope').size().to_dict()
    r = {'stratum': label, 'fit_id': fit_id,
         'inst_direct': n_inst.get(m.SCOPES[0], 0), 'inst_enzyme': n_inst.get(m.SCOPES[1], 0), **r}
    return r

rows = []
mvp = a['study'].str.startswith('MVP_')
strata = {
    'pooled_reference': a,
    'EUR_only': a[a.ancestry == 'EUR'],
    'nonEUR_only': a[a.ancestry != 'EUR'],
    'MVP_only_AFR_AMR_EAS': a[mvp],
    'AMR_only': a[a.ancestry == 'AMR'],
    'AFR_only': a[a.ancestry == 'AFR'],
    'EAS_only': a[a.ancestry == 'EAS'],
    'direct_EUR_vs_enzyme_nonEUR': a[((a.trait_scope == m.SCOPES[0]) & (a.ancestry == 'EUR')) |
                                     ((a.trait_scope == m.SCOPES[1]) & (a.ancestry != 'EUR'))],
    'direct_nonEUR_vs_enzyme_EUR': a[((a.trait_scope == m.SCOPES[0]) & (a.ancestry != 'EUR')) |
                                     ((a.trait_scope == m.SCOPES[1]) & (a.ancestry == 'EUR'))],
}
for k, v in strata.items():
    for fid in ('coding_mass', f'ancestry_strat:{k}'):
        rows.append(run(v.copy(), k, fid))
res = pd.DataFrame(rows)
cols = ['stratum', 'fit_id', 'inst_direct', 'inst_enzyme', 'n_clusters_direct', 'n_clusters_enzyme',
        'n_physical_clusters', 'n_shared_clusters', 'fraction_direct', 'fraction_enzyme', 'difference',
        'diff_boot_lo', 'diff_boot_hi', 'ratio', 'ratio_boot_lo', 'ratio_boot_hi', 'n_boot_finite']
res = res.reindex(columns=cols)
out = Path(__file__).parent / 'stratified_coding_contrast.tsv'
res.to_csv(out, sep='\t', index=False)
pd.set_option('display.width', 250)
print(res.to_string(index=False, float_format=lambda x: f'{x:.4f}'))

# Direct-cluster ancestry composition in the pooled clustering
p = m.assign_clusters(a.copy())
d = p[p.trait_scope == m.SCOPES[0]]
comp = d.groupby('cluster_id').agg(has_eur=('ancestry', lambda s: (s == 'EUR').any()),
                                   only_eur=('ancestry', lambda s: (s == 'EUR').all()),
                                   n=('ancestry', 'size'),
                                   anc=('ancestry', lambda s: ','.join(sorted(set(s)))),
                                   coding=('coding_mass', 'mean'), den=('arch_total', 'mean'))
print(comp.to_string())
print('direct clusters', len(comp), 'with any EUR instance', int(comp.has_eur.sum()),
      'EUR-only', int(comp.only_eur.sum()))
e = p[p.trait_scope == m.SCOPES[1]]
ce = e.groupby('cluster_id').ancestry.agg(lambda s: (s == 'EUR').any())
print('enzyme clusters', len(ce), 'with any EUR instance', int(ce.sum()))
# Class fractions by ancestry, separately, using within-class clustering
for scope in m.SCOPES:
    for anc_lab, mask in [('EUR', a.ancestry == 'EUR'), ('nonEUR', a.ancestry != 'EUR')]:
        sub = m.assign_clusters(a[(a.trait_scope == scope) & mask].copy())
        t = m.cluster_table(sub, 'coding_mass', 'arch_total')
        print(scope, anc_lab, 'clusters', len(t), 'fraction', round(t.num.sum() / t.den.sum(), 5))
