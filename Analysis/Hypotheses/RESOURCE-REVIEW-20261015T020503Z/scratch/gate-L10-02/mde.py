"""Gate L10-02: SE and minimum detectable difference for the Fig 2D EUR contrast at the
projected cluster counts, using the observed per-cluster (num, den) distributions.
Seed fixed. Normal approximation; arms treated as independent (shared-cluster covariance
ignored, which overstates SE if coding loci carry coding mass in both arms)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
S = Path(__file__).parent
ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
RUN = ROOT / 'figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816'
sys.path.insert(0, str(RUN / 'snapshot/scripts/figures'))
import fig2_multiplicity_correction as m
DIRECT, ENZ = m.SCOPES
rng = np.random.default_rng(20261015)
a = pd.read_csv(RUN / 'multiplicity/corrected_locus_coordinates.tsv', sep='\t', dtype={'locus': str})
a['coding_mass'] = a[m.CODING_COLS].sum(axis=1); a['arch_total'] = a[m.ARCH_COLS].sum(axis=1)

def tab(sub):
    sub = m.assign_clusters(sub.copy())
    return m.cluster_table(sub, 'coding_mass', 'arch_total')
te = tab(a[a.ancestry == 'EUR'])
d16 = te[te.trait_scope == DIRECT][['num', 'den']].to_numpy()
e14 = te[te.trait_scope == ENZ][['num', 'den']].to_numpy()
tp = tab(a)
e228 = tp[tp.trait_scope == ENZ][['num', 'den']].to_numpy()
frac = lambda x: x[:, 0].sum() / x[:, 1].sum()
print(f'EUR direct clusters {len(d16)} fraction {frac(d16):.4f}; per-cluster coding share >0.5 in {(d16[:,0]/d16[:,1] > 0.5).sum()}')
print('EUR direct per-cluster coding share:', np.round(np.sort(d16[:, 0] / d16[:, 1])[::-1], 3).tolist())
print(f'EUR enzyme clusters {len(e14)} fraction {frac(e14):.4f}; pooled enzyme clusters {len(e228)} fraction {frac(e228):.4f}')

B = 20000
def se_arm(pool, n):
    idx = rng.integers(0, len(pool), size=(B, n))
    v = pool[idx]
    f = v[..., 0].sum(1) / v[..., 1].sum(1)
    return f.mean(), f.std(ddof=1)

z = 1.959964 + 0.841621
print('\nscenario\tn_direct\tn_enzyme\tE_direct\tSE_direct\tE_enzyme\tSE_enzyme\tSE_diff\tMDE80\texpected_diff\tpower_normal')
from math import erf, sqrt
Phi = lambda x: 0.5 * (1 + erf(x / sqrt(2)))
rows = []
for lab, dpool_new, n_new in [('current_EUR_16v14', None, 0),
                              ('28 direct: 12 new like current direct', d16, 12),
                              ('28 direct: 12 new like enzyme clusters', e228, 12),
                              ('16 direct (no MVP NAFLD added)', None, 0)]:
    n_enz = 14 if lab.startswith('current') else 577
    epool = e14 if lab.startswith('current') else e228
    if n_new:
        idx_o = rng.integers(0, len(d16), size=(B, 16)); idx_n = rng.integers(0, len(dpool_new), size=(B, n_new))
        v = np.concatenate([d16[idx_o], dpool_new[idx_n]], axis=1)
        fd = v[..., 0].sum(1) / v[..., 1].sum(1); Ed, SEd = fd.mean(), fd.std(ddof=1); nd = 16 + n_new
    else:
        Ed, SEd = se_arm(d16, 16); nd = 16
    Ee, SEe = se_arm(epool, n_enz)
    sed = sqrt(SEd ** 2 + SEe ** 2); diff = Ed - Ee
    pw = Phi(diff / sed - 1.959964)
    print(f'{lab}\t{nd}\t{n_enz}\t{Ed:.4f}\t{SEd:.4f}\t{Ee:.4f}\t{SEe:.4f}\t{sed:.4f}\t{z*sed:.4f}\t{diff:.4f}\t{pw:.3f}')
