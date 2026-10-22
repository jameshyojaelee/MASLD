"""Gate L10-02: projected Fig 2D cluster census if six EUR GWAS were fine-mapped genome-wide.
Approximation: each new locus is represented by its distance-clumped (+-500 kb) GWS lead,
which is what 01_identify_lead_snps would centre a 0.5 Mb window on. Clusters are formed
with the Fig 2D rule (1 Mb single linkage on pooled positions of both classes)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
S = Path(__file__).parent
ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
RUN = ROOT / 'figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816'
sys.path.insert(0, str(RUN / 'snapshot/scripts/figures'))
import fig2_multiplicity_correction as m
DIRECT, ENZ = m.SCOPES

def clump(f):
    d = pd.read_csv(S / f'gws_{f}.tsv', sep='\t', header=None,
                    names=['chrom', 'pos', 'p', 'a1', 'a2', 'beta', 'se'])
    d['z2'] = (d.beta / d.se) ** 2          # rank by |z| so p == 0 underflow ties break
    d = d.sort_values('z2', ascending=False)
    taken, leads = {}, []
    for r in d.itertuples():
        L = taken.setdefault(r.chrom, [])
        if all(abs(r.pos - x) > 500000 for x in L):
            L.append(r.pos); leads.append((r.chrom, r.pos))
    n_p0 = int((d.p == 0).sum())
    return pd.DataFrame(leads, columns=['chrom', 'pos']), len(d), n_p0

new = {}
for f in ['MVP_NAFLD_EUR', 'MVP_ALT_EUR', 'MVP_AST_EUR', 'UKBB_ALT', 'UKBB_AST', 'UKBB_GGT']:
    leads, nv, n0 = clump(f)
    new[f] = leads
    print(f'{f}\tgws_variants {nv}\tp_eq_0 {n0}\tclumped_loci {len(leads)}')

a = pd.read_csv(RUN / 'multiplicity/corrected_locus_coordinates.tsv', sep='\t', dtype={'locus': str})
eur = a[a.ancestry == 'EUR'][['study', 'trait_scope', 'chrom', 'pos']].copy()
cur_dir = eur[eur.trait_scope == DIRECT]
cur_enz = eur[eur.trait_scope == ENZ]

def census(direct_pos, enz_pos, label):
    fr = pd.concat([direct_pos.assign(trait_scope=DIRECT), enz_pos.assign(trait_scope=ENZ)])
    fr = m.assign_clusters(fr[['chrom', 'pos', 'trait_scope']].reset_index(drop=True))
    g = fr.groupby('cluster_id').trait_scope.agg(lambda s: set(s))
    nd = sum(DIRECT in s for s in g); ne = sum(ENZ in s for s in g)
    sh = sum(len(s) == 2 for s in g)
    span = fr.groupby('cluster_id').pos.agg(lambda s: s.max() - s.min())
    big = int((span > 5e6).sum())
    print(f'{label}\tphysical {len(g)}\tdirect {nd}\tenzyme {ne}\tshared {sh}\t'
          f'direct_exclusive {nd - sh}\tclusters_spanning_gt5Mb {big}\tmax_span_Mb {span.max()/1e6:.1f}')
    return fr

census(cur_dir[['chrom', 'pos']], cur_enz[['chrom', 'pos']], 'current_EUR')
enz_new = pd.concat([cur_enz[['chrom', 'pos']]] + [new[f] for f in ['MVP_ALT_EUR', 'MVP_AST_EUR', 'UKBB_ALT', 'UKBB_AST', 'UKBB_GGT']])
dir_new = pd.concat([cur_dir[['chrom', 'pos']], new['MVP_NAFLD_EUR']])
census(dir_new, cur_enz[['chrom', 'pos']], 'add_MVP_NAFLD_EUR_only')
fr = census(dir_new, enz_new, 'projected_genome_wide')
# direct clusters alone (no enzyme chaining), to show how many direct clusters MVP NAFLD EUR adds
d_only = m.assign_clusters(dir_new.assign(trait_scope=DIRECT).reset_index(drop=True))
d_cur = m.assign_clusters(cur_dir[['chrom', 'pos']].assign(trait_scope=DIRECT).reset_index(drop=True))
print('direct-only clustering: current', d_cur.cluster_id.nunique(), 'with MVP_NAFLD_EUR', d_only.cluster_id.nunique())
# how many current direct clusters get merged together by enzyme chaining
cur_ids = m.assign_clusters(cur_dir[['chrom','pos']].assign(trait_scope=DIRECT).reset_index(drop=True))
fr2 = m.assign_clusters(pd.concat([cur_dir[['chrom','pos']].assign(trait_scope=DIRECT, tag='cur'),
                                   enz_new.assign(trait_scope=ENZ, tag='enz')]).reset_index(drop=True))
k = fr2[fr2.tag == 'cur'].groupby('cluster_id').size()
print('current 16 direct clusters after genome-wide enzyme chaining ->', fr2[fr2.tag=='cur'].cluster_id.nunique(), 'physical clusters')
