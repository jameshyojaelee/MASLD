"""Feasibility gate for candidate L03-01. Reads saved tables only (login-node light).

Does NOT compute the primary contrast (state vs matched non-state eGene OR).
Computes: universe and class n, eGene status under the three eQTL SuSiE fit sets,
matching-stratum occupancy, placebo-pool size, LOEUF coverage and the LOEUF
positive control (unmatched and stratum-adjusted), MDE.
"""
import glob
import numpy as np
import pandas as pd

R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S = R + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L03-01/'
P = R + 'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
rng = np.random.default_rng(20261015)


def base(s):
    return s.astype(str).str.split('.').str[0]


g = pd.read_parquet(P + 'passport_gene_index.parquet')
e = pd.read_parquet(P + 'passport_evidence_long.parquet')
g['ens'] = base(g.ensembl_id)
e['ens'] = base(e.ensembl_id)
tt = set(e[(e.evidence_domain == 'transcriptomics') & (e.call_state != 'untestable')].ens)
gt = set(e[(e.evidence_domain == 'genetics') & (e.call_state != 'untestable')].ens)
joint = tt & gt
print('joint universe', len(joint), '| transcriptomics testable', len(tt), '| genetics testable', len(gt))
gj = g[g.ens.isin(joint)].set_index('ens')
print('classes in joint:', gj.primary_evidence_class.value_counts().to_dict())


def load_eq(d):
    fs = sorted(glob.glob(R + f'GWAS/finemapping/results/{d}/chr*_summary.tsv'))
    x = pd.concat([pd.read_csv(f, sep='\t') for f in fs])
    x['ens'] = base(x.ensg)
    return x.drop_duplicates('ens', keep='last').set_index('ens')


old = load_eq('eqtl_susie')
pf = load_eq('eqtl_susie_polyfun')
rf = load_eq('eqtl_susie_refit')
# COLOC r2 fit set: polyfun, with the 204 non-converged genes replaced by the refit
cur = pf.copy()
cur.loc[rf.index.intersection(cur.index), ['n_snps', 'converged', 'n_cs', 'max_pip', 'status']] = \
    rf.loc[rf.index.intersection(cur.index), ['n_snps', 'converged', 'n_cs', 'max_pip', 'status']]
for nm, x in [('old eqtl_susie (Apr)', old), ('polyfun', pf), ('polyfun+refit', cur)]:
    xj = x[x.index.isin(joint)]
    eg = xj.n_cs.gt(0)
    print(f'{nm}: joint in panel {len(xj)}, status {xj.status.value_counts().to_dict()}, '
          f'eGenes {int(eg.sum())} ({eg.mean():.4f})')

gj['egene_old'] = old.n_cs.reindex(gj.index).gt(0)
gj['egene_cur'] = cur.n_cs.reindex(gj.index).gt(0)
gj['fit_status'] = cur.status.reindex(gj.index).fillna('absent')
gj['n_snps'] = cur.n_snps.reindex(gj.index)
print('genes whose eGene label differs old vs current, joint:', int((gj.egene_old != gj.egene_cur).sum()))
print('fit status in joint (current):', gj.fit_status.value_counts().to_dict())
for c in ['genetically_anchored', 'concordant']:
    x = gj[gj.primary_evidence_class == c]
    print(f'{c}: n {len(x)}, eGene (current) {int(x.egene_cur.sum())}')

# covariates
co = pd.read_csv(R + 'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/'
                 'primary_source_controls/coefficient_table.tsv.gz', sep='\t')
co['ens'] = base(co.gene)
co = co.drop_duplicates('ens').set_index('ens')
gl = pd.read_csv(S + 'gene_length_v49.tsv', sep='\t', header=None,
                 names=['ens', 'chr', 'start', 'end', 'len']).drop_duplicates('ens').set_index('ens')
gj['AveExpr'] = co.AveExpr.reindex(gj.index)
gj['logFC'] = co.logFC.reindex(gj.index)
gj['padj'] = co.padj.reindex(gj.index)
gj['len'] = gl.len.reindex(gj.index)
print('covariate missing: AveExpr', int(gj.AveExpr.isna().sum()), 'len', int(gj.len.isna().sum()),
      'n_snps', int(gj.n_snps.isna().sum()), 'n_snps==0', int((gj.n_snps == 0).sum()))

gj['case'] = gj.primary_evidence_class.isin(['established_state_associated', 'concordant'])
cc = gj.dropna(subset=['AveExpr', 'len', 'n_snps']).copy()
cc = cc[cc.fit_status == 'ok']
print('complete-case genes with ok fit:', len(cc), 'cases', int(cc.case.sum()))
for v in ['AveExpr', 'len', 'n_snps']:
    cc['d_' + v] = pd.qcut(cc[v].rank(method='first'), 10, labels=False)
cc['stratum'] = cc.d_AveExpr * 100 + cc.d_len * 10 + cc.d_n_snps

# Occupancy for the candidate's control pool (unresolved) and for all non-state genes
for pool_name, pool_mask in [('unresolved only (as written)', cc.primary_evidence_class == 'unresolved'),
                             ('all non-state genes', ~cc.case)]:
    t = pd.DataFrame({'case': cc[cc.case].groupby('stratum').size(),
                      'ctrl': cc[pool_mask].groupby('stratum').size()}).fillna(0)
    need = 5 * t.case
    short = t[t.ctrl < need]
    print(f'[{pool_name}] strata with cases {int((t.case > 0).sum())}; cases in strata with zero controls '
          f'{int(t.loc[t.ctrl == 0, "case"].sum())}; cases in strata with < 5 controls per case '
          f'{int(short.case.sum())}; max case:ctrl stratum {t[t.case > 0].eval("case/(ctrl+1e-9)").max():.2f}; '
          f'controls available for 5:1 without replacement {int(np.minimum(t.ctrl, need).sum())} of {int(need.sum())}')

print('case AveExpr deciles:', cc[cc.case].d_AveExpr.value_counts().sort_index().to_dict())

# placebo pool
pl = cc[(cc.primary_evidence_class == 'unresolved') & (cc.AveExpr >= cc.AveExpr.median()) &
        (cc.logFC.abs() < 0.05) & (cc.padj > 0.5)]
print('placebo pool (unresolved, AveExpr top half, |logFC|<0.05, padj>0.5):', len(pl), '(candidate needs 1,155)')

# LOEUF positive control
gn = pd.read_csv(R + 'data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv', sep='\t',
                 usecols=['gene', 'gene_id', 'mane_select', 'lof.oe_ci.upper'])
gn = gn[gn.mane_select.astype(str).str.lower().isin(['true', '1'])]
gn['ens'] = base(gn.gene_id)
gn = gn.drop_duplicates('ens').set_index('ens')
cc['LOEUF'] = gn['lof.oe_ci.upper'].reindex(cc.index)
lo = cc.dropna(subset=['LOEUF']).copy()
print('joint complete-case genes with MANE LOEUF:', len(lo))
lo['dL'] = pd.qcut(lo.LOEUF.rank(method='first'), 10, labels=False)
bot, top = lo[lo.dL == 0], lo[lo.dL == 9]
print(f'unmatched eGene rate bottom LOEUF decile {bot.egene_cur.mean():.3f} (n {len(bot)}), '
      f'top decile {top.egene_cur.mean():.3f} (n {len(top)})')
print('eGene rate by LOEUF decile:', lo.groupby('dL').egene_cur.mean().round(3).to_dict())


# Stratum-adjusted (CMH) OR, bottom vs top LOEUF decile, same 10x10x10 strata
def cmh_arr(code, grp, out):
    n = np.bincount(code).astype(float)
    a = np.bincount(code, weights=(grp & out)); b = np.bincount(code, weights=(grp & ~out))
    c = np.bincount(code, weights=(~grp & out)); d = np.bincount(code, weights=(~grp & ~out))
    m = n > 0
    return (a[m] * d[m] / n[m]).sum() / (b[m] * c[m] / n[m]).sum()


bt = lo[lo.dL.isin([0, 9])].sort_values('stratum').copy()
code = pd.factorize(bt.stratum)[0]
grp = (bt.dL == 0).values
out = bt.egene_cur.values.astype(bool)
or_obs = cmh_arr(code, grp, out)
nulls = np.empty(10000)
for i in range(10000):
    perm = np.lexsort((rng.random(len(code)), code))
    nulls[i] = cmh_arr(code, grp[perm], out)
p_lo = (1 + (nulls <= or_obs).sum()) / (len(nulls) + 1)
inf = pd.Series(grp).groupby(code).nunique().gt(1).sum()
print(f'LOEUF PC stratum-adjusted CMH OR bottom vs top decile {or_obs:.3f}; informative strata {inf}; '
      f'genes in informative strata {int(pd.Series(grp).groupby(code).transform("nunique").gt(1).sum())}; '
      f'null OR median {np.median(nulls):.3f}; one-sided stratified-permutation p (10,000 draws, seed 20261015) {p_lo:.5f}')

# Case-gene spot checks
for sym in ['PNPLA3', 'TM6SF2', 'MTARC1', 'HSD17B13', 'APOE', 'GCKR', 'MBOAT7']:
    x = g[g.symbol == sym]
    if len(x) == 0:
        print(sym, 'absent from catalog'); continue
    ens = x.ens.iloc[0]
    print(sym, ens, x.primary_evidence_class.iloc[0], '| in joint', ens in joint,
          '| fit', cur.status.get(ens, 'absent'), 'n_cs', cur.n_cs.get(ens, np.nan),
          '| old n_cs', old.n_cs.get(ens, np.nan))

# MDE, two-proportion normal approximation, independent genes
from scipy.stats import norm
z = norm.ppf(0.975) + norm.ppf(0.8)
n1 = int(cc.case.sum())
for p0 in [cc[~cc.case].egene_cur.mean()]:
    for k in [5]:
        n2 = k * n1
        se = np.sqrt(p0 * (1 - p0) * (1 / n1 + 1 / n2))
        d = z * se
        for sgn in [-1, 1]:
            p1 = p0 + sgn * d
            orr = (p1 / (1 - p1)) / (p0 / (1 - p0))
            print(f'MDE cases {n1} vs {n2} controls, p0 {p0:.4f}: abs diff {d:.4f}, OR {orr:.3f}')
cc[['primary_evidence_class', 'case', 'egene_cur', 'egene_old', 'fit_status', 'AveExpr', 'len', 'n_snps',
    'stratum', 'LOEUF']].to_csv(S + 'gate_gene_table.tsv', sep='\t')
