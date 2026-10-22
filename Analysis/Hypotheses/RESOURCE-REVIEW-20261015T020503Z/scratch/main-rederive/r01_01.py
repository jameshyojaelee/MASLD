# Independent re-derivation of R01-01 (main loop). Reads saved tables only.
import pandas as pd
from scipy.stats import fisher_exact
c = pd.read_csv('GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
d = pd.read_csv('RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv')
d['ens'] = d.gene.str.split('.').str[0]
print('coloc rows', len(c), 'unique ens', c.ensembl.nunique(), '| deg rows', len(d), 'unique ens', d.ens.nunique())
d = d.drop_duplicates('ens')
d['deg'] = (d.padj < 0.05) & (d.logFC.abs() > 0.5)
print('DEGs genome-wide', int(d.deg.sum()))
m = c.merge(d[['ens', 'deg']], left_on='ensembl', right_on='ens', how='inner')
m['gen'] = m.susie_state_t12.eq('positive')
print(m.susie_state_t12.value_counts().to_dict())
def tab(x, name):
    n = len(x); g = int(x.gen.sum()); e = int(x.deg.sum()); b = int((x.gen & x.deg).sum())
    exp = g * e / n
    t = [[b, g - b], [e - b, n - g - e + b]]
    orr, p = fisher_exact(t)
    print(f'{name}: n={n} genetic={g} deg={e} both={b} expected={exp:.3f} '
          f'exp%_colocrow={100*exp/g:.2f} exp%_degrow={100*exp/e:.2f} obs%_degrow={100*b/e:.2f} OR={orr:.3f} p={p:.3f}')
tab(m, 'joint universe')
print(m.groupby('susie_state_t12').deg.agg(['size', 'sum']))
tab(m[m.susie_state_t12.isin(['positive', 'tested_le_0.5'])], 'SuSiE-evaluated only')
