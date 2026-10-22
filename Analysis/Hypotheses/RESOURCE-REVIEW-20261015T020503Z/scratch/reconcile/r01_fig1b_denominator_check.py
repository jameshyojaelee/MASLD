# Re-derive Fig 1B interface and recompute within SuSiE-evaluated genes.
import pandas as pd
from scipy.stats import fisher_exact
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
c=pd.read_csv(R+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
d=pd.read_csv(R+'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv')
d['base']=d.gene.str.split('.').str[0]
j=c.merge(d[['base','logFC','padj']],left_on='ensembl',right_on='base')
j['deg']=(j.padj<0.05)&(j.logFC.abs()>0.5)
j['gen']=j.coloc_best_susie_pp4>0.5
def rep(x,label):
    t=pd.crosstab(x.gen,x.deg)
    or_,p=fisher_exact(t.values)
    n=len(x);g=x.gen.sum();e=x.deg.sum();b=(x.gen&x.deg).sum()
    print(f'{label}: n={n} genetic={g} deg={e} both={b} expected={g*e/n:.3f} '
          f'exp%_coloc_row={100*e/n:.2f} obs%_coloc_row={100*b/g:.2f} '
          f'exp%_deg_row={100*g/n:.2f} obs%_deg_row={100*b/e:.2f} OR={or_:.4f} p={p:.4f}')
print(pd.crosstab(j.susie_state_t12,j.deg,margins=True))
rep(j,'Fig1B universe (all joint)')
rep(j[j.susie_state_t12.isin(['positive','tested_le_0.5'])],'SuSiE-evaluated only')
