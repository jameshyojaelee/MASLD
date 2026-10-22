import pandas as pd, numpy as np
from scipy import stats
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z'
L=R+'/scratch/gate-L02-01'
m=pd.read_csv(L+'/liver_rna_meta.tsv',sep='\t')
g=pd.read_csv(L+'/note_grades.tsv',sep='\t')
f=pd.read_csv(R+'/jobs/3111/morph_fat.tsv',sep='\t'); f=f[f.status=='ok']
d=m.merge(g[['donor','grade','ord','fib_or_nash']],on='donor').merge(f[['donor','fat_frac','white_frac']],on='donor',how='left')
print('n meta',len(d),'n morph',d.fat_frac.notna().sum())
print('center x Hardy (all 226)'); print(pd.crosstab(d.SMCENTER,d.DTHHRDY.fillna(-1)))
print('platform', d.SMGEBTCHT.value_counts().to_dict())
k=d[d.fat_frac.notna()].copy(); k['rf']=k.fat_frac.rank()
for c in ['SMTSISCH','SMRIN','SMEXNCRT','SMRRNART','SMMPPD','SMTSPAX']:
    r=stats.spearmanr(k.fat_frac,k[c],nan_policy='omit'); print(f'morph fat vs {c}: n {k[c].notna().sum()} rho {r[0]:.3f} p {r[1]:.2g}')
print('morph fat by Hardy median/n:'); print(k.groupby('DTHHRDY').fat_frac.agg(['median','count']).round(4))
print('morph fat by center:'); print(k.groupby('SMCENTER').fat_frac.agg(['median','count']).round(4))
print('KW Hardy',stats.kruskal(*[x.fat_frac.values for _,x in k.groupby('DTHHRDY')]))
# note-ordinal on all 223 vs ischemic time and Hardy
d2=d[d.DTHHRDY.notna()].copy(); x=d2.ord.copy(); x[d2.grade=='no_mention']=0; x[d2.grade=='unspecified']=2; d2['x']=x
print('note ordinal vs SMTSISCH all',stats.spearmanr(d2.x,d2.SMTSISCH))
print('note ordinal by Hardy mean:'); print(d2.groupby('DTHHRDY').x.agg(['mean','count']).round(2))
print('graded>=moderate by Hardy:'); print(pd.crosstab(d2.DTHHRDY, d2.grade.isin(['moderate','severe'])))
print('no_mention by center:'); print(pd.crosstab(d2.SMCENTER, d2.grade=='no_mention'))
