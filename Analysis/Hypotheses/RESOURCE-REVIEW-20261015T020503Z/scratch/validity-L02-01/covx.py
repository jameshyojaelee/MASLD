import pandas as pd, numpy as np
from scipy import stats
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z'
L=R+'/scratch/gate-L02-01'
m=pd.read_csv(L+'/liver_rna_meta.tsv',sep='\t'); g=pd.read_csv(L+'/note_grades.tsv',sep='\t')
d=m.merge(g[['donor','grade','ord']],on='donor'); d=d[d.DTHHRDY.notna()].copy()
x=d.ord.copy(); x[d.grade=='no_mention']=0; x[d.grade=='unspecified']=2; d['x']=x.astype(float)
d['age']=d.AGE.str[:2].astype(int)
print('n',len(d))
print('note ordinal vs age band', stats.spearmanr(d.x,d.age))
print('note ordinal by sex', d.groupby('SEX').x.agg(['mean','count']).round(2).to_dict())
print('note ordinal vs RIN', stats.spearmanr(d.x,d.SMRIN))
print('note ordinal by center', d.groupby('SMCENTER').x.agg(['mean','count']).round(2).to_dict())
print('MW center', stats.mannwhitneyu(d.x[d.SMCENTER=='B1'], d.x[d.SMCENTER=='C1']))
# within-Hardy association of exposure with covariates (what raw permutation breaks)
import numpy.linalg as la
Z=pd.get_dummies(d.DTHHRDY.astype(int).astype(str),drop_first=True).astype(float); Z.insert(0,'i',1.0)
def res(v): 
    Q,_=la.qr(Z.values); return v-Q@(Q.T@v)
for c in ['age','SMTSISCH','SMRIN']:
    print('within-Hardy partial r note ordinal vs',c, round(np.corrcoef(res(d.x.values),res(d[c].values.astype(float)))[0,1],3))
print('within-Hardy partial r note ordinal vs sex', round(np.corrcoef(res(d.x.values),res((d.SEX==2).values.astype(float)))[0,1],3))
print('within-Hardy partial r note ordinal vs centerC1', round(np.corrcoef(res(d.x.values),res((d.SMCENTER=='C1').values.astype(float)))[0,1],3))
