import pandas as pd, numpy as np
from scipy import stats
def tie_ceiling(y):
    mr=stats.rankdata(y,method='average'); return float(np.corrcoef(np.sort(mr),np.arange(1,len(y)+1))[0,1])
B='Analysis/MASLD_Model_Benchmark/executions/'
sc=pd.read_csv(B+'external-validation-expansion-20260901T190000Z/GSE276114_scores.tsv',sep='\t')
ds=pd.read_csv(B+'etiology-transfer-fibrosis-distinctness-20260901T135543Z/gse276114_design.tsv',sep='\t')
m=sc.merge(ds,left_on='sample_id',right_on='matrix_column_name',how='inner')
print('joined',len(m),'of',len(sc))
mp={"F0-2":0,"F3":1,"F4":2}; m['y']=m.disease_group.map(mp)
def rep(d,lab):
    r=stats.spearmanr(d.latent_severity,d.y).statistic; c=tie_ceiling(d.y)
    print(f'{lab}: n={len(d)} marg={d.y.value_counts().sort_index().to_dict()} rho={r:.4f} ceil={c:.4f} ratio={r/c:.4f}')
rep(m,'all (as published)')
for e,d in m.groupby('disease'): rep(d,e)
rep(m[m.disease!='MASLD'],'non-MASLD')
print(pd.crosstab(m.disease,m.disease_group))
