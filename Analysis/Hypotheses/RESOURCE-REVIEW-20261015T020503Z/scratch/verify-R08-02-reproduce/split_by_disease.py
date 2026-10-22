import numpy as np, pandas as pd
from scipy import stats
E='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/'
d=pd.read_csv(E+'etiology-transfer-fibrosis-distinctness-20260901T135543Z/gse276114_design.tsv',sep='\t')
s=pd.read_csv(E+'external-validation-expansion-20260901T190000Z/GSE276114_scores.tsv',sep='\t')
m=s.merge(d,left_on='sample_id',right_on='matrix_column_name',how='inner')
print('joined',len(m),'scored',len(s),'unique gsm',m.geo_accession.nunique())
mp={'F0-2':0,'F3':1,'F4':2}
m['y']=m.disease_group.map(mp)
def tie_ceiling(y):
    y=np.asarray(y,float)
    if len(np.unique(y))<2: return float('nan')
    mr=stats.rankdata(y,method='average')
    return float(np.corrcoef(np.sort(mr),np.arange(1,len(y)+1))[0,1])
def row(lbl,sub):
    rho=stats.spearmanr(sub.latent_severity,sub.y).statistic
    c=tie_ceiling(sub.y)
    print(f"{lbl:12s} n={len(sub):3d} counts={dict(sub.y.value_counts().sort_index())} rho={rho:.4f} ceiling={c:.4f} ratio={rho/c:.5f} rho6={rho:.6f} c6={c:.6f}")
row('ALL',m)
for k,g in m.groupby('disease'): row(k,g)
row('nonMASLD',m[m.disease!='MASLD'])
row('MASLD+ARLD',m[m.disease!='CVH'])
# MASLD vs CVH mean score at F4
print(m.groupby(['disease','disease_group']).latent_severity.agg(['count','mean','median']).round(3))
