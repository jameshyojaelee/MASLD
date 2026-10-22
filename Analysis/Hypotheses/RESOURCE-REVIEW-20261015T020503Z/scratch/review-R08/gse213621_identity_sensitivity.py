import pandas as pd, numpy as np
from scipy import stats
def tie_ceiling(y):
    mr = stats.rankdata(y, method="average")
    return float(np.corrcoef(np.sort(mr), np.arange(1, len(y)+1))[0,1])
S='Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R08/'
m=pd.read_csv(S+'gse213621_scored_with_clusters.tsv',sep='\t')
meta=pd.read_csv('Analysis/Deconvolution/bulk/GSE213621/GSE213621_metadata.tsv',sep='\t')
m=m.merge(meta[['sample_id','fibrotic_stage']],on='sample_id')
arms={'3level':{"F0F1":0,"F2":1,"F3F4":2},'4level':{"Control":0,"F0F1":1,"F2":2,"F3F4":3}}
# conflicting clusters: what stage pairs
conf=m[m.exclusion_reason.fillna('').str.startswith('conflicting')]
print('conflicting clusters stage pairs:'); print(conf.groupby('screen_component').fibrotic_stage.agg(lambda s:'/'.join(sorted(s))).value_counts().to_string())
conc=m[m.screen_component.map(m.screen_component.value_counts())==2]
for arm,mp in arms.items():
    a=m.assign(y=m.fibrotic_stage.map(mp)).dropna(subset=['y'])
    rho=stats.spearmanr(a.latent_severity,a.y).statistic; c=tie_ceiling(a.y)
    q=a[a.exclusion_reason.isna()]
    rq=stats.spearmanr(q.latent_severity,q.y).statistic; cq=tie_ceiling(q.y)
    print(f'{arm}: published-style n={len(a)} clusters={a.screen_component.nunique()} rho={rho:.6f} ceil={c:.6f} ratio={rho/c:.4f}')
    print(f'{arm}: identity-primary n={len(q)} clusters={q.screen_component.nunique()} rho={rq:.6f} ceil={cq:.6f} ratio={rq/cq:.4f}')
    # pairs within arm
    dup=a.screen_component.value_counts(); print(f'  {arm}: clusters with 2 libraries inside arm = {(dup==2).sum()}')
    # concordance of duplicate pairs' scores
    pp=a[a.screen_component.map(dup)==2].sort_values('screen_component')
    if len(pp):
        g=pp.groupby('screen_component').latent_severity.agg(['first','last'])
        print('  within-pair score corr', round(np.corrcoef(g['first'],g['last'])[0,1],4), 'n pairs', len(g))
