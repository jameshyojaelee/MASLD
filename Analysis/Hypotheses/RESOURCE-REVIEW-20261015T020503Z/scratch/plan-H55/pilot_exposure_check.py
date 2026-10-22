# H55 planning: re-derive the pilot exposure-validity numbers (exposure vs pathology notes only; no expression).
import pandas as pd
from scipy import stats
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/'
f=pd.read_csv(R+'jobs/3111/morph_fat.tsv',sep='\t'); f=f[f.status=='ok']
g=pd.read_csv(R+'scratch/gate-L02-01/note_grades.tsv',sep='\t')
m=pd.read_csv(R+'scratch/plan-H55/design_meta_counts_only.tsv',sep='\t')
d=f.merge(g[['donor','ord','pct']],on='donor')
k=d[d.ord.notna()]; print('graded with morph: n',len(k),'spearman %.3f p %.2g'%stats.spearmanr(k.fat_frac,k.ord))
k=d[d.pct.notna()]; print('stated pct with morph: n',len(k),'spearman %.3f p %.2g'%stats.spearmanr(k.fat_frac,k.pct))
d=d.merge(m[['donor','SMCENTER','H']],on='donor',how='left')
for c,x in d[d.ord.notna()].groupby('SMCENTER'):
    print('center',c,'graded n',len(x),'spearman %.3f'%stats.spearmanr(x.fat_frac,x.ord)[0], 'median fat %.4f'%x.fat_frac.median())
print('median fat by grade', d.groupby(d.ord).fat_frac.median().round(4).to_dict())
print('tissue_px vs fat spearman (pilot 119): %.3f'%stats.spearmanr(f.tissue_px,f.fat_frac)[0])
