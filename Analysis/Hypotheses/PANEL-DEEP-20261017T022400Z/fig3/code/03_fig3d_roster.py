"""Recompute the Fig 3D endpoint roster from the 774-unit metadata using the legend rule:
advanced = NAS>=5 or fibrosis>=F3; strict control = source control, no observed NAS>1, no observed fibrosis>F0;
cohort kept if >=5 units in each arm."""
import pandas as pd, numpy as np
O="Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
m=pd.read_csv(O+"meta_774.tsv",sep="\t",keep_default_na=False,na_values=["NA",""])
binned=m.dataset.eq('GSE213621')
adv=(m.group_binary=='Disease')&((m.nas_score>=5)|(m.fibrosis_stage>=3))
ctl=(m.group_binary=='Control')&~(m.nas_score>1)&~(m.fibrosis_stage>0)
m['arm']=np.where(adv,'advanced',np.where(ctl,'strict_control','other'))
t=pd.crosstab(m.dataset,m.arm); print(t)
keep=t.index[(t.get('advanced',0)>=5)&(t.get('strict_control',0)>=5)]
r=m[m.dataset.isin(keep)&m.arm.isin(['advanced','strict_control'])].copy()
print("roster",len(r), r.arm.value_counts().to_dict())
print(pd.crosstab([r.dataset,r.arm],r.fibrosis_stage.fillna(-1)))
# advanced units whose fibrosis value is a GSE213621 bin
a=r[r.arm=='advanced']
print("advanced with binned fibrosis (GSE213621):",(a.dataset=='GSE213621').sum(), "of",len(a))
print("advanced fibrosis NA:",a.fibrosis_stage.isna().sum())
print("advanced qualifying by NAS only (fib<3):",((a.nas_score>=5)&~(a.fibrosis_stage>=3)).sum())
r[['sample_id','dataset','arm','group_binary','fibrosis_stage','nas_score','inferred_sex','condition']].to_csv(O+"fig3d_recomputed_roster.tsv",sep="\t",index=False)
