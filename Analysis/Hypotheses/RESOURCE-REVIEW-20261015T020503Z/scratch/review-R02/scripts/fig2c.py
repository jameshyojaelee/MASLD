import pandas as pd, numpy as np
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=B+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R02/'
m=pd.read_pickle(S+'r2_main.pkl')
m=m[m.gene.notna()&(m.gene!='')].copy()
m['pp4_best']=m[['PP.H4.susie','PP.H4.abf']].max(axis=1)
m=m.sort_values(['gene','pp4_best','gwas_name'],ascending=[True,False,True],na_position='last')
best=m.groupby('gene').head(1)
best=best[best.pp4_best>0.5]
print('genes',len(best))
vc=pd.read_csv(B+'figures/candidates/fig2-coloc-r2-20260924T053244Z/variant_class/variant_classification.csv')
vc=vc.drop_duplicates('variant_key')
b=best.merge(vc,left_on='top_snp',right_on='variant_key',how='left')
print('unclassified',b.fine_class.isna().sum())
b['fc']=b.fine_class.replace({'spliceSite':'coding','fiveUTR':'UTR','threeUTR':'UTR'})
print(b.fc.value_counts())
print(b.fine_class.value_counts())
print('coding_consequence among coding:'); print(b[b.fc=='coding'].coding_consequence.value_counts(dropna=False))
print(b[b.fc=='coding'][['gene','gwas_name','top_snp','PP.H4.susie','PP.H4.abf','fine_class','coding_consequence']].to_string())
b.to_csv(S+'fig2c_rederived.csv',index=False)
