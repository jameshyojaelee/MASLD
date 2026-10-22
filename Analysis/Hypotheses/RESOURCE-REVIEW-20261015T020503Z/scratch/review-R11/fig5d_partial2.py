import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',20)
P='Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv'
d=pd.read_csv(P,sep='\t',keep_default_na=False,na_values=[''])
t=d.groupby(['trait_class','lineage','evidence_state']).size().unstack(fill_value=0)
print(t)
# reclassify partial with any mass < 0.01 as no overlap
d['state_1pct']=np.where((d.evidence_state=='partial')&(d.any_accessible_mass<0.01),'indeterminate(<1% mass)',d.evidence_state)
print(d.groupby(['trait_class','lineage','state_1pct']).size().unstack(fill_value=0))
print('min any mass among partial:', d.loc[d.evidence_state=='partial','any_accessible_mass'].min())
print('partial rows with exactly 0 any mass:', ((d.evidence_state=='partial')&(d.any_accessible_mass==0)).sum())
print('mean any_accessible_mass by lineage, all rows:'); print(d.groupby('lineage').any_accessible_mass.mean().round(4))
