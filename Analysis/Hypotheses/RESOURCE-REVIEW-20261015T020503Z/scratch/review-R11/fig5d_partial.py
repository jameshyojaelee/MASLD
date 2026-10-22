import pandas as pd, numpy as np
P='Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv'
d=pd.read_csv(P,sep='\t',keep_default_na=False,na_values=[''])
print(d.shape, d.primary_signal_pair.astype(str).value_counts().to_dict())
d=d[d.primary_signal_pair.astype(str).str.upper()=='TRUE']
print(d.shape)
print(d.groupby(['trait_class','lineage','evidence_state']).size().unstack(fill_value=0))
p=d[d.evidence_state=='partial']
for thr in [1e-3,1e-2,0.05,0.1,0.25]:
    print('partial rows with any_accessible_mass <',thr, (p.any_accessible_mass<thr).sum(),'/',len(p), round((p.any_accessible_mass<thr).mean(),3))
print(p.groupby('lineage').any_accessible_mass.describe(percentiles=[.25,.5,.75]).round(4))
