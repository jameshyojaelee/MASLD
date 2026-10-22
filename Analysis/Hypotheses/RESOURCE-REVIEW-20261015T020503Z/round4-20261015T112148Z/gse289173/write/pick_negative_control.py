# Pick the negative-control hepatocyte program from discovery stage means only (no GSE289173 expression read).
import pandas as pd, numpy as np
G='Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gse289173/gate/'
m=pd.read_csv(G+'discovery_std_stage_means.tsv',sep='\t')
c=pd.read_csv(G+'program_gene_coverage_10x2020A.tsv',sep='\t')
r=pd.read_csv('figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv',sep='\t')
d=m.merge(c,on=['cell_type','module'])
d=d[(d.cell_type=='hepatocytes')&(d.weight_present>=0.80)&(d.present>=8)].copy()
d['maxabs']=d[['m1','m2']].abs().max(axis=1)
d=d.sort_values('maxabs')
nm=[x for x in r.columns if 'name' in x.lower()]
d=d.merge(r[['cell_type','module']+nm],on=['cell_type','module'],how='left')
print(d.head(6).to_string(index=False))
d.to_csv('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gse289173/write/negative_control_candidates.tsv',sep='\t',index=False)
print(d[d.module.isin([20,8])].to_string(index=False))
