import pandas as pd, numpy as np
f='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv'
studies=['UKBB_ALT','MVP_ALT_EUR','UKBB_AST','MVP_AST_EUR']
cols=['gwas_name','gene','ensembl','chr','PP.H3.susie','PP.H4.susie','PP.H4.abf','n_cs_pairs','method','top_snp','top_snp_PP','ld_panel','eligibility_source','n_cs_pairs_eligible','locus_lambda_s','ld_reliability']
d=pd.read_csv(f,usecols=cols,low_memory=False)
print('all rows',len(d), 'studies', d.gwas_name.nunique())
d=d[d.gwas_name.isin(studies)].copy()
print(d.groupby('gwas_name').size())
print(d.groupby(['gwas_name','method']).size())
print(d.groupby(['gwas_name','eligibility_source']).size())
# duplicates gene x study
dup=d.groupby(['gwas_name','ensembl']).size()
print('rows per gene x study distribution'); print(dup.value_counts().head())
d.to_pickle('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L04-01/four.pkl')
