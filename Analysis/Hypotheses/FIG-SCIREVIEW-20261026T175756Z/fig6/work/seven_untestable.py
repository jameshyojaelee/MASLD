import pandas as pd
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
tiers=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
t12=set(tiers[tiers.tier.isin([1,2])].study_name)
g=['PSMC5','ASAP3','RITA1','ADAM15','CKB','LILRB3','GSTM2']
rows=[]
for ch in pd.read_csv(C+'susie_coloc_all_gwas.csv',usecols=['gwas_name','gene','PP.H4.susie','method','eligibility_source'],chunksize=200000):
    rows.append(ch[ch.gene.isin(g)&ch.gwas_name.isin(t12)])
a=pd.concat(rows)
for s in g:
    x=a[a.gene==s]
    print(s,'finite SuSiE studies',x['PP.H4.susie'].notna().sum(),'max',round(x['PP.H4.susie'].max(),4),'methods',x.method.value_counts().to_dict())
