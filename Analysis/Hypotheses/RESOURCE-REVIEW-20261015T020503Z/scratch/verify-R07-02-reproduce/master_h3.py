import pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
m=pd.read_csv(R+'susie_coloc_all_gwas.csv',usecols=['gwas_name','gene','ensembl','PP.H3.susie','PP.H4.susie','method','n_cs_pairs','n_cs_pairs_eligible'])
print(m.shape, m.method.value_counts().to_dict())
su=m[m['PP.H4.susie'].notna()]
print('gene-study rows with SuSiE PP4:',len(su),' with PP.H3.susie non-NA:',su['PP.H3.susie'].notna().sum(),' genes:',su.ensembl.nunique())
s=pd.read_csv(R+'susie_coloc_signal_pairs.tsv.gz',sep='\t')
print('signal-pair file genes:',s.ensembl.nunique(),'gene-study combos:',s[['ensembl','gwas']].drop_duplicates().shape[0])
g=su.groupby('ensembl').agg(maxH3=('PP.H3.susie','max'),maxH4=('PP.H4.susie','max')).reset_index()
a=g[(g.maxH3>0.8)&(g.maxH4<=0.5)]; b=g[(g.maxH3>0.9)&(g.maxH4<0.1)]
print('master gene-study level (all studies): H3>0.8&maxH4<=0.5 genes:',len(a),' H3>0.9&maxH4<0.1:',len(b))
e=pd.read_parquet(C+'passport_evidence_long.parquet'); gen=e[e.evidence_domain=='genetics'][['ensembl_id','call_state']]
print(a.merge(gen,left_on='ensembl',right_on='ensembl_id',how='left').call_state.value_counts(dropna=False).to_dict())
