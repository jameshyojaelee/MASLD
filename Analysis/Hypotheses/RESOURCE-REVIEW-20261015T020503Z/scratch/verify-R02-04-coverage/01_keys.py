import pandas as pd
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=B+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-04-coverage/'
main=pd.read_csv(B+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t').query('placement=="main"').study_name.tolist()
c=pd.read_csv(B+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv',
  usecols=['gwas_name','gene','PP.H4.susie','PP.H4.abf','top_snp','top_snp_PP'],dtype={'top_snp':str,'gene':str})
print('all rows',len(c),'main studies',len(main))
c=c[c.gwas_name.isin(main)]
c=c[c.gene.notna() & (c.gene!='')]
c=c[c.top_snp.notna() & (c.top_snp!='')]
print('main rows with gene+top_snp',len(c))
c.to_csv(S+'coloc_main.csv',index=False)
k=c.top_snp.str.split(':',expand=True)
keys=pd.DataFrame({'s':c.gwas_name,'chr':k[0],'pos':k[1]}).drop_duplicates()
print('unique keys',len(keys))
keys.to_csv(S+'keys.csv',index=False,header=False)
