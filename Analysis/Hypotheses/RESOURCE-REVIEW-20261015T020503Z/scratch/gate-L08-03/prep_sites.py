import pandas as pd
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/'
c=pd.read_csv(D+'credible_sets.csv',low_memory=False)
cs=c[c.susie_cs>0]
u=cs[['chromosome','position']].drop_duplicates()
u['chr']='chr'+u.chromosome.astype(str)
u[['chr','position']].to_csv('sites_hg19.tsv',sep='\t',index=False)
print('unique sites',len(u), 'chromosomes',sorted(u.chromosome.astype(str).unique()))
