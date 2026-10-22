import pandas as pd
GW='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/'
G='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/gtex_v10/'
COMP=str.maketrans('ACGT','TGCA')
c=pd.read_csv(GW+'credible_sets.csv',low_memory=False); c=c[c.susie_cs>0]
lift=pd.read_csv('sites_hg19_to_hg38.tsv',sep='\t'); lift['chromosome']=lift.chr19.str.replace('chr','').astype(int)
c=c.merge(lift,left_on=['chromosome','position'],right_on=['chromosome','pos19'])
a1=c.allele1.str.upper(); a2=c.allele2.str.upper(); m=c.strand38=='-'
a1=a1.where(~m,a1.str.translate(COMP)); a2=a2.where(~m,a2.str.translate(COMP))
c['al']=['|'.join(sorted(x)) for x in zip(a1,a2)]
v=c[['chr38','pos38','al']].drop_duplicates()
q=pd.concat([pd.read_parquet(G+'Liver.v10.%s.SuSiE_summary.parquet'%k,columns=['variant_id']) for k in ['eQTLs','sQTLs','apaQTLs']]).drop_duplicates()
p=q.variant_id.str.split('_',expand=True); q['chr38']=p[0]; q['pos38']=p[1].astype(int); q['al']=['|'.join(sorted(x)) for x in zip(p[2],p[3])]
pm=v.merge(q[['chr38','pos38']].drop_duplicates(),on=['chr38','pos38'])
am=v.merge(q[['chr38','pos38','al']].drop_duplicates(),on=['chr38','pos38','al'])
print('GWAS CS variants (hg38 site+alleles)',len(v),'| site matches in any GTEx liver QTL CS',len(pm),'| site+allele matches',len(am))
x=pm.merge(q[['chr38','pos38','al']],on=['chr38','pos38'],suffixes=('_gwas','_gtex')); x=x[x.al_gwas!=x.al_gtex]
print('site-only mismatches (examples):'); print(x.head(8).to_string())
