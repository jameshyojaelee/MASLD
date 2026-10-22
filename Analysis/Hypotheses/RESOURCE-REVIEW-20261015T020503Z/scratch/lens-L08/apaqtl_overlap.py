import pandas as pd
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
a=pd.read_parquet(D+'/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/gtex_v10/Liver.v10.apaQTLs.SuSiE_summary.parquet')
p=a.variant_id.str.split('_',expand=True)
a['chr']=p[0];a['pos']=p[1].astype(int);a['al']=[frozenset((r,t)) for r,t in zip(p[2],p[3])]
s=pd.read_csv(D+'/GWAS/finemapping/results/seqfunc/variant_substrate_hg38.tsv',sep='\t',low_memory=False)
s['al']=[frozenset((str(r),str(t))) for r,t in zip(s.ref_hg38,s.alt_hg38)]
print('substrate variants',len(s),'unique sites',s[['chr','pos_hg38']].drop_duplicates().shape[0],'genes',s.gene.nunique())
m=a.merge(s,left_on=['chr','pos','al'],right_on=['chr','pos_hg38','al'],how='inner')
print('matched rows',len(m))
m['clpp']=m.pip*m.max_pip
cs=m.groupby(['phenotype_id','gene_name','cs_id']).agg(clpp=('clpp','sum'),nvar=('variant_id','nunique'),gwas_gene=('gene',lambda x:';'.join(sorted(set(map(str,x))))),coloc=('colocalizes','max')).reset_index()
print('apaQTL CS sharing >=1 variant with a GWAS CS variant:',len(cs),'apa genes',cs.gene_name.nunique())
for t in [0.01,0.05,0.1,0.2]:
    print('CLPP>=%s'%t,(cs.clpp>=t).sum(),'genes',cs[cs.clpp>=t].gene_name.nunique())
print(cs.sort_values('clpp',ascending=False).head(25).to_string())
cs.to_csv(D+'/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/apaqtl_gwas_cs_overlap.tsv',sep='\t',index=False)
