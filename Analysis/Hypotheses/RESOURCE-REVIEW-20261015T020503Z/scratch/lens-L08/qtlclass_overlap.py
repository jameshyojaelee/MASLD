import pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
G=R+'/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/gtex_v10/'
s=pd.read_csv(R+'/GWAS/finemapping/results/seqfunc/variant_substrate_hg38.tsv',sep='\t',low_memory=False)
s['al']=[frozenset((str(r),str(t))) for r,t in zip(s.ref_hg38,s.alt_hg38)]
s=s.drop_duplicates(['chr','pos_hg38','al'])
def key(df,col):
    p=df[col].str.split('_',expand=True); df['chr']=p[0]; df['pos']=p[1].astype(int); df['al']=[frozenset((r,t)) for r,t in zip(p[2],p[3])]; return df
out=[]
for cls in ['eQTLs','sQTLs','apaQTLs']:
    a=key(pd.read_parquet(G+'Liver.v10.%s.SuSiE_summary.parquet'%cls),'variant_id')
    ncs=a.groupby(['phenotype_id','cs_id']).ngroups
    m=a.merge(s,left_on=['chr','pos','al'],right_on=['chr','pos_hg38','al'])
    m['clpp']=m.pip*m.max_pip
    cs=m.groupby(['phenotype_id','gene_name','cs_id']).agg(clpp=('clpp','sum')).reset_index()
    print(cls,'liver CS',ncs,'genes',a.gene_name.nunique(),'| CS touching GWAS CS',len(cs),'| CLPP>=0.1 CS',(cs.clpp>=0.1).sum(),'genes',cs[cs.clpp>=0.1].gene_name.nunique())
    cs['class']=cls; out.append(cs)
o=pd.concat(out); o.to_csv(G+'../qtlclass_cs_overlap.tsv',sep='\t',index=False)
g01=o[o.clpp>=0.1].groupby('class').gene_name.apply(set)
e=g01.get('eQTLs',set()); 
for c in ['sQTLs','apaQTLs']:
    print(c,'CLPP>=0.1 genes not in eQTL CLPP>=0.1:',sorted(g01.get(c,set())-e))
# small RNA
eg=pd.read_csv(G+'GTEx_Analysis_v10_smallRNA_eQTL/Liver.v10.smallRNA.cis_egenes.txt.gz',sep='\t')
print('smallRNA egenes cols',list(eg.columns)[:12]); 
q=eg[eg.qval<=0.05] if 'qval' in eg else eg
print('smallRNA tested',len(eg),'eGenes q<=0.05',len(q))
sp=pd.read_parquet(G+'GTEx_Analysis_v10_smallRNA_eQTL/Liver.v10.smallRNA.signif_pairs.parquet')
print('signif pairs',sp.shape, list(sp.columns)[:8])
sp=key(sp,'variant_id')
mm=sp.merge(s,left_on=['chr','pos','al'],right_on=['chr','pos_hg38','al'])
print('smallRNA signif pairs at GWAS CS variants',len(mm),'small RNAs',mm.phenotype_id.nunique(),'GWAS genes',mm.gene.nunique())
print(mm.groupby('phenotype_id').agg(maxpip=('max_pip','max'),genes=('gene',lambda x:';'.join(sorted(set(map(str,x)))))).sort_values('maxpip',ascending=False).head(15).to_string())
