import pandas as pd
G='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/gtex_v10/'
for cls in ['eQTLs','sQTLs','apaQTLs']:
    a=pd.read_parquet(G+'Liver.v10.%s.SuSiE_summary.parquet'%cls)
    if cls=='eQTLs': print(list(a.columns))
    ncs=a.groupby(['phenotype_id','cs_id']).ngroups
    pipcs=a.groupby(['phenotype_id','cs_id']).pip.sum()
    print(cls,'rows',len(a),'unique variants',a.variant_id.nunique(),'phenotypes',a.phenotype_id.nunique(),'CS',ncs,'gene_name',a.gene_name.nunique(),
          'CS size median',a.groupby(['phenotype_id','cs_id']).size().median(),'pip-sum per CS q05/med',round(pipcs.quantile(.05),3),round(pipcs.median(),3))
eg=pd.read_csv(G+'GTEx_Analysis_v10_smallRNA_eQTL/Liver.v10.smallRNA.cis_egenes.txt.gz',sep='\t')
print('smallRNA egenes cols',list(eg.columns))
print('smallRNA tested',len(eg),'q<=0.05',(eg.qval<=0.05).sum(),'biotypes of q<=0.05',eg[eg.qval<=0.05].biotype.value_counts().head(8).to_dict())
print('ma_samples max (proxy lower bound for n)',eg.ma_samples.max())
sp=pd.read_parquet(G+'GTEx_Analysis_v10_smallRNA_eQTL/Liver.v10.smallRNA.signif_pairs.parquet')
print('signif pairs',sp.shape,list(sp.columns))
