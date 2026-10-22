import pandas as pd
S='Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/critic-cand-3/'
g=pd.read_csv(S+'catalog_annotation_preview.tsv',sep='\t')
L=pd.read_csv('Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv',sep='\t')
m=L.assign(ensembl_id=L.ensembl.str.split(',')).explode('ensembl_id')[['locus_id','ensembl_id','mhc']]
x=g[g.gen_sup==True].merge(m,on='ensembl_id',how='left')
x['dc']=x.dir.isin(['call_up','call_down']); x['f4']=x.f4_fdr<0.05; x['anyfib']=x.stage_any_fib_vs_F0==True
print('genes',len(x),'loci',x.locus_id.nunique(),'loci_dir',x[x.dc].locus_id.nunique(),'loci_dir_F4',x[x.dc&x.f4].locus_id.nunique(),'loci_dir_anyfib',x[x.dc&x.anyfib].locus_id.nunique())
