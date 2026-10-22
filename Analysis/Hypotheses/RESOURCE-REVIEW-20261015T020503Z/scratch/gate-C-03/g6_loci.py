# Counts eligible loci for the concordance readout WITHOUT computing agreement (outcome-blind).
import pandas as pd, numpy as np
L=pd.read_csv('Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv',sep='\t')
print('loci',len(L),'mhc',int(L.mhc.sum()),'genes listed',L.n_genes.sum())
m=L.assign(ensembl_id=L.ensembl.str.split(',')).explode('ensembl_id')[['locus_id','ensembl_id','mhc']]
print('gene-locus rows',len(m),'unique genes',m.ensembl_id.nunique(),'genes in >1 locus',int((m.groupby('ensembl_id').size()>1).sum()))
E='figures/candidates/fig3e-s3-774-20261011T151616Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv'
st=pd.read_csv(E,sep='\t',usecols=['gene_id_base','axis','contrast','logFC','FDR'])
f4=st[(st.axis=='fibrosis')&(st.contrast=='F4_vs_F0')].set_index('gene_id_base')
anyf=st[st.axis=='fibrosis'].groupby('gene_id_base').FDR.min()
print('F4 genes tested',len(f4),'F4 FDR<0.05',int((f4.FDR<0.05).sum()))
for rule,fn in [('noGTEx','Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/critic-cand-3/r2_gene_direction_collapse.tsv'),('H26_GTEx','Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-C-03/gene_dir_h26rule.tsv')]:
    d=pd.read_csv(fn,sep='\t'); d.columns=['ensembl_id','dir']
    x=m.merge(d,on='ensembl_id',how='left')
    x['call']=x.dir.isin(['call_up','call_down'])
    x['f4sig']=x.ensembl_id.map(f4.FDR).lt(0.05)
    x['f4tested']=x.ensembl_id.isin(f4.index)
    x['anysig']=x.ensembl_id.map(anyf).lt(0.05)
    e=x[x.call&x.f4sig]; ea=x[x.call&x.anysig]
    print(rule,'genes with call',x[x.call].ensembl_id.nunique(),'call & F4 tested',x[x.call&x.f4tested].ensembl_id.nunique(),'call & F4 FDR<0.05 genes',e.ensembl_id.nunique())
    print(' loci with call',x[x.call].locus_id.nunique(),'loci call&F4',e.locus_id.nunique(),'of which MHC',e[e.mhc].locus_id.nunique(),'loci call&anyF',ea.locus_id.nunique(),'MHC',ea[ea.mhc].locus_id.nunique())
    k=e[~e.mhc].groupby('locus_id').ensembl_id.nunique()
    print(' non-MHC eligible loci',len(k),'eligible genes per locus',k.value_counts().sort_index().to_dict(),'loci with even count (tie possible)',int((k%2==0).sum()))
    ka=ea[~ea.mhc].groupby('locus_id').ensembl_id.nunique()
    print(' anyF non-MHC loci',len(ka),'even',int((ka%2==0).sum()))
    # direction-call balance (outcome-blind): share up among eligible loci genes
    print(' dir among eligible genes',e.drop_duplicates('ensembl_id').dir.value_counts().to_dict())
