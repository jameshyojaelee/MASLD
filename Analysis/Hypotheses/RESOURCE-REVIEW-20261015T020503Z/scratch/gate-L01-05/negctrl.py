import pandas as pd, numpy as np
S='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-05/'
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/'
tss=pd.read_csv(S+'tss_v49.tsv',sep='\t',names=['chr','tss','ens']).drop_duplicates('ens')
st=pd.read_csv(B+'stage_extension_all_gene_results.tsv',sep='\t',usecols=['gene_id_base','axis','contrast'])
ids=set(st.gene_id_base)
h=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv',sep='\t')
c=pd.read_csv(S+'classes_rederived.tsv',sep='\t')
cg=set(c.ensembl)
lc=c[~c.locus.isin(['L090','L091','L092'])].groupby('locus').cls.agg(lambda s:'|'.join(sorted(set(s))))
pure=lc[lc.isin(['GGT_only','ALTAST_only'])]
tss=tss[tss.ens.isin(ids)&~tss.ens.isin(cg)]
tss['chr']=tss.chr.str.replace('chr','')
out=[]
for lid,cl in pure.items():
    r=h[h.locus_id==lid].iloc[0]
    sub=tss[(tss.chr==str(r.chr))&(tss.tss>=r.start-250000)&(tss.tss<=r.end+250000)]
    out.append((lid,cl,len(sub)))
o=pd.DataFrame(out,columns=['locus','cls','n_noncoloc_genes_in_table'])
print(o.groupby('cls').n_noncoloc_genes_in_table.describe())
print('loci with >=1 non-coloc gene:', o[o.n_noncoloc_genes_in_table>0].groupby('cls').size().to_dict())
