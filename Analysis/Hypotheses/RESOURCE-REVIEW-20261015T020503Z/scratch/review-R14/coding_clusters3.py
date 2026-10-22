import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np
d=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv',sep='\t')
d['coding']=d.protein_altering_pip_mass+d.canonical_splice_pip_mass
d=d.sort_values(['chrom','pos']).reset_index(drop=True)
cid=[];c=-1;prev=None
for ch,p in zip(d.chrom,d.pos):
    if prev is None or ch!=prev[0] or p-prev[1]>1_000_000: c+=1
    cid.append(c); prev=(ch,p)
d['cluster']=cid
dr=d[d.trait_scope=='tier1_direct_masld_pdff']
print(dr.groupby('trait').agg(n_cs=('coding','size'),studies=('study','nunique'),cs_mean=('coding','mean')).to_string())
print(dr.groupby('study').agg(n_cs=('coding','size'),trait=('trait','first'),cs_mean=('coding','mean')).to_string())
enz=d[d.trait_scope=='tier2_liver_enzyme'].groupby('cluster').coding.mean().mean()
for lab,sub in [('PDFF only',dr[dr.trait=='PDFF']),('NAFLD+NASH',dr[dr.trait!='PDFF'])]:
    cm=sub.groupby('cluster').coding.mean()
    print(lab,'clusters',len(cm),'cluster-weighted coding',round(cm.mean()*100,3),'ratio vs enzyme',round(cm.mean()/enz,3))
