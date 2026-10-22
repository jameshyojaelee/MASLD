import pandas as pd, numpy as np
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
tier=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
t12=set(tier.study_name[tier.tier.isin([1,2])])
cols=['gwas_name','ensembl','PP.H3.susie','PP.H4.susie','PP.H3.abf','PP.H4.abf']
parts=[]
for ch in pd.read_csv(R+'susie_coloc_all_gwas.csv',usecols=cols,chunksize=500000,low_memory=False):
    parts.append(ch[ch.gwas_name.isin(t12)])
d=pd.concat(parts)
for c in cols[2:]: d[c]=pd.to_numeric(d[c],errors='coerce')
print('rows t12',len(d),'genes',d.ensembl.nunique())
g=d.groupby('ensembl').agg(h3s=('PP.H3.susie','max'),h4s=('PP.H4.susie','max'),h3a=('PP.H3.abf','max'),h4a=('PP.H4.abf','max'))
print('genes with finite susie',g.h4s.notna().sum())
for thr in (0.8,0.9):
    m=(g.h3s>thr)&(g.h4s<0.2)
    ma=(g.h3a>thr)&(g.h4a<0.2)
    print(f'H3.susie>{thr} & H4.susie<0.2:',int(m.sum()),f' H3.abf>{thr} & H4.abf<0.2:',int(ma.sum()))
