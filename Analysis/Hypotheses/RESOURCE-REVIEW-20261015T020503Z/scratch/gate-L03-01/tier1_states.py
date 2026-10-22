"""Recompute tier1 direct-trait locus states from COLOC r2 under both rules (max over studies; >=2 studies)."""
import pandas as pd, numpy as np
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
tier=pd.read_csv(R+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
t1=set(tier.loc[tier.tier==1,'study_name'])
cols=['gwas_name','gene','ensembl','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf','PP.H4.susie','n_snps','method']
d=pd.read_csv(R+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv',usecols=cols)
d=d[d.gwas_name.isin(t1)].copy()
print('tier1 studies in r2:',d.gwas_name.nunique(),'of',len(t1),'rows',len(d),'genes',d.ensembl.nunique())
d['gs']=d['PP.H1.abf']+d['PP.H3.abf']+d['PP.H4.abf']
d['es']=d['PP.H2.abf']+d['PP.H3.abf']+d['PP.H4.abf']
print('gene-study rows with 0.2<=gs<0.8:',round(((d.gs>=0.2)&(d.gs<0.8)).mean(),3),' gs<0.2:',round((d.gs<0.2).mean(),3),' gs>=0.8:',round((d.gs>=0.8).mean(),3))
print('per study fraction gs<0.2:'); print(d.groupby('gwas_name').gs.apply(lambda x:(x<0.2).mean()).round(3).to_string())
def rule(x,k):
    a=(x['PP.H4.susie']>0.5).sum(); b=(x['PP.H4.abf']>0.5).sum()
    c=((x.gs>=0.8)&(x.es>=0.8)).sum(); dd=(x.gs>=0.8).sum()
    if a>=k: return 'A'
    if b>=k: return 'B'
    if c>=k: return 'C'
    if dd>=k: return 'D'
    if (x.gs<0.2).all(): return 'F'
    return 'G'
for k in [1,2]:
    st=d.groupby('ensembl').apply(lambda x: rule(x,k))
    print(f'k={k} (studies at threshold):',st.value_counts().sort_index().to_dict())
    if k==2: st2=st
    else: st1=st
sym=d.drop_duplicates('ensembl').set_index('ensembl').gene
for s in ['PNPLA3','TM6SF2','MTARC1','HSD17B13','APOE','GCKR','MBOAT7']:
    e=sym[sym==s].index
    print(s, [ (st1.get(i), st2.get(i)) for i in e])
# does the F/G boundary track cis-SNP count?
mx=d.groupby('ensembl').agg(gs=('gs','max'),n=('n_snps','median'))
mx['F']=d.groupby('ensembl').gs.max()<0.2
print('median n_snps F genes',mx[mx.F].n.median(),' non-F',mx[~mx.F].n.median())
print('Spearman(max gs, n_snps):',round(mx.gs.corr(mx.n,method='spearman'),3))
