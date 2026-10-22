import pandas as pd, numpy as np
S='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L04-01'
d=pd.read_pickle(S+'/four.pkl'); d=d[d.ensembl.notna()].copy()
def st(h4,h3):
    if pd.isna(h4): return 'no_susie'
    if h4>0.5: return 'H4'
    if h3>0.8: return 'H3'
    return 'other'
d['state']=[st(a,b) for a,b in zip(d['PP.H4.susie'],d['PP.H3.susie'])]
W=d.pivot(index='ensembl',columns='gwas_name',values='state')
P4=d.pivot(index='ensembl',columns='gwas_name',values='PP.H4.susie')
T=d.pivot(index='ensembl',columns='gwas_name',values='top_snp')
G=d.drop_duplicates('ensembl').set_index('ensembl')['gene']
def cluster(pos_df):
    # pos_df: index gene, cols chr,pos ; single linkage 1Mb
    pos_df=pos_df.sort_values(['chr','pos']); lab=[]; k=-1; pc=None; pp=None
    for c,p in zip(pos_df.chr,pos_df.pos):
        if c!=pc or p-pp>1_000_000: k+=1
        lab.append(k); pc,pp=c,p
    return pd.Series(lab,index=pos_df.index)
obs=[]; out={}
for trait,(a,b) in {'ALT':('UKBB_ALT','MVP_ALT_EUR'),'AST':('UKBB_AST','MVP_AST_EUR')}.items():
    either=W[(W[a]=='H4')|(W[b]=='H4')]
    both=either[(either[a]!='no_susie')&(either[b]!='no_susie')]
    ts=[T.loc[g,a] if W.loc[g,a]=='H4' else T.loc[g,b] for g in both.index]
    pdf=pd.DataFrame({'chr':[s.split(':')[0] for s in ts],'pos':[int(s.split(':')[1]) for s in ts]},index=both.index)
    pdf['chr']=pdf.chr.str.replace('chr','')
    lab=cluster(pdf)
    print(trait,'positive in either',len(either),'tested in both',len(both),'loci',lab.nunique())
    print('  genes per locus', lab.value_counts().value_counts().sort_index().to_dict())
    for A,B in [(a,b),(b,a)]:
        for g in both.index:
            if W.loc[g,A]=='H4':
                obs.append(dict(trait=trait,disc=A,rep=B,ensembl=g,gene=G[g],locus=f'{trait}_{lab[g]}',x=P4.loc[g,A],y=int(W.loc[g,B]=='H4')))
O=pd.DataFrame(obs); O.to_csv(S+'/pooled_obs.tsv',sep='\t',index=False)
for t,o in O.groupby('trait'):
    cs=o.groupby('locus').size()
    print(t,'pooled obs',len(o),'unique genes',o.ensembl.nunique(),'loci',o.locus.nunique(),'obs per locus mean',round(cs.mean(),2),'size-weighted mean', round((cs**2).sum()/cs.sum(),2), 'max', cs.max(), 'top loci', o.groupby('locus').gene.apply(lambda s: ','.join(sorted(set(s)))[:80]).loc[cs.sort_values().index[-3:]].tolist())
    print('  x mean',round(o.x.mean(),3),'sd',round(o.x.std(),3),' y mean',round(o.y.mean(),3))
    # loci with all genes replicated / none
    lr=o.groupby('locus').y.agg(['mean','size']); mult=lr[lr['size']>1]
    print('  multi-obs loci',len(mult),' all-or-none outcome',int(((mult['mean']==0)|(mult['mean']==1)).sum()))
# coverage of 436
t12=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
pos436=set(t12.loc[t12.susie_state_t12=='positive','ensembl']); print('436 check', len(pos436))
alst=['UKBB_ALT','MVP_ALT_EUR','UKBB_AST','MVP_AST_EUR']
anyH4=set(W.index[(W[alst]=='H4').any(axis=1)])
print('of 436: H4 in any of 4 ALT/AST EUR studies', len(pos436&anyH4))
testable=set(O.ensembl); print('of 436: positive in a 4-study and tested in partner (has replication state)', len(pos436&testable))
rep=set(O.loc[O.y==1,'ensembl']); print('of 436: replicated (H4 in both biobanks, same trait)', len(pos436&rep))
print('436 genes not H4 in any of the 4 studies', len(pos436-anyH4))
