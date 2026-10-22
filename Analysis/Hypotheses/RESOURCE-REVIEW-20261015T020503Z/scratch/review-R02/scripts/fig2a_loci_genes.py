import pandas as pd, numpy as np
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=B+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R02/'
m=pd.read_pickle(S+'r2_main.pkl')
tier=pd.read_csv(B+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t'); main=set(tier.study_name[tier.placement=='main'])
cs=pd.read_csv(B+'GWAS/finemapping/results/credible_sets.csv',usecols=['chromosome','position','study','locus','ancestry'])
cs=cs[cs.study.isin(main)]
print('main studies with credible sets:',cs.study.nunique(),'of',len(main)); print('missing:',sorted(main-set(cs.study)))
w=cs.groupby(['study','locus']).agg(ch=('chromosome','first'),lo=('position','min'),hi=('position','max'),anc=('ancestry','first')).reset_index()
w['lo']-=250000; w['hi']+=250000
w['ch']=w.ch.astype(str)
print('EUR loci windows', (w.anc=='EUR').sum(),'nonEUR',(w.anc!='EUR').sum())
pos=m[m['PP.H4.susie']>0.5].copy()
pos['c']=pos.top_snp.str.split(':').str[0]; pos['p']=pos.top_snp.str.split(':').str[1].astype(int)
def inwin(r, same_study=False):
    ww=w[(w.ch==r.c)&(w.lo<=r.p)&(w.hi>=r.p)]
    if same_study: ww=ww[ww.study==r.gwas_name]
    return len(ww)>0
pos['in_any']=pos.apply(inwin,axis=1)
pos['in_same']=pos.apply(lambda r: inwin(r,True),axis=1)
g=pos.groupby('gene').agg(any_in=('in_any','any'),same_in=('in_same','any'),studies=('gwas_name',lambda x: set(x)))
print('436 genes: any positive pair whose lead lies in any fine-mapped window:',g.any_in.sum(),'; in a window of the SAME study:',g.same_in.sum(),'of',len(g))
mvpeur={'MVP_ALT_EUR','MVP_AST_EUR','MVP_NAFLD_EUR'}
print('genes supported ONLY by MVP EUR studies (not fine-mapped):',g.studies.apply(lambda s: s<=mvpeur).sum())
print('genes supported only by studies absent from credible_sets.csv or with <=15 loci (MVP EUR + UKBB):',g.studies.apply(lambda s: s<=mvpeur|{'UKBB_ALT','UKBB_AST','UKBB_GGT'}).sum())
