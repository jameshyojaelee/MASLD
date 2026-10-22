import pandas as pd, numpy as np
S='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L04-01'
O=pd.read_csv(S+'/pooled_obs.tsv',sep='\t')
def anova_icc(df):
    g=df.groupby('locus').y; k=g.ngroups; N=len(df); ni=g.size().values
    grand=df.y.mean(); msb=(ni*(g.mean().values-grand)**2).sum()/(k-1)
    msw=((df.y-df.groupby('locus').y.transform('mean'))**2).sum()/(N-k)
    n0=(N-(ni**2).sum()/N)/(k-1)
    icc=(msb-msw)/(msb+(n0-1)*msw); 
    mw=(ni**2).sum()/N
    return icc, n0, mw
for t,o in O.groupby('trait'):
    icc,n0,mw=anova_icc(o); 
    print(t,'pooled: ICC(y within locus)',round(icc,3),'n0',round(n0,2),'size-weighted m',round(mw,2),'DE=1+(m_w-1)ICC',round(1+(mw-1)*icc,2),'n_eff',round(len(o)/(1+(mw-1)*icc),1))
    for dsc,od in o.groupby('disc'):
        icc,n0,mw=anova_icc(od); print('  ',dsc,'n',len(od),'loci',od.locus.nunique(),'ICC',round(icc,3),'DE',round(1+(mw-1)*icc,2))
