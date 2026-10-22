"""Independent re-derivation of the post hoc hepatocyte-specificity class difference:
pair-level and locus-clustered sign-flip, sign test, locus-mean sign test."""
import pandas as pd, numpy as np, sys
from scipy import stats
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
out=sys.argv[1]
p=pd.read_csv(root+'figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference/results/pairs_singlecell_hep_specificity.tsv',sep='\t')
L=pd.read_csv(root+'Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv',sep='\t')
g2l={}
for r in L.itertuples():
    for e in str(r.ensembl).split(','):
        g2l.setdefault(e.strip(),set()).add(r.locus_id)
multi=sum(len(v)>1 for v in g2l.values())
print('H10 loci',len(L),'genes mapped',len(g2l),'genes in >1 locus',multi, 'sum n_genes',L.n_genes.sum())
p['d']=p.y_genetic-p.y_state
p['locus']=p.g_genetic.map(lambda g: ','.join(sorted(g2l.get(g,{'NA'}))))
print('n pairs',len(p),'unique genetic',p.g_genetic.nunique(),'unique state',p.g_state.nunique(),'unmapped',(p.locus=='NA').sum(),'multi-locus',p.locus.str.contains(',').sum())
print('state genes also in genetic set',len(set(p.g_state)&set(p.g_genetic)))
d=p.d.to_numpy(); obs=d.mean()
npos=(d>0).sum(); nneg=(d<0).sum(); nzero=(d==0).sum()
print('mean diff',obs,'median',np.median(d),'pos',npos,'neg',nneg,'zero',nzero)
print('sign test p',stats.binomtest(int(npos),int(npos+nneg),0.5).pvalue)
loci=p.locus.unique(); idx=pd.Categorical(p.locus,categories=loci).codes
print('n loci',len(loci),'pairs per locus max',np.bincount(idx).max())
lm=p.groupby('locus').d.mean()
print('locus means pos',(lm>0).sum(),'neg',(lm<0).sum(),'locus sign test p',stats.binomtest(int((lm>0).sum()),int((lm!=0).sum()),0.5).pvalue)
res={}
for seed in [20260923,20260907]:
    N=100000
    rng=np.random.default_rng(seed)
    # pair-level
    s=rng.choice([-1.0,1.0],size=(N,len(d)))
    nullp=(s*d).mean(axis=1); pp=(1+np.sum(np.abs(nullp)>=abs(obs)))/(1+N)
    rng=np.random.default_rng(seed)
    sl=rng.choice([-1.0,1.0],size=(N,len(loci)))
    nulll=(sl[:,idx]*d).mean(axis=1); pl=(1+np.sum(np.abs(nulll)>=abs(obs)))/(1+N)
    res[seed]=(pp,pl,nullp.std(ddof=1),nulll.std(ddof=1))
    print('seed',seed,'pair signflip p',pp,'locus signflip p',pl,'null sd pair',nullp.std(ddof=1),'null sd locus',nulll.std(ddof=1))
# exact locus-level t-like check: locus sums
ls_=p.groupby('locus').d.sum()
print('locus-sum t-test (one-sample on locus sums, descriptive)',stats.ttest_1samp(ls_,0).pvalue)
p.to_csv(out+'/pairs_with_locus_recomputed.tsv',sep='\t',index=False)
