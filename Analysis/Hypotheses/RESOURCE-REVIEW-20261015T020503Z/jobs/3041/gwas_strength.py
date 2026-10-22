# For each pooled discovery observation (gene positive in study A, tested in partner B),
# read GWAS p at the discovery top SNP in A and B, and the minimum p in B within +/-250 kb.
import pandas as pd, numpy as np, sys
S='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L04-01'
J='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3041'
SS='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/data/sumstats/'
O=pd.read_csv(S+'/pooled_obs.tsv',sep='\t')
d=pd.read_pickle(S+'/four.pkl'); d=d[d.ensembl.notna()]
top=d.set_index(['gwas_name','ensembl'])['top_snp']
O['top_snp_A']=[top[(a,g)] for a,g in zip(O.disc,O.ensembl)]
O['chr']=O.top_snp_A.str.split(':').str[0].astype(int); O['pos']=O.top_snp_A.str.split(':').str[1].astype(int)
cache={}
def load(st):
    if st not in cache:
        x=pd.read_csv(SS+st+'_reformatted_hg19.tsv',sep='\t',usecols=['chromosome','position','pval'],dtype={'chromosome':str})
        x=x[x.chromosome.isin([str(c) for c in O.chr.unique()])]
        x['chromosome']=x.chromosome.astype(int)
        cache[st]={c:g.sort_values('position')[['position','pval']].to_numpy() for c,g in x.groupby('chromosome')}
        print('loaded',st,len(x),flush=True)
    return cache[st]
pA=[];pB_same=[];pB_min=[]
for r in O.itertuples():
    A=load(r.disc)[r.chr]; B=load(r.rep)[r.chr]
    ia=np.searchsorted(A[:,0],r.pos); pA.append(A[ia,1] if ia<len(A) and A[ia,0]==r.pos else np.nan)
    ib=np.searchsorted(B[:,0],r.pos); pB_same.append(B[ib,1] if ib<len(B) and B[ib,0]==r.pos else np.nan)
    lo=np.searchsorted(B[:,0],r.pos-250000); hi=np.searchsorted(B[:,0],r.pos+250000,side='right')
    pB_min.append(B[lo:hi,1].min() if hi>lo else np.nan)
O['p_A_top']=pA; O['p_B_same_snp']=pB_same; O['p_B_min_250kb']=pB_min
O.to_csv(J+'/pooled_obs_gwas_strength.tsv',sep='\t',index=False)
print(O[['p_A_top','p_B_same_snp','p_B_min_250kb']].isna().sum().to_dict())
print('done')
