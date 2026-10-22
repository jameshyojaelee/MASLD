#!/usr/bin/env python3
"""R14-10 re-derivation: hepatocyte-specificity class difference (281 pairs).
Reproduces the gene-level sign-flip, then a post hoc locus-clustered sign-flip in which
all pairs whose colocalized gene sits in the same H10 locus flip together.
Descriptive sensitivity only; not a replacement test."""
import numpy as np, pandas as pd
from scipy import stats
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
P=R+"figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference/results/pairs_singlecell_hep_specificity.tsv"
L=R+"Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv"
OUT=R+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/Claims/R14-10/"
SEED=20260907; N=100_000
p=pd.read_csv(P,sep="\t"); d=(p.y_genetic-p.y_state).to_numpy()
loci=pd.read_csv(L,sep="\t")
g2l={}
for _,r in loci.iterrows():
    for e in str(r.ensembl).split(","): g2l.setdefault(e,r.locus_id)
p["locus"]=p.g_genetic.map(g2l)
nmiss=int(p.locus.isna().sum())
p["locus"]=p.locus.fillna("unmapped_"+p.g_genetic)  # unmapped genes form their own cluster
p["d"]=d
rows=[]
rng=np.random.default_rng(SEED)
def gene_flip():
    s=0; cnt=0; obs=abs(d.mean())
    for _ in range(N//10000):
        sg=rng.choice([-1.0,1.0],size=(10000,len(d))); nm=(sg*d).mean(1); cnt+=np.sum(np.abs(nm)>=obs)
    return (1+cnt)/(1+N)
rng=np.random.default_rng(SEED); pg=gene_flip()
codes,uniq=pd.factorize(p.locus); K=len(uniq)
rng=np.random.default_rng(SEED); cnt=0; obs=abs(d.mean())
for _ in range(N//10000):
    sl=rng.choice([-1.0,1.0],size=(10000,K)); nm=(sl[:,codes]*d).mean(1); cnt+=np.sum(np.abs(nm)>=obs)
pl=(1+cnt)/(1+N)
# locus-level sign test on locus-mean differences
lm=p.groupby("locus").d.mean(); npos=int((lm>0).sum()); nneg=int((lm<0).sum())
npos_g=int((d>0).sum()); nneg_g=int((d<0).sum())
res=pd.DataFrame([dict(n_pairs=len(d),n_unique_genetic=p.g_genetic.nunique(),n_unique_state=p.g_state.nunique(),
  n_loci=K,n_unmapped_genetic=nmiss,mean_difference=d.mean(),
  gene_signflip_p=pg,locus_signflip_p=pl,gene_sign_pos=npos_g,gene_sign_neg=nneg_g,
  gene_sign_test_p=stats.binomtest(npos_g,npos_g+nneg_g).pvalue,
  locus_mean_pos=npos,locus_mean_neg=nneg,locus_sign_test_p=stats.binomtest(npos,npos+nneg).pvalue,
  seed=SEED,n_draws=N,status="post_hoc_descriptive")])
res.to_csv(OUT+"locus_signflip_result.tsv",sep="\t",index=False)
p.to_csv(OUT+"pairs_with_locus.tsv",sep="\t",index=False)
print(res.T.to_string())
print(p.locus.value_counts().head(8).to_string())
