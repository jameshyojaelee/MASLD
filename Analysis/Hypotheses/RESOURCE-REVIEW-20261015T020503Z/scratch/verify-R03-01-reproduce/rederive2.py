import pandas as pd, itertools
f="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv"
d=pd.read_csv(f,sep="\t",usecols=["gene_id_versioned","dataset","logFC","FDR"])
L=d.pivot(index="gene_id_versioned",columns="dataset",values="logFC")
print("all-gene Spearman of logFC:"); print(L.corr(method="spearman").round(3))
c=d[(d.FDR<0.05)&(d.logFC.abs()>0.5)].copy(); c["sgn"]=(c.logFC>0)*2-1
w=c.pivot(index="gene_id_versioned",columns="dataset",values="sgn")
for drop in [None]+list(w.columns):
    x=w if drop is None else w.drop(columns=drop)
    n=x.notna().sum(axis=1); up=(x==1).sum(axis=1); dn=(x==-1).sum(axis=1)
    m=n>=2; dd=m&(up>0)&(dn>0)
    print("drop",drop,"multi",m.sum(),"disc",dd.sum(),round(dd.sum()/m.sum(),4))
n=w.notna().sum(axis=1); up=(w==1).sum(axis=1); dn=(w==-1).sum(axis=1); disc=(up>0)&(dn>0)
print("discordant involving GSE135251 call", (disc & w.GSE135251.notna()).sum(), "involving both 126848 & 135251", (disc & w.GSE135251.notna() & w.GSE126848.notna()).sum())
# direction of 135251 vs 126848 shared
both=w[["GSE126848","GSE135251"]].dropna()
print(pd.crosstab(both.GSE126848,both.GSE135251))
