import pandas as pd, itertools
f="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv"
d=pd.read_csv(f,sep="\t",usecols=["gene_id_versioned","gene_id_base","gene_name","dataset","logFC","FDR"],keep_default_na=False,na_values=[""])
print("rows",len(d), "datasets", d.dataset.value_counts().to_dict())
print("dup gene x dataset", d.duplicated(["gene_id_versioned","dataset"]).sum())
c=d[(d.FDR<0.05)&(d.logFC.abs()>0.5)].copy()
c["sgn"]=(c.logFC>0).map({True:1,False:-1})
print("calls per cohort", c.dataset.value_counts().to_dict(), "total", len(c))
g=c.groupby("gene_id_versioned")
n=g.size(); nup=g.sgn.apply(lambda s:(s>0).sum()); ndn=n-nup
print("union", len(n))
tier=n.value_counts().sort_index(); print("tiers", tier.to_dict())
multi=n>=2
disc=(nup>0)&(ndn>0)
print("multi-cohort genes", multi.sum(), "discordant", (multi&disc).sum(), "frac", round((multi&disc).sum()/multi.sum(),4))
for k in range(2,6):
    m=n==k; print("tier",k,"n",m.sum(),"discordant",(m&disc).sum())
# GSE126848 involvement
w=c.pivot(index="gene_id_versioned",columns="dataset",values="sgn")
dd=w.loc[disc[disc].index]
print("discordant with GSE126848 call", dd["GSE126848"].notna().sum())
both=w[["GSE126848","GSE135251"]].dropna()
print("GSE126848&GSE135251 shared calls", len(both), "same sign", (both.GSE126848==both.GSE135251).sum(), round((both.GSE126848==both.GSE135251).mean(),4))
# same-direction counting: max(nup,ndn)
same=pd.concat([nup,ndn],axis=1).max(axis=1)
st=same.value_counts().sort_index(); print("same-direction max tier counts", st.to_dict(), "2+", (same>=2).sum())
# without GSE126848
w2=w.drop(columns="GSE126848")
n2=w2.notna().sum(axis=1); up2=(w2==1).sum(axis=1); dn2=(w2==-1).sum(axis=1)
m2=n2>=2; d2=(up2>0)&(dn2>0)
print("without GSE126848: multi", m2.sum(), "discordant", (m2&d2).sum(), round((m2&d2).sum()/m2.sum(),4))
# pairwise sign agreement for all pairs
cols=list(w.columns)
for a,b in itertools.combinations(cols,2):
    x=w[[a,b]].dropna(); print(a,b,"shared",len(x),"same",(x[a]==x[b]).sum(), round((x[a]==x[b]).mean(),3))
# per-cohort up/down
print(c.groupby("dataset").sgn.apply(lambda s:((s>0).sum(),(s<0).sum())).to_dict())
# 5-cohort genes: list discordant
five=n[n==5].index
print("5-cohort discordant genes", c[c.gene_id_versioned.isin(five[disc.loc[five].values])].groupby("gene_id_versioned").apply(lambda x: x.gene_name.iloc[0]+":"+",".join(f"{r.dataset[3:]}{'+' if r.sgn>0 else '-'}" for r in x.itertuples())).tolist())
