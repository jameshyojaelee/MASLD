import pandas as pd, itertools
f="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv"
d=pd.read_csv(f,sep="\t",usecols=["gene_id_versioned","dataset","logFC","FDR","bh_family_size"])
print("rows per dataset:",d.groupby("dataset").size().to_dict(), "family:",d.bh_family_size.unique())
h=d[(d.FDR<0.05)&(d.logFC.abs()>0.5)].copy()
h["sign"]=(h.logFC>0).map({True:1,False:-1})
print("calls per cohort:",h.groupby("dataset").size().to_dict(), "memberships:",len(h))
g=h.groupby("gene_id_versioned").agg(n=("dataset","nunique"),nup=("sign",lambda s:(s>0).sum()),ndn=("sign",lambda s:(s<0).sum()))
print("union:",len(g)); print("tiers:",g.n.value_counts().sort_index().to_dict())
g["disc"]=(g.nup>0)&(g.ndn>0)
m=g[g.n>=2]; print("multi-cohort genes:",len(m),"discordant:",m.disc.sum(), round(100*m.disc.mean(),1))
print("discordant by tier:",m.groupby("n").disc.agg(["sum","size"]).to_dict("index"))
# same-direction count: max(nup,ndn)
g["nsame"]=g[["nup","ndn"]].max(axis=1)
print("same-direction tiers (genes with >=2 cohorts same sign):",g.nsame.value_counts().sort_index().to_dict(), "total >=2:",(g.nsame>=2).sum())
# discordant involving GSE126848
hp=h.pivot_table(index="gene_id_versioned",columns="dataset",values="sign",aggfunc="first")
dg=g[g.disc].index
inv=hp.loc[dg,"GSE126848"].notna().sum(); print("discordant involving GSE126848:",inv,"of",len(dg))
# discordant where GSE126848 is the lone sign (removing it resolves)
# pairwise
for a,b in itertools.combinations(sorted(hp.columns),2):
    both=hp[[a,b]].dropna(); same=(both[a]==both[b]).sum()
    print(f"{a} x {b}: shared {len(both)}, same sign {same} ({100*same/len(both):.1f}%)")
# without GSE126848
h2=h[h.dataset!="GSE126848"]
g2=h2.groupby("gene_id_versioned").agg(n=("dataset","nunique"),nup=("sign",lambda s:(s>0).sum()),ndn=("sign",lambda s:(s<0).sum()))
m2=g2[g2.n>=2]; dd=((m2.nup>0)&(m2.ndn>0)); print("without GSE126848: multi",len(m2),"discordant",dd.sum(),round(100*dd.mean(),1))
# 5-cohort all-same
five=g[g.n==5]; print("5-cohort: all same sign",((five.nup==5)|(five.ndn==5)).sum(),"of",len(five))
