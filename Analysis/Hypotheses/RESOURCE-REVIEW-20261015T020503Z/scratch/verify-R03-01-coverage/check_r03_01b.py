import warnings; warnings.filterwarnings("ignore")
import pandas as pd
f="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv"
d=pd.read_csv(f,sep="\t",usecols=["gene_id_versioned","dataset","logFC","FDR","n_control","n_disease","covariates"])
print(d.groupby("dataset")[["n_control","n_disease","covariates"]].first())
w=d.pivot(index="gene_id_versioned",columns="dataset",values="logFC")
print("genome-wide Spearman of logFC (all 24,115 genes):"); print(w.corr(method="spearman").round(3))
h=d[(d.FDR<0.05)&(d.logFC.abs()>0.5)]
up=h.assign(s=h.logFC>0).groupby("dataset").s.mean().round(3); print("fraction up among calls:",up.to_dict())
# discordant multi-cohort genes: which pair drives them
hs=h.assign(sign=(h.logFC>0).astype(int)*2-1).pivot(index="gene_id_versioned",columns="dataset",values="sign")
m=hs[hs.notna().sum(axis=1)>=2]
disc=m[(m==1).any(axis=1)&(m==-1).any(axis=1)]
# genes where the discordance is resolved by dropping GSE126848
rest=disc.drop(columns="GSE126848")
resolved=((rest==1).any(axis=1)&(rest==-1).any(axis=1))==False
inv=disc.GSE126848.notna()
print("discordant:",len(disc),"involving GSE126848:",inv.sum(),"; resolved by dropping GSE126848:",(resolved&inv).sum())
# GSE126848 x GSE135251 discordant pairs: direction pattern
both=hs[["GSE126848","GSE135251"]].dropna()
print(pd.crosstab(both.GSE126848,both.GSE135251))
