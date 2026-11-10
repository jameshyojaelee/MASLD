import pandas as pd, numpy as np
from scipy.stats import beta
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
d=pd.read_csv(R+"figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_pair_rows_with_background.tsv",sep="\t")
d["mx"]=d[["gse244832_accessible_mass","gse281367_accessible_mass"]].max(axis=1)
for tc,g in d.groupby("trait_class"):
    gm=g.groupby("ensembl").mx.max(); gm2=g.groupby("ensembl").any_accessible_mass.max()
    print(tc,"genes",len(gm),"<5% everywhere (max src)",(gm<0.05).sum(),"(any_accessible_mass)",(gm2<0.05).sum())
def cp(k,n):
    lo=0 if k==0 else beta.ppf(0.025,k,n-k+1); hi=1 if k==n else beta.ppf(0.975,k+1,n-k); return round(100*lo,1),round(100*hi,1)
for k,n in [(0,22),(3,22),(4,22),(0,424),(15,424),(56,424)]: print("CP",k,n,cp(k,n))
E=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5E.tsv",sep="\t")
print("5E max n_open_genes",E.n_open_genes.max(), E.loc[E.n_open_genes.idxmax(),["cell_type","module","program_name"]].tolist())
C=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5C.tsv",sep="\t")
print("5C filled by feature:",C.groupby("feature").filled.sum().to_dict())
print("5C sign of rho by module x feature for filled:")
B=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5B.tsv",sep="\t")
disc=B[B["class"]=="discordant"]; print("5B discordant n",len(disc),"|mRNA log2FC| median",disc.bulk_logFC.abs().median().round(3),"max",disc.bulk_logFC.abs().max().round(3),"n with RNA padj<0.05",(disc.bulk_padj<0.05).sum())
