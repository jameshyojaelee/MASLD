import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
x=pd.read_pickle("t12.pkl")
enz={"ALT","AST","GGT"}
x["cls"]=np.where(x.trait.isin(enz),"enzyme","direct")
# per gene: max SuSiE PP.H4 over 35 studies (keyed by ensembl)
s=x[x["PP.H4.susie"].notna()].copy()
print("finite SuSiE pairs by class:",s.groupby("cls").size().to_dict(), " >0.5:", s[s["PP.H4.susie"]>0.5].groupby("cls").size().to_dict())
# gene-level positive
gmax=s.sort_values(["PP.H4.susie","gwas_name"],ascending=[False,True]).groupby("ensembl").head(1)
pos=gmax[gmax["PP.H4.susie"]>0.5]
print("genes SuSiE>0.5:",len(pos), pos.cls.value_counts().to_dict(), pos.trait.value_counts().to_dict())
# compare to gene_level file
gl=pd.read_csv(R+"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
gp=gl[gl.coloc_best_susie_pp4>0.5]
tr=pd.read_csv(R+"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t")
gp=gp.merge(tr[["study_name","trait"]],left_on="susie_driving_gwas",right_on="study_name",how="left")
print("gene_level file: ",len(gp), gp.trait.value_counts().to_dict())
# genes tested in both classes
t=s.groupby(["ensembl","cls"])["PP.H4.susie"].max().unstack()
both=t.dropna()
print("genes SuSiE-tested in both classes:",len(both)," enzyme-supported:",(both.enzyme>0.5).sum()," direct-supported:",(both.direct>0.5).sum())
