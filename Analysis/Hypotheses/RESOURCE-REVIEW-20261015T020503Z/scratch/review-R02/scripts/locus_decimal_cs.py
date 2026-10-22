import pandas as pd
p="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/credible_sets.csv"
cs=pd.read_csv(p, dtype={"locus":str}, low_memory=False)
print(cs.columns.tolist()[:20], len(cs))
cs["lchr"]=cs["locus"].str.split(".").str[0]
cs["lpos"]=pd.to_numeric(cs["locus"].str.split(".").str[1])
cs["d"]=(cs["position"]-cs["lpos"]).abs()
g=cs.groupby(["study","locus"]).agg(mind=("d","min"),pos=("position","median"),lpos=("lpos","first"))
print("loci",len(g),"loci whose nearest member is >1.5 Mb from decimal-parsed centre:",(g.mind>1.5e6).sum())
print(g[g.mind>1.5e6].head(10))
