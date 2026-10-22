import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/"
tier=pd.read_csv(R+"config/gwas_trait_tier.tsv",sep="\t")
t12=tier[tier.tier.isin([1,2])]
cs=pd.read_csv(R+"results/credible_sets.csv",usecols=["trait","locus","study","ancestry","chromosome","position","susie_cs","susie_pip"])
cs=cs[cs.study.isin(t12.study_name)]
loc=cs[["study","trait","locus","ancestry"]].drop_duplicates().copy()
loc["chr"]=loc.locus.astype(str).str.split(".").str[0]
def pos_of(l):
    s=str(l); a,b=s.split(".")
    return int(b)
loc["pos"]=loc.locus.map(pos_of)
def nclust(df,w):
    n=0
    for c,g in df.groupby("chr"):
        p=np.sort(g.pos.unique()); 
        n+=1+np.sum(np.diff(p)>w) if len(p) else 0
    return n
for w in [0,1e5,2.5e5,5e5,1e6,2e6,5e6]:
    per={t:nclust(g,w) for t,g in loc.groupby("trait")}
    print(int(w), per, sum(per.values()), nclust(loc,w))
