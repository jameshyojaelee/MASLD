import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/"
tier=pd.read_csv(R+"config/gwas_trait_tier.tsv",sep="\t")
t12=tier[tier.tier.isin([1,2])]
cs=pd.read_csv(R+"results/credible_sets.csv")
cs=cs[cs.study.isin(t12.study_name)]
print(cs.columns.tolist())
for col in ["susie_cs","carma_cs","both_in_cs","either_in_cs","in_joint_cs","susie_reliable"]:
    print(col, cs[col].value_counts(dropna=False).head(5).to_dict())
cs["chr"]=cs.locus.astype(str).str.split(".").str[0]
cs["lpos"]=cs.locus.astype(str).str.split(".").str[1].astype(int)
def nclust(pos_by_chr,w):
    n=0
    for c,p in pos_by_chr.items():
        p=np.sort(np.unique(p)); n+=1+np.sum(np.diff(p)>w)
    return n
filters={"all":cs,"susie_cs>0":cs[cs.susie_cs>0],"susie_cs>0&reliable":cs[(cs.susie_cs>0)&(cs.susie_reliable==True)],"either":cs[cs.either_in_cs==True]}
for name,d in filters.items():
    loc=d[["study","trait","locus","chr","lpos"]].drop_duplicates()
    for w in [0,5e5,1e6]:
        per={t:nclust({c:g.lpos.values for c,g in gg.groupby("chr")},w) for t,gg in loc.groupby("trait")}
        print(name,int(w),{k:int(v) for k,v in per.items()},sum(per.values()),nclust({c:g.lpos.values for c,g in loc.groupby("chr")},w), "study-loci",len(loc))
    # bins
    for b in [1e6]:
        loc2=loc.assign(bin=(loc.lpos//b).astype(int))
        per=loc2.groupby("trait").apply(lambda g: len(g[["chr","bin"]].drop_duplicates()))
        print(name,"bin1Mb",per.to_dict(),per.sum(),len(loc2[["chr","bin"]].drop_duplicates()))
