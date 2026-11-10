import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/"
tr=pd.read_csv(R+"config/gwas_trait_tier.tsv",sep="\t"); reg=pd.read_csv(R+"config/gwas_registry.tsv",sep="\t")
t12=tr[tr.tier.isin([1,2])].merge(reg[["study_name","ancestry"]],on="study_name")
cs=pd.read_csv(R+"results/credible_sets.csv",usecols=["chromosome","position","locus","study","trait"],dtype={"locus":str})
cs=cs[cs.study.isin(t12.study_name)].merge(t12[["study_name","trait"]].rename(columns={"study_name":"study","trait":"t12trait"}),on="study")
def nk(s): return s.astype(str).str.strip().str.rstrip("0").str.rstrip(".")
cs["k"]=nk(cs.locus)
w=cs.groupby(["t12trait","study","k","chromosome"]).position.agg(["min","max"]).reset_index()
w["lo"]=w["min"]-250_000; w["hi"]=w["max"]+250_000
# lead p
lp=[]
for s in w.study.unique():
    d=pd.read_csv(R+f"data/lead_snps/{s}_leadSNPs_tiered.tsv",sep="\t",dtype={"locus":str})
    d["k"]=nk(d.locus); d["study"]=s; lp.append(d[["study","k","index_p"]])
lp=pd.concat(lp).groupby(["study","k"]).index_p.min().reset_index()
w=w.merge(lp,on=["study","k"],how="left")
print("study-locus windows",len(w)," matched lead p:",w.index_p.notna().sum()," GWS:",(w.index_p<5e-8).sum()," nonGWS:",(w.index_p>=5e-8).sum()," unmatched:",w.index_p.isna().sum())
# merge within trait
out=[]
for (t,c),b in w.sort_values("lo").groupby(["t12trait","chromosome"]):
    cur=None
    for _,r in b.iterrows():
        if cur is None or r.lo>cur["hi"]:
            if cur: out.append(cur)
            cur=dict(trait=t,chr=c,lo=r.lo,hi=r.hi,ps=[r.index_p])
        else:
            cur["hi"]=max(cur["hi"],r.hi); cur["ps"].append(r.index_p)
    out.append(cur)
L=pd.DataFrame(out)
print("loci per trait:",L.groupby("trait").size().to_dict()," total",len(L))
def cat(ps):
    ps=np.array(ps,dtype=float)
    if np.any(ps<5e-8): return "has_GWS"
    if np.all(~np.isnan(ps)): return "all_nonGWS"
    return "unresolved"
L["cat"]=L.ps.apply(cat)
print(L.cat.value_counts().to_dict())
print(pd.crosstab(L.trait,L.cat))
w=w.merge(reg[["study_name","ancestry"]].rename(columns={"study_name":"study"}),on="study")
w["tierp"]=np.where(w.index_p.isna(),"unmatched",np.where(w.index_p<5e-8,"GWS","nonGWS"))
print(pd.crosstab(w.ancestry,w.tierp))
print(w[w.tierp=="unmatched"].study.value_counts().to_dict())
