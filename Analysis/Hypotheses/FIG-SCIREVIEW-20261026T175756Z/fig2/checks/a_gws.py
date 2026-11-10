# Fig 2A loci: how many of the 360 trait-specific loci contain a credible set at a genome-wide-significant lead
import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/"
tier=pd.read_csv(R+"config/gwas_trait_tier.tsv",sep="\t"); reg=pd.read_csv(R+"config/gwas_registry.tsv",sep="\t")
t12=tier[tier.tier.isin([1,2])].merge(reg[["study_name","ancestry"]],on="study_name")
trait=dict(zip(t12.study_name,t12.trait)); anc=dict(zip(t12.study_name,t12.ancestry))
cs=pd.read_csv(R+"results/credible_sets.csv",usecols=["chromosome","position","locus","study"],dtype={"locus":str})
cs=cs[cs.study.isin(t12.study_name)]
w=cs.groupby(["study","locus"]).agg(chrom=("chromosome","first"),mn=("position","min"),mx=("position","max")).reset_index()
L=[]
for s in w.study.unique():
    try: t=pd.read_csv(R+f"data/lead_snps/{s}_leadSNPs_tiered.tsv",sep="\t",dtype={"locus":str}); t["study"]=s; L.append(t)
    except FileNotFoundError: print("no lead file",s)
L=pd.concat(L)
w=w.merge(L[["study","locus","index_p"]],on=["study","locus"],how="left")
print("windows",len(w),"no lead match",w.index_p.isna().sum(), w[w.index_p.isna()].study.value_counts().to_dict())
w["gws"]=w.index_p<5e-8; w["trait"]=w.study.map(trait)
def merged(g,pad=250000):
    iv=sorted(zip(g.mn-pad,g.mx+pad,g.gws,g.index_p.isna())); out=[]
    for s,e,gw,un in iv:
        if out and s<=out[-1][0]: out[-1][0]=max(out[-1][0],e); out[-1][1]|=gw; out[-1][2]|=un
        else: out.append([e,gw,un])
    return out
res={}
for (t,c),g in w.groupby(["trait","chrom"]):
    for e,gw,un in merged(g):
        k=res.setdefault(t,[0,0,0]); k[0]+=1; k[1]+=int(gw); k[2]+=int((not gw) and un)
print("trait: [loci, loci with a GWS-lead window, non-GWS loci with unmatched window]",res)
print("total",sum(v[0] for v in res.values()),"with GWS",sum(v[1] for v in res.values()))
