"""Fig 2A fine-mapped loci column, recomputed from credible_sets.csv with the stated rule
(per study-locus credible-set member window, +/-250 kb padding, merged within trait;
each merged locus assigned to its strongest-PIP ancestry), plus the alternative that counts
merged loci within each trait x ancestry (no winner-take-all)."""
import pandas as pd, numpy as np, json, os
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/"
tier=pd.read_csv(R+"config/gwas_trait_tier.tsv",sep="\t"); reg=pd.read_csv(R+"config/gwas_registry.tsv",sep="\t")
t12=tier[tier.tier.isin([1,2])].merge(reg[["study_name","ancestry"]],on="study_name")
anc=dict(zip(t12.study_name,t12.ancestry.replace({"CSA":"SAS"}))); trait=dict(zip(t12.study_name,t12.trait))
cs=pd.read_csv(R+"results/credible_sets.csv",usecols=["chromosome","position","locus","study","recommended_pip","max_pip"])
cs=cs[cs.study.isin(t12.study_name)].copy()
cs["pip"]=cs.recommended_pip.fillna(cs.max_pip).fillna(0)
w=cs.groupby(["study","locus"]).agg(chrom=("chromosome","first"),mn=("position","min"),mx=("position","max"),pip=("pip","max")).reset_index()
w["anc"]=w.study.map(anc); w["trait"]=w.study.map(trait)
def merge(df,pad=250_000):
    out=[]
    df=df.assign(ps=df.mn-pad,pe=df.mx+pad).sort_values(["ps","pe","pip"])
    cur=None
    for r in df.itertuples():
        if cur is not None and r.ps<=cur["end"]:
            cur["end"]=max(cur["end"],r.pe)
            if r.pip>cur["pip"]: cur["pip"],cur["anc"]=r.pip,r.anc
            cur["ancs"].add(r.anc)
        else:
            if cur: out.append(cur)
            cur={"end":r.pe,"pip":r.pip,"anc":r.anc,"ancs":{r.anc}}
    if cur: out.append(cur)
    return out
at={}; at_any={}; per_trait={}
for (t,c),g in w.groupby(["trait","chrom"]):
    for m in merge(g):
        at[(m["anc"],t)]=at.get((m["anc"],t),0)+1
        for a in m["ancs"]: at_any[(a,t)]=at_any.get((a,t),0)+1
        per_trait[t]=per_trait.get(t,0)+1
glob=sum(len(merge(g)) for c,g in w.groupby("chrom"))
print("per trait", per_trait, sum(per_trait.values()), "global", glob)
print("winner-take-all by ancestry x trait", dict(sorted(at.items())))
print("loci with any member from ancestry x trait", dict(sorted(at_any.items())))
json.dump({"per_trait":per_trait,"total":sum(per_trait.values()),"global":glob,
           "wta":{f"{a}|{t}":v for (a,t),v in at.items()},"any":{f"{a}|{t}":v for (a,t),v in at_any.items()}},
          open("out_ABCEF/A_loci.json" if os.path.isdir("out_ABCEF") else "A_loci.json","w"),indent=1)
