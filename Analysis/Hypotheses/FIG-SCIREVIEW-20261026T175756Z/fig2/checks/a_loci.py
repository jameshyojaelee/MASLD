# Independent recount of Fig 2A loci column: per-study-locus windows +/-250 kb merged within trait (and globally)
import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/"
tier=pd.read_csv(R+"config/gwas_trait_tier.tsv",sep="\t"); reg=pd.read_csv(R+"config/gwas_registry.tsv",sep="\t")
t12=tier[tier.tier.isin([1,2])].merge(reg[["study_name","ancestry"]],on="study_name")
print("tier1/2 studies", t12.study_name.nunique()); print(t12.groupby(["trait"]).study_name.nunique().to_dict())
anc=dict(zip(t12.study_name,t12.ancestry.replace({"CSA":"SAS"}))); trait=dict(zip(t12.study_name,t12.trait))
cs=pd.read_csv(R+"results/credible_sets.csv",usecols=["chromosome","position","locus","study","recommended_pip","max_pip"],dtype={"locus":str})
cs=cs[cs.study.isin(t12.study_name)].copy(); cs["pip"]=cs.recommended_pip.fillna(cs.max_pip).fillna(0)
w=cs.groupby(["study","locus"]).agg(chrom=("chromosome","first"),mn=("position","min"),mx=("position","max"),pip=("pip","max")).reset_index()
w["anc"]=w.study.map(anc); w["trait"]=w.study.map(trait)
print("study-locus windows", len(w), "studies with any CS", w.study.nunique())
def merged(g,pad=250000):
    iv=sorted(zip(g.mn-pad,g.mx+pad,g.pip,g.anc)); out=[]
    for s,e,p,a in iv:
        if out and s<=out[-1]["e"]:
            out[-1]["e"]=max(out[-1]["e"],e); out[-1]["m"].append((p,a))
        else: out.append({"e":e,"m":[(p,a)]})
    return out
per={}; tied=0; wta={}
for (t,c),g in w.groupby(["trait","chrom"]):
    for m in merged(g):
        per[t]=per.get(t,0)+1
        top=max(p for p,a in m["m"]); win={a for p,a in m["m"] if p==top}
        if len(win)>1: tied+=1; k=(t,"tied")
        else: k=(t,win.pop())
        wta[k]=wta.get(k,0)+1
glob=sum(len(merged(g)) for c,g in w.groupby("chrom"))
print("per trait",per,"total",sum(per.values()),"global unique",glob,"tied",tied)
print(sorted(wta.items()))
