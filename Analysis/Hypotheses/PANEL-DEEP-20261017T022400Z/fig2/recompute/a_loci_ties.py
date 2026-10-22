"""Check how the ancestry allocation of Fig 2A merged loci depends on ties in the strongest PIP."""
import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/"
tier=pd.read_csv(R+"config/gwas_trait_tier.tsv",sep="\t"); reg=pd.read_csv(R+"config/gwas_registry.tsv",sep="\t")
t12=tier[tier.tier.isin([1,2])].merge(reg[["study_name","ancestry"]],on="study_name")
anc=dict(zip(t12.study_name,t12.ancestry.replace({"CSA":"SAS"}))); trait=dict(zip(t12.study_name,t12.trait))
cs=pd.read_csv(R+"results/credible_sets.csv",usecols=["chromosome","position","locus","study","recommended_pip","max_pip","susie_pip","mesusie_pip"],dtype={"locus":str})
cs=cs[cs.study.isin(t12.study_name)].copy()
cs["pip"]=cs.recommended_pip.fillna(cs.max_pip).fillna(0)
cs["rec_is_mesusie"]=np.isclose(cs.recommended_pip,cs.mesusie_pip)
print("rows where recommended_pip equals mesusie_pip:", int(cs.rec_is_mesusie.sum()), "of", len(cs))
w=cs.groupby(["study","locus"]).agg(chrom=("chromosome","first"),mn=("position","min"),mx=("position","max"),pip=("pip","max")).reset_index()
w["anc"]=w.study.map(anc); w["trait"]=w.study.map(trait)
def merge_tuple(ivs,pad=250_000):
    lst=sorted((mn-pad,mx+pad,pip,a) for mn,mx,pip,a in ivs)
    merged=[]
    for ps,pe,pip,a in lst:
        if merged and ps<=merged[-1][0]:
            if pe>merged[-1][0]: merged[-1][0]=pe
            if pip>merged[-1][1]: merged[-1][1],merged[-1][2]=pip,a
            merged[-1][3].append((pip,a))
        else: merged.append([pe,pip,a,[(pip,a)]])
    return merged
res={}; ties=0; tie_rows=[]
for (t,c),g in w.groupby(["trait","chrom"]):
    for m in merge_tuple(list(zip(g.mn,g.mx,g.pip,g.anc))):
        res[(m[2],t)]=res.get((m[2],t),0)+1
        top=max(p for p,a in m[3]); winners={a for p,a in m[3] if p==top}
        if len(winners)>1: ties+=1; tie_rows.append((t,c,top,sorted(winners),m[2]))
print("tuple-sort allocation:", dict(sorted(res.items())))
print("merged loci whose max PIP is tied across >1 ancestry:", ties)
for r in tie_rows[:15]: print(r)
import collections
print(collections.Counter((r[0],tuple(r[3])) for r in tie_rows).most_common(12))
