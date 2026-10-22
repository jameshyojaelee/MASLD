import pandas as pd, numpy as np
f="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv"
d=pd.read_csv(f,sep="\t",keep_default_na=False,na_values=[""])
print("rows",len(d),"pairs",d[["gwas_name","ensembl"]].drop_duplicates().shape[0])
print(d.evidence_state.value_counts())
d["maxsrc"]=d[["gse244832_accessible_mass","gse281367_accessible_mass"]].max(axis=1)
p=d[d.evidence_state=="partial"]
print("partial rows",len(p))
print("min partial maxsrc",p.maxsrc.min(),"min partial any",p.any_accessible_mass.min())
print("partial maxsrc quantiles",p.maxsrc.quantile([0.1,0.25,0.5,0.75,0.9]).to_dict())
# hepatocyte liver enzyme partial rows
hp=p[(p.lineage=="hepatocyte")&(p.trait_class=="liver_enzyme")]
print("hep LE partial",len(hp),"maxsrc<0.01",(hp.maxsrc<0.01).sum(),"any<0.01",(hp.any_accessible_mass<0.01).sum())
# per pair across lineages
for col in ["any_accessible_mass","maxsrc"]:
    g=d.groupby(["gwas_name","ensembl","trait_class"])[col].max().reset_index()
    print("==",col,"pairs",len(g))
    for thr in [0.01,0.5]:
        lt=g[g[col]<thr]
        print(f" max over lineages < {thr}: {len(lt)} of {len(g)}; by class", lt.trait_class.value_counts().to_dict(), "class totals", g.trait_class.value_counts().to_dict())
    print(" reach>=0.5 any lineage:",(g[col]>=0.5).sum())
# mean any-lineage mass per lineage/class
print(d.groupby(["trait_class","lineage"]).any_accessible_mass.mean().round(4))
print("mean any by lineage", d.groupby("lineage").any_accessible_mass.mean().round(4).to_dict())
# per-state counts per class x lineage, check panel
print(d.pivot_table(index=["trait_class","lineage"],columns="evidence_state",values="ensembl",aggfunc="count",fill_value=0))
# partial rows <0.01 overall, and by class
print("partial rows maxsrc<0.01:",(p.maxsrc<0.01).sum(),"of",len(p), p[p.maxsrc<0.01].trait_class.value_counts().to_dict())
print("partial rows any<0.01:",(p.any_accessible_mass<0.01).sum())
