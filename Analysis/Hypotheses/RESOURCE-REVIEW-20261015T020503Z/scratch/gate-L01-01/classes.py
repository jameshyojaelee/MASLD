import pandas as pd, numpy as np
from pathlib import Path
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
c=pd.read_csv(B/"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gwas_name","gene","ensembl","PP.H4.susie","PP.H4.abf","method"])
print("rows",len(c),"gwas",c.gwas_name.nunique())
print(c.groupby("gwas_name").size().to_string())
t=pd.read_csv(B/"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t")
c=c.merge(t[["study_name","trait","tier"]],left_on="gwas_name",right_on="study_name",how="left")
print("unmatched gwas names:",c[c.tier.isna()].gwas_name.unique())
h=c[(c["PP.H4.susie"]>0.5)&(c.tier.isin([1,2]))]
print("hits",len(h),"genes",h.ensembl.nunique())
print(h.groupby(["trait"]).ensembl.nunique())
print(h.groupby(["gwas_name"]).ensembl.nunique().to_string())
def cls(g):
    tr=set(g.trait); 
    if tr & {"NAFLD","PDFF","NASH"}: return "direct_any"
    a=bool(tr&{"ALT","AST"}); gg="GGT" in tr
    return "GGT_plus_ALTAST" if a and gg else ("GGT_only" if gg else "ALTAST_only")
k=h.groupby("ensembl").apply(cls).rename("cls_recomputed")
L=pd.read_csv(B/"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L01/gene_trait_class_locus.tsv",sep="\t").set_index("ensembl")
m=L.join(k,how="outer")
print("lens genes",len(L),"recomputed genes",len(k),"agree",(m.cls==m.cls_recomputed).sum())
print(pd.crosstab(m.cls.fillna("NA"),m.cls_recomputed.fillna("NA")))
# UKBB-only classes
u=h[h.gwas_name.isin(["UKBB_ALT","UKBB_AST","UKBB_GGT"])]
uk=u.groupby("ensembl").apply(cls).rename("cls_ukbb")
m=m.join(uk,how="left")
m["cls_ukbb"]=m.cls_ukbb.fillna("no_UKBB_enzyme_hit")
print(pd.crosstab(m.cls,m.cls_ukbb))
# GGT evidence by study for GGT-containing genes
gg=h[h.trait=="GGT"].groupby("ensembl").gwas_name.apply(lambda s: ",".join(sorted(set(s))))
print(gg.value_counts().to_string())
m.to_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-01/class_check.tsv",sep="\t")
