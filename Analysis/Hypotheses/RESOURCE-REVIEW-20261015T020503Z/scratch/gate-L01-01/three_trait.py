# Classes restricted to GWAS sources that measured ALT, AST and GGT (drop MVP, which has no GGT)
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
S=B/"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-01"
c=pd.read_csv(B/"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gwas_name","ensembl","PP.H4.susie"])
t=pd.read_csv(B/"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t")
c=c.merge(t[["study_name","trait","tier"]],left_on="gwas_name",right_on="study_name")
h=c[(c["PP.H4.susie"]>0.5)&c.tier.isin([1,2])]
m=pd.read_csv(S/"class_check.tsv",sep="\t",index_col=0)
three=["UKBB_ALT","UKBB_AST","UKBB_GGT","BBJ_ALT","BBJ_AST","BBJ_GGT","PanUKBB_AFR_ALT","PanUKBB_AFR_AST","PanUKBB_AFR_GGT","PanUKBB_CSA_ALT","PanUKBB_CSA_AST","PanUKBB_CSA_GGT"]
direct=set(h[h.tier==1].ensembl)
e=h[h.gwas_name.isin(three)]
def cls(g):
    tr=set(g.trait); a=bool(tr&{"ALT","AST"}); gg="GGT" in tr
    return "GGT_plus_ALTAST" if a and gg else ("GGT_only" if gg else "ALTAST_only")
k=e.groupby("ensembl").apply(cls)
m["cls3"]=k.reindex(m.index)
m.loc[m.index.isin(direct),"cls3"]="direct_any"
m["cls3"]=m.cls3.fillna("no_threetrait_enzyme_hit")
print(pd.crosstab(m.cls,m.cls3))
g=m[m.cls3!="no_threetrait_enzyme_hit"]
lc=g.groupby("locus").cls3.nunique(); p=lc[lc==1].index
pc=g[g.locus.isin(p)].groupby("cls3").locus.nunique(); print(pc.to_string())
zt=stats.norm.ppf(0.975)+stats.norm.ppf(0.8)
nG,nA=pc["GGT_only"],pc["ALTAST_only"]; print(f"MDE {zt*np.sqrt(1/nG+1/nA):.3f} SD for {nG} vs {nA}")
# MVP-only contribution to ALTAST_only in primary
mvp=h[h.gwas_name.str.startswith("MVP_")&h.trait.isin(["ALT","AST"])].ensembl
print("primary ALTAST_only genes whose only enzyme hits are MVP:",((m.cls=="ALTAST_only")&(m.cls3=="no_threetrait_enzyme_hit")).sum())
