# Class sizes when "only" requires the other enzyme to be tested and not colocalized
# by either SuSiE or ABF. Outcome-free (no expression index).
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
G=B/"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-01"
m=pd.read_csv(G/"class_check.tsv",sep="\t",index_col=0)       # 436 genes, locus ids
c=pd.read_csv(B/"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
   usecols=["gwas_name","ensembl","PP.H4.abf","PP.H4.susie"])
t=pd.read_csv(B/"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t")
c=c.merge(t[["study_name","trait","tier"]],left_on="gwas_name",right_on="study_name",how="left")
c=c[c.tier.isin([1,2])&c.ensembl.isin(m.index)]
c["h4any"]=c[["PP.H4.abf","PP.H4.susie"]].max(axis=1)
THREE=lambda s: s.split("_")[0] in ("UKBB","BBJ") or s.startswith("PanUKBB_AFR") or s.startswith("PanUKBB_CSA")
zt=stats.norm.ppf(0.975)+stats.norm.ppf(0.8)
def run(cc,label):
    rows={}
    for e,g in cc.groupby("ensembl"):
        tr=set(g.trait)
        if tr&{"NAFLD","PDFF","NASH"} and (g[g.trait.isin(["NAFLD","PDFF","NASH"])]["PP.H4.susie"]>0.5).any():
            rows[e]="direct_any"; continue
        sA=(g[g.trait.isin(["ALT","AST"])]["PP.H4.susie"]>0.5).any()
        sG=(g[g.trait=="GGT"]["PP.H4.susie"]>0.5).any()
        anyA=(g[g.trait.isin(["ALT","AST"])]["h4any"]>0.5).any()
        anyG=(g[g.trait=="GGT"]["h4any"]>0.5).any()
        testedG=(g.trait=="GGT").any(); testedA=g.trait.isin(["ALT","AST"]).any()
        if sA and sG: rows[e]="GGT_plus_ALTAST"
        elif sA: rows[e]="ALTAST_only" if (testedG and not anyG) else "ambiguous"
        elif sG: rows[e]="GGT_only" if (testedA and not anyA) else "ambiguous"
    k=pd.Series(rows).rename("cls").to_frame().join(m[["locus"]])
    lc=k.groupby("locus").cls.nunique(); pure=lc[lc==1].index
    pl=k[k.locus.isin(pure)].groupby("cls").locus.nunique()
    nG,nA=pl.get("GGT_only",0),pl.get("ALTAST_only",0)
    print(label,"genes:",k.cls.value_counts().to_dict(),"| pure loci:",pl.to_dict(),
          "| MDE %.3f SD"%(zt*np.sqrt(1/nG+1/nA)))
run(c,"23-study, strict absence:")
run(c[c.gwas_name.map(THREE)],"three-enzyme studies, strict absence:")
