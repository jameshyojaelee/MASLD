# v2: exclude 23-study direct_any genes from all runs; add ALT-containing reference group.
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
G=B/"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-01"
m=pd.read_csv(G/"class_check.tsv",sep="\t",index_col=0)
m=m[m.cls!="direct_any"]
c=pd.read_csv(B/"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
   usecols=["gwas_name","ensembl","PP.H4.abf","PP.H4.susie"])
t=pd.read_csv(B/"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t")
c=c.merge(t[["study_name","trait","tier"]],left_on="gwas_name",right_on="study_name",how="left")
c=c[c.tier.isin([1,2])&c.ensembl.isin(m.index)&c.trait.isin(["ALT","AST","GGT"])]
c["h4any"]=c[["PP.H4.abf","PP.H4.susie"]].max(axis=1)
THREE=lambda s: s.startswith(("UKBB_","BBJ_","PanUKBB_AFR_","PanUKBB_CSA_"))
zt=stats.norm.ppf(0.975)+stats.norm.ppf(0.8)
def classify(cc, alt_only_ref=False):
    rows={}
    for e,g in cc.groupby("ensembl"):
        S=lambda tr: (g[g.trait.isin(tr)]["PP.H4.susie"]>0.5).any()
        A=lambda tr: (g[g.trait.isin(tr)]["h4any"]>0.5).any()
        T=lambda tr: g.trait.isin(tr).any()
        sA,sG=S(["ALT","AST"]),S(["GGT"])
        if sA and sG: rows[e]="GGT_plus_ALTAST"
        elif sA:
            if not (T(["GGT"]) and not A(["GGT"])): rows[e]="ambiguous"
            elif alt_only_ref: rows[e]="ALT_containing" if S(["ALT"]) else "AST_only"
            else: rows[e]="ALTAST_only"
        elif sG: rows[e]="GGT_only" if (T(["ALT","AST"]) and not A(["ALT","AST"])) else "ambiguous"
    k=pd.Series(rows).rename("cls").to_frame().join(m[["locus"]])
    lc=k.groupby("locus").cls.nunique(); pure=lc[lc==1].index
    return k, k[k.locus.isin(pure)].groupby("cls").locus.nunique()
for lab,cc in [("23-study",c),("three-enzyme",c[c.gwas_name.map(THREE)])]:
    for ref in (False,True):
        k,pl=classify(cc,ref)
        ra="ALT_containing" if ref else "ALTAST_only"
        nG,nA=pl.get("GGT_only",0),pl.get(ra,0)
        print(f"{lab} strict, ref={ra}: genes {k.cls.value_counts().to_dict()} | pure loci {pl.to_dict()} | MDE GGT vs {ra} {zt*np.sqrt(1/nG+1/nA):.3f} SD")
