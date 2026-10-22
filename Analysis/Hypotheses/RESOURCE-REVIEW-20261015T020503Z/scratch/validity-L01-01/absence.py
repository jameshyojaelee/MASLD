# Is "only" a tested negative or an absent test? Outcome-free.
import pandas as pd, numpy as np
from pathlib import Path
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
G=B/"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-01"
m=pd.read_csv(G/"class_check.tsv",sep="\t",index_col=0)
lc=m.groupby("locus").cls.nunique(); pure=lc[lc==1].index
p=m[m.locus.isin(pure)&m.cls.isin(["GGT_only","ALTAST_only"])].copy()
c=pd.read_csv(B/"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
   usecols=["gwas_name","ensembl","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","PP.H4.susie","method","eligibility_source"])
t=pd.read_csv(B/"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t")
c=c.merge(t[["study_name","trait","tier"]],left_on="gwas_name",right_on="study_name",how="left")
c=c[c.tier.isin([1,2])&c.ensembl.isin(p.index)]
print("GGT studies:", sorted(c[c.trait=="GGT"].gwas_name.unique()))
print("ALT/AST studies:", sorted(c[c.trait.isin(["ALT","AST"])].gwas_name.unique()))
def summarize(genes, other, label):
    x=c[c.ensembl.isin(genes)&c.trait.isin(other)]
    rows=[]
    for gname in genes:
        y=x[x.ensembl==gname]
        if len(y)==0: rows.append("no_row_any_study"); continue
        hasS=y["PP.H4.susie"].notna().any()
        gw=(y["PP.H3.abf"]+y["PP.H4.abf"]).max()   # ABF posterior that the other trait has an association with a causal variant here
        if hasS: rows.append("susie_tested")
        elif gw>0.5: rows.append("no_susie_but_ABF_says_other_trait_signal")
        else: rows.append("no_susie_and_ABF_H3+H4<=0.5")
    s=pd.Series(rows).value_counts()
    print(label, s.to_dict())
    print("  methods seen:", x.method.value_counts().to_dict())
g=p[p.cls=="GGT_only"].index; a=p[p.cls=="ALTAST_only"].index
summarize(g,["ALT","AST"],"GGT_only genes, ALT/AST rows:")
summarize(a,["GGT"],"ALTAST_only genes, GGT rows:")
