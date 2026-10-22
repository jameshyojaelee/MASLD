"""Among Catalog genes with a Tier-1/2 SuSiE result and best PP.H4 <= 0.5, how many have PP.H3 > 0.8
(distinct-variant pattern) in the same Tier-1/2 study set, and which 6B cell they occupy."""
import pandas as pd, numpy as np
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
t=pd.read_csv(ROOT+"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t"); t12=set(t[t.tier.isin([1,2])].study_name)
m=pd.read_csv(ROOT+"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gwas_name","ensembl","PP.H3.susie","PP.H4.susie"],low_memory=False)
m=m[m.gwas_name.isin(t12)]
for c in ("PP.H3.susie","PP.H4.susie"): m[c]=pd.to_numeric(m[c],errors="coerce")
g=m.groupby("ensembl").agg(h3=("PP.H3.susie","max"),h4=("PP.H4.susie","max"))
c=pd.read_csv(ROOT+"Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute/catalog_calls_recomputed.tsv.gz",sep="\t")
x=c.join(g,on="gen_key")
sel=x[(x.h4<=0.5)&(x.h3>0.8)]
print("Catalog genes with SuSiE best PP.H4<=0.5 and max PP.H3>0.8:",len(sel))
print(pd.crosstab(sel.gen_call,sel.rna_call,margins=True))
print(sel[sel.symbol.isin(["PNPLA3","TM6SF2","MTARC1","HSD17B13","MBOAT7"])][["symbol","h3","h4"]])
