# Outcome-free checks on the class labels (no expression index used).
import pandas as pd, numpy as np
from pathlib import Path
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
G=B/"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-01"
m=pd.read_csv(G/"class_check.tsv",sep="\t",index_col=0)
lc=m.groupby("locus").cls.nunique(); pure=lc[lc==1].index
p=m[m.locus.isin(pure)&m.cls.isin(["GGT_only","ALTAST_only"])].copy()
print("pure genes by class", p.cls.value_counts().to_dict(), "loci", p.groupby("cls").locus.nunique().to_dict())
c=pd.read_csv(B/"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gwas_name","gene","ensembl","PP.H4.susie","PP.H3.susie","PP.H4.abf"])
t=pd.read_csv(B/"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t")
c=c.merge(t[["study_name","trait","tier"]],left_on="gwas_name",right_on="study_name",how="left")
c=c[c.tier.isin([1,2])&c.ensembl.isin(p.index)]
sym=c.drop_duplicates("ensembl").set_index("ensembl").gene
p["symbol"]=sym.reindex(p.index).values
# 1. how close are GGT-only genes to ALT/AST colocalization (absence-of-evidence labels)
def mx(tr):
    s=c[c.trait.isin(tr)].groupby("ensembl")["PP.H4.susie"].max()
    return s.reindex(p.index)
p["maxH4_altast"]=mx(["ALT","AST"]); p["maxH4_ggt"]=mx(["GGT"])
p["maxH4abf_altast"]=c[c.trait.isin(["ALT","AST"])].groupby("ensembl")["PP.H4.abf"].max().reindex(p.index)
p["maxH4abf_ggt"]=c[c.trait.isin(["GGT"])].groupby("ensembl")["PP.H4.abf"].max().reindex(p.index)
g=p[p.cls=="GGT_only"]; a=p[p.cls=="ALTAST_only"]
bins=[-0.01,0.1,0.25,0.5]
print("GGT_only genes: max SuSiE PP.H4 for ALT/AST (NaN = no SuSiE pair):", pd.cut(g.maxH4_altast,bins).value_counts(dropna=False).sort_index().to_dict())
print("ALTAST_only genes: max SuSiE PP.H4 for GGT:", pd.cut(a.maxH4_ggt,bins).value_counts(dropna=False).sort_index().to_dict())
print("GGT_only genes with ALT/AST ABF PP.H4>0.5:", int((g.maxH4abf_altast>0.5).sum()), "of", len(g))
print("ALTAST_only genes with GGT ABF PP.H4>0.5:", int((a.maxH4abf_ggt>0.5).sum()), "of", len(a))
# 2. ALT vs AST composition within ALTAST_only
h=c[(c["PP.H4.susie"]>0.5)]
tr=h.groupby("ensembl").trait.apply(lambda s: "+".join(sorted(set(s)&{"ALT","AST"})))
print("ALTAST_only genes by enzyme:", tr.reindex(a.index).value_counts().to_dict())
# 3. enzyme-encoding genes
enc=["GPT","GPT2","GOT1","GOT2","GGT1","GGT2","GGT5","GGT6","GGT7","GGTLC1"]
print("enzyme-encoding genes in pure classes:", p[p.symbol.isin(enc)][["symbol","cls","locus"]].to_dict("records"))
# 4. gene length / biotype from hg19 GENCODE table
gt=pd.read_csv(B/"GWAS/ldsc/program_heritability_20260830/gene_sets/gencode_v46lift37_genes_hg19.tsv",sep="\t")
print("gencode cols", list(gt.columns))
