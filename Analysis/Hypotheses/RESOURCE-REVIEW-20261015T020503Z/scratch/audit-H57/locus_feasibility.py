# Audit H57: can genetic genes be grouped into GWAS signal clusters without the decimal locus column?
# Reads genetic-side columns only; no DEG status is read.
import pandas as pd, numpy as np
C="GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
A="GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
c=pd.read_csv(C); c["g"]=c.ensembl.astype(str).str.split(".").str[0]
pos=c[c.susie_state_t12=="positive"].copy()
print("positive genes in tier12 table",len(pos))
a=pd.read_csv(A,usecols=["gwas_name","ensembl","chr","PP.H4.susie","method","top_snp","eligibility_source"],dtype={"chr":str,"top_snp":str})
a["g"]=a.ensembl.astype(str).str.split(".").str[0]
print("method counts:",a.method.value_counts(dropna=False).to_dict())
print("eligibility_source counts:",a.eligibility_source.value_counts(dropna=False).to_dict())
m=pos.merge(a,left_on=["g","susie_driving_gwas"],right_on=["g","gwas_name"],how="left")
print("rows after merge",len(m),"missing top_snp",int(m.top_snp.isna().sum()),"dup genes",int(m.g.duplicated().sum()))
m=m.drop_duplicates("g")
print("PP.H4.susie agrees with coloc_best_susie_pp4:",np.allclose(m["PP.H4.susie"].astype(float),m.coloc_best_susie_pp4.astype(float)))
print("top_snp examples",m.top_snp.head(5).tolist())
sp=m.top_snp.str.split(":",expand=True); m["c"]=sp[0]; m["p"]=pd.to_numeric(sp[1],errors="coerce")
print("chr mismatch top_snp vs chr col",int((m.c.astype(str)!=m.chr.astype(str)).sum()),"nonnumeric pos",int(m.p.isna().sum()))
for w in (250_000,500_000,1_000_000):
    k=0
    for ch,grp in m.sort_values(["c","p"]).groupby("c"):
        pp=grp.p.to_numpy(); k+=1+int((np.diff(pp)>w).sum())
    print(f"single-linkage clusters at {w} bp: {k}")
print("distinct (driving gwas, top_snp):",m[["susie_driving_gwas","top_snp"]].drop_duplicates().shape[0],"distinct top_snp:",m.top_snp.nunique())
