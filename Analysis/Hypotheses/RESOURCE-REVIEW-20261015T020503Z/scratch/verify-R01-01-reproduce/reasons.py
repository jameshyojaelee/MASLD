import pandas as pd, numpy as np
FM="GWAS/finemapping"
tiers=pd.read_csv(f"{FM}/config/gwas_trait_tier.tsv",sep="\t")
t12=set(tiers.loc[tiers.tier.isin([1,2]),"study_name"])
print("tier1/2 studies",len(t12))
m=pd.read_csv(f"{FM}/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
   usecols=["gwas_name","ensembl","method","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","PP.H4.susie"],low_memory=False)
m=m[m.gwas_name.isin(t12)].copy()
m["g"]=m.ensembl.astype(str).str.split(".").str[0]
print("tier12 rows",len(m)); print(m.method.value_counts(dropna=False).to_string())
C="GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
B="RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"
d=pd.read_csv(C); d["g"]=d.ensembl.astype(str).str.split(".").str[0]
deg=pd.read_csv(B); deg["g"]=deg.gene.astype(str).str.split(".").str[0]; deg=deg.drop_duplicates("g")
jt=deg.merge(d[["g","susie_state_t12","coloc_best_abf_pp4"]],on="g")
jt["is_deg"]=(jt.padj<0.05)&(jt.logFC.abs()>0.5)
ns=set(jt.g[jt.susie_state_t12=="no_susie_result"])
mm=m[m.g.isin(ns)].copy()
mm["gwas_sig"]=mm[["PP.H1.abf","PP.H3.abf","PP.H4.abf"]].sum(axis=1)
mm["eqtl_sig"]=mm[["PP.H2.abf","PP.H3.abf","PP.H4.abf"]].sum(axis=1)
per=mm.groupby("g").agg(any_fallback=("method",lambda s:(s=="abf_fallback").any()),
    all_abf_only=("method",lambda s:(s=="abf_only").all()),
    max_gwas=("gwas_sig","max"),max_eqtl=("eqtl_sig","max"),max_h4=("PP.H4.abf","max"))
print("\nno_susie_result genes in joint universe:",len(ns),"found in master:",len(per))
print("any abf_fallback row:",int(per.any_fallback.sum()),"; all rows abf_only (no eQTL SuSiE fit):",int(per.all_abf_only.sum()))
print("max ABF H1+H3+H4 >0.5 in some tier1/2 study:",int((per.max_gwas>0.5).sum()))
print("max ABF H2+H3+H4 >0.5 :",int((per.max_eqtl>0.5).sum()))
print("max ABF PP.H4 >0.5 :",int((per.max_h4>0.5).sum()))
x=jt.set_index("g").join(per,how="inner")
print("\nDEG rate by reason among no_susie_result:")
x["reason"]=np.where(x.all_abf_only,"abf_only_all","has_fallback")
print(x.groupby("reason").is_deg.agg(["sum","count","mean"]).to_string())
print("\nDEG rate by SuSiE state (joint):")
print(jt.groupby("susie_state_t12").is_deg.agg(["sum","count","mean"]).to_string())
print("\nmedian AveExpr by state:"); print(jt.groupby("susie_state_t12").AveExpr.median().to_string())
