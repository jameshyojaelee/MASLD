import pandas as pd, numpy as np, random
pd.set_option("display.width",220); pd.set_option("display.max_columns",40)
C="RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/"
g=pd.read_parquet(C+"passport_gene_index.parquet"); e=pd.read_parquet(C+"passport_evidence_long.parquet")
deg=pd.read_csv("RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"); deg["ens"]=deg.gene.str.replace(r"\..*","",regex=True)
col=pd.read_csv("GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
dd=deg.set_index("ens"); cc=col.set_index("ensembl")
def cls(ens):
    bt=ens in dd.index; gt=ens in cc.index
    bcall=bool(bt and dd.at[ens,"padj"]<0.05 and abs(dd.at[ens,"logFC"])>0.5)
    sp=cc.at[ens,"coloc_best_susie_pp4"] if gt else np.nan; st=cc.at[ens,"susie_state_t12"] if gt else ""
    gcall=bool(gt and sp>0.5)
    if not(bt and gt): return "untested"
    if st=="untestable": return "untested"
    if gcall and bcall: return "concordant"
    if gcall: return "genetically_anchored"
    if bcall: return "established_state_associated"
    return "unresolved"
g["re"]=g.ensembl_id.map(cls)
print("full recompute agree",(g.re==g.primary_evidence_class).sum(),"/",len(g)); print(pd.crosstab(g.primary_evidence_class,g.re).to_string())
bad=g[g.re!=g.primary_evidence_class]; print(bad[["symbol","ensembl_id","primary_evidence_class","re"]].head(20).to_string())
print("susie pp4 >0.5:",(col.coloc_best_susie_pp4>0.5).sum(),"positive state:",(col.susie_state_t12=="positive").sum(),"min positive pp4",col.loc[col.susie_state_t12=="positive","coloc_best_susie_pp4"].min(),"max tested_le",col.loc[col.susie_state_t12=="tested_le_0.5","coloc_best_susie_pp4"].max())
print("untestable state rows:",col[col.susie_state_t12=="untestable"][["ensembl","coloc_best_susie_pp4","coloc_best_abf_pp4","susie_untestable_t12_studies"]].to_string())
random.seed(20260923); rows=[]
for c in ["genetically_anchored","concordant","established_state_associated","unresolved","untested"]:
    rows+=random.sample(sorted(g.loc[g.primary_evidence_class==c,"ensembl_id"]),4)
s=g[g.ensembl_id.isin(rows)].set_index("ensembl_id").loc[rows].reset_index()
s["logFC"]=s.ensembl_id.map(lambda x: dd.at[x,"logFC"] if x in dd.index else np.nan); s["padj"]=s.ensembl_id.map(lambda x: dd.at[x,"padj"] if x in dd.index else np.nan)
s["susie_pp4"]=s.ensembl_id.map(lambda x: cc.at[x,"coloc_best_susie_pp4"] if x in cc.index else np.nan); s["in_coloc"]=s.ensembl_id.isin(cc.index)
s["ok"]=s.re==s.primary_evidence_class
print(s[["symbol","ensembl_id","primary_evidence_class","re","logFC","padj","in_coloc","susie_pp4","ok"]].to_string()); print("spot agree",s.ok.sum(),"/",len(s))
s[["symbol","ensembl_id","primary_evidence_class","re","logFC","padj","in_coloc","susie_pp4","ok"]].to_csv("Analysis/Hypotheses/r8-review-20261011T150915Z/spot_check_20.tsv",sep="\t",index=False)
# missing genes
md=set(deg.ens)-set(g.ensembl_id); mc=set(col.ensembl)-set(g.ensembl_id)
print("deg genes absent from catalog",len(md),"of which DEG calls",int(((deg.padj<0.05)&(deg.logFC.abs()>0.5)&deg.ens.isin(md)).sum()),"| coloc genes absent",len(mc),"of which susie positive",int(((col.susie_state_t12=="positive")&col.ensembl.isin(mc)).sum()))
q=pd.read_csv(C+"preflight/gen_identity_adjudication_v1/gen_identity_quarantine.tsv",sep="\t"); print("quarantine",q.shape, q.reason_code.value_counts().to_dict()); 
qe=set(";".join(q[[c for c in q.columns if "ensembl" in c]].astype(str).agg(";".join,axis=1)).replace(",",";").split(";"))
print("absent deg genes explained by quarantined IDs",len(md&qe),"/",len(md),"; absent coloc explained",len(mc&qe),"/",len(mc))
w=e.pivot(index="ensembl_id",columns="evidence_domain",values=["call_state","estimate"]); w.columns=[a+"_"+b for a,b in w.columns]
ind=w[w.call_state_genetics=="indeterminate"]; est=ind.estimate_genetics.astype(float)
print("indeterminate genetics rows with PP.H4 >0.5:",(est>0.5).sum(),">=0.8:",(est>=0.8).sum(),"of",len(ind))
x=w.reset_index().merge(col,left_on="ensembl_id",right_on="ensembl"); hi=x[(x.call_state_genetics=="indeterminate")&(x.estimate_genetics.astype(float)>0.5)]
print("those rows: estimate equals ABF",np.isclose(hi.estimate_genetics.astype(float),hi.coloc_best_abf_pp4).sum(),"susie state",hi.susie_state_t12.value_counts().to_dict())
