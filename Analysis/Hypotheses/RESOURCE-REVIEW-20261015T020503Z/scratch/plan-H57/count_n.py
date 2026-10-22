# H57 prespec: count set sizes (margins only). Does not compute overlap expectations,
# per-bin margins, or any test statistic.
import pandas as pd, numpy as np
C="GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
B="RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
P="figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/"
c=pd.read_csv(C); c["g"]=c.ensembl.astype(str).str.split(".").str[0]
print("coloc rows",len(c),"unique g",c.g.nunique())
c=c.sort_values(["coloc_best_susie_pp4","coloc_best_abf_pp4"],ascending=False).drop_duplicates("g")
d=pd.read_csv(B); d["g"]=d.gene.astype(str).str.split(".").str[0]
print("bulk rows",len(d),"unique g",d.g.nunique(),"PAR_Y-like ids",d.gene.astype(str).str.contains("PAR").sum())
d=d.drop_duplicates("g")
jt=d.merge(c[["g","coloc_best_susie_pp4","susie_state_t12"]],on="g",how="inner")
jt["susie"]=pd.to_numeric(jt.coloc_best_susie_pp4,errors="coerce")
jt["deg"]=(jt.padj<0.05)&(jt.logFC.abs()>0.5)
print("U_full",len(jt),"G(>=0.5)",int((jt.susie>=0.5).sum()),"G(>0.5)",int((jt.susie>0.5).sum()),"D",int(jt.deg.sum()))
print("AveExpr NA",int(jt.AveExpr.isna().sum()),"SE NA",int(jt.SE.isna().sum()))
print(pd.crosstab(jt.susie_state_t12.fillna("NA"),jt.susie.notna(),margins=True).to_string())
print("D by state:",jt.groupby(jt.susie_state_t12.fillna("NA")).deg.sum().to_dict())
ue=jt[jt.susie_state_t12.isin(["positive","tested_le_0.5"])]
print("U_eval",len(ue),"G",int((ue.susie>=0.5).sum()),"D",int(ue.deg.sum()),"susie NA",int(ue.susie.isna().sum()))
ue7=jt[jt.susie.notna()]
print("U_eval7",len(ue7),"G",int((ue7.susie>=0.5).sum()),"D",int(ue7.deg.sum()))
pos=jt[jt.susie_state_t12=="positive"]
print("positive state min PP.H4",pos.susie.min(),"; tested_le_0.5 max PP.H4",jt[jt.susie_state_t12=="tested_le_0.5"].susie.max())
print("untestable PP.H4 range",jt[jt.susie_state_t12=="untestable"].susie.describe()[["min","max"]].to_dict())
# AveExpr ranges by universe (distribution only)
for nm,u in [("U_full",jt),("U_eval",ue)]:
    q=u.AveExpr.quantile([0,.1,.5,.9,1]).round(3).to_dict(); print(nm,"AveExpr quantiles",q)
# control-set sizes (membership only, no DEG status read)
pi=pd.read_parquet(P+"passport_program_index.parquet"); m=pd.read_parquet(P+"passport_program_membership.parquet")
m["g"]=m.ensembl_id.astype(str).str.split(".").str[0]
fib=pi[(pi.cell_type=="fibroblasts")&pi.program_label.str.contains("ECM|collagen",regex=True)].program_uid
s=set(m.g[m.program_uid.isin(fib)])
print("fibroblast ECM/collagen programs",len(fib),"union genes",len(s),"in U_full",int(jt.g.isin(s).sum()),"in U_eval",int(ue.g.isin(s).sum()))
