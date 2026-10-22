import pandas as pd, numpy as np
from scipy import stats
C="GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
B="RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"
c=pd.read_csv(C); c["g"]=c["ensembl"].astype(str).str.split(".").str[0]
c=c.sort_values(["coloc_best_susie_pp4","coloc_best_abf_pp4"],ascending=False).drop_duplicates("g")
d=pd.read_csv(B); d["g"]=d["gene"].astype(str).str.split(".").str[0]; d=d.drop_duplicates("g")
print("coloc genes",len(c),"bulk genes",len(d))
print("coloc susie_state_t12 counts (all coloc genes):"); print(c["susie_state_t12"].value_counts(dropna=False))
jt=d.merge(c[["g","coloc_best_susie_pp4","coloc_best_abf_pp4","susie_state_t12"]],on="g")
jt["susie"]=pd.to_numeric(jt.coloc_best_susie_pp4,errors="coerce")
jt["deg"]=(jt.padj<0.05)&(jt.logFC.abs()>0.5)
jt["gen"]=jt.susie>0.5
print("joint",len(jt),"gen",jt.gen.sum(),"deg",jt.deg.sum(),"both",(jt.gen&jt.deg).sum())
print(pd.crosstab(jt.susie_state_t12.fillna("NA"),jt.susie.notna(),margins=True))
print("DEGs by state"); print(pd.crosstab(jt.susie_state_t12.fillna("NA"),jt.deg,margins=True))
def ov(sub,lab):
    n=len(sub);g=sub.gen;e=sub.deg;a=int((g&e).sum());b=int((g&~e).sum());cc=int((~g&e).sum());dd=n-a-b-cc
    orr,p=stats.fisher_exact([[a,b],[cc,dd]])
    exp=(a+b)*(a+cc)/n
    print(f"{lab}: n={n} gen={a+b} deg={a+cc} both={a} exp={exp:.4f} exp_frac_gen_row={exp/(a+b):.4f} exp_frac_deg_row={exp/(a+cc):.4f} obs_frac_gen={a/(a+b):.4f} obs_frac_deg={a/(a+cc):.4f} OR={orr:.4f} p={p:.4g} nonoverlap_exp={1-(a+cc)/n:.4f}")
ov(jt,"full interface")
fin=jt[jt.susie.notna()]
ov(fin,"finite SuSiE PP.H4")
st=jt[jt.susie_state_t12!="no_susie_result"]
ov(st,"state != no_susie_result")
