import pandas as pd, numpy as np
from scipy import stats
C="GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
B="RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"
d=pd.read_csv(C); d["g"]=d["ensembl"].astype(str).str.split(".").str[0]
print("coloc rows",len(d),"unique g",d.g.nunique())
print(d["susie_state_t12"].value_counts(dropna=False).to_string())
d=d.sort_values(["coloc_best_susie_pp4","coloc_best_abf_pp4"],ascending=False).drop_duplicates("g",keep="first")
deg=pd.read_csv(B); deg["g"]=deg["gene"].astype(str).str.split(".").str[0]; deg=deg.drop_duplicates("g",keep="first")
print("bulk genes",len(deg))
jt=deg.merge(d[["g","coloc_best_susie_pp4","coloc_best_abf_pp4","susie_state_t12"]],on="g",how="inner")
jt["susie"]=pd.to_numeric(jt.coloc_best_susie_pp4,errors="coerce")
jt["is_deg"]=(jt.padj<0.05)&(jt.logFC.abs()>0.5)
gen=(jt.susie>=0.5).fillna(False)
print("joint",len(jt),"genetic(>=0.5)",gen.sum(),"genetic(>0.5)",(jt.susie>0.5).sum(),"deg",jt.is_deg.sum(),"overlap",(gen&jt.is_deg).sum())
print("expected",gen.sum()*jt.is_deg.sum()/len(jt))
print("\nstate x finite susie within joint:")
print(pd.crosstab(jt.susie_state_t12.fillna("NA"),jt.susie.notna()).to_string())
print("\nstate x is_deg:")
print(pd.crosstab(jt.susie_state_t12.fillna("NA"),jt.is_deg).to_string())
print("\nstate x genetic:")
print(pd.crosstab(jt.susie_state_t12.fillna("NA"),gen).to_string())
# restricted universe: finite SuSiE PP.H4
f=jt[jt.susie.notna()]
g2=f.susie>=0.5
a=(g2&f.is_deg).sum(); b=(g2&~f.is_deg).sum(); c=(~g2&f.is_deg).sum(); dd=(~g2&~f.is_deg).sum()
orr,p=stats.fisher_exact([[a,b],[c,dd]])
print(f"\nFinite-SuSiE universe n={len(f)} genetic={a+b} deg={a+c} overlap={a} expected={(a+b)*(a+c)/len(f):.4f} OR={orr:.4f} p={p:.4g}")
print(f"  expected frac of genetic row={(a+c)/len(f):.4f}; expected frac of DEG row={(a+b)/len(f):.4f}; observed {a/(a+b):.4f} / {a/(a+c):.4f}")
# full universe
a=(gen&jt.is_deg).sum(); b=(gen&~jt.is_deg).sum(); c=(~gen&jt.is_deg).sum(); dd=(~gen&~jt.is_deg).sum()
orr,p=stats.fisher_exact([[a,b],[c,dd]])
print(f"Full universe n={len(jt)} OR={orr:.4f} p={p:.4g}; exp frac genetic row {(a+c)/len(jt):.4f}; DEG row {(a+b)/len(jt):.4f}")
# finite SuSiE OR untestable (exclude only no_susie_result)
h=jt[jt.susie_state_t12!="no_susie_result"]
g3=(h.susie>=0.5).fillna(False)
a=(g3&h.is_deg).sum(); b=(g3&~h.is_deg).sum(); c=(~g3&h.is_deg).sum(); dd=(~g3&~h.is_deg).sum()
orr,p=stats.fisher_exact([[a,b],[c,dd]])
print(f"Excl no_susie_result n={len(h)} genetic={a+b} deg={a+c} overlap={a} expected={(a+b)*(a+c)/len(h):.4f} OR={orr:.4f} p={p:.4g}")
