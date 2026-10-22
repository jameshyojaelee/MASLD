import pandas as pd, numpy as np
from scipy import stats
C="GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
B="RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"
c=pd.read_csv(C); c["g"]=c["ensembl"].astype(str).str.split(".").str[0]
c=c.sort_values(["coloc_best_susie_pp4","coloc_best_abf_pp4"],ascending=False).drop_duplicates("g")
d=pd.read_csv(B); d["g"]=d["gene"].astype(str).str.split(".").str[0]; d=d.drop_duplicates("g")
jt=d.merge(c[["g","coloc_best_susie_pp4","coloc_best_abf_pp4","susie_state_t12"]],on="g")
jt["deg"]=(jt.padj<0.05)&(jt.logFC.abs()>0.5)
print("DEG rate by SuSiE evaluability:")
print(jt.groupby(jt.susie_state_t12.eq("no_susie_result")).deg.agg(["sum","count","mean"]))
# matched controls: which SuSiE state?
m=pd.read_csv("RNA-seq/results/genetics_state_interface/20261009T142027Z_coloc_r2_corrected_controls_r8/matched_pairs_by_arm.tsv",sep="\t")
m=m[(m.bulk_arm=="fragment_23370")&(m.genetic_definition=="susie_pp4_0.5")]
st=c.set_index("g").susie_state_t12
print("matched pairs",len(m)); print("control gene SuSiE state:"); print(m.matched_control_gene.map(st).value_counts(dropna=False))
# ABF-tested check for SuSiE-or-ABF
jt["abf"]=pd.to_numeric(jt.coloc_best_abf_pp4,errors="coerce")
print("finite ABF in interface:",jt.abf.notna().sum(),"of",len(jt))
