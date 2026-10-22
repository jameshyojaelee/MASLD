"""Recompute Fig 3G: 774-unit disease-state log2FC vs tier-1/2 SuSiE PP.H4 (COLOC r2 20260924T053244Z)."""
import pandas as pd, numpy as np
from scipy.stats import spearmanr
O="Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
deg=pd.read_csv("RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
deg['ens']=deg.gene.str.split('.').str[0]
print("tested",len(deg),"canonical",((deg.padj<0.05)&(deg.logFC.abs()>0.5)).sum(), "up",((deg.padj<0.05)&(deg.logFC>0.5)).sum())
t12=pd.read_csv("GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
print(t12.susie_state_t12.value_counts(dropna=False).to_dict())
print("finite susie pp4:",t12.coloc_best_susie_pp4.notna().sum(), " >0.5:",(t12.coloc_best_susie_pp4>0.5).sum())
unt=t12[t12.susie_state_t12=='untestable']; print("untestable pp4 values:",unt[['ensembl','coloc_best_susie_pp4','coloc_best_abf_pp4']].to_dict('records'))
j=deg.merge(t12[['ensembl','coloc_best_susie_pp4','susie_state_t12','susie_driving_gwas']],left_on='ens',right_on='ensembl',how='inner')
f=j[j.coloc_best_susie_pp4.notna()].copy()
print("RNA x finite PP4 pairs:",len(f), "dup ens:",f.ens.duplicated().sum())
print("states among finite:",f.susie_state_t12.value_counts().to_dict())
f['tx']=(f.padj<0.05)&(f.logFC.abs()>0.5); f['gen']=f.coloc_best_susie_pp4>0.5
f['cls']=np.select([f.tx&f.gen,f.tx&~f.gen,~f.tx&f.gen],['Both','Transcriptomics','Genetics'],'Neither')
print(f.cls.value_counts().to_dict())
rho,p=spearmanr(f.coloc_best_susie_pp4,f.logFC.abs()); print("spearman |lfc| vs pp4",round(rho,4),p)
print("gen >0.5 up/down:",((f.gen)&(f.logFC>0)).sum(),((f.gen)&(f.logFC<0)).sum())
print("min PP4 in gen:",f[f.gen].coloc_best_susie_pp4.min(), "max PP4 neither:", f[~f.gen].coloc_best_susie_pp4.max())
print("lfc range:",f.logFC.min(),f.logFC.max(), "pp4 range", f.coloc_best_susie_pp4.min(), f.coloc_best_susie_pp4.max())
# primary SuSiE genes not in RNA
prim=t12[t12.coloc_best_susie_pp4>0.5]; print("primary 436 in RNA:",prim.ensembl.isin(deg.ens).sum())
f[['gene','symbol','logFC','padj','coloc_best_susie_pp4','susie_state_t12','susie_driving_gwas','cls']].to_csv(O+"fig3g_recomputed.tsv",sep="\t",index=False)
# trait composition of genetics class
print("genetics-class driving gwas:", f[f.gen].susie_driving_gwas.str.extract(r'_(ALT|AST|GGT|NAFLD|MASLD|PDFF|ALP|[A-Za-z]+)_')[0].value_counts().to_dict())
