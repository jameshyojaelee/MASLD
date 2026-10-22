import pandas as pd, numpy as np
A="figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_all_gene_results.tsv"
B="figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/stage_all_gene_results.tsv"
a=pd.read_csv(A,sep="\t",usecols=["gene_id_base","gene_id_versioned","logFC","P.Value","FDR","contrast","axis","n_samples","bh_family_size"])
b=pd.read_csv(B,sep="\t",usecols=["gene_id_base","logFC","P.Value","FDR","contrast","axis","n_samples"])
print("A contrasts", a.groupby(["axis","contrast"]).size().to_dict())
print("B contrasts", b.groupby(["axis","contrast"]).size().to_dict())
for t,name in ((a,"A"),(b,"B")):
    print(name, t[t.FDR<0.05].groupby(["axis","contrast"]).size().to_dict())
m=a.merge(b,on=["gene_id_base","contrast","axis"],suffixes=("_a","_b"))
print("max |dlogFC| A vs B", (m.logFC_a-m.logFC_b).abs().max(), len(m))
g=pd.read_csv("GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
pos=set(g.loc[g.coloc_best_susie_pp4>0.5,"ensembl"]); print("positive",len(pos))
tr=["F0_to_F1","F1_to_F2","F2_to_F3","F3_to_F4"]
s=a[a.contrast.isin(tr)&a.gene_id_base.isin(pos)]
sel=set(s.loc[s.FDR<0.05,"gene_id_base"]); print("3F genes",len(sel))
ss=s[s.gene_id_base.isin(sel)]; print("effects",len(ss),"marks",(ss.FDR<0.05).sum(), ss[ss.FDR<0.05].groupby("contrast").size().to_dict())
src=pd.read_csv("figures/main/fig3_bulk_transcriptomics/source_tables/fig3f_stage_genetics_matrix.tsv",sep="\t")
print("installed source rows",len(src), src.gene_id_base.nunique(), (src.FDR<0.05).sum(), set(src.gene_id_base)==sel)
