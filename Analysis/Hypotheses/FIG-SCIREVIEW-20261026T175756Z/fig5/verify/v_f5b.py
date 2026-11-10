import pandas as pd, numpy as np
from scipy.stats import spearmanr
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
B=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5B.tsv",sep="\t")
print("n",len(B),"rho",round(spearmanr(B.bulk_logFC,B.protein_logFC)[0],4),"same sign",(np.sign(B.bulk_logFC)==np.sign(B.protein_logFC)).mean().round(4))
print(B.rna_call.value_counts().to_dict(), B["class"].value_counts().to_dict())
ps=B[B.protein_padj<0.05]; ss=np.sign(ps.bulk_logFC)==np.sign(ps.protein_logFC)
print("protein sig",len(ps),"same",ss.sum())
for lab,g in ps.groupby(ps.bulk_padj<0.05): print(" RNA padj<0.05" if lab else " RNA padj>=0.05",len(g),(np.sign(g.bulk_logFC)==np.sign(g.protein_logFC)).sum())
call=(ps.bulk_padj<0.05)&(ps.bulk_logFC.abs()>0.5)
for lab,g in ps.groupby(call): print(" 5C-style call" if lab else " no 5C-style call",len(g),(np.sign(g.bulk_logFC)==np.sign(g.protein_logFC)).sum())
