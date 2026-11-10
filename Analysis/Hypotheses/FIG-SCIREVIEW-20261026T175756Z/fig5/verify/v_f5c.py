import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
de=pd.read_csv(R+"Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv",sep="\t")
m=len(de); p=de.pvalue.values; o=np.argsort(p); q=np.minimum.accumulate((p[o]*m/np.arange(1,m+1))[::-1])[::-1]; Q=np.empty(m); Q[o]=np.minimum(q,1)
print("n proteins",m,"max|BH-padj|",np.abs(Q-de.padj).max())
print(de[de.gene.isin(["GNMT","MAT1A","CYP2C19"])][["gene","logFC","padj"]].to_string())
C=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5C.tsv",sep="\t")
print(C.columns.tolist()); print(C.head(3).to_string())
print("---- 5C post hoc trio")
t=C[C.gene.isin(["GNMT","MAT1A","CYP2C19"])&C.feature.isin(["NAS","Fibrosis"])][["gene","feature","n","rho","pvalue"]].copy()
p=t.pvalue.values; m=len(p); o=np.argsort(p); q=np.minimum.accumulate((p[o]*m/np.arange(1,m+1))[::-1])[::-1]; Q=np.empty(m); Q[o]=q; t["q_own_6"]=Q
print(t.to_string())
print("125 BH:",C.padj_bh125.lt(0.05).sum(),"filled",C.filled.sum(),"n tests",len(C))
print("n values:",C.groupby("gene").n.min().loc[lambda s:s<46].to_dict())
print("covariates:",C.covariates.unique())
