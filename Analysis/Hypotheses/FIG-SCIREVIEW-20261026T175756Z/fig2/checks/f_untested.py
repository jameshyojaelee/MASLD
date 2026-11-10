# Fig 2F: SuSiE test status of each drawn cell
import pandas as pd, numpy as np
x=pd.read_csv("t12_coloc_r2.tsv.gz",sep="\t",low_memory=False)
x["key"]=x.ensembl.fillna(x.gene); x["s"]=pd.to_numeric(x["PP.H4.susie"],errors="coerce"); x["a"]=pd.to_numeric(x["PP.H4.abf"],errors="coerce")
genes=["GGT1","PANX1","EPHA2","C2orf16","ACTG1","SHROOM3","MLIP","MSL2","EFHD1","GOT2","GLDC","HSCB","CHEK2"]
y=x[x.gene.isin(genes)]
c=y.groupby(["gene","reg_anc"]).agg(abf=("a","max"),s=("s","max"),s_tested=("s",lambda v: v.notna().sum()),n=("s","size")).reset_index()
c["pp"]=np.fmax(c.abf.fillna(-1),c.s.fillna(-1))
c["state"]=np.where(c.s>0.5,"dark: SuSiE>0.5",np.where(c.pp>0.5,"light: ABF>0.5 only","gray: neither"))
c["susie"]=np.where(c.s_tested>0,"SuSiE tested","SuSiE untested")
print(pd.crosstab(c.state,c.susie))
print("cells:",len(c))
