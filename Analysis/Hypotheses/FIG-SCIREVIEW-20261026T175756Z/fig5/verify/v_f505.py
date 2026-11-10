import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
E=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5E.tsv",sep="\t")
T=E[E.testable==True].copy(); m=len(T); print("testable",m)
# own BH
p=T.pvalue.values; o=np.argsort(p); r=np.arange(1,m+1)
q=np.minimum.accumulate((p[o]*m/r)[::-1])[::-1]; Q=np.empty(m); Q[o]=np.minimum(q,1)
T["q_own"]=Q; print("max|q_own-padj|",np.abs(T.q_own-T.padj).max())
print("q<0.05:",(T.q_own<0.05).sum(),"robust:",(T.robust==True).sum(),"neg/pos robust:",((T.robust==True)&(T.rho<0)).sum(),((T.robust==True)&(T.rho>0)).sum())
Ts=T.sort_values("pvalue").reset_index(drop=True)
Ts["bh_bar"]=0.05*(Ts.index+1)/m
print(Ts[["cell_type","module","program_name","pvalue","padj","bh_bar","robust"]].iloc[18:27].to_string())
# Monte Carlo SE of p at p=0.0121 with 9999 perms
pp=0.0121; print("MC SE",np.sqrt(pp*(1-pp)/9999))
# Count k: p=(k+1)/(9999+1)
print("k for p=0.0121:",0.0121*10000-1)
# Text/legend mentions of boundary?
