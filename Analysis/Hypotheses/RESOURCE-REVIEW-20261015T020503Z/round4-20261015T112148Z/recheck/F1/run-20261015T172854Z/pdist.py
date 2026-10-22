import sys, numpy as np, pandas as pd
from scipy.stats import hypergeom
sys.argv=[sys.argv[0], sys.argv[1]]
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
c=pd.read_csv(R+"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
b=pd.read_csv(R+"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
c["g"]=c.ensembl.str.split(".").str[0]; b["g"]=b.gene.str.split(".").str[0]
b["deg"]=(b.padj<0.05)&(b.logFC.abs()>0.5)
u=b.merge(c,on="g")
nr=u[u.susie_state_t12=="no_susie_result"]
print("no_susie_result genes",len(nr),"finite ABF",int(nr.coloc_best_abf_pp4.notna().sum()),"finite SuSiE",int(nr.coloc_best_susie_pp4.notna().sum()))
def pd_(df,gcol,k=10):
    df=df.sort_values("g").reset_index(drop=True)
    df["bin"]=pd.qcut(df.AveExpr.rank(method="first"),k,labels=False)
    pmf=np.array([1.0]);E=0
    for _,s in df.groupby("bin"):
        n=len(s);G=int(s[gcol].sum());Dd=int(s.deg.sum());E+=G*Dd/n
        pmf=np.convolve(pmf,hypergeom.pmf(np.arange(min(G,Dd)+1),n,G,Dd))
    O=int((df[gcol]&df.deg).sum());x=np.arange(len(pmf))
    return O,E,float(pmf[np.abs(x-E)>=abs(O-E)-1e-9].sum())
u["gen"]=u.coloc_best_susie_pp4>=0.5; u["abf"]=u.coloc_best_abf_pp4>=0.5
ev=u[u.susie_state_t12.isin(["positive","tested_le_0.5"])]
print("T1 distance-p",pd_(ev,"gen")); print("D8 distance-p",pd_(u,"abf"))
