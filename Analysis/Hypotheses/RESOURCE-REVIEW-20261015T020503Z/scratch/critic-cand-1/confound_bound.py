# Bound on the expression-level route: rho(DE t, AveExpr) over the joined genes. Does not use TWAS z jointly with t.
import pandas as pd, glob
from scipy.stats import spearmanr
B="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
de=pd.read_csv(B+"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
de["sid"]=de.gene.str.replace(r"\.\d+$","",regex=True)
ids=pd.concat([pd.read_csv(f,sep="\t",usecols=["TargetID"]) for f in glob.glob(B+"RNA-seq/results/causal_inference/otters_broadaway/ukbb_alt/chr*/P0.05.txt")]).TargetID
m=de[de.sid.isin(set(ids))]
print("n",len(m),"rho(t,AveExpr) %.4f"%spearmanr(m.t,m.AveExpr)[0])
for g in ["GSE126848","GSE130970","GSE135251","GSE162694","GSE213621"]:
    f=glob.glob(B+f"RNA-seq/results/manuscript_release/candidates/bulk-loco-{g}-*/deg_results.csv")[0]
    l=pd.read_csv(f,nrows=2); print(g,"LOCO columns ok" if "t" in l.columns else list(l.columns))
