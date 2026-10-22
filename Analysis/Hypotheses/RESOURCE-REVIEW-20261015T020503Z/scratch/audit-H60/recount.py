# Audit recount for H60. Never pairs TWAS z with DE t (no joint statistic of z and t).
import glob, numpy as np, pandas as pd
B="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
OT=B+"RNA-seq/results/causal_inference/otters_broadaway/"
de=pd.read_csv(B+"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
de["sid"]=de.gene.str.replace(r"\.\d+$","",regex=True)
de["canon"]=(de.padj<0.05)&(de.logFC.abs()>0.5)
print("DE rows",len(de),"canonical",int(de.canon.sum()),"dup sid",int(de.sid.duplicated().sum()))
def load(t):
    fs=sorted(glob.glob(OT+f"{t}/chr*/P0.05.txt"))
    d=pd.concat([pd.read_csv(f,sep="\t") for f in fs])
    return d,len(fs)
for t in ["ukbb_alt","ghodsian"]:
    d,nf=load(t)
    print(t,"files",nf,"rows",len(d),"NaN z",int(d.FUSION_Z.isna().sum()),"dup",int(d.TargetID.duplicated().sum()))
    m=d.merge(de[["sid","AveExpr","canon"]],left_on="TargetID",right_on="sid")
    s=m.groupby("CHROM").size()
    print(" join",len(m),"chroms",len(s),"min",s.min(),"max",s.max(),"canon",int(m.canon.sum()),"|z|>3",int((m.FUSION_Z.abs()>3).sum()),"|z|>3&canon",int(((m.FUSION_Z.abs()>3)&m.canon).sum()))
    dec=pd.qcut(m.AveExpr.rank(method="first"),10,labels=False)
    print(" decile sizes",sorted(set(dec.value_counts().tolist())))
    if t=="ukbb_alt":
        for g,nm in [("ENSG00000108821","COL1A1"),("ENSG00000164692","COL1A2"),("ENSG00000167701","GPT")]:
            print(" in ALT join",nm,g,bool((m.TargetID==g).any()),"in DE",bool((de.sid==g).any()))
        # canonical share overall vs among |z|>3 (magnitude only, no sign)
        print(" canonical share all %.4f, among |z|>3 %.4f"%(m.canon.mean(), m.loc[m.FUSION_Z.abs()>3,"canon"].mean()))
        print(" mean ALT z in join %.4f, frac z>0 %.4f"%(m.FUSION_Z.mean(), (m.FUSION_Z>0).mean()))
c=pd.read_csv(OT+"concordance/otters_vs_spredixcan_concordance.csv")
gh,_=load("ghodsian")
pc=c.merge(gh[["TargetID","CHROM"]],left_on="ensembl_id",right_on="TargetID")
pc=pc[np.isfinite(pc.spx_z)]
s=pc.groupby("CHROM").size()
print("PC1 n",len(pc),"chroms",len(s),"min chr n",s.min(),"(chr",s.idxmin(),")")
