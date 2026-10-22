import pandas as pd, glob, numpy as np
B="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
de=pd.read_csv(B+"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
de["sid"]=de.gene.str.replace(r"\.\d+(_PAR_Y)?$","",regex=True)
print("DE rows",len(de),"unique stable",de.sid.nunique(),"PAR_Y rows",de.gene.str.contains("PAR_Y").sum())
de["canon"]=(de.padj<0.05)&(de.logFC.abs()>0.5)
print("canonical",de.canon.sum(),"up",((de.canon)&(de.logFC>0)).sum(),"down",((de.canon)&(de.logFC<0)).sum())
de=de[~de.gene.str.contains("PAR_Y")]
for t in ["ukbb_alt","ghodsian","pdff","ukbb_ast","ukbb_ggt"]:
    fs=sorted(glob.glob(B+f"RNA-seq/results/causal_inference/otters_broadaway/{t}/chr*/P0.05.txt"))
    tw=pd.concat([pd.read_csv(f,sep="\t") for f in fs])
    tw=tw.dropna(subset=["FUSION_Z"])
    dup=tw.TargetID.duplicated().sum()
    m=tw.merge(de,left_on="TargetID",right_on="sid")
    print(t,"twas genes",tw.TargetID.nunique(),"dupTarget",dup,"joined",len(m),"canonical in join",int(m.canon.sum()),
          "|z|>3",int((m.FUSION_Z.abs()>3).sum()),"|z|>3&canon",int(((m.FUSION_Z.abs()>3)&m.canon).sum()))
