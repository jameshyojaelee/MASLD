# Blinded counts for the H60 prespec. Reads TWAS z, positions, DE gene membership and AveExpr.
# Never pairs TWAS z with DE t or logFC. The observed statistic is not computed.
import glob, numpy as np, pandas as pd
B="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
OT=B+"RNA-seq/results/causal_inference/otters_broadaway/"
de=pd.read_csv(B+"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv", usecols=["gene","padj","logFC","AveExpr"])
de["sid"]=de.gene.str.replace(r"\.\d+$","",regex=True)
de["canon"]=(de.padj<0.05)&(de.logFC.abs()>0.5)
de=de[["sid","canon","AveExpr"]]
col=pd.read_csv(B+"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
gen=set(col.loc[col.coloc_best_susie_pp4>0.5,"ensembl"].astype(str).str.replace(r"\.\d+$","",regex=True))
print("SuSiE PP.H4>0.5 genes", len(gen))
def load(t):
    d=pd.concat([pd.read_csv(f,sep="\t") for f in sorted(glob.glob(OT+f"{t}/chr*/P0.05.txt"))]).dropna(subset=["FUSION_Z"])
    d["mid"]=(d.GeneStart+d.GeneEnd)/2
    return d.drop_duplicates("TargetID")
for t in ["ukbb_alt","ghodsian","ukbb_ast","ukbb_ggt"]:
    m=load(t).merge(de,left_on="TargetID",right_on="sid").sort_values(["CHROM","mid"]).reset_index(drop=True)
    sz=m.groupby("CHROM").size()
    ing=m.TargetID.isin(gen)
    keep=np.ones(len(m),bool)
    for ch,g in m[ing].groupby("CHROM"):
        idx=np.where((m.CHROM==ch).values)[0]
        d=np.min(np.abs(m.loc[idx,"mid"].values[:,None]-g.mid.values[None,:]),axis=1)
        keep[idx[d<=1e6]]=False
    szk=m[keep].groupby("CHROM").size()
    mhc=(m.CHROM==6)&(m.mid>=28477797)&(m.mid<=33448354)
    s=m[m.FUSION_Z.abs()>3].copy()
    newblk=(s.CHROM.diff().fillna(1)!=0)|(s.mid.diff().fillna(1e9)>1e6)
    nblk=int(newblk.sum())
    blk_sizes=newblk.cumsum().value_counts()
    print(f"{t}: joined {len(m)} chroms {len(sz)} min {sz.min()} max {sz.max()} | coloc in join {int(ing.sum())} | after coloc+1Mb exclusion {int(keep.sum())} chroms {len(szk)} min chr n {szk.min()} | MHC (hg19 chr6:28,477,797-33,448,354) genes {int(mhc.sum())} | |z|>3 {len(s)} in {nblk} locus blocks (1 Mb gap), largest block {blk_sizes.max()} median {int(blk_sizes.median())} | |z|>3 & canonical {int((s.canon).sum())} | |z|>3 & colocalized {int(s.TargetID.isin(gen).sum())}")
    dec=pd.qcut(m.AveExpr.rank(method='first'),10,labels=False)
    print(f"   decile sizes {dec.value_counts().sort_index().tolist()}")
