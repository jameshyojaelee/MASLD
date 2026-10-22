import pandas as pd, numpy as np
from pathlib import Path
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
S=B/"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-01"
loci=pd.read_csv(B/"Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv",sep="\t")
gt=pd.read_csv(B/"GWAS/ldsc/program_heritability_20260830/gene_sets/gencode_v46lift37_genes_hg19.tsv",sep="\t",header=None,names=["chr","start","end","strand","ens","sym","type"])
gt=gt[gt.type=="protein_coding"].copy(); gt["chr"]=gt.chr.str.replace("chr",""); gt["tss"]=np.where(gt.strand=="+",gt.start,gt.end)
gt["ens"]=gt.ens.str.split(".").str[0]
m=pd.read_csv(S/"class_check.tsv",sep="\t",index_col=0)
idx=pd.read_csv(S/"index_all_genes.tsv",sep="\t",index_col=0)
lc=m.groupby("locus").cls.nunique(); pure=lc[lc==1].index
pl=m[m.locus.isin(pure)&m.cls.isin(["GGT_only","ALTAST_only"])].groupby("locus").cls.first()
coloc=set(m.index)
for flank in (100_000,250_000,500_000):
    rows=[]
    for _,r in loci[loci.locus_id.isin(pl.index)].iterrows():
        g=gt[(gt.chr==str(r.chr))&(gt.tss>=r.start-flank)&(gt.tss<=r.end+flank)&(~gt.ens.isin(coloc))]
        g=g[g.sym.isin(idx.index)]
        rows.append((r.locus_id,pl[r.locus_id],len(g)))
    d=pd.DataFrame(rows,columns=["locus","cls","n_noncoloc_pc"])
    print(f"flank {flank//1000} kb: loci with >=1 non-COLOC protein-coding gene in pseudobulk by class:",d[d.n_noncoloc_pc>0].groupby("cls").size().to_dict(),"median genes/locus",d.n_noncoloc_pc.median())
# GGT1 profile and membership (label-free)
print("GGT1 in COLOC 436:", "ENSG00000100031" in coloc)
