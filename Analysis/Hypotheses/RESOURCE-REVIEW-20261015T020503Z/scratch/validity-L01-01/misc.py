import pandas as pd, numpy as np
from pathlib import Path
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
D=B/"Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
meta=pd.read_csv(D/"sample_celltype_metadata.tsv",sep="\t")
srr2donor={}
for ds in ("GSE244832","GSE185477","GSE202379","GSE136103"):
    dp=pd.read_csv(B/f"data/{ds}/metadata/donor_pairing.csv",dtype=str)
    for _,r in dp.iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s and s!="nan": srr2donor[s.strip()]=f"{ds}_{r['donor_id']}"
L=["Hepatocytes","Cholangiocytes","Fibroblasts","Endothelial cells","Macrophages","T cells"]
meta["donor"]=meta["sample"].map(srr2donor).fillna(meta["sample"])
x=meta[meta.cell_type.isin(L)]
k=x.groupby("donor").cell_type.nunique()
print("donors with all six lineages:",int((k==6).sum()),"| with Hep+Chol:",int(x[x.cell_type.isin(["Hepatocytes","Cholangiocytes"])].groupby("donor").cell_type.nunique().eq(2).sum()))
# cholangiocyte cells per donor
cc=x[x.cell_type=="Cholangiocytes"].groupby("donor").n_cells.sum()
print("cholangiocyte cells per donor: min",cc.min(),"median",cc.median(),"n<50:",int((cc<50).sum()),"of",len(cc))
# H10 locus spacing
lo=pd.read_csv(B/"Analysis/Hypotheses/H10-20261007T210859Z/loci.tsv",sep="\t").sort_values(["chr","start"])
gap=lo.groupby("chr").apply(lambda g: (g.start.values[1:]-g.end.values[:-1])).explode().astype(float)
print("H10 loci",len(lo),"min gap between adjacent loci (bp):",int(gap.min()),"n gaps <500kb:",int((gap<5e5).sum()), "n gaps <1Mb:",int((gap<1e6).sum()))
# canalicular / biliary-transport genes among the 436 and their class
m=pd.read_csv(B/"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-01/class_check.tsv",sep="\t",index_col=0)
c=pd.read_csv(B/"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gene","ensembl"]).drop_duplicates("ensembl").set_index("ensembl").gene
m["symbol"]=c.reindex(m.index).values
can=["ABCB4","ABCB11","ATP8B1","ABCC2","ABCG5","ABCG8","SLC10A1","SLCO1B1","SLCO1B3","NR1H4","CYP7A1","HNF1B","GGT1","ALPL","SERPINA1","PNPLA3","TM6SF2","GPT","GOT1"]
print(m[m.symbol.isin(can)][["symbol","cls","locus"]].to_string())
