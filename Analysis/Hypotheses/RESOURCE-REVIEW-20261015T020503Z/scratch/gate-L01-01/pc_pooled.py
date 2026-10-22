# Label-free precheck of positive control (b): ALL 436 COLOC genes (289 loci, classes ignored)
# hepatocyte z vs expression-decile-matched protein-coding background. Does not use GGT/ALTAST labels.
import pandas as pd, numpy as np
from pathlib import Path
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
S=B/"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-01"
D=B/"Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
srr2donor={}
for ds in ("GSE244832","GSE185477","GSE202379","GSE136103"):
    dp=pd.read_csv(B/f"data/{ds}/metadata/donor_pairing.csv",dtype=str)
    for _,r in dp.iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s and s!="nan": srr2donor[s.strip()]=f"{ds}_{r['donor_id']}"
LIN=["Hepatocytes","Fibroblasts","Cholangiocytes","Macrophages","T_cells","Endothelial_cells"]
prof={}
for ln in LIN:
    x=pd.read_csv(D/f"pseudobulk_{ln}.tsv.gz",sep="\t",index_col=0).groupby(level=0).sum()
    cols=pd.Series(x.columns).map(srr2donor).fillna(pd.Series(x.columns)).values
    x=x.T.groupby(cols).sum().T
    cpm=x.div(x.sum(axis=0).replace(0,np.nan),axis=1)*1e6
    prof[ln]=np.log2(cpm+1).mean(axis=1)
Pm=pd.DataFrame(prof)
Z=Pm.sub(Pm.mean(axis=1),axis=0).div(Pm.std(axis=1).replace(0,np.nan),axis=0)
me=Pm.mean(axis=1)
gt=pd.read_csv(B/"GWAS/ldsc/program_heritability_20260830/gene_sets/gencode_v46lift37_genes_hg19.tsv",sep="\t",header=None,names=["chr","start","end","strand","ens","sym","type"])
pcg=set(gt[gt.type=="protein_coding"].sym)
m=pd.read_csv(S/"class_check.tsv",sep="\t",index_col=0)
c=pd.read_csv(B/"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gene","ensembl"]).drop_duplicates()
sym=c.drop_duplicates("ensembl").set_index("ensembl").gene
m["symbol"]=sym.reindex(m.index).values
coloc=set(m.symbol)
print("COLOC genes in pseudobulk:",len(coloc & set(Z.index)),"of",len(coloc))
bg=[g for g in Z.index if g in pcg and g not in coloc and me[g]>0.05 and np.isfinite(Z.loc[g,"Hepatocytes"])]
print("background protein-coding genes expressed:",len(bg))
edges=np.quantile(me[bg],np.linspace(0,1,11)); edges[0]=-np.inf; edges[-1]=np.inf
dec=lambda v: np.clip(np.searchsorted(edges,v,side="right")-1,0,9)
bgd=pd.DataFrame({"z":Z.loc[bg,"Hepatocytes"],"dec":dec(me[bg].values)})
bgmean=bgd.groupby("dec").z.mean()
m=m[m.symbol.isin(Z.index)].copy()
m["zhep"]=Z.loc[m.symbol,"Hepatocytes"].values; m["me"]=me[m.symbol].values; m["dec"]=dec(m.me.values)
m["zhep_adj"]=m.zhep-bgmean.reindex(m.dec).values
L=m.groupby("locus").agg(adj=("zhep_adj","mean"))
obs=L.adj.mean()
print("COLOC decile distribution:",m.groupby("dec").size().reindex(range(10),fill_value=0).to_dict())
print(f"pooled 289-locus mean hepatocyte z minus decile-matched background: {obs:.4f} (z units); locus SD {L.adj.std():.4f}; n loci {len(L)}")
# null: random background genes, one per locus gene in same decile, locus-averaged
rng=np.random.default_rng(20261015); B_=2000; nul=np.empty(B_)
pools={d:bgd[bgd.dec==d].z.values for d in range(10)}
locs=m.groupby("locus").dec.apply(list)
for b in range(B_):
    nul[b]=np.mean([np.mean([rng.choice(pools[d]) for d in ds]) - np.mean([bgmean[d] for d in ds]) for ds in locs])
p=(1+np.sum(np.abs(nul)>=abs(obs)))/(B_+1)
print(f"null mean {nul.mean():.4f} sd {nul.std():.4f}; two-sided p {p:.4g} (B={B_}, seed 20261015, floor {1/(B_+1):.4g})")
# hepatocyte z of the 436 genes vs cholangiocyte z, label-free
m["zchol"]=Z.loc[m.symbol,"Cholangiocytes"].values
bgd["zc"]=Z.loc[bg,"Cholangiocytes"].values
m["zchol_adj"]=m.zchol-bgd.groupby("dec").zc.mean().reindex(m.dec).values
print("pooled locus mean cholangiocyte z minus matched bg:",round(m.groupby("locus").zchol_adj.mean().mean(),4))
