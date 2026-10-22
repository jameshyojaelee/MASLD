import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats
B=Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
S=B/"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-01"
D=B/"Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
m=pd.read_csv(S/"class_check.tsv",sep="\t",index_col=0)
# --- pure loci, primary and UKBB-only
def pure(col):
    g=m.dropna(subset=[col]); lc=g.groupby("locus")[col].nunique(); p=lc[lc==1].index
    return g[g.locus.isin(p)].groupby(col).locus.nunique(), g[g.locus.isin(p)]
pc,pg=pure("cls"); print("PRIMARY pure loci by class:\n",pc.to_string())
print("primary pure-loci genes by class:\n",pg.groupby("cls").size().to_string())
m2=m[m.cls_ukbb!="no_UKBB_enzyme_hit"].copy()
lc=m2.groupby("locus").cls_ukbb.nunique(); p=lc[lc==1].index
print("UKBB-only sensitivity pure loci (genes w/o UKBB enzyme hit dropped):\n",m2[m2.locus.isin(p)].groupby("cls_ukbb").locus.nunique().to_string())
# --- symbol mapping
c=pd.read_csv(B/"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gene","ensembl"]).drop_duplicates()
sym=c.set_index("ensembl").gene
sym=sym[~sym.index.duplicated()]
m["symbol"]=sym.reindex(m.index).values
# --- index exactly as label_vs_content.py
srr2donor={}
for ds in ("GSE244832","GSE185477","GSE202379","GSE136103"):
    dp=pd.read_csv(B/f"data/{ds}/metadata/donor_pairing.csv",dtype=str)
    for _,r in dp.iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s and s!="nan": srr2donor[s.strip()]=f"{ds}_{r['donor_id']}"
LIN=["Hepatocytes","Fibroblasts","Cholangiocytes","Macrophages","T_cells","Endothelial_cells"]
prof={}
for ln in LIN:
    x=pd.read_csv(D/f"pseudobulk_{ln}.tsv.gz",sep="\t",index_col=0)
    x=x.groupby(level=0).sum()
    cols=pd.Series(x.columns).map(srr2donor).fillna(pd.Series(x.columns)).values
    x=x.T.groupby(cols).sum().T
    cpm=x.div(x.sum(axis=0).replace(0,np.nan),axis=1)*1e6
    prof[ln]=np.log2(cpm+1).mean(axis=1)
P=pd.DataFrame(prof)
Z=P.sub(P.mean(axis=1),axis=0).div(P.std(axis=1).replace(0,np.nan),axis=0)
idx=(Z.Cholangiocytes-Z.Hepatocytes).rename("chol_minus_hep")
meanexp=P.mean(axis=1).rename("mean_log2cpm")
# --- gene-level marker check (label-free validation of the per-gene index)
mk={"hep":["ALB","APOA1","APOC3","CYP3A4","CYP2E1","SERPINA1","TF","HP"],
    "chol":["KRT19","KRT7","SOX9","EPCAM","CFTR","SPP1","ANXA4","BICC1"],
    "other":["COL1A1","PECAM1","CD68","CD3E","GGT1"]}
for k,v in mk.items():
    print(k, {g:(round(float(idx.get(g,np.nan)),2), round(float(meanexp.get(g,np.nan)),2)) for g in v})
# --- coverage of COLOC genes
m["in_pb"]=m.symbol.isin(P.index)
m["mean_log2cpm"]=meanexp.reindex(m.symbol).values
m["idx"]=idx.reindex(m.symbol).values
m["expressed"]=m.mean_log2cpm>0.1
pp=m[m.locus.isin(pg.locus.unique())&m.cls.isin(["GGT_only","ALTAST_only"])]
print("pure GGT/ALTAST genes:",len(pp),"symbol mapped:",pp.symbol.notna().sum(),"in pseudobulk:",pp.in_pb.sum(),"mean log2CPM>0.1:",pp.expressed.sum(),"idx finite:",np.isfinite(pp.idx).sum())
print("by class genes in pb / finite idx:",pp.groupby("cls").apply(lambda g: (int(g.in_pb.sum()), int(np.isfinite(g.idx).sum()), len(g))).to_dict())
loc=pp[np.isfinite(pp.idx)].groupby(["locus","cls"]).agg(idx=("idx","mean"),mexp=("mean_log2cpm","mean"),ngene=("idx","size")).reset_index()
print("loci with >=1 finite-index gene by class:",loc.groupby("cls").size().to_dict())
print("genes per locus distribution:",loc.ngene.value_counts().sort_index().to_dict())
# label-free pooled SD
sd=loc.idx.std(ddof=1); print("pooled locus-level SD of chol-minus-hep index (labels not used):",round(sd,4))
print("genome-wide gene-level SD of index (genes mean log2CPM>0.1):",round(idx[meanexp>0.1].std(),4), "n genes",int((meanexp>0.1).sum()))
# expression decile residual SD (label-free)
loc["dec"]=pd.qcut(loc.mexp,10,labels=False)
res=loc.idx-loc.groupby("dec").idx.transform("mean")
sdr=np.sqrt((res**2).sum()/(len(loc)-10)); print("within-decile pooled SD (label-free):",round(sdr,4))
nG=(loc.cls=="GGT_only").sum(); nA=(loc.cls=="ALTAST_only").sum()
for a in (0.05,):
    zt=stats.norm.ppf(1-a/2)+stats.norm.ppf(0.8)
    d=zt*np.sqrt(1/nG+1/nA)
    print(f"MDE (normal approx) nG={nG} nA={nA}: {d:.3f} SD = {d*sd:.3f} index units (raw SD), {d*sdr:.3f} (within-decile SD)")
# UKBB sensitivity MDE
loc_counts={}
pc2=m2[m2.locus.isin(p)].groupby("cls_ukbb").locus.nunique()
nG2,nA2=pc2.get("GGT_only",0),pc2.get("ALTAST_only",0)
d2=zt*np.sqrt(1/nG2+1/nA2); print(f"UKBB sensitivity MDE nG={nG2} nA={nA2}: {d2:.3f} SD")
# per-decile class counts are labels - only report decile sizes, not class means
print("GGT_only loci per expression decile:",loc[loc.cls=="GGT_only"].groupby("dec").size().reindex(range(10),fill_value=0).to_dict())
print("ALTAST_only loci per expression decile:",loc[loc.cls=="ALTAST_only"].groupby("dec").size().reindex(range(10),fill_value=0).to_dict())
# reference: empirical index difference between marker classes, in pooled-SD units
hep=[g for g in mk["hep"] if g in idx.index]; ch=[g for g in mk["chol"] if g in idx.index]
print("marker contrast chol-markers minus hep-markers, index units:",round(idx[ch].mean()-idx[hep].mean(),3))
pd.DataFrame({"idx":idx,"mean_log2cpm":meanexp}).to_csv(S/"index_all_genes.tsv",sep="\t")
