# Blinded gate checks for the signed TWAS-vs-disease-state candidate.
# The observed primary statistic (TWAS z vs disease t, unshifted) is NOT computed.
import pandas as pd, numpy as np, glob, json
from scipy.stats import rankdata, spearmanr
B="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
OT=B+"RNA-seq/results/causal_inference/otters_broadaway/"
SEED=20261015; NDRAW=2000
de=pd.read_csv(B+"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
de["sid"]=de.gene.str.replace(r"\.\d+$","",regex=True)
de["canon"]=(de.padj<0.05)&(de.logFC.abs()>0.5)
def load(t):
    tw=pd.concat([pd.read_csv(f,sep="\t") for f in glob.glob(OT+f"{t}/chr*/P0.05.txt")]).dropna(subset=["FUSION_Z"])
    tw["mid"]=(tw.GeneStart+tw.GeneEnd)/2
    return tw.drop_duplicates("TargetID")
tw={t:load(t) for t in ["ukbb_alt","ukbb_ast","ghodsian"]}
z=lambda t,g: tw[t].set_index("TargetID").FUSION_Z.get(g,np.nan)
# A. orientation
for g,nm in [("ENSG00000167701","GPT"),("ENSG00000120053","GOT1"),("ENSG00000125166","GOT2")]:
    print("orient",nm,g,"ALT z %.2f AST z %.2f NAFLD z %.2f"%(z("ukbb_alt",g),z("ukbb_ast",g),z("ghodsian",g)))
for a,b in [("ukbb_alt","ghodsian"),("ukbb_alt","ukbb_ast")]:
    j=tw[a].set_index("TargetID").FUSION_Z.to_frame("a").join(tw[b].set_index("TargetID").FUSION_Z.rename("b"),how="inner")
    print("cross-trait",a,b,"n",len(j),"spearman %.4f"%spearmanr(j.a,j.b)[0])
for t in tw: print("mean z^2",t,"%.3f"%np.mean(tw[t].FUSION_Z**2), "n",len(tw[t]))
rng=np.random.default_rng(SEED)
def std_rank(x):
    r=rankdata(x); return (r-r.mean())/r.std()
def make_shifts(chrom_sizes):
    out=[]
    for _ in range(NDRAW):
        s=[]
        for n in chrom_sizes:
            m=min(500,n//4); s.append(int(rng.integers(m,n-m+1)))
        out.append(s)
    return out
def within_chr_null(df,xcol,ycol,shifts=None,observed=False):
    df=df.sort_values(["CHROM","mid"]).reset_index(drop=True)
    xr=std_rank(df[xcol].values); yr=std_rank(df[ycol].values)
    bounds=df.groupby("CHROM").indices; keys=sorted(bounds)
    sizes=[len(bounds[k]) for k in keys]
    if shifts is None: shifts=make_shifts(sizes)
    null=[]
    for s in shifts:
        xs=np.concatenate([np.roll(xr[bounds[k]],sh) for k,sh in zip(keys,s)])
        null.append(np.mean(xs*yr))
    null=np.array(null)
    obs=np.mean(xr*yr) if observed else None
    return null,obs,sizes
res={}
# B. blinded null for primary (no observed)
for t in ["ukbb_alt","ghodsian"]:
    m=tw[t].merge(de[["sid","t","AveExpr","canon"]],left_on="TargetID",right_on="sid")
    null,_,sizes=within_chr_null(m,"FUSION_Z","t")
    q=np.percentile(null,[2.5,97.5,1.25,98.75])
    sd=null.std()
    print(f"PRIMARY-NULL {t} n={len(m)} chroms={len(sizes)} minchr={min(sizes)} null mean {null.mean():.4f} SD {sd:.4f} 2.5/97.5 {q[0]:.4f} {q[1]:.4f} 1.25/98.75 {q[2]:.4f} {q[3]:.4f} indepSD {1/np.sqrt(len(m)):.4f}")
    print(f"  MDE (two-sided alpha 0.025 Holm-first, 80% power) ~ {(2.241+0.842)*sd:.4f}; alpha 0.05 ~ {(1.96+0.842)*sd:.4f}")
    # D. negative controls through same null: TWAS z vs AveExpr; TWAS z vs log n_snps
    for nc in ["AveExpr"]:
        nn,obs,_=within_chr_null(m,"FUSION_Z",nc,observed=True)
        p=(1+np.sum(np.abs(nn)>=abs(obs)))/(NDRAW+1)
        print(f"  NEG-CONTROL {t} TWAS z vs {nc}: rho {obs:.4f} null SD {nn.std():.4f} two-sided p {p:.4f}")
    m["lognsnp"]=np.log(m.n_snps)
    nn,obs,_=within_chr_null(m,"FUSION_Z","lognsnp",observed=True)
    p=(1+np.sum(np.abs(nn)>=abs(obs)))/(NDRAW+1)
    print(f"  NEG-CONTROL {t} TWAS z vs log n_snps: rho {obs:.4f} null SD {nn.std():.4f} two-sided p {p:.4f}")
    res[t]=dict(n=len(m),null_sd=float(sd))
# C. positive controls through same null
j=tw["ukbb_alt"][["TargetID","CHROM","mid","FUSION_Z"]].merge(tw["ghodsian"][["TargetID","FUSION_Z"]].rename(columns={"FUSION_Z":"z2"}),on="TargetID")
nn,obs,_=within_chr_null(j,"FUSION_Z","z2",observed=True)
print(f"POS-CONTROL ALT vs NAFLD TWAS z: n {len(j)} rho(rank-product) {obs:.4f} null SD {nn.std():.4f} null max {np.abs(nn).max():.4f} p {(1+np.sum(np.abs(nn)>=abs(obs)))/(NDRAW+1):.5f}")
c=pd.read_csv(OT+"concordance/otters_vs_spredixcan_concordance.csv")
c=c.merge(tw["ghodsian"][["TargetID","CHROM","mid"]],left_on="ensembl_id",right_on="TargetID")
nn,obs,_=within_chr_null(c,"otters_z","spx_z",observed=True)
print(f"POS-CONTROL OTTERS vs S-PrediXcan (ghodsian): n {len(c)} rho {obs:.4f} null SD {nn.std():.4f} null max {np.abs(nn).max():.4f} p {(1+np.sum(np.abs(nn)>=abs(obs)))/(NDRAW+1):.5f}")
# E. sizes for sensitivities
col=pd.read_csv(B+"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
gen=set(col.loc[col.coloc_best_susie_pp4>0.5,"ensembl"])
for t in ["ukbb_alt","ghodsian"]:
    m=tw[t].merge(de[["sid","t","AveExpr","canon"]],left_on="TargetID",right_on="sid")
    ing=m.TargetID.isin(gen)
    keep=np.ones(len(m),bool)
    for ch,g in m[ing].groupby("CHROM"):
        idx=(m.CHROM==ch).values
        d=np.min(np.abs(m.loc[idx,"mid"].values[:,None]-g.mid.values[None,:]),axis=1)
        k=keep[idx]; k[d<=1e6]=False; keep[idx]=k
    print(f"SENS {t}: colocalized (SuSiE PP.H4>0.5) in join {int(ing.sum())}; genes left after removing them and genes within 1 Mb (window midpoints) {int(keep.sum())}")
    print(f"  |z|>3 & canonical {int(((m.FUSION_Z.abs()>3)&m.canon).sum())}; |z|>3 {int((m.FUSION_Z.abs()>3).sum())}")
json.dump(res,open("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/critic-cand-1/gate_checks.json","w"),indent=1)
