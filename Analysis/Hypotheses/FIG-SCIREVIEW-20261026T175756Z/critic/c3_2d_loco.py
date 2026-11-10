# Fig 2D independent recompute: GWS-lead credible sets, 1-Mb single-linkage clusters, cluster-mean fractions,
# physical-cluster bootstrap (own seed), ancestry strata with within-stratum clusters.
import pandas as pd, numpy as np
A="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
cs=pd.read_csv(A+"figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv",sep="\t",dtype={"locus":str})
print("sets after PIP-sum gate",len(cs), cs.trait_scope.value_counts().to_dict())
L=[]
for s in cs.study.unique():
    t=pd.read_csv(A+f"GWAS/finemapping/data/lead_snps/{s}_leadSNPs_tiered.tsv",sep="\t"); t["study"]=s; L.append(t)
L=pd.concat(L); L["chrom"]=L.CHR.astype(str); L["pos"]=L.BP.astype(int)
cs["chrom"]=cs.chrom.astype(str)
m=cs.merge(L[["study","chrom","pos","index_p"]],on=["study","chrom","pos"],how="left")
sent=pd.read_csv(A+"figures/candidates/fig2d-gws-leads-20261021T124551Z/values/unmatched_sentinel_p.tsv",sep="\t").set_index("credible_set_uid").sentinel_p
m["p"]=m.index_p.fillna(m.credible_set_uid.map(sent))
print("unmatched",m.index_p.isna().sum(),"unmatched with sentinel p",m.credible_set_uid.map(sent).notna().sum(),"unmatched non-EUR",(m.index_p.isna()&(m.ancestry!="EUR")).sum())
g=m[m.p<5e-8].copy(); print("GWS sets",len(g),g.trait_scope.value_counts().to_dict())
print(pd.crosstab(g.trait_scope,g.ancestry))
cod=["protein_altering_pip_mass","canonical_splice_pip_mass"]; arch=cod+["synonymous_or_utr_pip_mass","other_noncoding_pip_mass"]
g["coding"]=g[cod].sum(1); g["tot"]=g[arch].sum(1)
def clusters(d):
    d=d.copy(); d["ci"]=d.chrom.astype(int); d=d.sort_values(["ci","pos"])
    cid=[];c=-1;prev=None
    for ch,p in zip(d.ci,d.pos):
        if prev is None or ch!=prev[0] or p-prev[1]>1_000_000: c+=1
        cid.append(c); prev=(ch,p)
    d["cl"]=cid; return d
D="tier1_direct_masld_pdff"
def est(d, B=10000, seed=1, out=True):
    d=clusters(d)
    t=d.groupby(["trait_scope","cl"])[arch+["coding","tot"]].mean()
    res={}
    for sc in d.trait_scope.unique():
        x=t.loc[sc]; res[sc]={c:round(100*x[c].sum()/x.tot.sum(),2) for c in arch}; res[sc]["coding"]=100*x.coding.sum()/x.tot.sum(); res[sc]["nclus"]=len(x); res[sc]["nsets"]=(d.trait_scope==sc).sum()
    w=t[["coding","tot"]].unstack(0).fillna(0)
    nd=w[("coding",D)].values; dd=w[("tot",D)].values; ne=w[("coding","tier2_liver_enzyme")].values; de=w[("tot","tier2_liver_enzyme")].values
    rng=np.random.default_rng(seed); idx=rng.integers(0,len(w),(B,len(w)))
    with np.errstate(all="ignore"):
        diff=nd[idx].sum(1)/dd[idx].sum(1)-ne[idx].sum(1)/de[idx].sum(1)
    diff=diff[np.isfinite(diff)]
    pt=res[D]["coding"]-res["tier2_liver_enzyme"]["coding"]
    return res,len(w),pt,np.percentile(diff,[2.5,97.5])*100,len(diff)
# Critic addition: leave-one-direct-cluster-out sensitivity of the 2D pooled difference.
d = clusters(g)
dc = d[d.trait_scope == D].groupby("cl").agg(coding=("coding", "mean"), tot=("tot", "mean"), chrom=("chrom", "first"), pos=("pos", "min"), n=("coding", "size"),
                                            trait=("trait", lambda v: ",".join(sorted(set(v)))), anc=("ancestry", lambda v: ",".join(sorted(set(v)))))
dc["frac"] = dc.coding / dc.tot
dc = dc.sort_values("coding", ascending=False)
print(dc.round(3).to_string())
res, nphys, pt, ci, nf = est(g)
print("all: diff", round(pt, 2), "CI", ci.round(1))
for k in [1, 2, 3]:
    drop = dc.index[:k]
    keep = d[~((d.trait_scope == D) & d.cl.isin(drop))]
    keep = keep.drop(columns=["cl", "ci"])
    r = est(keep)
    print(f"drop top {k} direct clusters {[(dc.loc[c,'chrom'], int(dc.loc[c,'pos'])) for c in drop]}: direct {r[0][D]['coding']:.2f} enzyme {r[0]['tier2_liver_enzyme']['coding']:.2f} diff {r[2]:.2f} CI {r[3].round(1)}")
