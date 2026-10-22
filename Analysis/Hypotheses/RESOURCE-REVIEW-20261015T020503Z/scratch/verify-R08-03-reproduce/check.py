import pandas as pd, numpy as np
from scipy import stats
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
E=ROOT+"Analysis/MASLD_Model_Benchmark/executions/external-validation-expansion-20260901T190000Z/"
sc=pd.read_csv(E+"GSE213621_scores.tsv",sep="\t")
meta=pd.read_csv(ROOT+"Analysis/Deconvolution/bulk/GSE213621/GSE213621_metadata.tsv",sep="\t")
print("meta cols first 3:",list(meta.columns[:3]),"index name sample?", meta.index[:2].tolist())
print("fibrotic_stage counts:",meta["fibrotic_stage"].value_counts().to_dict())
m=sc.merge(meta[["sample_id","fibrotic_stage"]],on="sample_id",how="inner")
print("scored",len(sc),"joined",len(m))
mp={"F0F1":0,"F2":1,"F3F4":2}
m["y"]=m["fibrotic_stage"].map(mp)
d=m[m.y.notna()].copy()
def tie_ceiling(y):
    y=np.asarray(y,float); mr=stats.rankdata(y,method="average")
    return float(np.corrcoef(np.sort(mr),np.arange(1,len(y)+1))[0,1])
rho=stats.spearmanr(d.latent_severity,d.y).statistic; c=tie_ceiling(d.y)
print(f"published reproduction: n={len(d)} rho={rho:.6f} ceiling={c:.6f} ratio={rho/c:.6f}")
ro=pd.read_csv(ROOT+"figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv",sep="\t",low_memory=False)
g=ro[ro.dataset=="GSE213621"]
print("roster GSE213621 rows",len(g),"unique sample_id",g.sample_id.nunique(),"screen_components",g.screen_component.nunique())
# cross-dataset components
cc=ro.groupby("screen_component").dataset.nunique()
print("components spanning >1 dataset (all roster):",int((cc>1).sum()))
print("GSE213621 components spanning >1 dataset:",int((cc.loc[g.screen_component.unique()]>1).sum()))
print("GSE213621 component sizes:",g.groupby("screen_component").size().value_counts().sort_index().to_dict())
print("filtered_component unique",g.filtered_component.nunique(),"a1_component unique",g.a1_component.nunique())
print("identity_qualified:",g.identity_qualified.value_counts().to_dict(),"primary_selected:",g.primary_selected.value_counts().to_dict())
print("exclusion_reason:",g.exclusion_reason.value_counts(dropna=False).to_dict())
d=d.merge(g[["sample_id","screen_component","source_fibrosis_category","primary_selected","identity_qualified"]],on="sample_id",how="left")
print("299 missing roster:",int(d.screen_component.isna().sum()))
print("stage agreement roster vs meta:",(d.source_fibrosis_category.astype(str)==d.fibrotic_stage.astype(str)).sum(), "of", len(d))
print("clusters in 299:",d.screen_component.nunique())
sz=d.groupby("screen_component").size()
print("cluster size dist:",sz.value_counts().sort_index().to_dict())
multi=sz[sz>1].index
conf=[]
for k in multi:
    s=d[d.screen_component==k]
    st=sorted(s.fibrotic_stage.unique())
    if len(st)>1: conf.append("/".join(st))
print("conflicting clusters:",len(conf),pd.Series(conf).value_counts().to_dict())
# also check against full 367 with controls
m2=m.merge(g[["sample_id","screen_component"]],on="sample_id")
sz2=m2.groupby("screen_component").size()
print("367 set: clusters",m2.screen_component.nunique(),"multi",int((sz2>1).sum()))
cf2=m2.groupby("screen_component").fibrotic_stage.nunique()
mix=m2[m2.screen_component.isin(cf2[cf2>1].index)].groupby("screen_component").fibrotic_stage.apply(lambda x:"/".join(sorted(x.unique())))
print("367 conflicting:",mix.value_counts().to_dict())
# dedup analyses
nonconf=[k for k in sz.index if d[d.screen_component==k].fibrotic_stage.nunique()==1]
dd=d[d.screen_component.isin(nonconf)]
for rule in ["first","mean"]:
    if rule=="first":
        x=dd.sort_values("sample_id").groupby("screen_component").first()
        sv=x.latent_severity; yy=x.y
    else:
        x=dd.groupby("screen_component").agg(latent_severity=("latent_severity","mean"),y=("y","first"))
        sv=x.latent_severity; yy=x.y
    r=stats.spearmanr(sv,yy).statistic; c=tie_ceiling(yy)
    print(f"drop conflicting, one per cluster ({rule}): n={len(x)} rho={r:.6f} ceil={c:.6f} ratio={r/c:.6f}")
# variants: last record, random seeds
x=dd.sort_values("sample_id").groupby("screen_component").last()
r=stats.spearmanr(x.latent_severity,x.y).statistic;c=tie_ceiling(x.y)
print(f"drop conflicting, last: n={len(x)} rho={r:.6f} ceil={c:.6f} ratio={r/c:.6f}")
rs=[]
rng=np.random.default_rng(20261015)
for i in range(1000):
    x=dd.sample(frac=1,random_state=int(rng.integers(1e9))).groupby("screen_component").first()
    rs.append(stats.spearmanr(x.latent_severity,x.y).statistic/tie_ceiling(x.y))
print("random one-per-cluster ratio over 1000 draws: min %.4f median %.4f max %.4f"%(min(rs),np.median(rs),max(rs)))
# just drop conflicting, keep all records
r=stats.spearmanr(dd.latent_severity,dd.y).statistic;c=tie_ceiling(dd.y)
print(f"drop conflicting only, keep duplicates: n={len(dd)} rho={r:.6f} ceil={c:.6f} ratio={r/c:.6f}")
# keep all clusters one per cluster (including conflicting, first record)
x=d.sort_values("sample_id").groupby("screen_component").first()
r=stats.spearmanr(x.latent_severity,x.y).statistic;c=tie_ceiling(x.y)
print(f"all 274 clusters, first record: n={len(x)} rho={r:.6f} ceil={c:.6f} ratio={r/c:.6f}")
d.to_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R08-03-reproduce/gse213621_299_with_components.tsv",sep="\t",index=False)
