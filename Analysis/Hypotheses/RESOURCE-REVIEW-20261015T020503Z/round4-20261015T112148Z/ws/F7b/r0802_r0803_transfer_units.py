"""R08-02 / R08-03: re-derive Fig 7 transfer counts and tie-ceiling ratios by etiology and by identity cluster.
Inputs are the saved released-scorer outputs; no model refit. Seed 20260923."""
import numpy as np, pandas as pd
from scipy import stats
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
E=f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/external-validation-expansion-20260901T190000Z"
OUT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F7b"
rng=np.random.default_rng(20260923)
def tie_ceiling(y):  # verbatim from external_validation.py:49-53
    y=np.asarray(y,float)
    if len(np.unique(y))<2: return float("nan")
    mr=stats.rankdata(y,method="average")
    return float(np.corrcoef(np.sort(mr),np.arange(1,len(y)+1))[0,1])
def row(cohort,arm,d,unit,nunits,boot=None):
    rho=stats.spearmanr(d.score,d.y).statistic if d.y.nunique()>1 else np.nan
    c=tie_ceiling(d.y)
    r=dict(cohort=cohort,subset=arm,n_records=len(d),n_units=nunits,unit=unit,spearman=rho,tie_ceiling=c,ratio=rho/c if c==c else np.nan)
    if boot is not None: r.update(boot)
    return r
def cluster_boot(d,cl,B=4000):
    g=d.groupby(cl).indices; keys=list(g); out=[]
    for _ in range(B):
        idx=np.concatenate([g[k] for k in rng.choice(keys,len(keys))])
        s=d.iloc[idx]; out.append(stats.spearmanr(s.score,s.y).statistic)
    return dict(boot_lo=np.quantile(out,.025),boot_hi=np.quantile(out,.975),boot_B=B,boot_unit=cl)
rows=[]
# GSE276114 by etiology
s=pd.read_csv(f"{E}/GSE276114_scores.tsv",sep="\t")
dz=pd.read_csv(f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/etiology-transfer-fibrosis-distinctness-20260901T135543Z/gse276114_design.tsv",sep="\t")
d=s.merge(dz,left_on="sample_id",right_on="matrix_column_name",validate="1:1")
assert len(d)==177
d["y"]=d.disease_group.map({"F0-2":0,"F3":1,"F4":2}); d["score"]=d.latent_severity
print(pd.crosstab(d.disease,d.disease_group))
rows.append(row("GSE276114","all records (published)",d,"record",177))
for dis in ["MASLD","CVH","ARLD"]:
    sub=d[d.disease==dis]; rows.append(row("GSE276114",dis,sub,"record",len(sub)))
sub=d[d.disease!="MASLD"]; rows.append(row("GSE276114","non-MASLD (CVH+ARLD)",sub,"record",len(sub)))
# GSE213621 identity
s=pd.read_csv(f"{E}/GSE213621_scores.tsv",sep="\t")
m=pd.read_csv(f"{ROOT}/Analysis/Deconvolution/bulk/GSE213621/GSE213621_metadata.tsv",sep="\t")  # ragged N+1 header: pandas puts the extra first field in the index; named columns align
ro=pd.read_csv(f"{ROOT}/figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv",sep="\t",low_memory=False)
ro=ro[ro.dataset=="GSE213621"][["sample_id","screen_component","primary_selected"]]
print("roster GSE213621 rows",len(ro),"components",ro.screen_component.nunique())
d=s.merge(m[["sample_id","fibrotic_stage"]],on="sample_id",validate="1:1").merge(ro,on="sample_id",how="left",validate="1:1")
assert d.screen_component.notna().all()
d["score"]=d.latent_severity
for arm,mp in [("3level",{"F0F1":0,"F2":1,"F3F4":2}),("4level",{"Control":0,"F0F1":1,"F2":2,"F3F4":3})]:
    a=d[d.fibrotic_stage.isin(mp)].copy(); a["y"]=a.fibrotic_stage.map(mp)
    ncl=a.screen_component.nunique()
    pairs=a.groupby("screen_component").filter(lambda g:len(g)>1)
    conf=pairs.groupby("screen_component").y.nunique().gt(1)
    print(arm,"records",len(a),"clusters",ncl,"multi-library clusters",pairs.screen_component.nunique(),"conflicting",int(conf.sum()))
    if arm=="3level":
        print(pairs[pairs.screen_component.isin(conf[conf].index)].groupby("screen_component").fibrotic_stage.apply(lambda x:"/".join(sorted(x))).value_counts())
    rows.append(row("GSE213621",f"{arm} all records (published)",a,"record",len(a),cluster_boot(a.reset_index(drop=True),"screen_component")))
    bad=set(conf[conf].index)
    p=a[(a.primary_selected==True)&~a.screen_component.isin(bad)].reset_index(drop=True)
    assert p.screen_component.is_unique
    rows.append(row("GSE213621",f"{arm} one record per identity, conflicting-stage identities dropped",p,"identity",len(p),cluster_boot(p,"screen_component")))
    f=a.drop_duplicates("screen_component").reset_index(drop=True)
    rows.append(row("GSE213621",f"{arm} one record per identity (first), all identities",f,"identity",len(f)))
rows.append(dict(cohort="GSE268273",subset="all records (published, release expected_ood.json)",n_records=109,n_units=109,unit="record (not identity-screened)",spearman=0.6631190836675596,tie_ceiling=0.9523223203684403,ratio=0.6631190836675596/0.9523223203684403))
t=pd.DataFrame(rows); t.to_csv(f"{OUT}/fig7_transfer_units_by_etiology_and_identity.tsv",sep="\t",index=False,float_format="%.6f")
pd.set_option("display.width",250); print(t.round(4).to_string())
