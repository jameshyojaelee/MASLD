import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
S=R+"figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions/"
sm=pd.read_csv(S+"stage_extension_sample_manifest.tsv",sep="\t")
print(sm.columns.tolist())
fib=sm[sm.axis.isin(["fibrosis","fibrosis_adjacent"])]
print(pd.crosstab([fib.contrast],[fib.dataset]))
print("unique units in fibrosis contrasts:",fib.analysis_unit_id.nunique())
u=fib.drop_duplicates("analysis_unit_id")
print(pd.crosstab(u.dataset,u.group,margins=True))
ssm=pd.read_csv(S+"stage_sample_manifest.tsv",sep="\t")
print("stage_sample_manifest rows",len(ssm),"unique",ssm.analysis_unit_id.nunique())
print(pd.crosstab(ssm.dataset,ssm.fib_group,margins=True))
nas=sm[sm.axis.isin(["NAS","NAS_adjacent"])]
print("NAS datasets", nas.groupby("dataset").analysis_unit_id.nunique().to_dict())
f5=pd.read_csv(S+"five_cohort_sample_manifest.tsv",sep="\t")
print("five cohort rows",len(f5))
print(pd.crosstab(f5.dataset,f5.fibrosis_stage.fillna(-1),margins=True))
g=f5[f5.dataset=="GSE213621"]
print(pd.crosstab([g.group_binary,g.diagnosis_harmonized.fillna("NA")],g.fibrosis_stage.fillna(-1)))
g=f5[f5.dataset=="GSE162694"]
print(pd.crosstab(g.group_binary,g.fibrosis_stage.fillna(-1)))
# F0 composition in stage fits
f0=ssm[ssm.fib_group=="F0"].merge(f5[["analysis_unit_id","group_binary"]],on="analysis_unit_id",how="left")
print(pd.crosstab(f0.dataset,f0.group_binary,margins=True))
