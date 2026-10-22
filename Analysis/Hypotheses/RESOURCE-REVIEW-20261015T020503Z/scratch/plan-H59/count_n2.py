# Counts only (no score values): per program, donors per stage within each informative dataset, 60-donor universe.
import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
H=R+"Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"
M=R+"Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
A=R+"Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/results/donor_program_scores.tsv.gz"
s=pd.read_csv(H,sep="\t",keep_default_na=False,usecols=["cell_type","module","donor","dataset"])
m=pd.read_csv(M,sep="\t",keep_default_na=False).rename(columns={"sample":"donor"})
m=m[["donor","disease_stage_coarse","exclude_stage_analysis"]]
d=s.merge(m,on="donor",how="inner",validate="many_to_one")
excl=d.exclude_stage_analysis.astype(str).str.lower().isin(["true","t","1"])
d=d[(~excl)&d.disease_stage_coarse.isin(["Healthy","Steatosis","Steatohepatitis"])&(d.dataset!="GSE189600")]
t=d.groupby(["cell_type","module","dataset","disease_stage_coarse"]).size().unstack(fill_value=0)
t=t.loc[(slice(None),slice(None),["GSE202379","GSE244832"]),:]
t["n"]=t.sum(1); t["levels"]=(t[["Healthy","Steatosis","Steatohepatitis"]]>0).sum(1)
print(t.groupby(level="dataset").agg(["min","median","max"]).T.to_string())
low=t[t[["Healthy","Steatosis","Steatohepatitis"]].replace(0,99).min(1)<2]
print("\nprogram x dataset rows with a non-empty stage level holding a single donor:",len(low))
print("rows with Healthy==0:",(t.Healthy==0).sum(), " Healthy==1:",(t.Healthy==1).sum())
print("programs with <2 stage levels in GSE202379:",(t.xs("GSE202379",level="dataset").levels<2).sum(), " in GSE244832:",(t.xs("GSE244832",level="dataset").levels<2).sum())
a=pd.read_csv(A,sep="\t",nrows=3); print("\nambient columns:",list(a.columns))
a=pd.read_csv(A,sep="\t",keep_default_na=False)
print("ambient rows",len(a)); 
for c in a.columns:
    if a[c].dtype==object and a[c].nunique()<12: print(c,a[c].value_counts().to_dict())
