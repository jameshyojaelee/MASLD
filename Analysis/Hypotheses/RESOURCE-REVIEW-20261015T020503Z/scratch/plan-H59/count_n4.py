# Counts only: informative-fit n per program (GSE202379 + GSE244832, 60-donor universe), by lineage.
import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
H=R+"Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"
M=R+"Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
s=pd.read_csv(H,sep="\t",keep_default_na=False,usecols=["cell_type","module","donor","dataset"])
m=pd.read_csv(M,sep="\t",keep_default_na=False).rename(columns={"sample":"donor"})[["donor","disease_stage_coarse","exclude_stage_analysis"]]
d=s.merge(m,on="donor",how="inner",validate="many_to_one")
excl=d.exclude_stage_analysis.astype(str).str.lower().isin(["true","t","1"])
d=d[(~excl)&d.disease_stage_coarse.isin(["Healthy","Steatosis","Steatohepatitis"])]
print("62-universe donors per program: min/median/max", d.groupby(["cell_type","module"]).size().agg(["min","median","max"]).to_dict())
inf=d[d.dataset.isin(["GSE202379","GSE244832"])]
n=inf.groupby(["cell_type","module"]).size()
print("informative n per program: min/median/max", n.agg(["min","median","max"]).to_dict())
print(n.groupby(level="cell_type").agg(["min","median","max"]).to_string())
