# Identity checks only: ambient raw vs primary score; donor universes. No stage statistic computed.
import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
F=R+"figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/"
amb=pd.read_csv(F+"ambient_program_effects.tsv",sep="\t")
print("ambient effects rows",len(amb),amb.analysis_universe.value_counts().to_dict())
cc=amb[amb.analysis_universe=="complete_case_common_universe"]
print("cc rows",len(cc),"failed datasets",cc.failed_decontx_datasets.value_counts().to_dict())
v=pd.read_csv(R+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R04-01-reproduce/refit_drop_gse189600_family.tsv",sep="\t")
m=cc.merge(v,on=["cell_type","module"])
print("cc raw_q vs verifier q max abs diff",np.nanmax(np.abs(m.raw_qvalue-m.q)),"n diff",np.nanmax(np.abs(m.n_model_donors-m.n)))
print("Ductular raw_q",cc[(cc.cell_type=="hepatocytes")&(cc.module==20)].raw_qvalue.values)
s=pd.read_csv(R+"Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv",sep="\t",usecols=["cell_type","module","donor","score","dataset"])
meta=pd.read_csv(R+"Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv",sep="\t")
keep=meta[(meta.exclude_stage_analysis.astype(str).str.upper()=="FALSE")&meta.disease_stage_coarse.isin(["Healthy","Steatosis","Steatohepatitis"])&(meta.dataset!="GSE189600")]
s60=s[s.donor.isin(keep["sample"])]
a=pd.read_csv(R+"Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/results/donor_program_scores.tsv.gz",sep="\t")
print("ambient score rows",len(a),"donors",a.donor.nunique(),"programs",a.groupby(["cell_type","module"]).ngroups, "datasets",a.dataset_from_cells.value_counts().to_dict())
j=s60.merge(a,on=["cell_type","module","donor"],how="left",indicator=True)
print("60-universe rows",len(s60),"matched in ambient",(j._merge=="both").sum())
print("max |raw - primary score|",np.nanmax(np.abs(j.raw-j.score)),"n raw NA",j.raw.isna().sum(),"n corrected NA",j.corrected.isna().sum())
print("corrected NA by dataset:",j[j.corrected.isna()].dataset.value_counts().to_dict())
print("corrected NA by cell_type:",j[j.corrected.isna()].cell_type.value_counts().to_dict())
print("cc corrected_n vs raw_n differ:",(cc.corrected_n!=cc.raw_n).sum())
na_don=j[j.corrected.isna()].groupby("dataset").donor.unique()
print({k:list(v) for k,v in na_don.items()})
