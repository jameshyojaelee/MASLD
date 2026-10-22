# Counts only: gate G1 feasibility (levels with >=2 donors per informative dataset), stability_median of named programs.
import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
H=R+"Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"
M=R+"Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
G=R+"figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"
s=pd.read_csv(H,sep="\t",keep_default_na=False,usecols=["cell_type","module","donor","dataset"])
m=pd.read_csv(M,sep="\t",keep_default_na=False).rename(columns={"sample":"donor"})[["donor","disease_stage_coarse","exclude_stage_analysis"]]
d=s.merge(m,on="donor",how="inner",validate="many_to_one")
excl=d.exclude_stage_analysis.astype(str).str.lower().isin(["true","t","1"])
d=d[(~excl)&d.disease_stage_coarse.isin(["Healthy","Steatosis","Steatohepatitis"])&(d.dataset!="GSE189600")]
t=d[d.dataset.isin(["GSE202379","GSE244832"])].groupby(["cell_type","module","dataset","disease_stage_coarse"]).size().unstack(fill_value=0)
t["lev2"]=(t[["Healthy","Steatosis","Steatohepatitis"]]>=2).sum(1); t["n"]=t[["Healthy","Steatosis","Steatohepatitis"]].sum(1)
ok=t.groupby(level=["cell_type","module"]).apply(lambda x: ((x.lev2>=2)&(x.n>=8)).all() and len(x)==2)
print("programs meeting G1 (both informative datasets: >=2 stage levels with >=2 donors, n>=8):",int(ok.sum()),"of",len(ok))
print("program x dataset rows with exactly 2 levels having >=2 donors:",int((t.lev2==2).sum()))
r=pd.read_csv(G,sep="\t",keep_default_na=False)
r["stability_median"]=pd.to_numeric(r.stability_median)
print("stability_median >=0.5:",int((r.stability_median>=0.5).sum()),"of",len(r))
print(r[(r.cell_type=="hepatocytes")&(r.module.isin([8,20]))][["cell_type","module","module_name","stability_median","source_stability"]].to_string(index=False))
print("single-stage datasets carry 0 stage information; donors in GSE174748+GSE185477 per program: max",d[d.dataset.isin(["GSE174748","GSE185477"])].groupby(["cell_type","module"]).size().max())
