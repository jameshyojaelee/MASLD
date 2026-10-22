# Counts only: donors by dataset x stage per program. Reads no score values beyond presence.
import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
H=R+"Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"
M=R+"Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
s=pd.read_csv(H,sep="\t",keep_default_na=False,usecols=["cell_type","module","donor","dataset"])
m=pd.read_csv(M,sep="\t",keep_default_na=False).rename(columns={"sample":"donor"})
m=m[["donor","dataset","disease_stage_coarse","exclude_stage_analysis"]].rename(columns={"dataset":"ds_m"})
d=s.merge(m,on="donor",how="inner",validate="many_to_one")
excl=d.exclude_stage_analysis.astype(str).str.lower().isin(["true","t","1"])
d=d[(~excl)&d.disease_stage_coarse.isin(["Healthy","Steatosis","Steatohepatitis"])]
print("programs",d.groupby(["cell_type","module"]).ngroups,"donors",d.donor.nunique())
print(d.drop_duplicates("donor").groupby(["dataset","disease_stage_coarse"]).size().unstack(fill_value=0))
for ct,mod,name in [("hepatocytes",20,"Ductular injury"),("hepatocytes",8,"Stromal ECM")]:
    g=d[(d.cell_type==ct)&(d.module==mod)]
    print("\n",name,ct,mod,"n",len(g))
    print(g.groupby(["dataset","disease_stage_coarse"]).size().unstack(fill_value=0))
# informative datasets per program (>=2 stage levels within dataset) in 60-donor universe
d60=d[d.dataset!="GSE189600"]
inf=d60.groupby(["cell_type","module","dataset"]).disease_stage_coarse.nunique().reset_index()
inf=inf[inf.disease_stage_coarse>=2]
k=inf.groupby(["cell_type","module"]).dataset.apply(lambda x: ",".join(sorted(x)))
print("\nInformative-dataset sets across 117 programs (60-donor universe):")
print(k.value_counts())
n60=d60.groupby(["cell_type","module"]).size()
print("n donors per program, 60-universe: min",n60.min(),"median",n60.median(),"max",n60.max())
# per program n in informative datasets
dd=d60.merge(inf[["cell_type","module","dataset"]],on=["cell_type","module","dataset"])
ni=dd.groupby(["cell_type","module","dataset"]).size().unstack(fill_value=0)
print(ni.describe().loc[["min","50%","max"]])
print("cell types:",d.groupby("cell_type").module.nunique().to_dict())
