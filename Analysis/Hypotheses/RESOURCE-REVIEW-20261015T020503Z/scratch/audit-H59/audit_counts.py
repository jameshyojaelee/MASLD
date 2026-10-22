# Audit: independent n counts and score-file vs metadata stage agreement. Reads no score values except presence.
import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
H=R+"Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"
M=R+"Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
PE=R+"figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/primary_effects.tsv"
s=pd.read_csv(H,sep="\t",keep_default_na=False,usecols=["cell_type","module","donor","dataset","disease_stage_coarse","stage_ordinal","exclude_stage_analysis"])
m=pd.read_csv(M,sep="\t",keep_default_na=False)
print("metadata rows",len(m),"datasets",m.dataset.value_counts().to_dict())
print(m.groupby(["dataset","disease_stage_coarse","exclude_stage_analysis"]).size().to_string())
m=m.rename(columns={"sample":"donor","dataset":"ds_m","disease_stage_coarse":"stage_m","exclude_stage_analysis":"excl_m"})
d=s.merge(m[["donor","ds_m","stage_m","excl_m"]],on="donor",how="left",validate="many_to_one")
print("score rows",len(s),"unmatched in metadata",d.ds_m.isna().sum())
print("dataset mismatch",(d.dataset!=d.ds_m).sum())
dd=d.drop_duplicates("donor")
print("stage mismatch score vs meta (donor level):")
print(pd.crosstab(dd.disease_stage_coarse+"|"+dd.exclude_stage_analysis.astype(str), dd.stage_m+"|"+dd.excl_m.astype(str)).to_string())
keep=(d.excl_m.astype(str).str.upper()=="FALSE")&d.stage_m.isin(["Healthy","Steatosis","Steatohepatitis"])
u=d[keep]
print("62-universe donors",u.donor.nunique(),"programs",u.groupby(["cell_type","module"]).ngroups)
print(u.drop_duplicates("donor").groupby(["dataset","stage_m"]).size().unstack(fill_value=0).to_string())
print("GSE189600 donors:",sorted(u[u.dataset=="GSE189600"].drop_duplicates("donor")[["donor","stage_m"]].values.tolist()))
# compare per-program n to primary_effects
pe=pd.read_csv(PE,sep="\t")
n62=u.groupby(["cell_type","module"]).size().rename("n62").reset_index()
c=pe.merge(n62,on=["cell_type","module"])
print("n62 == primary_effects n_donors:",(c.n62==c.n_donors).all(), "rows",len(c))
u60=u[u.dataset!="GSE189600"]
print("60 universe donors",u60.donor.nunique(), "per program min/median/max",u60.groupby(["cell_type","module"]).size().agg(["min","median","max"]).to_dict())
inf=u60[u60.dataset.isin(["GSE202379","GSE244832"])]
n=inf.groupby(["cell_type","module"]).size()
print("fitted n by lineage:",n.groupby(level="cell_type").agg(["min","max"]).to_string())
pd_=inf.groupby(["cell_type","module","dataset"]).size().unstack()
print("per dataset min/max",pd_.min().to_dict(),pd_.max().to_dict())
# G1 check
cnt=inf.groupby(["cell_type","module","dataset","stage_m"]).size().unstack(fill_value=0)
ok=(cnt.sum(axis=1)>=8)&((cnt>=2).sum(axis=1)>=2)
g1=ok.groupby(level=["cell_type","module"]).agg(lambda x: x.all() and len(x)==2)
print("G1 pass",g1.sum(),"of",len(g1))
print("rows with a single-donor stage cell:",((cnt==1).any(axis=1)).sum())
print("rows with zero-donor stage cell:",((cnt==0).any(axis=1)).sum())
print("min cell by dataset/stage:",cnt.groupby(level="dataset").min().to_string())
# stage levels per dataset in fitted set
print("programs where a fitted dataset lacks a stage level:",((cnt==0).any(axis=1)).groupby(level=["cell_type","module"]).any().sum())
# smallest single-stage-level counts per dataset
for ct,mod in [("hepatocytes",20),("hepatocytes",8)]:
    print(ct,mod,cnt.loc[(ct,mod)].to_string())
# GSE244832 programs with only 2 levels etc
lv=(cnt>0).sum(axis=1)
print("stage levels per program-dataset:",lv.value_counts().to_dict())
