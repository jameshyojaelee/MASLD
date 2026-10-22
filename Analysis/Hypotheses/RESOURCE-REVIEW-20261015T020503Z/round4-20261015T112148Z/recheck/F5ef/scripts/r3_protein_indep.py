import pandas as pd, numpy as np
from scipy.stats import rankdata, pearsonr
from statsmodels.stats.multitest import multipletests
P="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
RC=P+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F5ef/"
ph=pd.read_csv(P+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef/out/protein_histology_117.tsv",sep="\t")
reg=pd.read_csv(P+"figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv",sep="\t")
print("registry file order == sorted(cell_type,module):", list(reg.program_uid)==list(reg.sort_values(["cell_type","module"]).program_uid))
t=ph[ph.testable].copy()
print("testable",len(t),"n values",t.n.unique())
q=multipletests(t.pvalue,method="fdr_bh")[1]
print("BH max diff",np.abs(q-t.padj).max(), "q<.05",(q<.05).sum())
rob=(q<.05)&t.sensitivity_sign_agree&t.batch_sign_agree
# recompute sign agree
sa=(np.sign(t.rho)==np.sign(t.equal_rho))&(np.sign(t.rho)==np.sign(t.leave_top_rho))
ba=(np.sign(t.rho)==np.sign(t.rho_2019))&(np.sign(t.rho)==np.sign(t.rho_2020))
rob2=(q<.05)&sa&ba
print("robust",rob2.sum(),"stored",t.robust.sum(),"neg",(rob2&(t.rho<0)).sum(),"pos",(rob2&(t.rho>0)).sum())
print(t.loc[rob2,["cell_type","module","module_name","rho","ci_low","ci_high","padj"]].sort_values("rho").to_string())
print("seeds: perm-20260923 == position:", (t.perm_seed-20260923==np.arange(1,len(t)+1)).all(), (t.boot_seed-30260923==np.arange(1,len(t)+1)).all())
print("bootstrap finite min",t.n_boot_finite.min())
# CI contains rho?
print("rho outside CI:",((t.rho<t.ci_low)|(t.rho>t.ci_high)).sum())

# independent partial Spearman from stored per-sample scores and stored metadata
sc=pd.read_csv(P+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef/out/protein_program_scores_117.tsv",sep="\t")
pm=pd.read_csv(P+"Analysis/Multimodal_Program_Projection/results/proteomics/panel4c_metadata.tsv",sep="\t")
comp=["Steatosis","Ballooning","Inflammation","Fibrosis"]
mm=pm[(pm.group=="MASLD")&pm[comp].notna().all(1)].copy()
print("MASLD complete histology n",len(mm),"of",len(pm),"samples; MASLD total",(pm.group=="MASLD").sum())
Z=(mm[comp]-mm[comp].mean())/mm[comp].std(ddof=1)
U,S,Vt=np.linalg.svd(Z.to_numpy(),full_matrices=False)
load=Vt[0]; pc=Z.to_numpy()@load
if load.sum()<0: load=-load; pc=-pc
print("PC1 loadings",dict(zip(comp,np.round(load,3))),"var expl",round(S[0]**2/(S**2).sum(),3))
mm["burden"]=pc
def nuis(d):
    M=pd.DataFrame({"i":1.0,"b":(d.acquisition_batch.astype(str)==str(sorted(d.acquisition_batch.astype(str).unique())[-1])).astype(float),"age":d.age,"bmi":d.bmi})
    sx=pd.get_dummies(d.sex.astype(str),drop_first=True).astype(float); return np.column_stack([M.to_numpy(),sx.to_numpy()])
def resid(y,M):
    r=rankdata(y); beta,*_=np.linalg.lstsq(M,r,rcond=None); return r-M@beta
out=[]
for pid,g in sc[sc.group=="MASLD"].groupby("program_id"):
    d=mm.merge(g[["sample_id","score"]],on="sample_id")
    M=nuis(d); rho=np.corrcoef(resid(d.score,M),resid(d.burden,M))[0,1]
    out.append((pid,len(d),rho))
o=pd.DataFrame(out,columns=["program_uid","n_mine","rho_mine"]).merge(t[["program_uid","n","rho"]],on="program_uid")
print("programs",len(o),"max|rho diff|",np.abs(o.rho_mine-o.rho).max(),"n",o.n_mine.unique())
