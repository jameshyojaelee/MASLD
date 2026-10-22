import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
F=R+"figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/"
reg=pd.read_csv(F+"program_registry_v2.tsv",sep="\t")
print("registry rows",len(reg))
for ct,mod in [("hepatocytes",20),("hepatocytes",8)]:
    r=reg[(reg.cell_type==ct)&(reg.module==mod)].iloc[0]
    print(ct,mod,r.module_name,"stab_med",r.stability_median,"beta",r.primary_beta,"q",r.primary_qvalue,"hc3q",r.primary_hc3_qvalue,"n",r.primary_n_donors,"maxlev",r.primary_max_leverage,"eq_beta",r.equal_run_beta)
    for c in ["primary_selected","robust_display","registry_state","source_stability"]:
        if c in reg.columns: print("  ",c,r[c])
print("stability_median>=0.5:",(reg.stability_median>=0.5).sum())
print("primary q<0.05:",(reg.primary_qvalue<0.05).sum(), reg.loc[reg.primary_qvalue<0.05,["cell_type","module","module_name"]].values.tolist())
print("lineage counts",reg.cell_type.value_counts().to_dict())
pe=pd.read_csv(F+"primary_effects.tsv",sep="\t")
print("pe Ductular", pe[(pe.cell_type=="hepatocytes")&(pe.module==20)][["beta","qvalue","hc3_qvalue","n_donors","max_leverage"]].to_dict("records"))
v=pd.read_csv(R+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R04-01-reproduce/refit_drop_gse189600_family.tsv",sep="\t")
print("refit_drop cols",list(v.columns)[:30],"rows",len(v))
qc=[c for c in v.columns if "q" in c.lower()]
print("q cols",qc)
for c in qc:
    print(c,"n<0.05",(v[c]<0.05).sum())
print(v[(v.cell_type=="hepatocytes")&(v.module.isin([20,8]))].T.to_string())
