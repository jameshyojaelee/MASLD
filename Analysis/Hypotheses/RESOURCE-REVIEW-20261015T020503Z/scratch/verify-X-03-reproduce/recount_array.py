import pandas as pd
D="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V4-20261008T015029Z"
a=pd.read_csv(f"{D}/per_array_results.tsv",sep="\t",dtype={"individual":str})
h=a[a.dataset=="HRA007511_HMSMA"].copy()
prim=h[~h.individual.isin(["27","28","29"])]
print("arrays primary:",len(prim),"unique ind",prim.individual.nunique())
for c in ["residual_I","centered_d","equal_centered","leave_top_centered"]:
    print(c,"positive:",(prim[c]>0).sum(),"/",len(prim))
print("residual_I > null_q95:",(prim.residual_I>prim.null_q95).sum())
print("spotperm observed > q95:",(prim.spotperm_observed>prim.spotperm_null_q95).sum())
print("matched_null_p <=0.05:",(prim.matched_null_p<=0.05).sum())
print(prim.sort_values("centered_d")[["individual","residual_I","null_mean","null_q95","centered_d","matched_null_p","equal_centered","leave_top_centered"]].head(4).to_string())
