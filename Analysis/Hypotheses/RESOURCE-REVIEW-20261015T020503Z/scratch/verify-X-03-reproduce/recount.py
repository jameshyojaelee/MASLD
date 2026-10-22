import pandas as pd, numpy as np
D="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V4-20261008T015029Z"
p=pd.read_csv(f"{D}/per_donor_results.tsv",sep="\t")
print(p.groupby(["dataset","mode","analysis_set"]).size())
for s in ["primary_32","all_35","primary_drop_top4"]:
    q=p[(p.dataset=="HRA007511_HMSMA")&(p.analysis_set==s)]
    print(s,"n",len(q),"n_unique_ind",q.individual.nunique(),
          "raw I>0:",(q.donor_residual_I>0).sum(),"d>0:",(q.donor_centered_d>0).sum(),
          "min raw I",round(q.donor_residual_I.min(),5),"min d",round(q.donor_centered_d.min(),5),
          "mean d",round(q.donor_centered_d.mean(),6))
q=p[(p.dataset=="HRA007511_HMSMA")&(p.analysis_set=="primary_32")].sort_values("donor_centered_d")
print(q[["individual","donor_residual_I","donor_null_mean","donor_centered_d"]].head(6).to_string())
print("individuals:",sorted(q.individual.tolist()))
n=len(q); print("exact one-sided p, all positive, n=32:",1/2**n, " Holm x2:",2/2**n)
