import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
S="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-G2"
s=pd.read_csv(S+"/dryrun/stage_roster_with_complexity.tsv",sep="\t")
rows=[]
for ref,cmp in [("F0","F1"),("F0","F2"),("F0","F3"),("F0","F4"),("F1","F2"),("F2","F3"),("F3","F4")]:
    for ds in sorted(s.dataset.unique()):
        x=s[(s.dataset==ds)&s.fib_group.isin([ref,cmp])]
        for v in ["z_loglib","z_gini","z_nf"]:
            a=x.loc[x.fib_group==ref,v]; b=x.loc[x.fib_group==cmp,v]
            p=mannwhitneyu(b,a).pvalue if len(a)>1 and len(b)>1 else np.nan
            rows.append(dict(contrast=f"{cmp}_vs_{ref}",dataset=ds,metric=v,n_ref=len(a),n_cmp=len(b),mean_diff_cmp_minus_ref=round(b.mean()-a.mean(),3),mwu_p=p))
r=pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
r["bh_q_all_63"]=multipletests(r.mwu_p,method="fdr_bh")[1]
r.to_csv(S+"/covariate_balance_by_contrast_cohort.tsv",sep="\t",index=False)
pd.set_option("display.width",200)
print(r[r.contrast.isin(["F1_vs_F0","F4_vs_F0"])].to_string(index=False))
print(r[r.bh_q_all_63<0.05].to_string(index=False))
