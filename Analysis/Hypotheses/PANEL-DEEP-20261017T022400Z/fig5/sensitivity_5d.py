"""5D display sensitivity: uniform-posterior (locus-matched) reference per bar state, and 1% floor.
Input: F5bcd/fig5d_pair_rows_with_background.tsv (pair x lineage rows with n_mapped, n_any_peak, exp_any).
Observed states recomputed with own code from per-source masses in that table and cross-checked
against the replay output used in recompute_5d.py.
"""
import numpy as np, pandas as pd
B="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
O=B+"Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5/recompute/"
d=pd.read_csv(B+"figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_pair_rows_with_background.tsv",sep="\t")
d["cls"]=np.where(d.trait_class.str.contains("direct"),"direct","enzyme")
mx=np.maximum(d.gse244832_accessible_mass,d.gse281367_accessible_mass)
d["obs_ge05"]=mx>=0.05; d["obs_ge50"]=mx>=0.5
# uniform-posterior reference: posterior spread evenly over mapped tested variants -> mass in any-source peaks = exp_any
d["ref_ge05"]=d.exp_any>=0.05; d["ref_ge50"]=d.exp_any>=0.5
g=d.groupby(["cls","lineage","ensembl"]).agg(obs05=("obs_ge05","max"),obs50=("obs_ge50","max"),ref05=("ref_ge05","max"),ref50=("ref_ge50","max"),
    obs_mass=("any_accessible_mass","max"),exp_mass=("exp_any","max")).reset_index()
s=g.groupby(["cls","lineage"]).agg(n=("ensembl","size"),obs_ge5pct=("obs05","sum"),ref_ge5pct=("ref05","sum"),obs_ge50pct=("obs50","sum"),ref_ge50pct=("ref50","sum"),
    median_obs_mass=("obs_mass","median"),median_exp_mass=("exp_mass","median"))
print(s.to_string())
s.to_csv(O+"5d_sensitivity_uniform_reference.tsv",sep="\t")
# 1% floor alternative
d["st01"]=np.select([np.minimum(d.gse244832_accessible_mass,d.gse281367_accessible_mass)>=0.5,mx>=0.5,mx>=0.01],[3,2,1],0)
g1=d.groupby(["cls","lineage","ensembl"]).st01.max().reset_index()
print("1% floor gene-level counts (0 under1,1 partial,2 one,3 both)\n",g1.groupby(["cls","lineage"]).st01.value_counts().unstack(fill_value=0).to_string())
