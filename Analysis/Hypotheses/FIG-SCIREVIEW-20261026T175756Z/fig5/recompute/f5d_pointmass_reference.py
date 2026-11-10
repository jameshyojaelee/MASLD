"""5D: observed % genes with >=50% posterior in peaks (any source / both) vs two references:
uniform posterior over tested variants (drawn on panel) and a point-mass posterior at one random tested variant
(probability = fraction of mapped tested variants in peaks). Gene = best over its study pairs."""
import pandas as pd, numpy as np
d=pd.read_csv("figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_pair_rows_with_background.tsv",sep="\t")
mx=np.maximum(d.gse244832_accessible_mass,d.gse281367_accessible_mass)
mn=np.minimum(d.gse244832_accessible_mass,d.gse281367_accessible_mass)
d["obs_any"]=mx>=0.5; d["obs_both"]=mn>=0.5
d["unif_any"]=d.exp_any>=0.5
d["pm_any"]=d.n_any_peak/d.n_mapped; d["pm_both"]=d.n_shared_peak/d.n_mapped
def agg(g):
    return pd.Series(dict(obs_any=g.obs_any.any(),obs_both=g.obs_both.any(),unif_any=g.unif_any.any(),
        pm_any=1-np.prod(1-g.pm_any),pm_both=1-np.prod(1-g.pm_both)))
G=d.groupby(["trait_class","lineage","ensembl"]).apply(agg).reset_index()
S=G.groupby(["trait_class","lineage"]).agg(n=("ensembl","size"),obs_any_pct=("obs_any","mean"),obs_both_pct=("obs_both","mean"),
   unif_any_pct=("unif_any","mean"),pointmass_any_pct=("pm_any","mean"),pointmass_both_pct=("pm_both","mean"))
S.iloc[:,1:]*=100
print(S.round(2).to_string())
S.to_csv("Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig5/recompute/f5d_pointmass_reference.tsv",sep="\t")
