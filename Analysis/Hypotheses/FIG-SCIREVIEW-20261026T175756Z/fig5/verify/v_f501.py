# Independent check of F5-01 (5D reference). Verifier code, not reviewer's.
import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
d=pd.read_csv(R+"figures/candidates/review-r3-fixes-20261015T112148Z/F5bcd/fig5d_pair_rows_with_background.tsv",sep="\t")
print("rows",len(d),"pairs",d[["gwas_name","ensembl","signal_pair_index"]].drop_duplicates().shape[0])
print("lineages",d.lineage.value_counts().to_dict())
# exp_any = n_any_peak/n_mapped ?
print("max |exp_any - n_any_peak/n_mapped|", (d.exp_any-d.n_any_peak/d.n_mapped).abs().max())
print("exp_any quantiles", d.exp_any.quantile([0,.5,.9,.99,1]).round(4).to_dict())
print("pairs with exp_any>=0.5:", (d.exp_any>=0.5).sum())
print("n_mapped quantiles", d.n_mapped.quantile([0,.05,.5,1]).to_dict())
# Observed any-source flag
d["o_any"]=d[["gse244832_accessible_mass","gse281367_accessible_mass"]].max(axis=1)>=0.5
d["o_both"]=d[["gse244832_accessible_mass","gse281367_accessible_mass"]].min(axis=1)>=0.5
d["p"]=d.n_any_peak/d.n_mapped
rows=[]
for (tc,lin),g in d.groupby(["trait_class","lineage"]):
    gg=g.groupby("ensembl")
    obs=gg.o_any.any()
    p_indep=gg.p.apply(lambda x:1-np.prod(1-x.values))
    p_max=gg.p.max()
    p_mean=gg.p.mean()
    rows.append(dict(trait_class=tc,lineage=lin,n=len(obs),obs_any=100*obs.mean(),obs_both=100*gg.o_both.any().mean(),
      ref_indep=100*p_indep.mean(),ref_max=100*p_max.mean(),ref_meanp=100*p_mean.mean(),
      ref_markov_upper=100*np.minimum(1,2*p_max).mean()))
S=pd.DataFrame(rows); print(S.round(2).to_string())
S.to_csv(R+"Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig5/verify/v_f501_refs.tsv",sep="\t",index=False)
# Panel table
P=pd.read_csv(R+"figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5D.tsv",sep="\t")
print(P[["trait_class","lineage","measure","n_genes_meeting","n_genes","pct","ci95_lo_pct","ci95_hi_pct"]].to_string())
