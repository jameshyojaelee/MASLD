import pandas as pd, numpy as np
from statsmodels.stats.multitest import multipletests
D="figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis"
sec="hotspot_hepatocytes_f05c535ae5bbc0b9"; duc="hotspot_hepatocytes_48f39dd4d817a10e"
for arm in ["new774_three_cohorts","new774_plus_GSE174478"]:
    me=pd.read_csv(f"{D}/{arm}/meta_effects.tsv",sep="\t")
    for mk,g in me.groupby("model_kind"):
        g=g[g.estimable==True]
        q=multipletests(g.p_value_meta,method="fdr_bh")[1]
        g=g.assign(qre=q)
        print(arm,mk,"family",len(g),"max|q-q_stored|",np.nanmax(np.abs(g.qre-g.q_value)))
        for ax in ["fibrosis","nas"]:
            h=g[g.axis==ax]
            print("  ",ax,"n_q<0.05=",int((h.qre<0.05).sum()),"min_q=%.4f"%h.qre.min(),"n_cohorts=",sorted(h.n_cohorts.unique()))
        for fid,name in [(sec,"StromalECM"),(duc,"Ductular")]:
            for ax in ["fibrosis","nas"]:
                r=g[(g.feature_id==fid)&(g.axis==ax)]
                print("    ",name,ax,"beta=%.4f p=%.4g q=%.4f"%(r.beta_meta.iloc[0],r.p_value_meta.iloc[0],r.qre.iloc[0]))
