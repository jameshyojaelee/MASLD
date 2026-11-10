# Fig 2G/2H: re-derive printed lead p, PIPs, gene fills, and compare gene fills with the adopted bulk fit
import pandas as pd, numpy as np, glob
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
S=R+"figures/candidates/panel-deep-fixes-20261017T035633Z/fig2/source_tables/"
for p in ["2G","2H"]:
    h=pd.read_csv(S+p+".tsv",sep="\t",low_memory=False)
    gw=h[h.track=="gwas"]; print(p,"gwas n",len(gw),"min p",gw.p.min(),"at",gw.loc[gw.p.idxmin(),"position"],"lead_pos",gw.lead_pos.iloc[0],"coloc_top",gw.coloc_top.iloc[0], "rsid",gw.lead_rsid.iloc[0])
    print(" n variants p<5e-8:",(gw.p<5e-8).sum())
    pip=h[h.track=="pip"]
    if len(pip): 
        t=pip.sort_values("value",ascending=False).head(4); print(" top PIP:",t[[c for c in ["position","value","in_cs"] if c in t]].values.tolist(), "n in_cs", pip.in_cs.sum() if 'in_cs' in pip else '')
    eq=h[h.track=="eqtl"]; print(" eqtl min p",eq.p.min(),"at",eq.loc[eq.p.idxmin(),"position"])
    print(h[h.track=="gene"][["label","lfc","padj","state","colour"]].to_string())
    print(h[h.track=="text"][["label"]].to_string())
