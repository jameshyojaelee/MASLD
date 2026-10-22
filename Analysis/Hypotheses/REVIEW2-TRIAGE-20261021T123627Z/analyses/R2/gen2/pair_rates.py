# Per-trait-class SuSiE pair rates before/after the AN1 symmetric check (descriptive).
import importlib.util, pandas as pd, os
ROOT=os.getcwd()
D="Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/gene_level_20261021T165252Z/"
spec=importlib.util.spec_from_file_location("rt","GWAS/finemapping/src/perf/rebuild_tier12.py")
rt=importlib.util.module_from_spec(spec); spec.loader.exec_module(rt)
tier,trait,placement,label=rt.load_tiers()
for nm in ["baseline","sensitivity"]:
    m=pd.read_csv(D+f"master_reduced_{nm}.csv",usecols=["gwas_name","ensembl","PP.H4.susie"],float_precision="round_trip")
    m["tier"]=m.gwas_name.map(tier)
    f=m[m["PP.H4.susie"].notna() & m.tier.isin([1,2])]
    for t in [1,2]:
        x=f[f.tier==t]; pos=(x["PP.H4.susie"]>0.5).sum()
        print(nm,"tier",t,"finite",len(x),"pos",pos,"rate %.4f"%(pos/len(x)))
    print(nm,"finite t12 rows",len(f),"genes",f.ensembl.nunique())
