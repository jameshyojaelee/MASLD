import pandas as pd, numpy as np
from scipy.stats import fisher_exact
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
tier=pd.read_csv(ROOT+"/GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t")
main=tier[tier.placement=="main"].set_index("study_name")
c=pd.read_csv(ROOT+"/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
  usecols=["gwas_name","gene","PP.H4.susie","n_cs_pairs","n_cs_pairs_eligible","eligibility_source","method"],low_memory=False)
c=c[c.gwas_name.isin(main.index)&c.gene.notna()].copy()
c["tier"]=c.gwas_name.map(main.tier)
print(c.groupby(["eligibility_source"])["PP.H4.susie"].describe().to_string())
print("rerun rows by tier & positive:\n", c[c.eligibility_source=="rerun_shared_posterior_check"].groupby("tier")["PP.H4.susie"].agg(n="size", pos=lambda x:(x>0.5).sum(), min="min").to_string())
t=c[c["PP.H4.susie"].notna()]
tab=t.groupby("tier").agg(rows=("PP.H4.susie","size"), pos=("PP.H4.susie",lambda x:(x>0.5).sum()), pairs=("n_cs_pairs","sum"))
tab["row_rate"]=tab.pos/tab.rows
print(tab.to_string())
a=tab.loc[1]; b=tab.loc[2]
print("Fisher rows tier1 vs tier2:", fisher_exact([[a.pos,a.rows-a.pos],[b.pos,b.rows-b.pos]]))
print("hit pairs / all tested CS pairs: tier1 30/%d=%.4f tier2 810/%d=%.4f"%(a.pairs,30/a.pairs,b.pairs,810/b.pairs))
print("Fisher pairs eligible 30/171 vs 810/3844:", fisher_exact([[30,141],[810,3034]]))
