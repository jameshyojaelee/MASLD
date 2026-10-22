import pandas as pd
L=pd.read_csv("GWAS/finemapping/results/gwas_atac/variant_overlap_summary.csv",dtype=str,keep_default_na=False)
print("max_pip min/quantiles", pd.to_numeric(L.max_pip).describe().round(4).to_dict())
print("best_method", L.best_method.value_counts().to_dict())
D="GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/"
m=pd.read_csv(D+"credible_set_members_annotated.tsv",sep="\t",dtype=str,keep_default_na=False)
g=m[(m.pip_sum_gate_passed=="true")&(m.study.str.startswith("BBJ"))].copy()
g["pip"]=g.susie_pip.astype(float)
print(g.groupby(g.promoter_context_status).pip.describe().round(4))
cs=pd.read_csv("GWAS/finemapping/results/credible_sets.csv",dtype=str,keep_default_na=False,usecols=["study","susie_pip","susie_cs"])
print("credible_sets.csv studies", cs.study.nunique(), sorted(cs.study.unique())[:60])
