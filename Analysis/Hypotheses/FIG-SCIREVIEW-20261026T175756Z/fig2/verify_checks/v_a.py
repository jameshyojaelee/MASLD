import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
tr=pd.read_csv(R+"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t"); reg=pd.read_csv(R+"GWAS/finemapping/config/gwas_registry.tsv",sep="\t")
t12=tr[tr.tier.isin([1,2])].merge(reg[["study_name","ancestry"]],on="study_name")
cs=pd.read_csv(R+"GWAS/finemapping/results/credible_sets.csv")
print(cs.columns.tolist()[:30]); print(len(cs))
cs12=cs[cs.study.isin(t12.study_name)]
nloc=cs12.groupby("study").locus.nunique()
print("t12 studies with any CS row:",cs12.study.nunique())
print("missing:",sorted(set(t12.study_name)-set(cs12.study)))
print(nloc.reindex(t12.study_name).fillna(0).astype(int).to_string())
