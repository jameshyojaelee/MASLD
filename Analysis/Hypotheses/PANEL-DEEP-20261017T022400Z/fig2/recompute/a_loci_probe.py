import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/"
tier=pd.read_csv(R+"config/gwas_trait_tier.tsv",sep="\t")
reg=pd.read_csv(R+"config/gwas_registry.tsv",sep="\t")
t12=tier[tier.tier.isin([1,2])].merge(reg[["study_name","ancestry"]],on="study_name")
print(len(t12)); print(t12.groupby("ancestry").size()); print(t12.groupby("trait").size())
cs=pd.read_csv(R+"results/credible_sets.csv",usecols=["trait","locus","study","ancestry","susie_cs","susie_pip","chromosome","position"])
cs=cs[cs.study.isin(t12.study_name)]
print("studies in cs", cs.study.nunique(), sorted(set(t12.study_name)-set(cs.study)))
loc=cs[["study","trait","locus","ancestry"]].drop_duplicates()
print("study-locus", len(loc))
print(loc.groupby("trait").locus.nunique(), loc.groupby("trait").locus.nunique().sum(), loc.locus.nunique())
print(loc.groupby(["trait","ancestry"]).locus.nunique())
