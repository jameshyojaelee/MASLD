# H80: sampling-frame counts for the no-result stratum (genes with converged eQTL fit, 0 CS).
import pandas as pd
FM = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
R2 = f"{FM}/results/susie_coloc_r2_20260924T053244Z"
W = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H80"
tier = pd.read_csv(f"{FM}/config/gwas_trait_tier.tsv", sep="\t")
t12 = set(tier.loc[tier.tier.isin([1, 2]), "study_name"])
reg = pd.read_csv(f"{FM}/config/gwas_registry.tsv", sep="\t")
g = pd.read_csv(f"{W}/gene_state_by_polyfun_fit.tsv", sep="\t")
nores = set(g.loc[(g.susie_state_t12 == "no_susie_result") & (g.fit == "cs0"), "ensembl"])
m = pd.read_csv(f"{R2}/susie_coloc_all_gwas.csv", usecols=["gwas_name", "ensembl", "method"], low_memory=False)
m = m[m.gwas_name.isin(t12) & m.ensembl.isin(nores)]
m = m.merge(reg[["study_name", "ancestry", "trait_type"]], left_on="gwas_name", right_on="study_name")
print("no-result cs0 genes", len(nores), "with a tier-1/2 pair", m.ensembl.nunique())
print("pairs by method\n", m.method.value_counts().to_string())
w = m[["ensembl", "ancestry"]].drop_duplicates()
print("gene x ancestry windows", len(w)); print(w.ancestry.value_counts().to_string())
print("tier-1/2 studies by ancestry and type\n", reg[reg.study_name.isin(t12)].groupby(["ancestry", "trait_type"]).size().to_string())
