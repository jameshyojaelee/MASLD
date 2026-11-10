# Extract the 35 tier-1/2 strata rows (few columns) from COLOC r2 into a small local file
import pandas as pd
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
tier=pd.read_csv(R+"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t"); reg=pd.read_csv(R+"GWAS/finemapping/config/gwas_registry.tsv",sep="\t")
t12=tier[tier.tier.isin([1,2])].merge(reg[["study_name","ancestry"]],on="study_name")
x=pd.read_csv(R+"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
  usecols=["gwas_name","gene","ensembl","PP.H4.abf","PP.H4.susie","method","top_snp","ancestry","eligibility_source","n_cs_pairs","n_cs_pairs_eligible"],low_memory=False)
x=x[x.gwas_name.isin(t12.study_name)]
x=x.merge(t12[["study_name","trait","tier","ancestry"]].rename(columns={"study_name":"gwas_name","ancestry":"reg_anc"}),on="gwas_name")
x.to_csv("t12_coloc_r2.tsv.gz",sep="\t",index=False)
print(x.shape, x.gwas_name.nunique()); print(t12[["study_name","trait","ancestry"]].to_string())
