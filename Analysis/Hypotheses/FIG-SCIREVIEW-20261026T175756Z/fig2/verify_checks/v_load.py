import pandas as pd, numpy as np, time
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
tier=pd.read_csv(R+"GWAS/finemapping/config/gwas_trait_tier.tsv",sep="\t")
print(tier.columns.tolist()); 
reg=pd.read_csv(R+"GWAS/finemapping/config/gwas_registry.tsv",sep="\t")
t12=tier[tier.tier.isin([1,2])].merge(reg[["study_name","ancestry"]],on="study_name")
print(len(t12)); print(t12.groupby(["trait","ancestry"]).size().to_dict())
t=time.time()
x=pd.read_csv(R+"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gwas_name","gene","ensembl","chr","PP.H4.abf","PP.H4.susie","method","top_snp","top_snp_PP","ancestry","n_cs_pairs","n_cs_pairs_eligible","eligibility_source"])
print("read",time.time()-t,len(x))
x=x[x.gwas_name.isin(t12.study_name)].merge(t12[["study_name","trait","tier","ancestry"]].rename(columns={"study_name":"gwas_name","ancestry":"reg_anc"}),on="gwas_name")
print(len(x), x.gwas_name.nunique())
print(x.method.value_counts(dropna=False).to_dict())
print(x.eligibility_source.value_counts(dropna=False).to_dict())
x.to_pickle("t12.pkl")
