import pandas as pd, numpy as np, time
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
t=time.time()
ct=pd.read_csv(R+"GWAS/finemapping/results/combined_finemapping_coloc_targeted.csv",usecols=["chromosome","position","study","locus","susie_pip","susie_cs","recommended_pip"])
print("read",time.time()-t, len(ct))
print(ct.study.value_counts().to_dict())
ct.to_pickle("ct.pkl")
b=pd.read_csv(R+"figures/candidates/panel-deep-fixes-20261017T035633Z/fig2/source_tables/2B.tsv",sep="\t")
print(b.gwas.value_counts().to_dict())
