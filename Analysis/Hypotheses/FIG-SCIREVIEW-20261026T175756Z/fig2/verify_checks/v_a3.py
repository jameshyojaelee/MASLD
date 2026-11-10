import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
gl=pd.read_csv(R+"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
ms=gl[gl.coloc_best_susie_pp4>0.5]; ss=gl[(~(gl.coloc_best_susie_pp4>0.5))&(gl.coloc_best_abf_pp4>0.5)]
band=lambda v: pd.cut(v,[0.5,0.8,0.9,1.0001],right=True,labels=["0.5-0.8","0.8-0.9",">0.9"]).value_counts().to_dict()
print("multi",len(ms),band(ms.coloc_best_susie_pp4)); print("single",len(ss),band(ss.coloc_best_abf_pp4)); print(ss.susie_state_t12.value_counts().to_dict())
# if multi banded by max(susie,abf)
print("multi banded by max(SuSiE,ABF):",band(ms[["coloc_best_susie_pp4","coloc_best_abf_pp4"]].max(axis=1)))
