import pandas as pd
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
e=pd.read_parquet(C+'passport_evidence_long.parquet')
gen=e[e.evidence_domain=='genetics']
t=pd.read_csv(R+'gene_level_coloc_tier12.csv')
m=gen.merge(t[['ensembl','susie_state_t12','coloc_best_susie_pp4','coloc_best_abf_pp4']],left_on='ensembl_id',right_on='ensembl',how='left')
print(pd.crosstab(m.call_state,m.susie_state_t12.fillna('not_in_tier12')))
ind=m[m.call_state=='indeterminate']
print('indeterminate with max(SuSiE,ABF) PP4<0.1:', ((ind.coloc_best_susie_pp4.fillna(0)<0.1)&(ind.coloc_best_abf_pp4.fillna(0)<0.1)).sum(), 'of', len(ind))
