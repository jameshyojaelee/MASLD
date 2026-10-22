import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_columns',50)
e=pd.read_parquet('figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_evidence_long.parquet')
c=e[e.assay=='SuSiE_COLOC_or_ABF_fallback'].copy()
c['est']=pd.to_numeric(c['estimate'],errors='coerce').astype(float)
c['call_state']=c['call_state'].astype(str); c['symbol']=c['symbol'].astype(str); c['ensembl_id']=c['ensembl_id'].astype(str)
t=pd.read_csv('GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv')
m=c.merge(t,left_on='ensembl_id',right_on='ensembl',how='left',indicator=True)
s=m.coloc_best_susie_pp4.astype(float); a=m.coloc_best_abf_pp4.astype(float); est=m.est
m['src']=np.select([est.isna(), (est-s).abs()<1e-9, (est-a).abs()<1e-9],['na','susie','abf'],'other')
print(pd.crosstab(m.call_state,m.src))
ind=m[(m.call_state=='indeterminate')&(m.est>0.5)]
print('indet>0.5 n',len(ind)); print(ind.src.value_counts())
si=ind.coloc_best_susie_pp4
print('of these, finite susie <=0.5:',(si.notna()&(si<=0.5)).sum(), ' susie NaN:',si.isna().sum(), ' susie>0.5:', (si>0.5).sum())
print(ind.susie_state_t12.value_counts())
print('indet>=0.8 n',(ind.est>=0.8).sum())
sup=m[m.call_state=='supported']; print('supported src', sup.src.value_counts().to_dict(), 'susie>0.5', (sup.coloc_best_susie_pp4>0.5).sum())
allind=m[m.call_state=='indeterminate']
print('all indet src',allind.src.value_counts().to_dict())
print('all indet with finite susie but abf shown:',((allind.src=='abf')&allind.coloc_best_susie_pp4.notna()).sum())
print('abf>susie among finite-susie indet:', (allind.coloc_best_abf_pp4>allind.coloc_best_susie_pp4).sum())
# dataset_id: does it match susie_driving_gwas or abf_driving_gwas?
m['ds']=m.dataset_id.astype(str)
print('dataset matches susie_driving_gwas:', (m.ds==m.susie_driving_gwas).sum(), ' abf_driving_gwas:', (m.ds==m.abf_driving_gwas).sum(), ' driving_gwas:',(m.ds==m.driving_gwas).sum())
for g in ['GNMT','DGAT2','CFLAR','EZR']:
    print(m[m.symbol==g][['symbol','dataset_id','est','coloc_best_susie_pp4','coloc_best_abf_pp4','susie_driving_gwas','abf_driving_gwas','susie_state_t12','call_state','src']].to_string(header=(g=='GNMT')))
# rule: indeterminate rows with finite susie <=0.5 -- which value shown overall?
x=allind[allind.susie_state_t12=='tested_le_0.5']
print('tested_le_0.5 indet n',len(x),'src',x.src.value_counts().to_dict())
