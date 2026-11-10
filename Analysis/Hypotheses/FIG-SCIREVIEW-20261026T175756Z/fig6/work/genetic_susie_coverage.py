import pandas as pd
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/'
G='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/gen-corrected-r8-20261009T142027Z/genetics_context/frozen_evidence_classes.tsv'
f=pd.read_csv(G,sep='\t',low_memory=False)
print('raw rows',len(f)); print(f.static_class.value_counts())
print('primary_genetic true',(f.primary_genetic.astype(str).str.lower()=='true').sum())
pg=f[f.primary_genetic.astype(str).str.lower()=='true']
print('primary_genetic: susie>0.5',(pg.coloc_best_susie_pp4>0.5).sum(),' susie NA',pg.coloc_best_susie_pp4.isna().sum(), ' abf>0.5 only', ((pg.coloc_best_susie_pp4.isna())&(pg.coloc_best_abf_pp4>0.5)).sum())
e=pd.read_parquet(B+'passport_evidence_long.parquet')
g=pd.read_parquet(B+'passport_gene_index.parquet')
gen=e[e.evidence_domain=='genetics'][['ensembl_id','symbol','call_state','testability_reason','estimate','dataset_id']]
rna=e[e.evidence_domain=='transcriptomics'][['ensembl_id','call_state','estimate','q_value']].rename(columns={'call_state':'rna_call','estimate':'lfc','q_value':'q'})
m=gen.merge(rna,on='ensembl_id').merge(g[['ensembl_id','primary_evidence_class']],on='ensembl_id')
# join susie/abf from GEN_RAW on ensembl_genetic or bulk
f['key']=f.ensembl_genetic.fillna(f.ensembl_bulk)
fk=f.drop_duplicates('key').set_index('key')
for col in ['coloc_best_susie_pp4','coloc_best_abf_pp4','genetic_row_count','primary_genetic_map_tested']:
    m[col]=m.ensembl_id.map(fk[col])
# second try via bulk id
fb=f.dropna(subset=['ensembl_bulk']).drop_duplicates('ensembl_bulk').set_index('ensembl_bulk')
miss=m.coloc_best_susie_pp4.isna()&m.coloc_best_abf_pp4.isna()
for col in ['coloc_best_susie_pp4','coloc_best_abf_pp4','genetic_row_count']:
    m.loc[miss,col]=m.loc[miss,'ensembl_id'].map(fb[col])
gi=m[m.call_state=='indeterminate']
print('genetics indeterminate',len(gi))
print(' with finite SuSiE',gi.coloc_best_susie_pp4.notna().sum(),' no SuSiE',gi.coloc_best_susie_pp4.isna().sum())
print(' no SuSiE and ABF>0.5',(gi.coloc_best_susie_pp4.isna()&(gi.coloc_best_abf_pp4>0.5)).sum())
print(' no SuSiE and no ABF', (gi.coloc_best_susie_pp4.isna()&gi.coloc_best_abf_pp4.isna()).sum())
print(pd.crosstab(gi.coloc_best_susie_pp4.isna(), gi.primary_evidence_class))
x=gi[gi.coloc_best_susie_pp4.isna()&(gi.coloc_best_abf_pp4>0.5)]
print('ABF>0.5 no SuSiE by class'); print(x.primary_evidence_class.value_counts())
gs=m[m.call_state=='supported']; print('supported',len(gs),'susie>0.5',(gs.coloc_best_susie_pp4>0.5).sum(), 'min susie',gs.coloc_best_susie_pp4.min())
gu=m[m.call_state=='untestable']; print('untestable',len(gu), gu.testability_reason.value_counts().to_dict(), 'finite susie', gu.coloc_best_susie_pp4.notna().sum(), 'finite abf', gu.coloc_best_abf_pp4.notna().sum())
print(gu[gu.testability_reason=='insufficient_shared_posterior'][['symbol','coloc_best_susie_pp4','coloc_best_abf_pp4','rna_call','primary_evidence_class']])
m.to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig6/work/catalog_gene_calls_with_coloc.tsv',sep='\t',index=False)
