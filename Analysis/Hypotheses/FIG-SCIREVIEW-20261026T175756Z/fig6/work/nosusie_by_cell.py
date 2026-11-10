import pandas as pd
W='/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig6/work/'
m=pd.read_csv(W+'catalog_gene_calls_with_coloc.tsv',sep='\t')
d=pd.read_parquet('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet',columns=['ensembl_id','coding_cs_state'])
m=m.merge(d,on='ensembl_id')
m['susie_state']=m.coloc_best_susie_pp4.notna().map({True:'finite SuSiE',False:'no SuSiE result'})
m.loc[m.call_state=='untestable','susie_state']='untestable'
print(pd.crosstab([m.call_state,m.susie_state],m.rna_call,margins=True))
c=m[m.coding_cs_state=='fine_mapped_coding_pip_ge_0.5']
print(pd.crosstab([c.call_state,c.susie_state],c.rna_call))
print(c[(c.call_state=='indeterminate')][['symbol','susie_state','coloc_best_susie_pp4','coloc_best_abf_pp4','rna_call']].sort_values('susie_state').to_string())
