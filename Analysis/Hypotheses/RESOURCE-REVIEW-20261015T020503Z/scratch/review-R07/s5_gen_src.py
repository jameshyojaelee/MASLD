import pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_rows',300)
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/gen-corrected-r8-20261009T142027Z/genetics_context/frozen_evidence_classes.tsv'
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
s=pd.read_csv(F,sep='\t',dtype=str,keep_default_na=False)
g=pd.read_parquet(D+'passport_gene_index.parquet')
print(len(s), s.gene_symbol.nunique())
print(s.static_class.value_counts())
print('primary_genetic true',(s.primary_genetic=='true').sum(),'state true',(s.established_state_associated=='true').sum())
print('genetic_untestable_shared_posterior true',(s.genetic_untestable_shared_posterior=='true').sum())
bulk_ids=set(x for v in s.ensembl_bulk for x in v.split(';') if x)
gen_ids=set(x for v in s.ensembl_genetic for x in v.split(';') if x)
cat=set(g.ensembl_id)
print('bulk ids',len(bulk_ids),'missing from catalog',len(bulk_ids-cat))
print('genetic ids',len(gen_ids),'missing from catalog',len(gen_ids-cat))
multi_gen=s[s.ensembl_genetic.str.contains(';')]; print('symbols with multi genetic IDs',len(multi_gen))
# genetic primary in source vs catalog
src_pg=s[s.primary_genetic=='true']; print('source primary genetic symbols',len(src_pg), 'in catalog by symbol', src_pg.gene_symbol.isin(g.symbol).sum())
src_st=s[s.established_state_associated=='true']; print('source state symbols',len(src_st),'in catalog', src_st.gene_symbol.isin(g.symbol).sum())
print(src_st[~src_st.gene_symbol.isin(g.symbol)][['gene_symbol','ensembl_bulk','bulk_logFC','bulk_padj','bulk_row_count']])
# multi-ID bulk symbols kept in catalog: which bulk value?
mb=s[s.ensembl_bulk.str.contains(';') & s.gene_symbol.isin(g.symbol)]
print('multi-bulk-ID symbols kept in catalog', len(mb)); print(mb[['gene_symbol','ensembl_bulk','ensembl_genetic','bulk_padj','static_class']])
# genetic tested but not bulk -> catalog?
print('genetic-tested IDs missing from catalog examples', list(gen_ids-cat)[:20])
