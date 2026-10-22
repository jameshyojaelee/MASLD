import pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_rows',200)
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/'
g=pd.read_parquet(D+'passport_gene_index.parquet')
c=pd.read_csv(B+'coefficient_table.tsv.gz',sep='\t'); c['ens']=c.gene.str.split('.').str[0]
c['canon']=(c.padj<0.05)&(c.logFC.abs()>0.5)
miss=c[~c.ens.isin(g.ensembl_id)]
print(len(miss)); print(miss[miss.canon][['gene','symbol','logFC','padj']])
print(miss.symbol.value_counts().head(10))
print('missing with padj<0.05:', (miss.padj<0.05).sum())
q=pd.read_csv(D+'passport_identity_quarantine.tsv',sep='\t'); print(q.shape); print(q.columns.tolist()); print(q.head(3).T)
for col in q.columns:
    if q[col].astype(str).str.contains('ENSG').any(): print('col with ENSG',col)
