import pandas as pd
pd.set_option('display.width',250); pd.set_option('display.max_rows',300)
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/'
g=pd.read_parquet(D+'passport_gene_index.parquet')
c=pd.read_csv(B+'coefficient_table.tsv.gz',sep='\t'); c['ens']=c.gene.str.split('.').str[0]
c['canon']=(c.padj<0.05)&(c.logFC.abs()>0.5)
c['symdup']=c.symbol.duplicated(keep=False)
c['inCat']=c.ens.isin(g.ensembl_id)
print(pd.crosstab(c.symdup,c.inCat))
# for duplicated-symbol groups: how many copies in catalog
grp=c[c.symdup].groupby('symbol').agg(n=('ens','size'),nin=('inCat','sum'))
print(grp.nin.value_counts())
print(grp[grp.nin>0].head(20))
# symbol in catalog for duplicated symbols that kept one copy: does the catalog symbol equal the source symbol?
x=c[c.symdup & c.inCat].merge(g[['ensembl_id','symbol']],left_on='ens',right_on='ensembl_id')
print(x[['gene','symbol_x','symbol_y','logFC','padj']].head(20))
# missing non-dup
print(c[(~c.symdup)&(~c.inCat)][['gene','symbol','logFC','padj']])
# among dropped dups, signif
d=c[c.symdup&~c.inCat]; print('dropped dup rows',len(d),'padj<0.05',(d.padj<0.05).sum(),'canon',d.canon.sum())
print(d[d.padj<0.05][['gene','symbol','logFC','padj']])
