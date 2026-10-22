import pandas as pd
pd.set_option('display.width',250)
C='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
d=pd.read_parquet(C)
print(d.shape); print(d.columns.tolist())
for c in d.columns:
    if d[c].dtype==object and d[c].nunique()<20: print(c, d[c].value_counts(dropna=False).to_dict())
