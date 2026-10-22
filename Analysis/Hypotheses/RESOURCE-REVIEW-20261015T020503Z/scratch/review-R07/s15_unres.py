import pandas as pd
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
e=pd.read_parquet(D+'passport_evidence_long.parquet'); g=pd.read_parquet(D+'passport_gene_index.parquet')
r=e[e.evidence_domain=='transcriptomics'].set_index('ensembl_id')
m=g.set_index('ensembl_id').join(r[['q_value','estimate']])
for cl in ['unresolved','genetically_anchored','untested']:
    x=m[m.primary_evidence_class==cl]; print(cl,len(x),'RNA q<0.05:',int((x.q_value<0.05).sum()))
