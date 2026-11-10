import pandas as pd
W='/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig6/work/'
m=pd.read_csv(W+'catalog_gene_calls_with_coloc.tsv',sep='\t')
t=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/tier12/gene_level_coloc_tier12.csv')
u=m[m.testability_reason=='underpowered_source']
ids=set(t.ensembl)
print('9312 genes with a direct-ID tier12 record:',u.ensembl_id.isin(ids).sum())
sym=set(t.ensembl)  # tier12 has no symbol; check full gene_level
f=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc.csv',usecols=['gene','ensembl'])
hit=u[u.symbol.isin(set(f.gene.dropna()))]
print('9312 genes whose symbol has a COLOC gene-level record under another ID:',len(hit))
print(hit[['symbol','ensembl_id','rna_call','primary_evidence_class']].head(20).to_string())
print(hit.rna_call.value_counts())
