import pandas as pd
pd.set_option('display.width',400); pd.set_option('display.max_columns',30); pd.set_option('display.max_rows',300); pd.set_option('display.max_colwidth',80)
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
F='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/gen-corrected-r8-20261009T142027Z/genetics_context/frozen_evidence_classes.tsv'
g=pd.read_parquet(D+'passport_gene_index.parquet'); e=pd.read_parquet(D+'passport_evidence_long.parquet')
s=pd.read_csv(F,sep='\t',dtype=str,keep_default_na=False)
genes=['GLP1R','HKDC1','MTARC1','THRB','CAST','LPAL2','PNPLA3','TM6SF2','GCKR','HSD17B13','APOE','MBOAT7','GNMT','PDE8B']
print(g[g.symbol.isin(genes)][['ensembl_id','symbol','primary_evidence_class','role_hypothesis','next_experiment_rule_id']])
print(e[e.symbol.isin(genes)][['symbol','evidence_domain','dataset_id','phenotype','estimate','q_value','direction','testability_state','call_state']].sort_values(['symbol','evidence_domain']))
print(s[s.gene_symbol.isin(genes)][['gene_symbol','ensembl_bulk','ensembl_genetic','static_class','coloc_best_susie_pp4','coloc_best_abf_pp4','driving_gwas','driving_trait','bulk_logFC','bulk_padj']])
