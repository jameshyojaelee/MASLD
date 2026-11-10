import pandas as pd
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/'
g=pd.read_parquet(B+'passport_gene_index.parquet')
e=pd.read_parquet(B+'passport_evidence_long.parquet')
print('genes',len(g),'unique ens',g.ensembl_id.nunique(),'evidence rows',len(e))
w=e.pivot_table(index='ensembl_id',columns='evidence_domain',values='call_state',aggfunc='first')
print(pd.crosstab(w['genetics'],w['transcriptomics'],margins=True))
print(g.primary_evidence_class.value_counts())
print(g.next_experiment_rule_id.value_counts())
m=g.set_index('ensembl_id').join(w)
print(pd.crosstab([m.genetics,m.transcriptomics],m.primary_evidence_class))
print(e.groupby(['evidence_domain','call_state','testability_reason']).size())
print('tested_negative rows', (e.call_state=='tested_negative').sum(), 'negative_call_passed nonnull', e.negative_call_passed.notna().sum())
