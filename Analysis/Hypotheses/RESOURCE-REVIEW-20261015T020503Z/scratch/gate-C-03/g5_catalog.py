import pandas as pd
C='figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
g=pd.read_parquet(C+'passport_gene_index.parquet'); ev=pd.read_parquet(C+'passport_evidence_long.parquet')
print('gene index',g.shape,'class',g.primary_evidence_class.value_counts().to_dict())
print('ev cols',list(ev.columns))
print('domains',ev.evidence_domain.value_counts().to_dict())
ge=ev[ev.evidence_domain=='genetics']
for c in ['direction','call_state']:
    if c in ge.columns: print(c,ge[c].value_counts(dropna=False).to_dict())
if 'next_experiment_rule' in g.columns: print('next_exp',g.next_experiment_rule.value_counts().head(8).to_dict())
sup=set(ge[ge.call_state=='supported'].ensembl_id)
d=pd.read_csv('Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-C-03/gene_dir_h26rule.tsv',sep='\t')
print('direction genes in supported',d.ensembl.isin(sup).sum(),'of',len(d))
x=g.assign(dir=g.ensembl_id.map(d.set_index('ensembl').dir))
print(pd.crosstab(x.primary_evidence_class,x.dir.fillna('none')))
