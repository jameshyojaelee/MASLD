import pandas as pd, numpy as np
pd.set_option('display.width',300); pd.set_option('display.max_columns',30); pd.set_option('display.max_colwidth',120)
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/'
e=pd.read_parquet(D+'passport_evidence_long.parquet'); g=pd.read_parquet(D+'passport_gene_index.parquet')
r=e[e.evidence_domain=='transcriptomics']; ge=e[e.evidence_domain=='genetics']
print('RNA indeterminate with q<0.05:', ((r.call_state=='indeterminate')&(r.q_value<0.05)).sum())
c=pd.read_csv(B+'coefficient_table.tsv.gz',sep='\t'); print('bulk fit padj<0.05 & |lfc|<=0.5:', ((c.padj<0.05)&(c.logFC.abs()<=0.5)).sum())
print('genetics indeterminate with estimate>0.5:', ((ge.call_state=='indeterminate')&(ge.estimate>0.5)).sum())
print('genetics effect_unit', ge.effect_unit.value_counts().to_dict(), 'assay', ge.assay.value_counts().to_dict())
print('genetics allowed_wording by call', ge.groupby(['call_state','testability_state']).allowed_wording.unique().to_dict())
print('RNA allowed_wording by call', r.groupby(['call_state','testability_state']).allowed_wording.unique().to_dict())
print('source_dependent', e.source_dependent.value_counts().to_dict(), 'provenance', e.provenance_state.value_counts().to_dict())
print('direction RNA', r.direction.value_counts().to_dict())
# untested composition
w=e.pivot_table(index='ensembl_id',columns='evidence_domain',values='testability_state',aggfunc='first')
m=g.set_index('ensembl_id').join(w)
u=m[m.primary_evidence_class=='untested']
print(u.groupby(['genetics','transcriptomics']).size())
print(u.groupby(['genetics','transcriptomics']).allowed_wording.unique())
# 10 genetics-supported but RNA untestable
cs=e.pivot_table(index='ensembl_id',columns='evidence_domain',values='call_state',aggfunc='first')
x=cs[(cs.genetics=='supported')&(cs.transcriptomics=='untestable')].join(g.set_index('ensembl_id').symbol)
print(x)
# top RNA-supported genes classed untested by |lfc|
y=cs[(cs.genetics=='untestable')&(cs.transcriptomics=='supported')].index
ry=r[r.ensembl_id.isin(y)].copy(); ry['a']=ry.estimate.abs()
print(ry.sort_values('q_value').head(15)[['symbol','estimate','q_value']])
