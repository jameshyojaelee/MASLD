import pandas as pd, numpy as np
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',200); pd.set_option('display.max_rows',200)
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/'
e=pd.read_parquet(D+'passport_evidence_long.parquet')
c=pd.read_csv(B+'coefficient_table.tsv.gz',sep='\t')
c['ens']=c.gene.str.split('.').str[0]
c['canon']=(c.padj<0.05)&(c.logFC.abs()>0.5)
print('bulk genes',len(c),'unique ens',c.ens.nunique(),'canon',c.canon.sum())
dup=c[c.ens.duplicated(keep=False)]; print('dup ens rows',len(dup)); print(dup.head(20))
r=e[e.evidence_domain=='transcriptomics'].copy()
m=r.merge(c[['gene','ens','symbol','logFC','padj','canon']],left_on='ensembl_id',right_on='ens',how='outer',indicator=True)
print(m._merge.value_counts())
print(pd.crosstab(m.call_state.fillna('ABSENT'),m.canon.fillna('ABSENT')))
x=m[(m.canon==True)&(m.call_state!='supported')]; print(x[['ensembl_id','symbol_x','symbol_y','ens','logFC','padj','call_state','testability_state','estimate','q_value']])
y=m[(m.canon!=True)&(m.call_state=='supported')]; print(y[['ensembl_id','symbol_x','symbol_y','logFC','padj','estimate','q_value']].head())
# estimates vs bulk
mm=m[m._merge=='both']
print('max |estimate-logFC|',(mm.estimate.astype(float)-mm.logFC).abs().max(),'max |q-padj|',(mm.q_value.astype(float)-mm.padj).abs().max())
print('testable rows in catalog not in bulk:', ((m._merge=='left_only')&(m.testability_state=='testable')).sum())
print(m[m._merge=='left_only'].testability_state.value_counts())
print(m[m._merge=='right_only'][['gene','symbol_y','logFC','padj','canon']].head(30)); print('right_only canon', m[m._merge=='right_only'].canon.sum())
# symbol mismatches
sm=mm[mm.symbol_x!=mm.symbol_y]; print('symbol mismatches',len(sm)); print(sm[['ensembl_id','symbol_x','symbol_y','call_state']].head(20))
