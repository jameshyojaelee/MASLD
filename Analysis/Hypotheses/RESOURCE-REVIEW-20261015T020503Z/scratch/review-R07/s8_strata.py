import pandas as pd, numpy as np
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/'
e=pd.read_parquet(D+'passport_evidence_long.parquet')
w=e.pivot_table(index='ensembl_id',columns='evidence_domain',values='call_state',aggfunc='first')
J=w[w.genetics.isin(['supported','indeterminate'])&w.transcriptomics.isin(['supported','indeterminate'])].copy()
J['G']=J.genetics=='supported'; J['R']=J.transcriptomics=='supported'
c=pd.read_csv(B+'coefficient_table.tsv.gz',sep='\t'); c['ens']=c.gene.str.split('.').str[0]
J=J.join(c.set_index('ens')[['AveExpr','SE']])
print('mean AveExpr G %.3f notG %.3f ; R %.3f notR %.3f'%(J[J.G].AveExpr.mean(),J[~J.G].AveExpr.mean(),J[J.R].AveExpr.mean(),J[~J.R].AveExpr.mean()))
for k in [5,10,20,50]:
    J['q']=pd.qcut(J.AveExpr,k,labels=False)
    s=J.groupby('q').apply(lambda d: d.G.sum()*d.R.sum()/len(d)).sum()
    print('bins',k,'stratified expected %.2f'%s)
# stratify on SE (precision) too
J['q']=pd.qcut(J.SE,10,labels=False); print('SE deciles expected %.2f'%J.groupby('q').apply(lambda d: d.G.sum()*d.R.sum()/len(d)).sum())
