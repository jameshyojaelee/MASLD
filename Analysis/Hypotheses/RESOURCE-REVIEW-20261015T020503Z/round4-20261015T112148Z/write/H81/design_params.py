"""H81 design parameters from existing per-gene tables (no simulation outcome). Read-only."""
import numpy as np, pandas as pd
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z'
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv'
g=pd.read_csv(D+'/scratch/gate-L09-01/gate_per_gene.tsv.gz',sep='\t')
d=pd.read_csv(R)
print('deg rows',len(d),'calls',int(((d.padj<0.05)&(d.logFC.abs()>0.5)).sum()))
print('gate rows',len(g),'canonical',int(g.canonical.sum()),'tost05 q<0.05',int((g.q_tost05<0.05).sum()))
tau=np.sqrt(g.tau2.clip(lower=0))
for lab,m in [('all',np.ones(len(g),bool)),('canonical',g.canonical.astype(bool))]:
    print(lab,'tau quantiles 25/50/75/90:',np.round(np.quantile(tau[m],[.25,.5,.75,.9]),3).tolist())
d=d.sort_values('gene').reset_index(drop=True)
d['dec']=pd.qcut(d.AveExpr.rank(method='first'),10,labels=False)
print('decile n',d.dec.value_counts().sort_index().tolist())
print('decile AveExpr edges',np.round(d.groupby('dec').AveExpr.agg(['min','max']).values,3).tolist())
print('decile median SE',np.round(d.groupby('dec').SE.median().values,3).tolist())
