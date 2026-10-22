"""H81 audit: independent recount of numbers stated in PRESPEC_H81.md. Read-only."""
import numpy as np, pandas as pd
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
R=ROOT+'/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls'
A=ROOT+'/Analysis/Hypotheses/refit-adoption-review-20261011T151011Z/out'
F=ROOT+'/figures/candidates/bulk-context-774-20261014T140904Z_22236311'
d=pd.read_csv(A+'/05_loco_per_call_detail.tsv',sep='\t')
print('per-call rows',len(d),'genes',d.gene.nunique(),'cohorts',sorted(d.cohort.unique()))
k=d.groupby('gene').kept.all()
print('kept in all five',int(k.sum()),'of',len(k),'-> SD-any flags',int((~k).sum()))
print('kept per cohort',d.groupby('cohort').kept.sum().to_dict())
c=pd.read_csv(F+'/cohort_disease_all_gene_results.tsv',sep='\t')
print('cohort rows',len(c),'genes',c.gene_id_versioned.nunique(),'columns',list(c.columns)[:20])
c['SE']=c.logFC/c.t
print('median per-cohort SE',c.groupby('dataset').SE.median().round(3).to_dict())
print(open(F+'/cohort_disease_design_audit.tsv').read())
g=pd.read_csv(R+'/deg_results.csv'); print('deg cols',list(g.columns)); print('gene id example',g.gene.head(3).tolist(),'all ENSG.version',bool(g.gene.str.match(r'^ENSG\d+\.\d+').all()))
print('calls',int(((g.padj<0.05)&(g.logFC.abs()>0.5)).sum()),'frac',round(((g.padj<0.05)&(g.logFC.abs()>0.5)).mean(),4),'2800/24115',round(2800/24115,4),'3200/24115',round(3200/24115,4))
w=pd.read_csv(R+'/lvqw_sample_weights.tsv',sep='\t'); print('weights cols',list(w.columns),len(w))
import pickle
