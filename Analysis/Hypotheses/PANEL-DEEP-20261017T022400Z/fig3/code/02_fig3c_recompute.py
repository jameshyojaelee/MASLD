"""Recompute Fig 3C (per-cohort disease-vs-control call membership and overlap classes)
from the five per-cohort limma fits (cohort_disease_all_gene_results.tsv)."""
import pandas as pd, numpy as np
O="Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
r=pd.read_csv("figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv",sep="\t")
print(r.groupby('dataset').agg(n=('gene_id_versioned','size'),nc=('n_control','first'),nd=('n_disease','first'),fam=('bh_family_size','first'),cov=('covariates','first')))
c=r[(r.FDR<0.05)&(r.logFC.abs()>0.5)].copy()
c['sgn']=np.sign(c.logFC)
print("calls per cohort:",c.dataset.value_counts().to_dict(), "memberships",len(c))
print("up/down per cohort:",c.groupby(['dataset','sgn']).size().to_dict())
g=c.groupby('gene_id_versioned').agg(k=('dataset','size'),nup=('sgn',lambda s:(s>0).sum()),ndn=('sgn',lambda s:(s<0).sum()))
def cls(row):
    if row.k==1: return 'unique'
    if row.nup==0 or row.ndn==0: return f'same{row.k}'
    return 'mixed'
g['cls']=g.apply(cls,axis=1)
print("union",len(g)); print(g.cls.value_counts().to_dict())
c=c.merge(g[['cls','k']],left_on='gene_id_versioned',right_index=True)
mat=pd.crosstab(c.dataset,c.cls,margins=True); print(mat); mat.to_csv(O+"fig3c_recomputed_cohort_by_class.tsv",sep="\t")
mixed=g[g.cls=='mixed'].index
inc126=c[(c.gene_id_versioned.isin(mixed))&(c.dataset=='GSE126848')].gene_id_versioned.nunique()
print("mixed including GSE126848:",inc126)
# without GSE126848
c2=c[c.dataset!='GSE126848']
g2=c2.groupby('gene_id_versioned').agg(k=('dataset','size'),nup=('sgn',lambda s:(s>0).sum()),ndn=('sgn',lambda s:(s<0).sum()))
multi=g2[g2.k>=2]; mx=((multi.nup>0)&(multi.ndn>0)).sum()
print("without GSE126848: multi",len(multi),"mixed",mx, round(mx/len(multi),4))
print("mixed by k:",g[g.cls=='mixed'].k.value_counts().sort_index().to_dict())
g.to_csv(O+"fig3c_recomputed_gene_classes.tsv",sep="\t")
# zero-read GSE126848 calls: logFC exactly 0.554?
z=c[(c.dataset=='GSE126848')]
print("GSE126848 calls with logFC rounded 0.554:", (z.logFC.round(3)==0.554).sum())
