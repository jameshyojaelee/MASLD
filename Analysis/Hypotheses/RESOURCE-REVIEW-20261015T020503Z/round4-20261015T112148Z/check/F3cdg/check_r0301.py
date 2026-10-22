import pandas as pd, numpy as np
d=pd.read_csv('figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv',sep='\t',usecols=['gene_id_versioned','dataset','logFC','FDR','bh_family_size'])
print(d.groupby('dataset').size().to_dict(), d.bh_family_size.unique())
c=d[(d.FDR<0.05)&(d.logFC.abs()>0.5)]
print('calls per cohort',c.groupby('dataset').size().to_dict(),'memberships',len(c))
s=np.sign(c.logFC)
g=c.assign(s=s).groupby('gene_id_versioned').s.agg(['size',lambda x:(x>0).sum()])
g.columns=['n','up']; g['dn']=g.n-g.up
print('union',len(g),'by n',g.n.value_counts().sort_index().to_dict())
m=g[g.n>1]; mixed=(m.up>0)&(m.dn>0)
print('multi',len(m),'mixed',mixed.sum(),round(mixed.mean(),4))
print('same by n',m[~mixed].n.value_counts().sort_index().to_dict(),'mixed by n',m[mixed].n.value_counts().sort_index().to_dict())
mx=set(m[mixed].index)
print('mixed with 126848',c[(c.dataset=='GSE126848')&c.gene_id_versioned.isin(mx)].gene_id_versioned.nunique())
c2=c[c.dataset!='GSE126848']; g2=c2.assign(s=np.sign(c2.logFC)).groupby('gene_id_versioned').s.agg(['size',lambda x:(x>0).sum()]); g2.columns=['n','up']
m2=g2[g2.n>1]; mx2=((m2.up>0)&(m2.up<m2.n)).sum(); print('without126848 multi',len(m2),'mixed',mx2,round(mx2/len(m2),4))
g['mx']=g[['up','dn']].max(axis=1); print('max-sign tiers among multi', g[g.n>1].mx.value_counts().sort_index().to_dict())
# ribbons per cohort/tier
g['tier']=np.where(g.n==1,'unique',np.where((g.up>0)&(g.dn>0),'mixed',g.n.astype(str)))
print(c.merge(g[['tier']],left_on='gene_id_versioned',right_index=True).groupby(['dataset','tier']).size().unstack())
