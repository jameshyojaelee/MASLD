import pandas as pd, numpy as np
r=pd.read_csv('bg_rows_check.tsv',sep='\t')
r['mx']=r[['gse244832_accessible_mass','gse281367_accessible_mass']].max(axis=1)
P=r.drop_duplicates(['gwas_name','ensembl','signal_pair_index'])
print(P.trait_class.value_counts().to_dict(), 'genes',P.groupby('trait_class').ensembl.nunique().to_dict(),
  'both classes', len(set(P[P.trait_class=='direct_MASLD'].ensembl)&set(P[P.trait_class=='liver_enzyme'].ensembl)),'unique genes',P.ensembl.nunique())
print('pairs per gene in direct and enzyme max:',P.groupby(['trait_class','ensembl']).size().groupby('trait_class').max().to_dict())
# pair-level: <1% in every lineage
pm=r.groupby(['trait_class','gwas_name','ensembl','signal_pair_index']).agg(maxany=('any_accessible_mass','max'),maxsrc=('mx','max')).reset_index()
for th in [0.01,0.05]:
    print('pairs maxany<',th, (pm.maxany<th).sum(), pm[pm.maxany<th].trait_class.value_counts().to_dict(), '| maxsrc', (pm.maxsrc<th).sum(), pm[pm.maxsrc<th].trait_class.value_counts().to_dict())
pa=r[r.evidence_state=='partial']; print('partial rows',len(pa),'<0.01 any',(pa.any_accessible_mass<0.01).sum(),'<0.01 maxsrc',(pa.mx<0.01).sum(), 'median',pa.any_accessible_mass.median())
print(r.evidence_state.value_counts().to_dict())
# gene level: every lineage, every pair < 5%
gm=r.groupby(['trait_class','ensembl']).agg(maxsrc=('mx','max'),maxany=('any_accessible_mass','max')).reset_index()
for th in [0.05]:
    print('genes all <5% (maxsrc):',gm[gm.maxsrc<th].groupby('trait_class').size().to_dict(),' (any):',gm[gm.maxany<th].groupby('trait_class').size().to_dict())
# gene-level best state per lineage, floor 0.05 on maxsrc
def st(x):
    if x.evidence_state=='replicated_accessible': return 0
    if x.evidence_state=='source_dependent': return 1
    return 2 if x.mx>=0.05 else 3
r['rk']=r.apply(st,axis=1)
g=r.groupby(['trait_class','lineage','ensembl']).rk.min().reset_index()
print(g.pivot_table(index=['trait_class','lineage'],columns='rk',values='ensembl',aggfunc='count',fill_value=0).to_string())
# background
r['exp']=r.exp_any
gb=r.groupby(['trait_class','lineage','ensembl']).agg(o=('any_accessible_mass','mean'),e=('exp','mean')).reset_index()
rng=np.random.default_rng(20260923)
for (tc,l),x in gb.groupby(['trait_class','lineage']):
    o,e=x.o.values,x.e.values; n=len(x)
    B=np.array([ (o[i].mean()/e[i].mean()) for i in (rng.integers(0,n,n) for _ in range(2000))])
    print(tc,l,n,'obs%.4f exp%.4f ratio %.3f CI %.3f-%.3f'%(o.mean(),e.mean(),o.mean()/e.mean(),*np.percentile(B,[2.5,97.5])),
          'n obs>exp',int((o>e).sum()),'n obs<exp',int((o<e).sum()),'median diff %.5f'%np.median(o-e))
