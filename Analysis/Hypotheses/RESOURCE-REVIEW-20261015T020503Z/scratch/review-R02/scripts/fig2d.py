import pandas as pd, numpy as np
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
D=B+'figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/'
c=pd.read_csv(D+'corrected_locus_coordinates.tsv',sep='\t')
print(c.shape, c.pip_sum_gate_passed.value_counts().to_dict())
ex=pd.read_csv(D+'excluded_credible_sets.tsv',sep='\t'); print('excluded',len(ex))
c=c[~c.credible_set_uid.isin(ex.credible_set_uid)]
c=c[c.pip_sum_gate_passed.astype(str).str.lower()=='true']
print('gated',c.groupby('trait_scope').size().to_dict())
c=c.sort_values(['chrom','pos']).copy()
cl=[]
for ch,b in c.groupby('chrom',sort=False):
    new=np.concatenate([[True],np.diff(b.pos.values)>1_000_000])
    cl.append(pd.Series([f'{ch}:{i}' for i in np.cumsum(new)],index=b.index))
c['cl']=pd.concat(cl).reindex(c.index)
c['num']=c.protein_altering_pip_mass+c.canonical_splice_pip_mass
c['den']=c.annotated_pip_mass
t=c.groupby(['trait_scope','cl'])[['num','den']].mean().reset_index()
for s,g in t.groupby('trait_scope'): print(s,len(g),g.num.sum()/g.den.sum())
fr=t.groupby('trait_scope').apply(lambda g:g.num.sum()/g.den.sum())
print('ratio',fr.iloc[0]/fr.iloc[1] if fr.index[0].startswith('tier1') else fr.iloc[1]/fr.iloc[0])
print('n physical clusters',c.cl.nunique(),'shared',(t.groupby('cl').trait_scope.nunique()==2).sum())
span=c.groupby('cl').agg(lo=('pos','min'),hi=('pos','max'),n=('pos','size'),nscope=('trait_scope','nunique'))
span['span_mb']=(span.hi-span.lo)/1e6
print(span.span_mb.describe())
print(span.sort_values('span_mb',ascending=False).head(12))
print('clusters spanning >5Mb:',(span.span_mb>5).sum(),' >10Mb:',(span.span_mb>10).sum())
# which direct clusters are shared and their spans
dcl=set(t[t.trait_scope.str.startswith('tier1')].cl)
print(span.loc[sorted(dcl)].sort_values('span_mb',ascending=False))
# equal-weight instance (no clustering) fractions for comparison
for s,g in c.groupby('trait_scope'): print('instance-level',s,g.num.sum()/g.den.sum())
