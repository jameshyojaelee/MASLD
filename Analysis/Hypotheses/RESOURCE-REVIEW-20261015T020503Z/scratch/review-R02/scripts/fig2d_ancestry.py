import pandas as pd, numpy as np
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
D=B+'figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/'
c=pd.read_csv(D+'corrected_locus_coordinates.tsv',sep='\t')
ex=pd.read_csv(D+'excluded_credible_sets.tsv',sep='\t')
c=c[~c.credible_set_uid.isin(ex.credible_set_uid)]
c=c[c.pip_sum_gate_passed.astype(str).str.lower()=='true'].copy()
c['num']=c.protein_altering_pip_mass+c.canonical_splice_pip_mass
c['den']=c.annotated_pip_mass
c['cls']=np.where(c.trait_scope.str.startswith('tier1'),'direct','enzyme')
c['cohort']=np.select([c.study.str.startswith('MVP'),c.study.str.startswith('BBJ'),c.study.str.startswith('PanUKBB'),c.study.str.startswith('UKBB'),c.study.str.startswith('FinnGen')],['MVP','BBJ','PanUKBB','UKBB','FinnGen'],'EUR_consortium')
print(pd.crosstab([c.ancestry],[c.cls]))
print(pd.crosstab([c.cohort,c.ancestry],[c.cls]))
def est(sub,seed=20260817,nb=10000,label=''):
    sub=sub.sort_values(['chrom','pos']).copy()
    cl=[]
    for ch,b in sub.groupby('chrom',sort=False):
        new=np.concatenate([[True],np.diff(b.pos.values)>1_000_000])
        cl.append(pd.Series([f'{ch}:{i}' for i in np.cumsum(new)],index=b.index))
    sub['cl']=pd.concat(cl).reindex(sub.index)
    t=sub.groupby(['cls','cl'])[['num','den']].mean().reset_index()
    w=t.pivot(index='cl',columns='cls',values=['num','den']).fillna(0)
    na=(t.cls=='direct').sum(); nb_=(t.cls=='enzyme').sum()
    fa=w[('num','direct')].sum()/w[('den','direct')].sum(); fb=w[('num','enzyme')].sum()/w[('den','enzyme')].sum()
    rng=np.random.default_rng(seed)
    idx=rng.integers(0,len(w),size=(nb,len(w)))
    with np.errstate(divide='ignore',invalid='ignore'):
        da=w[('num','direct')].values[idx].sum(1)/w[('den','direct')].values[idx].sum(1)
        db=w[('num','enzyme')].values[idx].sum(1)/w[('den','enzyme')].values[idx].sum(1)
    ok=np.isfinite(da)&np.isfinite(db)
    diff=(da-db)[ok]
    r=np.where(db[ok]>0,da[ok]/db[ok],np.nan)
    print(f"{label}: n_inst direct={int((sub.cls=='direct').sum())} enzyme={int((sub.cls=='enzyme').sum())}; clusters direct={na} enzyme={nb_}; frac direct={fa:.4f} enzyme={fb:.4f}; ratio={fa/fb if fb>0 else np.nan:.3f}; diff={fa-fb:.4f} [{np.percentile(diff,2.5):.4f},{np.percentile(diff,97.5):.4f}]; ratio CI [{np.nanpercentile(r,2.5):.3f},{np.nanpercentile(r,97.5):.3f}]; finite={ok.sum()}")
est(c,label='ALL (reproduce)')
est(c[c.ancestry=='EUR'],label='EUR only')
est(c[c.ancestry!='EUR'],label='non-EUR only')
est(c[c.cohort=='MVP'],label='MVP only (same cohort, AFR/AMR/EAS)')
est(c[(c.cohort=='MVP')&(c.ancestry=='AMR')],label='MVP AMR only')
est(c[(c.cohort=='MVP')&(c.ancestry=='AFR')],label='MVP AFR only')
# enzyme fraction by ancestry (enzyme only), cluster-mean within ancestry
for a,g in c[c.cls=='enzyme'].groupby('ancestry'):
    g=g.sort_values(['chrom','pos']).copy(); cl=[]
    for ch,b in g.groupby('chrom',sort=False):
        new=np.concatenate([[True],np.diff(b.pos.values)>1_000_000]); cl.append(pd.Series([f'{ch}:{i}' for i in np.cumsum(new)],index=b.index))
    g['cl']=pd.concat(cl).reindex(g.index); t=g.groupby('cl')[['num','den']].mean()
    print('enzyme',a,'inst',len(g),'clusters',len(t),'frac',round(t.num.sum()/t.den.sum(),4))
for a,g in c[c.cls=='direct'].groupby('ancestry'):
    g=g.sort_values(['chrom','pos']).copy(); cl=[]
    for ch,b in g.groupby('chrom',sort=False):
        new=np.concatenate([[True],np.diff(b.pos.values)>1_000_000]); cl.append(pd.Series([f'{ch}:{i}' for i in np.cumsum(new)],index=b.index))
    g['cl']=pd.concat(cl).reindex(g.index); t=g.groupby('cl')[['num','den']].mean()
    print('direct',a,'inst',len(g),'clusters',len(t),'frac',round(t.num.sum()/t.den.sum(),4))
# credible set size by ancestry/class
print(c.groupby(['cls','ancestry']).n_members.median())
