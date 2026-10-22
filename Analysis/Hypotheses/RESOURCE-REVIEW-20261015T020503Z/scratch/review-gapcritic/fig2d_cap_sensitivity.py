import pandas as pd, numpy as np
D='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/'
c=pd.read_csv(D+'corrected_locus_coordinates.tsv',sep='\t')
ex=pd.read_csv(D+'excluded_credible_sets.tsv',sep='\t')
c=c[~c.credible_set_uid.isin(ex.credible_set_uid)]
c=c[c.pip_sum_gate_passed.astype(str).str.lower()=='true'].copy()
# locus fit key on restored coordinates, not the decimal string
c['fit']=c.study+'|'+c.chrom.astype(str)+':'+c.pos.astype(str)
c['k']=c.groupby('fit').credible_set_uid.transform('size')
print('fits:',c.fit.nunique(),' fits at k==10:',c[c.k==10].fit.nunique(),' max k',c.k.max())
print('instances in k==10 fits by scope:',c[c.k==10].groupby('trait_scope').size().to_dict())
print('k==10 fits by study:',c[c.k==10].groupby('study').fit.nunique().to_dict())
def stat(c):
    c=c.sort_values(['chrom','pos']).copy(); cl=[]
    for ch,b in c.groupby('chrom',sort=False):
        new=np.concatenate([[True],np.diff(b.pos.values)>1_000_000])
        cl.append(pd.Series([f'{ch}:{i}' for i in np.cumsum(new)],index=b.index))
    c['cl']=pd.concat(cl).reindex(c.index)
    c['num']=c.protein_altering_pip_mass+c.canonical_splice_pip_mass
    c['den']=c.annotated_pip_mass
    t=c.groupby(['trait_scope','cl'])[['num','den']].mean().reset_index()
    fr={s:g.num.sum()/g.den.sum() for s,g in t.groupby('trait_scope')}
    n={s:len(g) for s,g in t.groupby('trait_scope')}
    d,e=fr['tier1_direct_masld_pdff'],fr['tier2_liver_enzyme']
    return round(d,4),round(e,4),round(d/e,3),n
print('all (reproduce 0.2334/0.0387/6.027):',stat(c))
print('drop fits at k==10:',stat(c[c.k<10]))
print('keep fits with k<=5:',stat(c[c.k<=5]))
print('keep fits with k<=3:',stat(c[c.k<=3]))
