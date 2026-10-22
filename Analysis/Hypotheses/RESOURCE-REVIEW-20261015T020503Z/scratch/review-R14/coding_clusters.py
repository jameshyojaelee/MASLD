import pandas as pd, numpy as np
d=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv',sep='\t')
print('rows',len(d), 'gate', d.pip_sum_gate_passed.value_counts().to_dict())
g=d[d.pip_sum_gate_passed==True].copy()
g['coding']=g.protein_altering_pip_mass+g.canonical_splice_pip_mass
print(g.trait_scope.value_counts().to_dict())
# single-linkage 1Mb clustering per chrom, pooled
g=g.sort_values(['chrom','pos']).reset_index(drop=True)
cid=[];c=-1;prev=None
for ch,p in zip(g.chrom,g.pos):
    if prev is None or ch!=prev[0] or p-prev[1]>1_000_000: c+=1
    cid.append(c); prev=(ch,p)
g['cluster']=cid
print('clusters',g.cluster.nunique())
res={}
for sc,sub in g.groupby('trait_scope'):
    cm=sub.groupby('cluster').coding.mean()
    res[sc]=cm
    print(sc,'nclust',len(cm),'frac',round(cm.mean()*100,4))
both=set(res['tier1_direct_masld_pdff'].index)&set(res['tier2_liver_enzyme'].index)
print('shared clusters',len(both))
dc=res['tier1_direct_masld_pdff'].sort_values(ascending=False)
info=g[g.trait_scope=='tier1_direct_masld_pdff'].groupby('cluster').agg(chrom=('chrom','first'),pos_min=('pos','min'),pos_max=('pos','max'),n_cs=('credible_set_uid','size'),studies=('study','nunique'))
out=info.join(dc.rename('mean_coding')).sort_values('mean_coding',ascending=False)
out['cum_share']=out.mean_coding.cumsum()/out.mean_coding.sum()
print(out.to_string())
print('ratio',dc.mean()/res['tier2_liver_enzyme'].mean())
# drop top-k direct clusters
for k in range(1,6):
    print('drop top',k,'direct frac',round(dc.iloc[k:].mean()*100,3),'ratio',round(dc.iloc[k:].mean()/res['tier2_liver_enzyme'].mean(),3))
ec=res['tier2_liver_enzyme'].sort_values(ascending=False)
print('enzyme clusters with coding mean>0.2:',(ec>0.2).sum(),'of',len(ec))
print('top-19 enzyme clusters by max |z| not available here')
