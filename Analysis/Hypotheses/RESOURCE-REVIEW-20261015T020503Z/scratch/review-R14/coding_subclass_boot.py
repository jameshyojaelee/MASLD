import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np
d=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv',sep='\t')
d['coding']=d.protein_altering_pip_mass+d.canonical_splice_pip_mass
d=d.sort_values(['chrom','pos']).reset_index(drop=True)
cid=[];c=-1;prev=None
for ch,p in zip(d.chrom,d.pos):
    if prev is None or ch!=prev[0] or p-prev[1]>1_000_000: c+=1
    cid.append(c); prev=(ch,p)
d['cluster']=cid
d['cls']=np.where(d.trait_scope=='tier2_liver_enzyme','ENZ',np.where(d.trait=='PDFF','PDFF','NAFLDNASH'))
tab=d.groupby(['cluster','cls']).coding.mean().unstack()  # cluster x class means (NaN if absent)
ncl=len(tab)
rng=np.random.default_rng(20261015)
B=5000
out={k:[] for k in ['PDFF','NAFLDNASH','DIRECT']}
# direct pooled: recompute cluster mean over PDFF+NAFLD/NASH sets
dirm=d[d.cls!='ENZ'].groupby('cluster').coding.mean().reindex(tab.index)
tab['DIRECT']=dirm
obs={k:tab[k].mean()/tab['ENZ'].mean() for k in ['PDFF','NAFLDNASH','DIRECT']}
A=tab.values; cols=list(tab.columns)
for b in range(B):
    idx=rng.integers(0,ncl,ncl)
    s=tab.iloc[idx]
    e=s['ENZ'].mean()
    for k in out: out[k].append(s[k].mean()/e)
print('clusters',ncl,'counts', tab.notna().sum().to_dict())
for k in out:
    a=np.array(out[k]); a=a[np.isfinite(a)]
    print(k,'obs ratio',round(obs[k],3),'95% cluster-bootstrap',np.round(np.percentile(a,[2.5,97.5]),3),'finite',len(a), 'frac',round(tab[k].mean()*100,3))
print('enzyme frac',round(tab['ENZ'].mean()*100,4),'seed 20261015 B',B)
