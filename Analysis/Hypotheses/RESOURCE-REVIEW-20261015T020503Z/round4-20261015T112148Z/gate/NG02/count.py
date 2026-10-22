import h5py, numpy as np, pandas as pd, json, sys
np.random.seed(20260923)
G='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG02'
V='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Spatial/data/public_expansion_2026/yakubovsky2026/v.mat'
SCH='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/acquisition/yakubovsky2026/yakubovsky_hdf5_join_schema.tsv'
CAT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
PRO='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/scdvp/scdvp_table3_healthy.tsv'
sch=pd.read_csv(SCH,sep='\t')
sch=sch[sch.field.isin(['gene_name','mat','zon_struct/zone_index'])]
f=h5py.File(V,'r')
def strs(ds):
    out=[]
    for r in ds[()].ravel():
        o=f[r]; out.append(''.join(chr(c) for c in o[()].ravel()))
    return out
cat=pd.read_parquet(CAT)
classes=cat.primary_evidence_class.astype(str)
cent=['CYP2E1','GLUL','CYP1A2','CYP3A4','OAT']; port=['ASS1','HAL','SDS','CPS1','PCK1']
res={}; slopes={}; det={}
for d in ['M1','M2','M3','M4','M5','M6','M7','M8']:
    s=sch[sch.patient==d].set_index('field').hdf5_path
    genes=strs(f[s['gene_name']])
    X=f[s['mat']][()]; z=f[s['zon_struct/zone_index']][()].ravel()
    if X.shape[0]!=len(z): X=X.T
    info={'shape':list(X.shape),'zone_unique':sorted(set(np.round(z,3).tolist()))[:12],'is_integer':bool(np.allclose(X,np.round(X))),'rowsum_median':float(np.median(X.sum(1))),'max':float(X.max())}
    Xn=X/X.sum(1,keepdims=True)*1e4
    zc=z-z.mean()
    sl=(zc@ (Xn-Xn.mean(0)))/(zc@zc)
    m=Xn.mean(0)
    sl_rel=sl/np.where(m>0,m,np.nan)   # slope per zone step relative to the gene's mean
    slopes[d]=pd.Series(sl_rel,index=genes); det[d]=pd.Series((X>0).mean(0),index=genes)
    info['markers_rel_slope']={g:(float(slopes[d].get(g,np.nan))) for g in cent+port}
    res[d]=info; print(d,json.dumps(info),flush=True)
    del X,Xn
S=pd.DataFrame(slopes); D=pd.DataFrame(det)
S=S[~S.index.duplicated()]; D=D[~D.index.duplicated()]
live=['M2','M3','M4','M5','M6','M7','M8']
# marker check: central markers should share one sign and portal the opposite, in each donor
mk=S.loc[[g for g in cent+port if g in S.index]]
mk.to_csv(G+'/marker_slopes.tsv',sep='\t')
sepd={d:bool(mk.loc[[g for g in cent if g in mk.index],d].max()<mk.loc[[g for g in port if g in mk.index],d].min() or mk.loc[[g for g in cent if g in mk.index],d].min()>mk.loc[[g for g in port if g in mk.index],d].max()) for d in S.columns}
print('marker perfect separation per donor',sepd)
# coverage
sym=set(cat.symbol.dropna().astype(str))
detect_all=(D[live]>=0.05).all(1)
cov={}
for cl in ['ALL']+sorted(classes.unique()):
    sub=cat if cl=='ALL' else cat[classes==cl]
    ss=set(sub.symbol.dropna().astype(str))
    inV=ss & set(D.index)
    cov[cl]={'n_catalog':len(sub),'n_symbol_in_visium':len(inV),'n_detected_5pct_all7':int(detect_all.reindex(list(inV)).fillna(False).sum())}
pro=pd.read_csv(PRO,sep='\t')
psym=set(pro.protein.astype(str))
for cl in cov:
    sub=cat if cl=='ALL' else cat[classes==cl]
    ss=set(sub.symbol.dropna().astype(str))
    cov[cl]['n_protein_coef']=len(ss&psym)
    cov[cl]['n_protein_zonated']=int(pro[pro.protein.isin(ss)]['zonated'].astype(str).str.lower().isin(['true','1','yes']).sum())
print(json.dumps(cov,indent=1))
# protein known-answer
pm=pro.set_index('protein')
print('protein markers', {g:(float(pm.coefficient.get(g,np.nan)), str(pm.q_value.get(g,np.nan))) for g in cent+port})
json.dump({'donors':res,'separation':sepd,'coverage':cov},open(G+'/count_summary.json','w'),indent=1,default=str)
