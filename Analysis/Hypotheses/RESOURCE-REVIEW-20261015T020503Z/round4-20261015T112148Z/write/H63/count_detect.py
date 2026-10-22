# H63 pre-outcome counts: detection per donor (zones 1-8 only) and set sizes by Catalog class.
# Computes NO zonation statistic (no slope, no center of mass, no agreement).
import h5py, numpy as np, pandas as pd, json, sys
W=sys.argv[1]
V='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Spatial/data/public_expansion_2026/yakubovsky2026/v.mat'
SCH='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/acquisition/yakubovsky2026/yakubovsky_hdf5_join_schema.tsv'
CAT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
PRO='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/scdvp/scdvp_table3_healthy.tsv'
LM=json.load(open(W+'/landmarks.json')); LM=set(LM['central_LM'])|set(LM['portal_LM'])
sch=pd.read_csv(SCH,sep='\t'); f=h5py.File(V,'r')
def strs(ds): return [''.join(chr(c) for c in f[r][()].ravel()) for r in ds[()].ravel()]
det={}; genes=None
for d in ['M2','M3','M4','M5','M6','M7','M8']:
    s=sch[sch.patient==d].set_index('field').hdf5_path
    g=strs(f[s['gene_name']]); genes=g
    z=f[s['zon_struct/zone_index']][()].ravel()
    X=f[s['mat']][()]
    if X.shape[0]!=len(z): X=X.T
    keep=(z>=1)&(z<=8)
    X=X[keep]
    det[d]=pd.Series((X>0).mean(0),index=g)
    print(d,'spots zones1-8',int(keep.sum()),'excluded',int((~keep).sum()),'min row sum',float(X.sum(1).min()),flush=True)
    del X
D=pd.DataFrame(det); D=D[~D.index.duplicated(keep=False)]   # drop the 10 duplicated symbols entirely
D.to_csv(W+'/detection_fraction.tsv.gz',sep='\t')
print('symbols after dropping duplicated',len(D))
cat=pd.read_parquet(CAT); pro=pd.read_csv(PRO,sep='\t'); pro=pro[~pro.protein.duplicated(keep=False)]
print('protein rows after dropping duplicated ids',len(pro))
psym=set(pro.protein); pzon=set(pro.protein[pro.zonated=='Yes']); pstr=set(pro.protein[pro.strong_zonated=='Yes'])
sets={'primary_M4_M8':['M4','M5','M6','M7','M8'],'sens_M2_M8':['M2','M3','M4','M5','M6','M7','M8']}
rows=[]
for sn,dl in sets.items():
    robust=set(D.index[(D[dl]>=0.05).all(axis=1)]); anyd=set(D.index[(D[dl]>0).any(axis=1)])
    for cl in ['ALL','genetically_anchored','established_state_associated','concordant','unresolved','untested']:
        ss=set((cat if cl=='ALL' else cat[cat.primary_evidence_class==cl]).symbol)
        inV=ss&set(D.index)
        r=dict(donor_set=sn,cls=cl,n_catalog=len(ss),in_visium=len(inV),detected_any=len(inV&anyd),robust=len(inV&robust),
               robust_nonLM=len((inV&robust)-LM),protein=len(ss&psym),robust_protein=len(inV&robust&psym),
               robust_protein_nonLM=len((inV&robust&psym)-LM),robust_protein_zonated_nonLM=len((inV&robust&pzon)-LM),
               robust_protein_strong_nonLM=len((inV&robust&pstr)-LM),robust_protein_LM=len(inV&robust&psym&LM))
        rows.append(r)
out=pd.DataFrame(rows); out.to_csv(W+'/count_detect.tsv',sep='\t',index=False); print(out.to_string())
# protein genes absent from robust set: why
pr=set(pro.protein)&set(cat.symbol)
print('protein&catalog',len(pr),'in visium',len(pr&set(D.index)))
