# Gate C-02 step 2 (login node, < 1 GB): PC-A variants. Same scoring as job 3921 (per-array z, proxy + depth residualization,
# V4-style islands, self-inclusive row-standardized 6-NN, Lee's L). Variant a: PP panel without CYP3A4.
# Variant b: same, residualized on depth only (no hepatocyte proxy).
import importlib.util, os, sys, numpy as np, pandas as pd, scipy.sparse as sp, h5py
J='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3921/gate_c02_feasibility.py'
src=open(J).read()
head=src.split("Lmat = np.zeros")[0]   # loads data, panels, islands(), lees_L()
os.environ['NARR']='32'
exec(head)
rows=[]
for arr in prim:
    d=data[arr]; C=d['C']; n=C.shape[0]; lib=d['tot']
    N=sp.diags((1e4/lib).astype(np.float32))@C; N=N.tocsr(); N.data=np.log1p(N.data)
    Z=N.toarray(); Z=(Z-Z.mean(0))/np.where(Z.std(0,ddof=1)>0,Z.std(0,ddof=1),1)
    prox=Z[:,[col[g] for g in PROXY if det[arr][col[g]]]].mean(1)
    isl=islands(d['xy'])
    pcs=Z[:,[col[g] for g in pan_use['PC']]].mean(1,keepdims=True)
    pps=Z[:,[col[g] for g in pan_use['PP'] if g!='CYP3A4']].mean(1,keepdims=True)
    out=dict(array_id=arr)
    for lab,cols in (('proxy_depth',[prox,np.log1p(lib),np.log1p(d['ngen'])]),('depth_only',[np.log1p(lib),np.log1p(d['ngen'])])):
        D=np.column_stack([np.ones(n)]+cols); D[:,1:]=(D[:,1:]-D[:,1:].mean(0))/D[:,1:].std(0,ddof=1)
        R=lambda V: V-D@(np.linalg.pinv(D)@V)
        out['L_'+lab]=lees_L(R(pcs),R(pps),isl)[0]
    rows.append(out)
r=pd.DataFrame(rows); r.to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-C-02/g2_pca_variants.tsv',sep='\t',index=False)
for c in ['L_proxy_depth','L_depth_only']:
    print(c,'| donors L<0:',int((r[c]<0).sum()),'/32 | median',round(r[c].median(),4))
