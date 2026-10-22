import glob, h5py, numpy as np, pandas as pd, scipy.sparse as sp, os
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
mem=pd.read_csv(R+'/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv',sep='\t')
progs={'hotspot_hepatocytes_f05c535ae5bbc0b9':'Stromal ECM','hotspot_hepatocytes_48f39dd4d817a10e':'Ductular injury'}
mem=mem[mem.program_uid.isin(progs)].copy()
files={}
for s in ['JBO014','JBO015','JBO018','JBO019','JBO022']:
    files[('GSE192741',s)]=glob.glob(R+f'/Analysis/Spatial/data/GSE192741/*filtered_feature_bc_matrix_{s}.h5')[0]
for s in ['VLP115_A','VLP115_D','VLP116_A','VLP116_D','VLP119_A','VLP119_D','VLP120_A','VLP120_D','VLP121_A','VLP121_D']:
    files[('Vu',s)]=R+f'/data/vu_et_al_2025/raw/Visium/{s}/filtered_feature_bc_matrix.h5'
c2l={'GSE192741':R+'/Analysis/Spatial/results/cell2location/spatial_model/spatial_deconvolved.h5ad','Vu':R+'/Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad'}
c2lgenes={}
for k,p in c2l.items():
    with h5py.File(p,'r') as f:
        v=f['var']; idx=v.attrs.get('_index','_index')
        c2lgenes[k]=set(x.decode() if isinstance(x,bytes) else x for x in v[idx][:])
rows=[]
for (ds,s),p in files.items():
    with h5py.File(p,'r') as f:
        m=f['matrix']; names=[x.decode() for x in m['features/name'][:]]
        X=sp.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]))
        X=X.tocsr()  # genes x spots
        nsp=X.shape[1]
        pos={}
        for i,n in enumerate(names): pos.setdefault(n,i)
        for _,r in mem.iterrows():
            g=r.mapped_symbol; i=pos.get(g)
            det = X[i].getnnz()/nsp if i is not None else np.nan
            tot = X[i].sum() if i is not None else np.nan
            rows.append(dict(dataset=ds,sample=s,n_spots=nsp,program=progs[r.program_uid],gene=g,l1=r.original_l1_weight,in_raw=i is not None,in_c2l=g in c2lgenes[ds],frac_spots_detected=det,umi_total=tot))
d=pd.DataFrame(rows); d.to_csv('gene_sample_detection.tsv',sep='\t',index=False)
a=d.groupby(['dataset','program','gene','l1','in_raw','in_c2l']).agg(med_frac=('frac_spots_detected','median'),min_frac=('frac_spots_detected','min')).reset_index()
a.to_csv('gene_detection_summary.tsv',sep='\t',index=False)
for (ds,pg),x in a.groupby(['dataset','program']):
    tot=x.l1.sum()
    for lab,mask in [('in_c2l',x.in_c2l),('in_raw',x.in_raw),('missing_c2l_in_raw',(~x.in_c2l)&x.in_raw),('missing_c2l_in_raw_det>=1pct',(~x.in_c2l)&x.in_raw&(x.med_frac>=0.01)),('missing_c2l_in_raw_det>=5pct',(~x.in_c2l)&x.in_raw&(x.med_frac>=0.05))]:
        print(ds,pg,lab,'n',int(mask.sum()),'of',len(x),'L1 %.3f'%(x.l1[mask].sum()/tot))
print(a[~a.in_c2l].sort_values(['dataset','program','l1']).to_string())
