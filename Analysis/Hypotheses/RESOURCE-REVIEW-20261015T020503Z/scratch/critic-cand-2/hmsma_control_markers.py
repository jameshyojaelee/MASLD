# Detection of positive-control marker panels (repo marker table) in the 32 primary HMSMA arrays; overlap with program members and V4 proxy.
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
mk=pd.read_csv(B+'data/external/liver_celltype_markers_combined.tsv',sep='\t')
panels=['SAMmac_scar_associated','HSC_activated_SAMmes','Hepatocyte_PC (pericentral)','Hepatocyte_PP (periportal)']
proxy={'ALB','APOA1','APOA2','APOC3','APOH','AHSG','FGA','FGB','FGG','HP','SERPINA1','TF','TTR','HPX','RBP4'}
mem=pd.read_csv(B+'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv',sep='\t')
members=set(mem.mapped_symbol.dropna())
dm=pd.read_csv(B+'Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv',sep='\t')
want=set(mk[mk.cell_type.isin(panels)].gene)
def read_sym(g):
    if isinstance(g,h5py.Group):
        cats=np.array([c.decode() if isinstance(c,bytes) else c for c in g['categories'][:]]); return cats[g['codes'][:]]
    return np.array([s.decode() if isinstance(s,bytes) else s for s in g[:]])
cnt={g:0 for g in want}
for arr in dm.array_id[dm.group=='primary']:
    f=h5py.File(B+f'Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled/{arr}.h5ad','r')
    sym=read_sym(f['var']['gene_symbol']); X=f['X']
    M=sp.csr_matrix((X['data'][:],X['indices'][:],X['indptr'][:]),shape=tuple(X.attrs['shape'])); f.close()
    keep=np.isin(sym,list(want)); s=sym[keep]; u,inv=np.unique(s,return_inverse=True)
    agg=sp.csr_matrix((np.ones(len(inv)),(np.arange(len(inv)),inv)),shape=(len(inv),len(u)))
    det=np.asarray(((M[:,np.where(keep)[0]]@agg)>0).mean(axis=0)).ravel()
    for g,v in zip(u,det):
        if v>=0.01: cnt[g]+=1
for p in panels:
    gs=list(mk[mk.cell_type==p].gene)
    all32=[g for g in gs if cnt.get(g,0)==32]
    print(p,'| panel',len(gs),'| detected>=1% in all 32:',len(all32),'| of those not in proxy:',len([g for g in all32 if g not in proxy]),'| of those in any program:',len([g for g in all32 if g in members]))
    print('   ',all32)
