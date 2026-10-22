# Per array and program: number of member symbols detected in >=1% of spots and weight coverage (all-member denominator).
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT=B+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/critic-cand-2/'
mem=pd.read_csv(B+'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv',sep='\t')
dmap=pd.read_csv(B+'Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv',sep='\t')
ok=mem.mapped_symbol_status.isin(['gencode_v49_unique_symbol_confirmed','gencode_v49_unambiguous_ensembl_to_symbol'])
want=set(mem.mapped_symbol[ok])
def read_sym(g):
    if isinstance(g,h5py.Group):
        cats=np.array([c.decode() if isinstance(c,bytes) else c for c in g['categories'][:]]); return cats[g['codes'][:]]
    return np.array([s.decode() if isinstance(s,bytes) else s for s in g[:]])
rows=[]
for arr in dmap.array_id:
    f=h5py.File(B+f'Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled/{arr}.h5ad','r')
    sym=read_sym(f['var']['gene_symbol']); X=f['X']
    M=sp.csr_matrix((X['data'][:],X['indices'][:],X['indptr'][:]),shape=tuple(X.attrs['shape'])); f.close()
    keep=np.isin(sym,list(want)); s_keep=sym[keep]; uniq,inv=np.unique(s_keep,return_inverse=True)
    agg=sp.csr_matrix((np.ones(len(inv)),(np.arange(len(inv)),inv)),shape=(len(inv),len(uniq)))
    det=np.asarray(((M[:,np.where(keep)[0]]@agg)>0).mean(axis=0)).ravel()
    d=set(uniq[det>=0.01])
    m=mem.assign(det=ok & mem.mapped_symbol.isin(d))
    g=m.groupby('program_uid').apply(lambda x: pd.Series({'n_detected':int(x.det.sum()),'n_members':len(x),'coverage':float((x.original_l1_weight*x.det).sum()/x.original_l1_weight.sum())}), include_groups=False)
    g['array_id']=arr; rows.append(g.reset_index())
pd.concat(rows).to_csv(OUT+'hmsma_program_detected_genes_per_array.tsv',sep='\t',index=False)
print('done')
