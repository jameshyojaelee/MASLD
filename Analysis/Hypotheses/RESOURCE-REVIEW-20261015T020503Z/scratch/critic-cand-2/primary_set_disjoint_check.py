# For the primary pair set (lineage-matched, cross-lineage, >=50% full-membership coverage in all 32 primary arrays),
# recheck coverage on disjoint membership (shared genes removed from both programs), all-member denominator.
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=B+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/critic-cand-2/'
mem=pd.read_csv(B+'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv',sep='\t')
ok=mem.mapped_symbol_status.isin(['gencode_v49_unique_symbol_confirmed','gencode_v49_unambiguous_ensembl_to_symbol'])
allw={p:dict(zip(g.mapped_symbol.fillna('NA_'+g.source_gene),g.original_l1_weight)) for p,g in mem.groupby('program_uid')}
mapped={p:set(g.mapped_symbol[ok.loc[g.index]]) for p,g in mem.groupby('program_uid')}
dm=pd.read_csv(B+'Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv',sep='\t')
ps=pd.read_csv(S+'primary_pair_set_cov50_lineage_matched.tsv',sep='\t')
ps=ps[ps.cls.isin(['coordinated','informative_null','global_only'])]
want=set(mem.mapped_symbol[ok])
def read_sym(g):
    if isinstance(g,h5py.Group):
        cats=np.array([c.decode() if isinstance(c,bytes) else c for c in g['categories'][:]]); return cats[g['codes'][:]]
    return np.array([s.decode() if isinstance(s,bytes) else s for s in g[:]])
def cov(genes,det,mp):
    tot=sum(genes.values()); return sum(w for g,w in genes.items() if g in det and g in mp)/tot if tot>0 else 0
okcount=np.zeros(len(ps),dtype=int)
for arr in dm.array_id[dm.group=='primary']:
    f=h5py.File(B+f'Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled/{arr}.h5ad','r')
    sym=read_sym(f['var']['gene_symbol']); X=f['X']
    M=sp.csr_matrix((X['data'][:],X['indices'][:],X['indptr'][:]),shape=tuple(X.attrs['shape'])); f.close()
    keep=np.isin(sym,list(want)); s=sym[keep]; u,inv=np.unique(s,return_inverse=True)
    agg=sp.csr_matrix((np.ones(len(inv)),(np.arange(len(inv)),inv)),shape=(len(inv),len(u)))
    det=set(u[np.asarray(((M[:,np.where(keep)[0]]@agg)>0).mean(axis=0)).ravel()>=0.01])
    res=[]
    for a,b in zip(ps.program_a,ps.program_b):
        ga,gb=allw[a],allw[b]; sh=set(ga)&set(gb)
        da={g:w for g,w in ga.items() if g not in sh}; db={g:w for g,w in gb.items() if g not in sh}
        res.append(cov(da,det,mapped[a])>=0.5 and cov(db,det,mapped[b])>=0.5)
    okcount+=np.array(res)
ps['n_arrays_disjoint_cov50']=okcount
ps.to_csv(S+'primary_pair_set_cov50_lineage_matched_disjoint.tsv',sep='\t',index=False)
x=ps[ps.n_arrays_disjoint_cov50==32]
print('pairs passing disjoint >=0.5 in all 32:', x.cls.value_counts().to_dict(), '| programs:', len(set(x.program_a)|set(x.program_b)))
print(pd.crosstab(['-'.join(sorted([i,j])) for i,j in zip(x.cell_type_a,x.cell_type_b)], x.cls).to_string())
