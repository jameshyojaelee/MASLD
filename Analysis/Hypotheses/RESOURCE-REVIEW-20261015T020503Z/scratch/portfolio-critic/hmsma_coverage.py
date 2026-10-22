# Coverage of frozen program weight among genes present (and detected in >=1% of spots) in one HMSMA array, and
# how many cross-lineage coordinated edges have both programs at >=80% retained weight. No program score is computed.
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
mem=pd.read_csv(B+'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv',sep='\t')
print(mem.columns.tolist()[:12], mem.shape)
f=h5py.File(B+'Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled/HRA_11.h5ad','r')
sym=pd.Series([s.decode() if isinstance(s,bytes) else s for s in f['var']['gene_symbol'][:]]) if 'categories' not in f['var']['gene_symbol'] else None
if sym is None:
    g=f['var']['gene_symbol']; cats=np.array([c.decode() for c in g['categories'][:]]); sym=pd.Series(cats[g['codes'][:]])
X=f['X']
if isinstance(X,h5py.Group):
    M=sp.csr_matrix((X['data'][:],X['indices'][:],X['indptr'][:]),shape=tuple(X.attrs['shape']))
else: M=sp.csr_matrix(X[:])
det=np.asarray((M>0).mean(axis=0)).ravel()
detected=set(sym[det>=0.01]); present=set(sym)
print('spots',M.shape[0],'genes',M.shape[1],'detected>=1%',len(detected))
gcol=[c for c in mem.columns if 'symbol' in c.lower() or c.lower()=='gene'][0]
wcol=[c for c in mem.columns if 'weight' in c.lower()][0]
pcol=[c for c in mem.columns if 'program' in c.lower() or 'feature' in c.lower()][0]
print('cols',gcol,wcol,pcol)
mem['w']=mem[wcol].abs()
cov=mem.assign(det=mem[gcol].isin(detected)).groupby(pcol).apply(lambda x:(x.w*x.det).sum()/x.w.sum())
T4='RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction/'
e=pd.read_csv(B+T4+'coordinated_edges.tsv',sep='\t')
x=e[e.same_cell_type==False]
ok=x.program_a.map(cov).ge(0.8)&x.program_b.map(cov).ge(0.8)
print('programs with weight coverage>=0.8 (detected>=1%):',int((cov>=0.8).sum()),'of',len(cov))
print('cross-lineage coordinated edges',len(x),'with both programs >=0.8 coverage',int(ok.sum()))
print('ids matched to membership:',x.program_a.isin(cov.index).mean(),x.program_b.isin(cov.index).mean())
