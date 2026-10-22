# H63 pre-outcome counts: landmark genes, zone counts, key checks. No gene slope or agreement statistic is computed here.
import h5py, numpy as np, pandas as pd, json
V='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Spatial/data/public_expansion_2026/yakubovsky2026/v.mat'
SCH='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/acquisition/yakubovsky2026/yakubovsky_hdf5_join_schema.tsv'
CAT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
PRO='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/lens-L08/scdvp/scdvp_table3_healthy.tsv'
sch=pd.read_csv(SCH,sep='\t'); f=h5py.File(V,'r')
def strs(ds): return [''.join(chr(c) for c in f[r][()].ravel()) for r in ds[()].ravel()]
genes0=None; lm={}
for d in ['M1','M2','M3','M4','M5','M6','M7','M8']:
    s=sch[sch.patient==d].set_index('field').hdf5_path
    root=s['mat'].rsplit('/',1)[0]
    g=strs(f[s['gene_name']])
    if genes0 is None: genes0=g
    z=f[s['zon_struct/zone_index']][()].ravel()
    pc=[g[int(i)-1] for i in f[root+'/LM_pc_ind'][()].ravel()]
    pp=[g[int(i)-1] for i in f[root+'/LM_pp_ind'][()].ravel()]
    lm[d]=(sorted(pc),sorted(pp))
    nm=''.join(chr(c) for c in f[s['main_feature']][()].ravel())
    print(d,nm,'same_gene_order',g==genes0,'n_spots',len(z),'zone_counts',{int(k):int(v) for k,v in zip(*np.unique(z,return_counts=True))},
          'n_dup_symbols',len(g)-len(set(g)),'has_mat_norm',(root+'/mat_norm') in f, 'has_bg_vector',(root+'/bg_vector') in f)
print('LM sets identical across donors', all(lm[d]==lm['M2'] for d in lm))
print('central LM',lm['M2'][0]); print('portal LM',lm['M2'][1])
LM=set(lm['M2'][0])|set(lm['M2'][1])
cent=['CYP2E1','GLUL','CYP1A2','CYP3A4','OAT']; port=['ASS1','HAL','SDS','CPS1','PCK1']
print('markers that are landmarks',[m for m in cent+port if m in LM])
cat=pd.read_parquet(CAT); print('catalog cols',list(cat.columns)[:30]); print('catalog rows',len(cat),'symbol NA',int(cat.symbol.isna().sum()),'dup symbols',int(cat.symbol.dropna().duplicated().sum()))
print(cat.primary_evidence_class.value_counts().to_dict())
pro=pd.read_csv(PRO,sep='\t'); print('protein rows',len(pro),'dup',int(pro.protein.duplicated().sum()),'zonated',pro.zonated.value_counts().to_dict(),'strong',pro.strong_zonated.value_counts().to_dict())
print('protein ids with ; or -',int(pro.protein.astype(str).str.contains(';').sum()))
print('LM in protein table',len(LM & set(pro.protein)),'LM in catalog',len(LM & set(cat.symbol.dropna())))
sym=set(cat.symbol.dropna()); print('protein in catalog',len(set(pro.protein)&sym),'protein in visium',len(set(pro.protein)&set(genes0)),'protein in both',len(set(pro.protein)&sym&set(genes0)))
json.dump({'central_LM':lm['M2'][0],'portal_LM':lm['M2'][1]},open('landmarks.json','w'),indent=1)
