import pandas as pd, numpy as np
P='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
pr=pd.read_csv(P+'Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv',sep='\t')
reg=pd.read_csv(P+'figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv',sep='\t')
mem=pd.read_csv(P+'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv',sep='\t')
mem=mem[mem.mapped_symbol_status=='gencode_v49_unique_symbol_confirmed']
LIN={'hepatocytes':'hepatocyte','cholangiocytes':'cholangiocyte','macrophages':'macrophage','fibroblasts':'stellate','tcells':'t_nk'}
op=pr[pr.evidence_state.isin(['replicated_accessible','source_dependent'])]
res=[]
for r in reg.itertuples():
    genes=set(mem.loc[mem.program_uid==r.program_uid,'mapped_symbol'].astype(str).str.strip())
    o=op[(op.lineage==LIN[r.cell_type])&op.gene.isin(genes)]
    res.append(dict(prog=r.module_name,ct=r.cell_type,n_open=o.gene.nunique(),direct=sorted(o[o.trait_class=='direct_MASLD'].gene.unique()),
                    open=sorted(o.gene.unique())))
res=pd.DataFrame(res)
print('programs',len(res),'with >=1 open gene',(res.n_open>0).sum())
print(res[res.direct.str.len()>0][['prog','ct','direct','open']].to_string())
print(res[res.open.apply(lambda x:'RORA' in x)][['prog','ct','direct','open']].to_string())
ro=pr[(pr.gene=='RORA')]
print(ro[['gwas_name','trait_class','lineage','evidence_state','gse244832_accessible_mass','gse281367_accessible_mass']].to_string())
print('gwas names',sorted(pr.gwas_name.unique()))
g=pd.read_csv(P+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef/out/genetic_context_117.tsv',sep='\t')
print('F5ef table programs with open',(g.n_open_genes>0).sum())
h=pd.read_csv(P+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef/out/protein_histology_117.tsv',sep='\t')
print('testable',h.testable.sum(),'n set',h.n.dropna().unique())
t=h[h.testable==True].copy()
from scipy import stats
q=stats.false_discovery_control(t.pvalue)
print('max |q diff|',np.nanmax(np.abs(q-t.padj)))
print('padj<0.05',(t.padj<0.05).sum(),'robust',t.robust.sum(), 'robust&padj',((t.padj<0.05)&t.robust).sum())
rb=t[t.robust==True]
print('neg',(rb.rho<0).sum(),'pos',(rb.rho>0).sum())
print(rb[['cell_type','module_name','rho','padj']].sort_values('rho').to_string())
