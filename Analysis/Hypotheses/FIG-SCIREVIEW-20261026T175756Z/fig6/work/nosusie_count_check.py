import pandas as pd
W='/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig6/work/'
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/'
m=pd.read_csv(W+'catalog_gene_calls_with_coloc.tsv',sep='\t')
t=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/tier12/gene_level_coloc_tier12.csv')
gi=m[m.call_state=='indeterminate']
# A: direct id match to tier12
st=t.set_index('ensembl').susie_state_t12
gi=gi.assign(state_direct=gi.ensembl_id.map(st))
print('direct-ID tier12 state among genetics-indeterminate:'); print(gi.state_direct.value_counts(dropna=False))
# GEN_ID adjudicated classes file
a=pd.read_csv(B+'preflight/gen_identity_adjudication_v1/gen_frozen_classes_ensembl_adjudicated.tsv',sep='\t',low_memory=False)
print(a.columns.tolist()[:40])
