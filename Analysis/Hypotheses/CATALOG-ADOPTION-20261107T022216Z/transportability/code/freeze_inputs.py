from pathlib import Path
import csv
import hashlib

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
OUT = ROOT/'Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/transportability'
INPUTS = {
 'cohort': 'figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv',
 'native126': 'Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/gse126848_filterByExpr_refit_genes.tsv.gz',
 'classes': 'Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/fig3c_gene_classes_filterByExpr_refit.tsv.gz',
 'pooled': 'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv',
 'dge': 'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds',
 'metadata': 'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds',
 'roster': 'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_unit_roster.tsv',
 'ancestry': 'Analysis/Hypotheses/H79-20261015T125125Z/results/run-20261015T132031Z/D1_gene_table.tsv.gz',
 'ss_sidecar': 'Analysis/Hypotheses/SS-20261025T135237Z/out_s3b-20261028T035515Z/steatosite_catalog_sidecar.tsv.gz',
 'ss_gene': 'Analysis/Hypotheses/SS-20261025T135237Z/out_s3a-20261026T234211Z/gene_bin.tsv.gz',
 'ss_family': 'Analysis/Hypotheses/SS-20261025T135237Z/out_s3a-20261026T234211Z/F_SS1_primary.tsv',
 'ss_units': 'Analysis/Hypotheses/SS-20261025T135237Z/out_s3a-20261026T234211Z/units_used.tsv',
 'config': 'config/human_datasets.yaml',
}
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''):h.update(block)
 return h.hexdigest()
target=OUT/'input_manifest.tsv'
assert not target.exists()
with target.open('w') as f:
 w=csv.DictWriter(f,fieldnames=['input','path','sha256'],delimiter='\t');w.writeheader()
 for k,v in INPUTS.items():
  p=ROOT/v;w.writerow(dict(input=k,path=str(p),sha256=sha(p)))
spec=OUT/'ANALYSIS_SPEC.md'
(OUT/'ANALYSIS_SPEC.sha256').write_text(f'{sha(spec)}  ANALYSIS_SPEC.md\n')
print('Frozen analysis specification and',len(INPUTS),'input files; no new result fit run.')
