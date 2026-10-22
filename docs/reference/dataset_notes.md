# Dataset-Specific Notes

Per-dataset quirks and special handling.

## eQTL Datasets

- **Broadaway et al. 2024 (AJHG)**: Meta-analysis of 4 liver eQTL studies (N=1,183 total, 6,564 eGenes). Component studies include GTEx v8 (N=208, post-mortem donors), METSIM (Finnish metabolic syndrome cohort, liver biopsies), and two additional surgical resection/transplant cohorts. Samples are **general liver tissue** (not MASLD-specific), which is appropriate for COLOC: standard eQTLs reflect constitutive germline regulation of gene expression, detectable regardless of disease state. Downloaded from Mohlke Lab FTP; hg19 coordinates; marginal summary statistics per chromosome. Primary eQTL source for all COLOC analyses.
- **GTEx v8 Liver**: N=208 post-mortem donors. Used for TWAS (pre-trained elastic net models via FUSION/PredictDB) and as independent replication for COLOC (smaller N = lower power, but independent cohort). hg38 coordinates.
- **Ghodsian sc-eQTL (2024)**: N=312 from MASLD patients. Cell-type-specific eQTLs (ieQTLs) across 4 cell types (hepatocyte, endothelial, cholangiocyte, stellate). Unlike Broadaway/GTEx bulk eQTLs, these are **disease-context and cell-type resolved** — the variant's effect on expression can differ between MASLD and control conditions, and between cell types. Used for sc-eQTL COLOC and sc-TWAS. hg38 coordinates.

## RNA-seq Cohorts

- **PRJNA512027**: 7 samples excluded due to persistent NCBI SRA corruption (see `RNA-seq/Human/Patient_Cohorts/README.md`). **Dropped from paper 2026-05-15** after L0/S0 library-prep batch was found perfectly confounded with diagnosis label (all 34 Controls in L0, all 102 NASH in S0). Was already excluded from the Disease-vs-Control mega-analysis (`include_in_mega: false`); retained in YAML for fibrosis-stage modeling (05d, where per-stage contrasts are identifiable across L0/S0 batches). PRJNA512027 is **not** counted in the "cohorts collected" headline.
- **GSE167523**: No healthy controls (NAFL vs NASH only); excluded from mega-analysis
- **GSE162876/GSE159911**: Shallow sequencing (1.2-4.8M reads); library size QC threshold lowered to 1M
- **GSE213621**: 367 samples, MASLD spectrum; sex not annotated (uses inferred_sex)
- **GSE162694**: 143 samples, single-end HiSeq 3000; 31 normal histology controls + 112 NASH F0-F4; bariatric surgery cohort (Bril et al. 2021)
- **GSE174478**: 94 samples, paired-end NovaSeq 6000; Japanese NAFLD cohort, no healthy controls; excluded from mega-analysis
- **GSE193066**: 164 samples; NAFLD cohort, HCC-naive; no healthy controls; excluded from mega-analysis
- **GSE240729**: 67 samples; MASLD fibrosis stages F0-F4
