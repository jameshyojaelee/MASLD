# Reference Data (Consolidated)

Sources: `docs/reference/pipeline_outputs.md`, `docs/reference/dataset_notes.md`

---

## Pipeline Outputs

### Human Integration (Scripts 00-26)

| Script | Output |
|--------|--------|
| 00 | `metadata/unified_metadata.csv` (~1,400+ samples, 10 cohorts) |
| 01 | `qc/sample_qc_report.csv` (`pass_technical` for non-sex; `pass_sex` for sex-stratified) |
| 02 | `results/per_study/{DATASET}_de_results.csv` |
| 03 | `results/integration/merged_counts_tmm.rds` |
| 04 | `results/integration/variance_partition.pdf` |
| 05 | `results/integration/dream_results.csv` (34,704 genes, 19,552 DEGs padj<0.1) |
| 08 | `results/integration/gsea_results.csv` |
| 11 | `results/integration/ortholog_mapped_degs.csv` |
| 12 | `results/integration/library_intersection.csv` |
| 13 | `results/disease_signatures/nafl_vs_nash_dream.csv` (uses `diagnosis_harmonized`) |
| 14-16 | `results/disease_signatures/fibrosis_dream.csv`, `unified_disease_signatures.csv` |
| 25 | `results/causal_inference/deconv_attribution_scores.csv` (2,924 hep-intrinsic) |
| 26 | `results/integration/sex_deg_classification.csv` (4,430 sex-dimorphic) |

### Novel Analyses (Scripts 19-48, in `RNA-seq/results/`)

| Script | Output |
|--------|--------|
| 19 | `causal_inference/{GWAS}/mr_results_all.csv`, `twas_spredixcan_liver.csv`, `coloc_results.csv` |
| 20 | `drug_repurposing/cgp_reversal_hits.csv` (65), `multi_layer_drug_targets.csv` (20) |
| 27a | `multi_evidence/multi_evidence_atlas.csv` (14,080 x 61 cols) |
| 27b | `multi_evidence/benchmark_results.csv`, `preset_configurations.json` |
| 30 | `gwas_spatial_convergence/gwas_convergence_evidence_cards.csv` (14 loci) |
| 31 | `drug_repurposing/lincs_annotated_compounds.csv` |
| 32 | `drug_repurposing/network_proximity/proximity_results.csv` |
| 33 | `audit_sensitivity/sex_divergent_permissive.csv`, `cgp_liver_metabolic_filtered.csv`, `conserved_core_enrichment.csv` |
| 34 | `audit_sensitivity/combat_seq/combat_seq_concordance.csv` |
| 40 | `drug_repurposing/clinical_drug_validation_table.csv` (10 drugs) |
| 41 | `pathway_programs/ferroptosis_degs.csv` (63 DEGs) |
| 42 | `convergence/regulon_drug_gwas_convergence.csv` |
| 43 | `zonation/deg_zonation_classification.csv` |
| 44 | `subtypes/nmf_assignments.csv` (1,254 samples, k=2) |
| 46 | `causal_inference/finngen_{nafld,nash,hcc}/coloc_results.csv` |
| 47 | `causal_inference/bbj_{alt,ast,ggt}/coloc_results.csv` |
| 48 | `causal_inference/cross_ancestry/cross_ancestry_comparison.csv` |

### Mouse Integration (Scripts M00-M07)

| Script | Output |
|--------|--------|
| M00 | `metadata/unified_mouse_metadata.csv` (444 samples) |
| M01 | `qc/sample_qc_report.csv` (443 QC-passing) |
| M02 | `results/per_diet/{DIET}_de_results.csv` |
| M03 | `results/meta_analysis/meta_per_diet.csv`, `dream_pooled_results.csv` |
| M04 | `results/mouse_consensus_degs.csv` |

---

## Dataset Notes

### eQTL Sources

| Dataset | N | Build | Notes |
|---------|---|-------|-------|
| Broadaway et al. 2024 (AJHG) | 1,183 (meta of 4 studies) | hg19 | Primary COLOC source; general liver tissue (not MASLD-specific); 6,564 eGenes |
| GTEx v8 Liver | 208 | hg38 | TWAS (FUSION/PredictDB); independent COLOC replication; lower power |
| Ghodsian sc-eQTL (2024) | 312 | hg38 | MASLD patients; cell-type-specific ieQTLs (hepatocyte, endothelial, cholangiocyte, stellate); disease-context resolved |

### RNA-seq Cohort Quirks

| Dataset | Key Notes |
|---------|-----------|
| PRJNA512027 | 7 samples excluded (NCBI SRA corruption). **Dropped from paper 2026-05-15** after L0/S0 library-prep batch was found perfectly confounded with diagnosis (all 34 Controls in L0, all 102 NASH in S0). Retained in YAML/05d fibrosis-stage modeling only — not counted in "cohorts collected" headline. |
| GSE167523 | No healthy controls (NAFL vs NASH only); excluded from mega-analysis |
| GSE162876/GSE159911 | Shallow sequencing (1.2-4.8M reads); libsize threshold lowered to 1M |
| GSE213621 | 367 samples; sex not annotated (uses inferred_sex) |
| GSE162694 | 143 samples, SE HiSeq 3000; 31 controls + 112 NASH F0-F4; bariatric (Bril 2021) |
| GSE174478 | 94 samples, PE NovaSeq 6000; Japanese NAFLD; no healthy controls; excluded from mega-analysis |
| GSE193066 | 164 samples; NAFLD HCC-naive; no healthy controls; excluded from mega-analysis |
| GSE240729 | 67 samples; MASLD F0-F4 |
