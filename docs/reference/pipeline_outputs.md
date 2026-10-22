# Pipeline Output Inventory

Extracted from CLAUDE.md for reference. See CLAUDE.md for pipeline commands.

## Human (Scripts 00-18, 25-26)
| Script | Key Output |
|--------|-----------|
| 00 | `metadata/unified_metadata.csv` (~1,400+ samples, 10 cohorts) |
| 01 | `qc/sample_qc_report.csv` (~1,100+ QC-passing projected); `pass_technical` for non-sex analyses, `pass_sex` for sex-stratified |
| 02 | `results/per_study/{DATASET}_de_results.csv` |
| 03 | `results/integration/merged_counts_tmm.rds` (merged DGEList) |
| 04 | `results/integration/variance_partition.pdf` |
| 05 | `results/integration/dream_results.csv` (34,704 genes, 19,552 DEGs at padj<0.1) |
| 08 | `results/integration/gsea_results.csv` |
| 11 | `results/integration/ortholog_mapped_degs.csv` |
| 12 | `results/integration/library_intersection.csv` |
| 13 | `results/disease_signatures/nafl_vs_nash_dream.csv`, `nafl_vs_nash_per_study.csv` — uses `diagnosis_harmonized` (NAS-based) |
| 14-16 | `results/disease_signatures/fibrosis_dream.csv`, `unified_disease_signatures.csv` |
| 25 | `results/causal_inference/deconv_attribution_scores.csv` (2,924 hepatocyte-intrinsic), `dream_results_deconvolution_adjusted.csv` |
| 26 | `results/integration/sex_deg_classification.csv` (4,430 sex-dimorphic), `sex_interaction_dream.csv` |

## Novel Analyses (Scripts 19-44, in `RNA-seq/results/`)
| Script | Key Output |
|--------|-----------|
| 19 | `causal_inference/{GWAS_NAME}/mr_results_all.csv`, `twas_spredixcan_liver.csv`, `coloc_results.csv` |
| 20 | `drug_repurposing/cgp_reversal_hits.csv` (65, 28 liver-specific), `top50_reversal_compounds.csv`, `multi_layer_drug_targets.csv` (20) |
| 27a | `multi_evidence/multi_evidence_atlas.csv` (14,080 genes x 61 cols, 8 layers + ferroptosis + zonation) |
| 27b | `multi_evidence/benchmark_results.csv`, `preset_configurations.json` (5 presets; Therapeutic recovery: 179 genes, F1=0.352) |
| 30 | `gwas_spatial_convergence/gwas_human_mapping.csv`, `gwas_convergence_evidence_cards.csv` (14 priority loci) |
| 31 | `drug_repurposing/lincs_annotated_compounds.csv` (compound names, MOA, DGIdb, composite scores) |
| 32 | `drug_repurposing/network_proximity/proximity_results.csv` (drug-disease proximity z-scores) |
| 33 | `audit_sensitivity/sex_divergent_permissive.csv`, `cgp_liver_metabolic_filtered.csv`, `conserved_core_enrichment.csv` |
| 34 | `audit_sensitivity/combat_seq/combat_seq_concordance.csv` (r=0.918, Jaccard=0.615, dir=100%) |
| 40 | `drug_repurposing/clinical_drug_validation_table.csv` (10 drugs: 4 Strong, 1 Moderate, 5 Absent) |
| 41 | `pathway_programs/ferroptosis_degs.csv` (63 DEGs), `ferroptosis_enrichment.csv`, `lipotoxicity_program.csv` |
| 42 | `convergence/regulon_drug_gwas_convergence.csv` (GWAS enrichment 3.5x, p=0.014), `convergence_network.csv` |
| 43 | `zonation/deg_zonation_classification.csv` (14,080 genes, PC/PP/non-zoned), `zonation_enrichment.csv` |
| 44 | `subtypes/nmf_assignments.csv` (1,254 samples, k=2), `subtype_markers.csv` (200 genes), `subtype_clinical_associations.csv` |
| 46 | `causal_inference/finngen_{nafld,nash,hcc}/coloc_results.csv` (FinnGen R12 × Broadaway COLOC, case-control) |
| 47 | `causal_inference/bbj_{alt,ast,ggt}/coloc_results.csv` (BBJ × Broadaway COLOC, cross-ancestry East Asian) |
| 48 | `causal_inference/cross_ancestry/cross_ancestry_comparison.csv`, `cross_ancestry_validated_targets.csv`, `ancestry_specific_loci.csv` |

## Mouse (Scripts M00-M07)
| Script | Key Output |
|--------|-----------|
| M00 | `metadata/unified_mouse_metadata.csv` (444 samples) |
| M01 | `qc/sample_qc_report.csv` (443 QC-passing) |
| M02 | `results/per_diet/{DIET}_de_results.csv` |
| M03 | `results/meta_analysis/meta_per_diet.csv`, `dream_pooled_results.csv` |
| M04 | `results/mouse_consensus_degs.csv` |

> See `paper_outline.md` Section 2 for full competitive landscape (10 competitors with gap analysis).
