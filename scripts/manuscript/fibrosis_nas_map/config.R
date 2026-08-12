project_root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

release_id <- "resource-f-five-coloc-v6-candidate-2026-08-10"
workstream_id <- "BULK-PROGRAM-MAP-v8"
candidate_root <- file.path(
  project_root, "RNA-seq/results/manuscript_release/candidates", release_id
)
workstream_root <- file.path(candidate_root, "workstreams", workstream_id)
precontract_nonholdout_dge <- file.path(
  candidate_root, "workstreams/BULK-PROGRAM-MAP-v2/sealed_inputs/nonholdout_dge.rds"
)
precontract_nonholdout_meta <- file.path(
  candidate_root, "workstreams/BULK-PROGRAM-MAP-v2/sealed_inputs/nonholdout_meta.rds"
)
precontract_testability_audit <- file.path(
  candidate_root, "workstreams/BULK-PROGRAM-MAP-v4/precontract/testability_by_program_cohort.tsv.gz"
)

input_dge <- file.path(
  candidate_root,
  "inputs/BG001-DECISION/arms/F_legacy/results/integration/merged_dge.rds"
)
input_meta <- file.path(
  candidate_root,
  "inputs/BG001-DECISION/arms/F_legacy/results/integration/meta_matched.rds"
)
gene_annotation <- file.path(
  candidate_root,
  "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"
)
stage_preflight <- file.path(
  candidate_root,
  "workstreams/BULK-STAGE/PREFLIGHT-v1"
)
gse193066_crosswalk <- file.path(
  stage_preflight, "audits/gse193066_participant_crosswalk.tsv"
)

program_root <- file.path(
  project_root,
  "Analysis/Multimodal_Program_Projection/candidates",
  "program-context-v2-candidate-2026-08-07/hotspot"
)
program_registry <- file.path(program_root, "program_registry_v2.tsv")
program_membership <- file.path(program_root, "program_membership_v2.tsv")
hallmark_gmt <- file.path(
  project_root,
  "Analysis/downstream_analysis/pathway_analysis/data/genesets/hallmark.gmt"
)
published_panel_root <- file.path(project_root, "data/published_gene_panels")
published_panels <- c(
  Govaere_2020 = "govaere_2020_panel.tsv",
  Pantano_2021 = "pantano_2021_panel.tsv",
  Moylan_2014 = "moylan_2014_panel.tsv",
  Arendt_2015 = "arendt_2015_panel.tsv"
)

contract_file <- file.path(workstream_root, "contract", "analysis_contract.json")
sealed_input_root <- file.path(workstream_root, "sealed_inputs")
nonholdout_dge <- file.path(sealed_input_root, "nonholdout_dge.rds")
nonholdout_meta <- file.path(sealed_input_root, "nonholdout_meta.rds")
holdout_dge <- file.path(sealed_input_root, "holdout_dge.rds")
holdout_meta <- file.path(sealed_input_root, "holdout_meta.rds")
partition_manifest <- file.path(sealed_input_root, "partition_manifest.tsv")

discovery_root <- file.path(workstream_root, "discovery")
holdout_root <- file.path(workstream_root, "holdout")
figure_root <- file.path(workstream_root, "figures")

composition_candidate <- file.path(
  project_root, "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv"
)
composition_testability <- file.path(
  program_root, "composition_sample_qc_testability.tsv"
)
composition_ready <- file.path(program_root, "COMPOSITION_QC_READY")
nonholdout_composition <- file.path(sealed_input_root, "nonholdout_composition.tsv.gz")
holdout_composition <- file.path(sealed_input_root, "holdout_composition.tsv.gz")

holdout_cohort <- "GSE193066"
discovery_cohorts <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478")
discovery_expected_n <- c(
  GSE130970 = 76L,
  GSE135251 = 214L,
  GSE162694 = 86L,
  GSE174478 = 93L
)
transport_cohorts <- c("GSE240729", "GSE126848", "GSE167523", "GSE213621")

seed <- 20260811L
n_permutations <- as.integer(Sys.getenv("MASLD_PROGRAM_MAP_PERMUTATIONS", "10000"))
n_bootstrap <- as.integer(Sys.getenv("MASLD_PROGRAM_MAP_BOOTSTRAPS", "10000"))
program_weight_coverage <- 0.80
expected_testable_programs <- 113L
gene_set_coverage <- 0.80
gene_set_min_genes <- 10L
composition_pseudocount <- 1e-6
composition_expected_available <- 1219L

required_input_files <- c(
  input_dge,
  input_meta,
  gene_annotation,
  gse193066_crosswalk,
  program_registry,
  program_membership,
  composition_candidate,
  composition_testability,
  composition_ready,
  hallmark_gmt,
  file.path(published_panel_root, unname(published_panels))
)
