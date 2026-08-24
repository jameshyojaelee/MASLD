#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "lib_continuum.R"))

pre <- read_prespec()
root <- project_root()
out <- out_root()
inputs_out <- file.path(out, "inputs")
ensure_new_dir(inputs_out)

defaults <- list(
  resource_dge = file.path(
    root, "RNA-seq/results/manuscript_release/candidates",
    "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION",
    "arms/F_five/results/integration/merged_dge.rds"
  ),
  resource_manifest = file.path(
    root, "figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis",
    "stage_extensions/five_cohort_sample_manifest.tsv"
  ),
  gene_annotation = file.path(
    root, "RNA-seq/results/manuscript_release/candidates",
    "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE",
    "frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"
  ),
  program_registry = file.path(
    root, "Analysis/Multimodal_Program_Projection/candidates",
    "program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv"
  ),
  program_membership = file.path(
    root, "Analysis/Multimodal_Program_Projection/candidates",
    "program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
  ),
  composition = file.path(
    root, "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv"
  ),
  paired_crosswalk = file.path(
    root, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/candidates",
    "program-context-v2-candidate-2026-08-07",
    "fibrosis-adjacent-true-kleiner-lvqw-v1/audits",
    "gse193066_donor_biopsy_crosswalk.tsv"
  ),
  paired_counts = file.path(
    root, "results/remediation/bg001",
    "bg001-fragment-v211-gencode49-20260807T194845Z/frozen_sets/read_counts",
    "GSE193066/gene_counts.txt"
  ),
  paired_metadata = file.path(
    root, "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066",
    "metadata/metadata.tsv"
  ),
  native_nmf_axis = file.path(
    root, "RNA-seq/results/continuous_axis_benchmark/20260814T184338Z",
    "arms/A3_axis.tsv"
  ),
  prespecification = file.path(script_dir, "00_prespecification.json")
)

env_names <- c(
  resource_dge = "HAC_DGE_PATH",
  resource_manifest = "HAC_MANIFEST_PATH",
  gene_annotation = "HAC_GENE_ANNOTATION_PATH",
  program_registry = "HAC_PROGRAM_REGISTRY_PATH",
  program_membership = "HAC_PROGRAM_MEMBERSHIP_PATH",
  composition = "HAC_COMPOSITION_PATH",
  paired_crosswalk = "HAC_PAIRED_CROSSWALK_PATH",
  paired_counts = "HAC_PAIRED_COUNTS_PATH",
  paired_metadata = "HAC_PAIRED_METADATA_PATH",
  native_nmf_axis = "HAC_NMF_AXIS_PATH"
)

paths <- defaults
for (key in names(env_names)) {
  value <- Sys.getenv(env_names[[key]], unset = "")
  if (nzchar(value)) paths[[key]] <- value
}
for (key in names(paths)) {
  assert_true(file.exists(paths[[key]]), paste0("Missing input ", key, ": ", paths[[key]]))
  paths[[key]] <- normalizePath(paths[[key]], mustWork = TRUE)
}

kamzolas_root <- read_env_path("KAMZOLAS_V1_ROOT")
git_commit <- system2("git", c("-C", kamzolas_root, "rev-parse", "HEAD"), stdout = TRUE)
assert_true(length(git_commit) == 1L && identical(git_commit, pre$competitor_release$commit),
            paste0("Kamzolas release drift: expected ", pre$competitor_release$commit,
                   ", observed ", paste(git_commit, collapse = "")))
git_dirty <- system2("git", c("-C", kamzolas_root, "status", "--porcelain"), stdout = TRUE)
assert_true(length(git_dirty) == 0L, "Kamzolas v1 checkout is dirty")

message("Reading Resource substrate")
dge <- readRDS(paths$resource_dge)
assert_true(!is.null(dge$counts), "Resource RDS is not a DGEList-like object")
assert_true(identical(dim(dge$counts), c(23370L, 844L)),
            paste0("Resource count dimensions drifted: ", paste(dim(dge$counts), collapse = " x ")))
assert_true(all(dge$counts == floor(dge$counts)), "Resource counts are not integer fragments")
assert_true(!anyDuplicated(colnames(dge$counts)), "Resource sample IDs are duplicated")

manifest <- fread(paths$resource_manifest, na.strings = c("", "NA"))
required_manifest <- c("sample_id", "analysis_unit_id", "dataset", "inferred_sex",
                       "fibrosis_stage", "nas_score")
assert_true(all(required_manifest %in% names(manifest)), "Resource manifest schema drift")
assert_true(nrow(manifest) == 844L && !anyDuplicated(manifest$sample_id),
            "Resource manifest donor census drift")
assert_true(setequal(colnames(dge$counts), manifest$sample_id),
            "Resource DGE and manifest participant sets differ")
assert_true(!anyDuplicated(manifest$analysis_unit_id),
            "Resource analysis unit is not one donor per row")

expected <- c(GSE126848 = 53L, GSE130970 = 76L, GSE135251 = 214L,
              GSE162694 = 140L, GSE213621 = 361L)
census <- manifest[, .(
  n_donors = .N,
  n_fibrosis = sum(!is.na(fibrosis_stage)),
  n_nas = sum(!is.na(nas_score))
), by = dataset][order(dataset)]
assert_true(identical(setNames(census$n_donors, census$dataset), expected),
            "Resource cohort census drift")

resource_genes <- unique(base_gene_id(rownames(dge$counts)))
assert_true(length(resource_genes) == 23370L, "Resource Ensembl base IDs are not unique")
write_tsv_once(data.table(gene_id_base = resource_genes),
               file.path(inputs_out, "resource_gene_ids.tsv"))
write_tsv_once(census, file.path(inputs_out, "cohort_census.tsv"))

registry <- fread(paths$program_registry)
membership <- fread(paths$program_membership)
assert_true(nrow(registry) == pre$programs$family_size &&
              uniqueN(registry$program_uid) == pre$programs$family_size,
            "Frozen program registry must contain exactly 117 unique programs")
assert_true(all(unique(membership$program_uid) %in% registry$program_uid),
            "Program membership contains an unregistered program")

paired <- fread(paths$paired_crosswalk)
paired <- paired[pass_technical == TRUE]
paired_n <- paired[, .N, by = donor_id][N == 2L, .N]
assert_true(paired_n == 54L, paste0("Expected 54 QC-valid paired donors; observed ", paired_n))

input_manifest <- rbindlist(lapply(names(paths), function(key) {
  info <- file.info(paths[[key]])
  data.table(
    input_id = key,
    path = paths[[key]],
    bytes = as.numeric(info$size),
    sha256 = sha256_file(paths[[key]])
  )
}))
input_manifest <- rbind(
  input_manifest,
  data.table(
    input_id = "kamzolas_v1_git_checkout",
    path = kamzolas_root,
    bytes = NA_real_,
    sha256 = pre$competitor_release$commit
  )
)
write_tsv_once(input_manifest, file.path(inputs_out, "input_manifest.tsv"))

# Freeze the current manuscript Figure 3 and prior Slingshot benchmark byte-for-
# byte. The candidate validator re-hashes these protected scopes after the run;
# pre-existing user edits are allowed, but this workflow may not change them.
protected_roots <- c(
  current_figure3 = file.path(root, "figures/main/fig3_bulk_transcriptomics"),
  prior_continuum_benchmark = file.path(
    root, "RNA-seq/results/continuous_axis_benchmark/20260814T184338Z"
  )
)
assert_true(all(dir.exists(protected_roots)),
            paste0("Protected scope missing: ",
                   paste(protected_roots[!dir.exists(protected_roots)], collapse = ", ")))
protected_snapshot <- rbindlist(lapply(names(protected_roots), function(scope) {
  scope_root <- protected_roots[[scope]]
  files <- list.files(scope_root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  data.table(
    protected_scope = scope,
    protected_root = normalizePath(scope_root),
    relative_path = substring(normalizePath(files), nchar(normalizePath(scope_root)) + 2L),
    bytes = as.numeric(file.size(files)),
    sha256 = vapply(files, sha256_file, character(1L))
  )
}))
write_tsv_once(protected_snapshot,
               file.path(inputs_out, "protected_scope_snapshot.tsv"))

preflight <- data.table(
  check_id = c(
    "kamzolas_commit", "kamzolas_clean_checkout", "resource_integer_counts",
    "resource_donor_unit", "resource_participant_set", "resource_gene_ids", "cohort_census",
    "program_family", "paired_donor_census"
  ),
  status = "PASS",
  observed = c(
    git_commit, "clean", "23370x844", "844 unique analysis_unit_id",
    "DGE/manifest exact set; downstream reindexed to manifest order",
    length(resource_genes), paste(names(expected), expected, sep = "=", collapse = ";"),
    nrow(registry), paired_n
  )
)
write_tsv_once(preflight, file.path(inputs_out, "preflight_checks.tsv"))

writeLines(capture.output(sessionInfo()), file.path(inputs_out, "sessionInfo.txt"))
message("INPUT_PREPARATION_COMPLETE: ", inputs_out)
