#!/usr/bin/env Rscript
# Candidate-only substrate gate for the corrected exact-Kleiner stage rebuild.
# This script deliberately performs no gene filtering, normalization, voom,
# coefficient fitting, hypothesis testing, or figure generation.

suppressPackageStartupMessages(library(data.table))

CANDIDATE_ID <- "resource-f-five-coloc-v6-candidate-2026-08-10"
WORKSTREAM_ID <- "bulk-stage-preflight-v1"
EXPECTED_HASHES <- c(
  merged_counts_raw = "16afc00bf1db706731d07df1dfeabdad0bc423b1caa1b822631685a55ec7f225",
  merged_dge = "56d1cf97f791a81404ddc6c5ea4c5df9e6a38e7cf011b39d0d7c3279c87afa4d",
  meta_matched = "90ba6aca68643c4680c08fc982a6c72b76edb740ab726874cb50603711756a83",
  sample_qc = "1c7fef991c599554369ecc79de2b2f441597cf9b6c1a7ff3ecbe05f1c5969560",
  unified_metadata = "b9f6afc1e916391951c638a60d8ec95f1b380d20311469c0830f638f16649699",
  gse193066_metadata = "b66158d20da2127e8c3cc7ebae4f1d85a2d6813b7dd6f2af4bdac60070ff2bc5"
)
TRUE_KLEINER_COHORTS <- c(
  "GSE130970", "GSE135251", "GSE162694",
  "GSE174478", "GSE193066", "GSE240729"
)
EXPECTED_STAGE_COUNTS <- c(`0` = 126L, `1` = 187L, `2` = 174L, `3` = 132L, `4` = 42L)
EXPECTED_TRANSITIONS <- data.table(
  transition = c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"),
  contrast = c("F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3"),
  low_stage = 0:3,
  high_stage = 1:4,
  n_low = c(126L, 187L, 174L, 107L),
  n_high = c(187L, 174L, 132L, 42L),
  n_samples = c(313L, 361L, 306L, 149L),
  n_cohorts = c(6L, 6L, 6L, 5L),
  design_rank = c(8L, 8L, 8L, 7L)
)

candidate_root <- normalizePath(
  Sys.getenv("MASLD_RESOURCE_CANDIDATE_ROOT", unset = ""),
  mustWork = TRUE
)
if (basename(candidate_root) != CANDIDATE_ID) {
  stop("Candidate identity drift: ", candidate_root, call. = FALSE)
}
sentinel <- file.path(candidate_root, ".resource_candidate_root")
if (!file.exists(sentinel) || trimws(readLines(sentinel, warn = FALSE)[[1L]]) != CANDIDATE_ID) {
  stop("Candidate sentinel is missing or invalid", call. = FALSE)
}

bg001_root <- file.path(candidate_root, "inputs/BG001-DECISION")
arm_root <- file.path(bg001_root, "arms/F_legacy")
stage_metadata_root <- file.path(candidate_root, "inputs/BULK-NINE-COHORT/source_metadata")
input_paths <- c(
  merged_counts_raw = file.path(arm_root, "results/integration/merged_counts_raw.rds"),
  merged_dge = file.path(arm_root, "results/integration/merged_dge.rds"),
  meta_matched = file.path(arm_root, "results/integration/meta_matched.rds"),
  sample_qc = file.path(arm_root, "qc/sample_qc_report.csv"),
  unified_metadata = file.path(
    bg001_root,
    "source_snapshot/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
  ),
  gse193066_metadata = file.path(stage_metadata_root, "GSE193066_metadata.tsv")
)
final_root <- file.path(candidate_root, "workstreams/BULK-STAGE/PREFLIGHT-v1")
parent_root <- dirname(final_root)
publisher_path <- file.path(candidate_root, "code/publish_noreplace.py")
snapshot_io_path <- file.path(candidate_root, "code/snapshot_io.py")
python_path <- normalizePath(Sys.getenv("MASLD_PYTHON_PATH"), mustWork = TRUE)
rscript_path <- normalizePath(Sys.getenv("MASLD_RSCRIPT_PATH"), mustWork = TRUE)
producer_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(producer_path) != 1L) stop("Cannot resolve stage producer path", call. = FALSE)
producer_path <- normalizePath(sub("^--file=", "", producer_path), mustWork = TRUE)
job_tag <- Sys.getenv("SLURM_JOB_ID", as.character(Sys.getpid()))
tmp_root <- file.path(parent_root, paste0(".PREFLIGHT-v1.tmp.", job_tag, ".", Sys.getpid()))

thread_variables <- c(
  "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
)
if (!identical(Sys.getenv("SLURM_CPUS_PER_TASK"), "4") ||
    any(vapply(thread_variables, function(variable) {
      !identical(Sys.getenv(variable), "4")
    }, logical(1)))) {
  stop("Thread-count contract drift; stage preflight requires four threads", call. = FALSE)
}

is_symlink <- function(path) {
  target <- Sys.readlink(path)
  !is.na(target) && nzchar(target)
}
sha256_file <- function(path) {
  normalized <- normalizePath(path, mustWork = TRUE)
  output <- system2("sha256sum", normalized, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L) {
    stop("sha256sum failed for ", normalized, call. = FALSE)
  }
  value <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", value)) stop("Invalid SHA256 for ", normalized, call. = FALSE)
  value
}
write_tsv_once <- function(value, path) {
  connection <- file(path, open = "wx")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  write.table(
    as.data.frame(value), connection, sep = "\t", row.names = FALSE,
    col.names = TRUE, quote = FALSE, na = "NA", eol = "\n"
  )
  close(connection)
}

if (file.exists(final_root) || is_symlink(final_root) || file.exists(tmp_root) || is_symlink(tmp_root)) {
  stop("Refusing existing stage-preflight output", call. = FALSE)
}
for (path in c(input_paths, publisher_path, snapshot_io_path, producer_path)) {
  if (!file.exists(path) || file.info(path)$isdir || is_symlink(path)) {
    stop("Missing, non-regular, or symlinked input: ", path, call. = FALSE)
  }
}
observed_hashes <- vapply(input_paths, sha256_file, character(1))
if (!identical(unname(observed_hashes[names(EXPECTED_HASHES)]), unname(EXPECTED_HASHES))) {
  stop("Frozen stage-input hash mismatch", call. = FALSE)
}

raw <- readRDS(input_paths[["merged_counts_raw"]])
dge <- readRDS(input_paths[["merged_dge"]])
meta <- as.data.table(readRDS(input_paths[["meta_matched"]]))
qc <- fread(input_paths[["sample_qc"]])
unified <- fread(input_paths[["unified_metadata"]])
gse_raw <- fread(input_paths[["gse193066_metadata"]], check.names = FALSE)

if (!identical(dim(raw), c(86369L, 1281L))) stop("Raw fragment matrix dimension drift", call. = FALSE)
if (!identical(dim(dge), c(24196L, 1257L))) stop("Corrected DGE dimension drift", call. = FALSE)
if (anyDuplicated(rownames(raw)) || anyDuplicated(colnames(raw)) ||
    anyDuplicated(rownames(dge)) || anyDuplicated(colnames(dge))) {
  stop("Duplicate raw/DGE gene or sample identifiers", call. = FALSE)
}
required_meta <- c("sample_id", "dataset", "fibrosis_stage", "inferred_sex")
required_qc <- c("sample_id", "pass_technical")
required_unified <- c("sample_id", "dataset", "fibrosis_stage")
required_gse <- c("Run", "sample_id", "!Sample_title", "biopsy", "fibrosis stage")
if (!all(required_meta %in% names(meta)) || !all(required_qc %in% names(qc)) ||
    !all(required_unified %in% names(unified)) || !all(required_gse %in% names(gse_raw))) {
  stop("Stage source schema drift", call. = FALSE)
}
if (nrow(meta) != 1281L || nrow(qc) != 1281L || nrow(unified) != 1469L || nrow(gse_raw) != 164L ||
    anyDuplicated(meta$sample_id) || anyDuplicated(qc$sample_id) ||
    anyDuplicated(unified$sample_id) || anyDuplicated(gse_raw$sample_id)) {
  stop("Stage source row-count or sample-key drift", call. = FALSE)
}
raw_ids <- colnames(raw)
meta_ids <- as.character(meta$sample_id)
qc_ids <- as.character(qc$sample_id)
if (!setequal(meta_ids, qc_ids) || !setequal(meta_ids, raw_ids)) {
  stop("Raw/meta/QC sample universes disagree", call. = FALSE)
}
if (!identical(raw_ids, qc_ids)) stop("Raw/QC sample order drift", call. = FALSE)
meta_order <- match(raw_ids, meta_ids)
if (anyNA(meta_order) || anyDuplicated(meta_order)) {
  stop("Metadata-to-raw sample alignment is not bijective", call. = FALSE)
}
meta <- meta[meta_order]
if (!identical(raw_ids, as.character(meta$sample_id))) {
  stop("Metadata realignment to raw/QC order failed", call. = FALSE)
}
pass_ids <- qc[pass_technical == TRUE, sample_id]
if (length(pass_ids) != 1257L || !identical(colnames(dge), pass_ids)) {
  stop("Corrected QC/DGE sample identity drift", call. = FALSE)
}
if (!setequal(meta$sample_id, qc$sample_id) || !setequal(meta$sample_id, colnames(raw))) {
  stop("Raw/meta/QC sample universes disagree", call. = FALSE)
}

# The deposited Hoshida title root is the only permitted participant key.
gse <- gse_raw[, .(
  run_id = as.character(Run),
  metadata_sample_id = as.character(sample_id),
  source_title = as.character(get("!Sample_title")),
  biopsy = as.character(biopsy),
  source_fibrosis_stage = suppressWarnings(as.integer(get("fibrosis stage")))
)]
if (anyNA(gse) || any(gse$run_id != gse$metadata_sample_id) ||
    anyDuplicated(gse$run_id) || anyDuplicated(gse$source_title) ||
    any(!grepl("^HUnafld[0-9]{3}(_[12])?$", gse$source_title)) ||
    any(!gse$source_fibrosis_stage %in% 0:4)) {
  stop("GSE193066 deposited-key contract failed", call. = FALSE)
}
gse[, participant_token := sub("_[12]$", "", source_title)]
gse[, analysis_unit_id := paste("GSE193066", participant_token, sep = "::")]
gse[, is_first_biopsy := biopsy == "1st biopsy"]
gse[, is_second_biopsy := biopsy == "2nd biopsy"]
if (sum(gse$is_first_biopsy) != 106L || sum(gse$is_second_biopsy) != 58L ||
    uniqueN(gse$analysis_unit_id) != 106L) {
  stop("GSE193066 participant/biopsy census drift", call. = FALSE)
}
participant_gate <- gse[, .(
  n_rows = .N,
  n_first = sum(is_first_biopsy),
  n_second = sum(is_second_biopsy)
), by = analysis_unit_id]
if (any(participant_gate$n_first != 1L) || any(!participant_gate$n_rows %in% c(1L, 2L)) ||
    any(participant_gate$n_second != participant_gate$n_rows - 1L) ||
    sum(participant_gate$n_rows == 2L) != 58L ||
    any(gse[is_first_biopsy == TRUE, grepl("_2$", source_title)]) ||
    any(gse[is_second_biopsy == TRUE, !grepl("_2$", source_title)])) {
  stop("GSE193066 first-biopsy participant gate failed", call. = FALSE)
}

gse_crosswalk <- merge(
  gse,
  meta[dataset == "GSE193066", .(
    sample_id,
    harmonized_fibrosis_stage = suppressWarnings(as.integer(fibrosis_stage)),
    inferred_sex
  )],
  by.x = "run_id", by.y = "sample_id", all.x = TRUE, sort = FALSE
)
gse_crosswalk <- merge(
  gse_crosswalk,
  qc[, .(sample_id, pass_technical)],
  by.x = "run_id", by.y = "sample_id", all.x = TRUE, sort = FALSE
)
gse_crosswalk[, `:=`(
  raw_present = run_id %in% colnames(raw),
  corrected_dge_present = run_id %in% colnames(dge)
)]
if (nrow(gse_crosswalk) != 164L || anyNA(gse_crosswalk$harmonized_fibrosis_stage) ||
    any(gse_crosswalk$source_fibrosis_stage != gse_crosswalk$harmonized_fibrosis_stage)) {
  stop("GSE193066 source/harmonized-stage join failed", call. = FALSE)
}
gse_crosswalk[, included_primary :=
  is_first_biopsy & pass_technical & raw_present & corrected_dge_present]
gse_crosswalk[, exclusion_reason := fcase(
  included_primary, "included_first_biopsy",
  is_second_biopsy, "excluded_second_biopsy",
  !pass_technical, "excluded_failed_corrected_technical_qc",
  !raw_present, "excluded_absent_from_raw_matrix",
  !corrected_dge_present, "excluded_absent_from_corrected_dge",
  default = "excluded_unclassified"
)]
if (sum(gse_crosswalk$included_primary) != 105L ||
    sum(gse_crosswalk$exclusion_reason == "excluded_second_biopsy") != 58L ||
    sum(gse_crosswalk$exclusion_reason == "excluded_failed_corrected_technical_qc") != 1L) {
  stop("GSE193066 corrected first-biopsy inclusion census failed", call. = FALSE)
}

# Reconstruct the stage universe from frozen unified metadata, corrected QC,
# and the authoritative GSE193066 first-biopsy key. meta_matched contributes
# inferred sex only after the source-stage set has been frozen.
stage_source <- unified[
  dataset %in% TRUE_KLEINER_COHORTS &
    suppressWarnings(as.integer(fibrosis_stage)) %in% 0:4
]
stage_source[, fibrosis_stage := suppressWarnings(as.integer(fibrosis_stage))]
stage_source <- merge(
  stage_source,
  qc[, .(sample_id, pass_technical)],
  by = "sample_id", all.x = TRUE, sort = FALSE
)
stage_source <- merge(
  stage_source,
  meta[, .(sample_id, inferred_sex, matched_dataset = dataset,
           matched_fibrosis_stage = suppressWarnings(as.integer(fibrosis_stage)))],
  by = "sample_id", all.x = TRUE, sort = FALSE
)
stage_source <- merge(
  stage_source,
  gse[, .(run_id, source_title, biopsy, analysis_unit_id, is_first_biopsy)],
  by.x = "sample_id", by.y = "run_id", all.x = TRUE, sort = FALSE
)
stage_source[dataset != "GSE193066", `:=`(
  source_title = sample_id,
  biopsy = "single deposited biopsy",
  analysis_unit_id = paste(dataset, sample_id, sep = "::"),
  is_first_biopsy = TRUE
)]
stage_source[, analysis_unit_key_source := fifelse(
  dataset == "GSE193066",
  "deposited_title_root",
  "single_deposited_sample_accession"
)]
stage_source[, `:=`(
  raw_present = sample_id %in% colnames(raw),
  corrected_dge_present = sample_id %in% colnames(dge)
)]
stage_source[, `:=`(
  qc_present = !is.na(pass_technical),
  meta_present = !is.na(matched_dataset)
)]
if (nrow(stage_source) != 731L || anyDuplicated(stage_source$sample_id) ||
    any(stage_source[meta_present == TRUE, dataset != matched_dataset]) ||
    any(stage_source[meta_present == TRUE, fibrosis_stage != matched_fibrosis_stage])) {
  stop("Frozen unified/meta/QC stage join failed", call. = FALSE)
}
if (anyNA(stage_source[meta_present == TRUE, inferred_sex]) ||
    any(!nzchar(trimws(stage_source[meta_present == TRUE, inferred_sex])))) {
  stop("Matched stage metadata lacks inferred sex", call. = FALSE)
}
stage_source[, included_primary :=
  is_first_biopsy & qc_present & pass_technical & meta_present &
    raw_present & corrected_dge_present]
stage_source[, exclusion_reason := fcase(
  included_primary, "included_primary",
  dataset == "GSE193066" & !is_first_biopsy, "excluded_gse193066_second_biopsy",
  !qc_present, "excluded_absent_from_corrected_qc",
  !meta_present, "excluded_absent_from_matched_metadata",
  !raw_present, "excluded_absent_from_raw_matrix",
  !pass_technical, "excluded_failed_corrected_technical_qc",
  !corrected_dge_present, "excluded_absent_from_corrected_dge",
  default = "excluded_unclassified"
)]
reason_counts <- stage_source[, .N, by = exclusion_reason]
expected_reasons <- data.table(
  exclusion_reason = c(
    "included_primary", "excluded_gse193066_second_biopsy",
    "excluded_absent_from_corrected_qc",
    "excluded_failed_corrected_technical_qc"
  ),
  N = c(661L, 58L, 2L, 10L)
)
setorder(reason_counts, exclusion_reason)
setorder(expected_reasons, exclusion_reason)
if (!identical(reason_counts, expected_reasons)) stop("Stage eligibility reason census drift", call. = FALSE)

stage_samples <- stage_source[included_primary == TRUE, .(
  sample_id,
  analysis_unit_id,
  analysis_unit_key_source,
  dataset,
  source_title,
  biopsy,
  fibrosis_stage,
  inferred_sex,
  pass_technical,
  raw_present,
  corrected_dge_present
)]
setorder(stage_samples, dataset, fibrosis_stage, sample_id)
if (nrow(stage_samples) != 661L || anyDuplicated(stage_samples$sample_id) ||
    anyDuplicated(stage_samples$analysis_unit_id) ||
    !setequal(unique(stage_samples$dataset), TRUE_KLEINER_COHORTS)) {
  stop("Corrected exact-stage biological-unit gate failed", call. = FALSE)
}
stage_census <- stage_samples[, .(n_samples = .N), by = fibrosis_stage][order(fibrosis_stage)]
if (!identical(stage_census$fibrosis_stage, 0:4) ||
    !identical(stage_census$n_samples, unname(EXPECTED_STAGE_COUNTS))) {
  stop("Corrected 661-sample stage census drift", call. = FALSE)
}

stage_grid <- CJ(dataset = sort(TRUE_KLEINER_COHORTS), fibrosis_stage = 0:4)
cohort_stage <- stage_samples[, .(n_samples = .N), by = .(dataset, fibrosis_stage)]
cohort_stage <- merge(stage_grid, cohort_stage, by = c("dataset", "fibrosis_stage"), all.x = TRUE)
cohort_stage[is.na(n_samples), n_samples := 0L]
setorder(cohort_stage, dataset, fibrosis_stage)

transition_samples <- list()
transition_census <- list()
transition_cohorts <- list()
design_audits <- list()
model_designs <- list()
for (index in seq_len(nrow(EXPECTED_TRANSITIONS))) {
  expected <- EXPECTED_TRANSITIONS[index]
  eligible_cohorts <- stage_samples[
    fibrosis_stage %in% c(expected$low_stage, expected$high_stage),
    .(n_stages = uniqueN(fibrosis_stage)), by = dataset
  ][n_stages == 2L, dataset]
  selected <- copy(stage_samples[
    dataset %in% eligible_cohorts &
      fibrosis_stage %in% c(expected$low_stage, expected$high_stage)
  ])
  selected[, `:=`(
    transition = expected$transition,
    contrast = expected$contrast,
    arm = fifelse(fibrosis_stage == expected$high_stage, "high", "low")
  )]
  observed <- data.table(
    transition = expected$transition,
    contrast = expected$contrast,
    low_stage = expected$low_stage,
    high_stage = expected$high_stage,
    n_low = sum(selected$arm == "low"),
    n_high = sum(selected$arm == "high"),
    n_samples = nrow(selected),
    n_cohorts = uniqueN(selected$dataset),
    design_rank = expected$design_rank
  )
  compare <- c("transition", "contrast", "low_stage", "high_stage", "n_low", "n_high", "n_samples", "n_cohorts")
  if (!identical(as.list(observed[, ..compare]), as.list(expected[, ..compare]))) {
    stop("Transition census drift for ", expected$transition, call. = FALSE)
  }
  expected_cohorts <- if (expected$transition == "F3_to_F4") {
    setdiff(TRUE_KLEINER_COHORTS, "GSE193066")
  } else TRUE_KLEINER_COHORTS
  if (!setequal(eligible_cohorts, expected_cohorts)) {
    stop("Transition cohort-set drift for ", expected$transition, call. = FALSE)
  }
  transition_samples[[expected$transition]] <- selected[, .(
    sample_id, analysis_unit_id, dataset, fibrosis_stage, inferred_sex,
    transition, contrast, arm
  )]
  transition_census[[expected$transition]] <- observed[, !"design_rank"]
  transition_cohorts[[expected$transition]] <- selected[, .(
    n_low = sum(arm == "low"),
    n_high = sum(arm == "high"),
    n_total = .N,
    both_arms_present = all(c("low", "high") %in% arm)
  ), by = .(transition, contrast, dataset)]

  info <- data.frame(
    dataset = factor(selected$dataset),
    inferred_sex = factor(selected$inferred_sex),
    fib_group = factor(selected$arm, levels = c("low", "high")),
    row.names = selected$sample_id
  )
  design <- model.matrix(~ dataset + inferred_sex + fib_group, data = info)
  rank <- qr(design)$rank
  if (rank != ncol(design) || rank != expected$design_rank ||
      !"fib_grouphigh" %in% colnames(design)) {
    stop("No-fit design estimability failed for ", expected$transition, call. = FALSE)
  }
  design_audits[[expected$transition]] <- data.table(
    transition = expected$transition,
    contrast = expected$contrast,
    formula = "~ dataset + inferred_sex + fib_group",
    coefficient = "fib_grouphigh",
    design_rank = rank,
    n_design_columns = ncol(design),
    design_columns = paste(colnames(design), collapse = ";"),
    n_samples = nrow(design),
    n_cohorts = uniqueN(selected$dataset),
    cohort_set = paste(sort(unique(selected$dataset)), collapse = ";"),
    gene_filter_status = "not_run",
    model_fit_status = "not_run"
  )
  model_dt <- as.data.table(design)
  model_dt[, sample_id := rownames(design)]
  model_dt[, `:=`(transition = expected$transition, contrast = expected$contrast)]
  setcolorder(model_dt, c("sample_id", "transition", "contrast", colnames(design)))
  model_designs[[expected$transition]] <- model_dt
}

transition_sample_manifest <- rbindlist(transition_samples, use.names = TRUE)
transition_summary <- rbindlist(transition_census, use.names = TRUE)
transition_cohort_audit <- rbindlist(transition_cohorts, use.names = TRUE)
design_preflight <- rbindlist(design_audits, use.names = TRUE)
model_design_preflight <- rbindlist(model_designs, use.names = TRUE, fill = TRUE)
if (nrow(transition_sample_manifest) != 1129L || any(!transition_cohort_audit$both_arms_present)) {
  stop("Transition sample/cohort manifest gate failed", call. = FALSE)
}

dir.create(parent_root, recursive = TRUE, showWarnings = FALSE, mode = "0750")
if (!dir.exists(parent_root) || is_symlink(parent_root)) {
  stop("Unable to create safe stage parent root", call. = FALSE)
}
if (!isTRUE(dir.create(
  tmp_root, recursive = FALSE, showWarnings = FALSE, mode = "0750"
))) {
  stop("Unable to create unique stage-preflight temporary root", call. = FALSE)
}
for (child in c("audits", "manifests")) {
  if (!isTRUE(dir.create(
    file.path(tmp_root, child), recursive = FALSE,
    showWarnings = FALSE, mode = "0750"
  ))) {
    stop("Unable to create stage-preflight child root: ", child, call. = FALSE)
  }
}

setorder(stage_source, dataset, fibrosis_stage, sample_id)
setorder(gse_crosswalk, analysis_unit_id, biopsy, run_id)
write_tsv_once(stage_source, file.path(tmp_root, "audits/stage_eligibility_audit.tsv"))
write_tsv_once(gse_crosswalk, file.path(tmp_root, "audits/gse193066_participant_crosswalk.tsv"))
write_tsv_once(stage_samples, file.path(tmp_root, "audits/stage_sample_manifest.tsv"))
write_tsv_once(stage_census, file.path(tmp_root, "audits/stage_census.tsv"))
write_tsv_once(cohort_stage, file.path(tmp_root, "audits/cohort_stage_census.tsv"))
write_tsv_once(transition_summary, file.path(tmp_root, "audits/transition_census.tsv"))
write_tsv_once(transition_cohort_audit, file.path(tmp_root, "audits/transition_cohort_census.tsv"))
write_tsv_once(transition_sample_manifest, file.path(tmp_root, "audits/transition_sample_manifest.tsv"))
write_tsv_once(design_preflight, file.path(tmp_root, "audits/design_preflight.tsv"))
write_tsv_once(model_design_preflight, file.path(tmp_root, "audits/model_design_preflight.tsv"))

input_manifest <- data.table(
  input_id = names(input_paths),
  candidate_path = substring(input_paths, nchar(candidate_root) + 2L),
  size_bytes = as.numeric(file.info(input_paths)$size),
  sha256 = unname(observed_hashes),
  role = c(
    "primary_stage_native_count_source", "full_resource_filter_reference_only",
    "corrected_sample_metadata", "corrected_fragment_qc",
    "frozen_harmonized_stage_source", "authoritative_gse193066_participant_key"
  )
)
write_tsv_once(input_manifest, file.path(tmp_root, "manifests/input_manifest.tsv"))
producer_code_manifest <- data.table(
  artifact_id = c(
    "producer", "atomic_publisher", "atomic_publisher_io", "Rscript", "Python"
  ),
  role = c(
    "resource_candidate_code", "resource_candidate_code",
    "resource_candidate_code", "analysis_runtime_rscript",
    "resource_candidate_python"
  ),
  path = c(
    producer_path,
    normalizePath(publisher_path, mustWork = TRUE),
    normalizePath(snapshot_io_path, mustWork = TRUE),
    rscript_path,
    python_path
  ),
  sha256 = vapply(
    c(
      producer_path, publisher_path, snapshot_io_path,
      rscript_path, python_path
    ),
    sha256_file,
    character(1)
  )
)
setcolorder(producer_code_manifest, c("artifact_id", "path", "role", "sha256"))
write_tsv_once(
  producer_code_manifest,
  file.path(tmp_root, "manifests/producer_code_manifest.tsv")
)
run_contract <- data.table(
  contract_key = c(
    "candidate_id", "workstream_id", "source_layer", "biological_unit",
    "allowlist", "gse193066_timepoint", "gse193066_participant_key",
    "n_stage_samples", "stage_counts", "gene_filter_status", "normalization_status",
    "model_fit_status", "hypothesis_test_status", "canonical_write",
    "figure_write", "manual_review_required"
  ),
  value = c(
    CANDIDATE_ID, WORKSTREAM_ID, "human_bulk_corrected_ninecohort_candidate",
    "one deposited biopsy per analysis unit; GSE193066 restricted to first biopsy per title-root participant",
    paste(TRUE_KLEINER_COHORTS, collapse = ";"), "1st biopsy only",
    "deposited !Sample_title with terminal _1/_2 removed", "661",
    "F0=126;F1=187;F2=174;F3=132;F4=42", "not_run", "not_run",
    "not_run", "not_run", "false", "false", "true"
  )
)
setnames(run_contract, "contract_key", "key")
write_tsv_once(run_contract, file.path(tmp_root, "manifests/run_contract.tsv"))
environment_manifest <- data.table(
  component = c("R", "data.table"),
  version = c(as.character(getRversion()), as.character(packageVersion("data.table"))),
  hostname = Sys.info()[["nodename"]],
  locale = Sys.getlocale(),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "not_slurm"),
  slurm_partition = Sys.getenv("SLURM_JOB_PARTITION", "not_slurm"),
  slurm_qos = Sys.getenv("SLURM_JOB_QOS", "not_slurm"),
  slurm_nodelist = Sys.getenv("SLURM_JOB_NODELIST", "not_slurm"),
  slurm_cpus_per_task = Sys.getenv("SLURM_CPUS_PER_TASK", "not_slurm"),
  omp_num_threads = Sys.getenv("OMP_NUM_THREADS", "not_set"),
  openblas_num_threads = Sys.getenv("OPENBLAS_NUM_THREADS", "not_set"),
  mkl_num_threads = Sys.getenv("MKL_NUM_THREADS", "not_set"),
  veclib_maximum_threads = Sys.getenv("VECLIB_MAXIMUM_THREADS", "not_set"),
  numexpr_num_threads = Sys.getenv("NUMEXPR_NUM_THREADS", "not_set")
)
write_tsv_once(environment_manifest, file.path(tmp_root, "manifests/environment_manifest.tsv"))
producer_status <- data.table(
  status = "PREFLIGHT_COMPLETE_PENDING_INDEPENDENT_VALIDATION",
  candidate_id = CANDIDATE_ID,
  workstream_id = WORKSTREAM_ID,
  n_stage_samples = nrow(stage_samples),
  n_transition_rows = nrow(transition_sample_manifest),
  gene_filter_status = "not_run",
  model_fit_status = "not_run",
  manual_review_required = TRUE
)
write_tsv_once(producer_status, file.path(tmp_root, "manifests/producer_status.tsv"))

payload_paths <- sort(list.files(tmp_root, recursive = TRUE, full.names = TRUE))
payload_paths <- payload_paths[!file.info(payload_paths)$isdir]
producer_manifest <- data.table(
  candidate_path = substring(payload_paths, nchar(tmp_root) + 2L),
  size_bytes = as.numeric(file.info(payload_paths)$size),
  sha256 = vapply(payload_paths, sha256_file, character(1))
)
write_tsv_once(producer_manifest, file.path(tmp_root, "manifests/producer_artifact_manifest.tsv"))

publication <- system2(
  python_path,
  c(
    "-S",
    shQuote(normalizePath(publisher_path, mustWork = TRUE)),
    "--marker",
    "manifests/producer_artifact_manifest.tsv",
    shQuote(normalizePath(tmp_root, mustWork = TRUE)),
    shQuote(final_root)
  ),
  stdout = TRUE,
  stderr = TRUE
)
publication_status <- attr(publication, "status")
source_marker <- file.path(tmp_root, "manifests/producer_artifact_manifest.tsv")
official_marker <- file.path(final_root, "manifests/producer_artifact_manifest.tsv")
if ((!is.null(publication_status) && publication_status != 0L) ||
    !dir.exists(final_root) || is_symlink(final_root) ||
    !dir.exists(tmp_root) || is_symlink(tmp_root) ||
    !file_test("-f", source_marker) || is_symlink(source_marker) ||
    !file_test("-f", official_marker) || is_symlink(official_marker) ||
    !identical(sha256_file(source_marker), sha256_file(official_marker))) {
  stop(
    "Staged no-replace publication of stage preflight failed: ",
    paste(publication, collapse = " | "),
    call. = FALSE
  )
}
cat("PREFLIGHT_COMPLETE_PENDING_INDEPENDENT_VALIDATION: ", final_root, "\n", sep = "")
