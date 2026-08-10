#!/usr/bin/env Rscript
# Candidate-only adjacent-fibrosis LVQW refit on six true-Kleiner cohorts.
#
# This producer deliberately does not read or write any canonical disease-
# signature result. It reuses the vetted de_engine_lvqw.R implementation and
# writes one immutable candidate bundle for later Plan-60 adjudication.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ashr)
})

PROJECT_ROOT <- normalizePath(
  Sys.getenv(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
  ),
  mustWork = TRUE
)
INT <- file.path(
  PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration"
)
SCRIPT_DIR <- file.path(INT, "scripts")
INPUT_DIR <- file.path(INT, "results/integration")
CANDIDATE_ID <- "program-context-v2-candidate-2026-08-07"
ANALYSIS_ID <- "fibrosis-adjacent-true-kleiner-lvqw-v1"
FINAL_ROOT <- file.path(
  INT, "results/candidates", CANDIDATE_ID, ANALYSIS_ID
)

PRODUCER <- file.path(SCRIPT_DIR, "14h_true_kleiner_adjacent_lvqw_candidate.R")
VALIDATOR <- file.path(SCRIPT_DIR, "validate_true_kleiner_adjacent_lvqw_candidate.py")
WRAPPER <- file.path(SCRIPT_DIR, "run_true_kleiner_adjacent_lvqw_candidate.sbatch")
ENGINE <- file.path(SCRIPT_DIR, "de_engine_lvqw.R")
DGE_PATH <- file.path(INPUT_DIR, "merged_dge.rds")
META_PATH <- file.path(INPUT_DIR, "meta_matched.rds")
QC_PATH <- file.path(INT, "qc/sample_qc_report.csv")
GSE193066_METADATA_PATH <- file.path(
  PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/metadata/metadata.tsv"
)

TRUE_KLEINER_COHORTS <- c(
  "GSE130970", "GSE135251", "GSE162694",
  "GSE174478", "GSE193066", "GSE240729"
)
HARD_EXCLUDED_COHORTS <- c("GSE213621", "PRJNA512027")
EXPECTED_N_GENES <- 27638L
SHUFFLE_SEED <- 20260808L
PREFLIGHT_ONLY <- identical(tolower(Sys.getenv("MASLD_PREFLIGHT_ONLY", "false")), "true")

EXPECTED_TRANSITIONS <- data.table(
  transition = c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"),
  contrast = c("F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3"),
  low_stage = 0:3,
  high_stage = 1:4,
  n_low = c(126L, 187L, 174L, 108L),
  n_high = c(187L, 174L, 133L, 44L),
  n_samples = c(313L, 361L, 307L, 152L),
  n_cohorts = c(6L, 6L, 6L, 5L)
)
EXPECTED_COHORTS <- list(
  F0_to_F1 = TRUE_KLEINER_COHORTS,
  F1_to_F2 = TRUE_KLEINER_COHORTS,
  F2_to_F3 = TRUE_KLEINER_COHORTS,
  F3_to_F4 = setdiff(TRUE_KLEINER_COHORTS, "GSE193066")
)

stopifnot(
  dir.exists(INT),
  all(file.exists(c(
    PRODUCER, VALIDATOR, WRAPPER, ENGINE, DGE_PATH, META_PATH, QC_PATH,
    GSE193066_METADATA_PATH
  )))
)
is_symlink <- function(path) {
  target <- Sys.readlink(path)
  !is.na(target) && nzchar(target)
}
if (file.exists(FINAL_ROOT) || is_symlink(FINAL_ROOT)) {
  stop("Refusing existing or symlinked candidate root: ", FINAL_ROOT)
}

source(ENGINE)

sha256_file <- function(path) {
  normalized <- normalizePath(path, mustWork = TRUE)
  value <- system2("sha256sum", shQuote(normalized), stdout = TRUE, stderr = TRUE)
  status <- attr(value, "status")
  if ((!is.null(status) && status != 0L) || length(value) != 1L) {
    stop("sha256sum failed for ", normalized)
  }
  digest <- strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", digest)) stop("Invalid SHA256 for ", normalized)
  digest
}

write_tsv <- function(object, path) {
  if (file.exists(path) || is_symlink(path)) stop("Refusing overwrite: ", path)
  fwrite(object, path, sep = "\t", quote = FALSE, na = "NA")
}

write_csv <- function(object, path) {
  if (file.exists(path) || is_symlink(path)) stop("Refusing overwrite: ", path)
  fwrite(object, path, quote = FALSE, na = "NA")
}

write_lines_once <- function(lines, path) {
  if (file.exists(path) || is_symlink(path)) stop("Refusing overwrite: ", path)
  writeLines(lines, path, useBytes = TRUE)
}

file_manifest <- function(paths, ids) {
  stopifnot(length(paths) == length(ids), all(file.exists(paths)))
  info <- file.info(paths)
  data.table(
    artifact_id = ids,
    path = vapply(paths, function(x) {
      normalized <- normalizePath(x, mustWork = TRUE)
      prefix <- paste0(PROJECT_ROOT, .Platform$file.sep)
      if (!startsWith(normalized, prefix)) stop("Manifest path escapes project root: ", normalized)
      substring(normalized, nchar(prefix) + 1L)
    }, character(1)),
    bytes = as.numeric(info$size),
    sha256 = vapply(paths, sha256_file, character(1))
  )
}

cat("=== Candidate true-Kleiner adjacent-fibrosis LVQW refit ===\n")
cat("Candidate root:", FINAL_ROOT, "\n")
cat("Preflight only:", PREFLIGHT_ONLY, "\n")

dge_all <- readRDS(DGE_PATH)
meta_all <- as.data.table(readRDS(META_PATH))
qc <- fread(QC_PATH)
gse193066_source_raw <- fread(GSE193066_METADATA_PATH, check.names = FALSE)

required_meta <- c("sample_id", "dataset", "fibrosis_stage", "inferred_sex")
required_qc <- c("sample_id", "pass_technical")
required_gse193066 <- c(
  "Run", "sample_id", "!Sample_title", "biopsy", "fibrosis stage"
)
if (!all(required_meta %in% names(meta_all))) stop("Metadata schema is incomplete")
if (!all(required_qc %in% names(qc))) stop("QC schema is incomplete")
if (!all(required_gse193066 %in% names(gse193066_source_raw))) {
  stop("GSE193066 deposited metadata schema is incomplete")
}
if (anyDuplicated(meta_all$sample_id)) stop("Metadata sample_id is not unique")
if (anyDuplicated(qc$sample_id)) stop("QC sample_id is not unique")
if (anyDuplicated(colnames(dge_all))) stop("DGE sample names are not unique")
if (nrow(dge_all) != EXPECTED_N_GENES) {
  stop("Current DGE gene universe drift: observed ", nrow(dge_all),
       ", expected ", EXPECTED_N_GENES)
}

# GSE193066 contains 106 first biopsies and 58 paired second biopsies.  The
# deposited title suffix (_1/_2) is the authoritative donor key.  Keep one
# first biopsy per donor so repeated biopsies never enter an independent-sample
# LVQW fit.
gse193066_source <- gse193066_source_raw[, .(
  run_id = as.character(Run),
  metadata_sample_id = as.character(sample_id),
  source_title = as.character(get("!Sample_title")),
  biopsy = as.character(biopsy),
  source_fibrosis_stage = as.integer(get("fibrosis stage"))
)]
if (nrow(gse193066_source) != 164L) stop("GSE193066 source row count drift")
if (anyNA(gse193066_source) || any(!nzchar(gse193066_source$source_title))) {
  stop("GSE193066 source metadata contains missing required values")
}
if (any(gse193066_source$run_id != gse193066_source$metadata_sample_id)) {
  stop("GSE193066 Run/sample_id disagreement")
}
if (anyDuplicated(gse193066_source$run_id) ||
    anyDuplicated(gse193066_source$source_title)) {
  stop("GSE193066 deposited run/title key is not unique")
}
if (any(!grepl("^HUnafld[0-9]{3}(_[12])?$", gse193066_source$source_title))) {
  stop("GSE193066 deposited title does not match the donor-key contract")
}
if (!setequal(unique(gse193066_source$biopsy), c("1st biopsy", "2nd biopsy"))) {
  stop("GSE193066 biopsy labels drifted")
}
gse193066_source[, donor_token := sub("_[12]$", "", source_title)]
gse193066_source[, donor_id := paste("GSE193066", donor_token, sep = "::")]
gse193066_source[, is_first_biopsy := biopsy == "1st biopsy"]
if (sum(gse193066_source$is_first_biopsy) != 106L ||
    sum(!gse193066_source$is_first_biopsy) != 58L ||
    uniqueN(gse193066_source$donor_id) != 106L) {
  stop("GSE193066 first/second-biopsy census drift")
}
gse193066_donor_gate <- gse193066_source[, .(
  n_rows = .N,
  n_first = sum(is_first_biopsy),
  n_second = sum(!is_first_biopsy)
), by = donor_id]
if (any(gse193066_donor_gate$n_first != 1L) ||
    any(!gse193066_donor_gate$n_rows %in% c(1L, 2L)) ||
    any(gse193066_donor_gate$n_second != gse193066_donor_gate$n_rows - 1L) ||
    sum(gse193066_donor_gate$n_rows == 2L) != 58L) {
  stop("GSE193066 donor pairing is not one first biopsy plus at most one second biopsy")
}
if (any(gse193066_source[is_first_biopsy == TRUE, grepl("_2$", source_title)]) ||
    any(gse193066_source[is_first_biopsy == FALSE, !grepl("_2$", source_title)])) {
  stop("GSE193066 biopsy label/title suffix disagreement")
}

gse193066_join <- gse193066_source[, .(
  sample_id = run_id,
  source_title,
  donor_id,
  biopsy,
  is_first_biopsy,
  source_fibrosis_stage
)]
if (!setequal(meta_all[dataset == "GSE193066", sample_id], gse193066_join$sample_id)) {
  stop("GSE193066 harmonized/source sample universe disagreement")
}
meta_all[, `:=`(
  donor_id = paste(dataset, sample_id, sep = "::"),
  donor_id_source = "single_biopsy_sample_id",
  source_title = sample_id,
  biopsy = "single deposited biopsy",
  is_first_biopsy = TRUE
)]
meta_all[gse193066_join, on = .(sample_id), `:=`(
  donor_id = i.donor_id,
  donor_id_source = "deposited_title_root",
  source_title = i.source_title,
  biopsy = i.biopsy,
  is_first_biopsy = i.is_first_biopsy
)]
gse_meta_check <- meta_all[dataset == "GSE193066"][
  gse193066_join, on = .(sample_id)
]
if (anyNA(gse_meta_check$source_fibrosis_stage) ||
    any(gse_meta_check$fibrosis_stage != gse_meta_check$source_fibrosis_stage)) {
  stop("GSE193066 harmonized/source fibrosis-stage disagreement")
}

pass_ids <- qc[pass_technical == TRUE, sample_id]
meta_pass <- meta_all[
  sample_id %in% pass_ids &
    sample_id %in% colnames(dge_all) &
    dataset %in% TRUE_KLEINER_COHORTS &
    !is.na(fibrosis_stage) &
    (dataset != "GSE193066" | is_first_biopsy)
]
meta_pass[, fibrosis_stage := as.integer(fibrosis_stage)]
setorder(meta_pass, dataset, fibrosis_stage, sample_id)

if (!setequal(unique(meta_pass$dataset), TRUE_KLEINER_COHORTS)) {
  stop("True-Kleiner allowlist did not resolve exactly")
}
if (any(meta_pass$dataset %in% HARD_EXCLUDED_COHORTS)) {
  stop("Hard-excluded cohort entered the candidate universe")
}
if (any(!meta_pass$fibrosis_stage %in% 0:4)) stop("Invalid exact fibrosis stage")
if (anyNA(meta_pass$inferred_sex)) stop("Missing inferred_sex in candidate universe")
if (anyDuplicated(meta_pass$sample_id)) stop("Candidate sample_id is not unique")
if (anyDuplicated(meta_pass$donor_id)) stop("A biological donor occurs more than once")
if (any(meta_pass$dataset == "GSE193066" & meta_pass$biopsy != "1st biopsy")) {
  stop("A GSE193066 second biopsy entered the candidate universe")
}

gse193066_crosswalk <- merge(
  gse193066_join,
  meta_all[dataset == "GSE193066", .(
    sample_id,
    harmonized_fibrosis_stage = fibrosis_stage,
    inferred_sex
  )],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
gse193066_crosswalk <- merge(
  gse193066_crosswalk,
  qc[, .(sample_id, pass_technical)],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
gse193066_crosswalk[, dge_present := sample_id %in% colnames(dge_all)]
gse193066_crosswalk[, included_primary :=
  is_first_biopsy & pass_technical & dge_present &
    !is.na(harmonized_fibrosis_stage)]
gse193066_crosswalk[, exclusion_reason := fcase(
  included_primary, "included_first_biopsy",
  !is_first_biopsy, "excluded_second_biopsy",
  !pass_technical, "excluded_failed_technical_qc",
  !dge_present, "excluded_absent_from_dge",
  is.na(harmonized_fibrosis_stage), "excluded_missing_exact_stage",
  default = "excluded_unclassified"
)]
setorder(gse193066_crosswalk, donor_id, biopsy, sample_id)
if (nrow(gse193066_crosswalk) != 164L ||
    sum(gse193066_crosswalk$included_primary) != 105L ||
    uniqueN(gse193066_crosswalk[included_primary == TRUE, donor_id]) != 105L ||
    any(gse193066_crosswalk[included_primary == TRUE, biopsy] != "1st biopsy")) {
  stop("GSE193066 primary donor crosswalk census failed")
}

subset_dge <- function(metadata) {
  idx <- colnames(dge_all) %in% metadata$sample_id
  dge <- dge_all[, idx]
  covariates <- metadata[match(colnames(dge), sample_id)]
  if (anyNA(covariates$sample_id)) stop("DGE-to-metadata join failed")
  dge$samples$dataset <- factor(covariates$dataset, levels = sort(TRUE_KLEINER_COHORTS))
  dge$samples$inferred_sex <- factor(covariates$inferred_sex)
  dge$samples$fibrosis_stage <- covariates$fibrosis_stage
  dge
}

preflight_rows <- list()
pair_metadata <- list()
for (i in seq_len(nrow(EXPECTED_TRANSITIONS))) {
  expected <- EXPECTED_TRANSITIONS[i]
  expected_cohorts <- EXPECTED_COHORTS[[expected$transition]]
  eligible <- meta_pass[fibrosis_stage %in% c(expected$low_stage, expected$high_stage)]
  cohort_gate <- eligible[, .(
    n_low = sum(fibrosis_stage == expected$low_stage),
    n_high = sum(fibrosis_stage == expected$high_stage),
    n_stages = uniqueN(fibrosis_stage)
  ), by = dataset]
  cohort_gate <- cohort_gate[n_low > 0L & n_high > 0L & n_stages == 2L]
  selected <- eligible[dataset %in% cohort_gate$dataset]
  observed <- data.table(
    transition = expected$transition,
    contrast = expected$contrast,
    low_stage = expected$low_stage,
    high_stage = expected$high_stage,
    n_low = sum(selected$fibrosis_stage == expected$low_stage),
    n_high = sum(selected$fibrosis_stage == expected$high_stage),
    n_samples = nrow(selected),
    n_cohorts = uniqueN(selected$dataset)
  )
  compare_cols <- c(
    "transition", "contrast", "low_stage", "high_stage",
    "n_low", "n_high", "n_samples", "n_cohorts"
  )
  if (!identical(as.list(observed[, ..compare_cols]), as.list(expected[, ..compare_cols]))) {
    stop("Predeclared transition census mismatch for ", expected$transition)
  }
  if (!setequal(unique(selected$dataset), expected_cohorts)) {
    stop("Transition-specific true-Kleiner cohort set drift for ", expected$transition)
  }
  if (any(selected$dataset %in% HARD_EXCLUDED_COHORTS)) {
    stop("Hard exclusion failed for ", expected$transition)
  }
  preflight_rows[[expected$transition]] <- observed
  pair_metadata[[expected$transition]] <- selected
}
preflight <- rbindlist(preflight_rows)
print(preflight)

if (PREFLIGHT_ONLY) {
  cat("PREFLIGHT_PASS: donor-deduplicated 6/6/6/5-cohort arm census verified; no files written.\n")
  quit(save = "no", status = 0L)
}

dir.create(dirname(FINAL_ROOT), recursive = TRUE, showWarnings = FALSE)
job_tag <- Sys.getenv("SLURM_JOB_ID", as.character(Sys.getpid()))
TMP_ROOT <- paste0(FINAL_ROOT, ".tmp.", job_tag)
if (file.exists(TMP_ROOT) || is_symlink(TMP_ROOT)) {
  stop("Refusing existing temporary root: ", TMP_ROOT)
}
dir.create(file.path(TMP_ROOT, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(TMP_ROOT, "audits"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(TMP_ROOT, "manifests"), recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(TMP_ROOT)) stop("Failed to create temporary candidate root")

main_results <- list()
null_results <- list()
summary_rows <- list()
null_summary_rows <- list()
sample_rows <- list()
cohort_rows <- list()
design_rows <- list()
model_design_rows <- list()
shuffle_rows <- list()

for (i in seq_len(nrow(EXPECTED_TRANSITIONS))) {
  expected <- EXPECTED_TRANSITIONS[i]
  transition <- expected$transition
  contrast <- expected$contrast
  expected_cohorts <- EXPECTED_COHORTS[[transition]]
  selected <- copy(pair_metadata[[transition]])
  dge <- subset_dge(selected)
  fib_group <- factor(
    ifelse(dge$samples$fibrosis_stage == expected$high_stage, "high", "low"),
    levels = c("low", "high")
  )
  info <- data.frame(
    dataset = factor(dge$samples$dataset, levels = sort(TRUE_KLEINER_COHORTS)),
    inferred_sex = factor(dge$samples$inferred_sex),
    fib_group = fib_group,
    row.names = colnames(dge)
  )
  guarded <- build_design_guarded(info, c("dataset", "inferred_sex", "fib_group"))
  coefficient <- "fib_grouphigh"
  if (length(guarded$dropped) != 0L) stop("Design term was dropped for ", transition)
  if (!identical(guarded$formula_used, "~ dataset + inferred_sex + fib_group")) {
    stop("Unexpected design formula for ", transition, ": ", guarded$formula_used)
  }
  if (!coefficient %in% colnames(guarded$design)) stop("Coefficient absent for ", transition)
  if (qr(guarded$design)$rank != ncol(guarded$design)) stop("Rank-deficient design for ", transition)

  result <- fit_lvqw(dge, guarded$design, coefficient, do_ashr = TRUE)
  if (nrow(result) != EXPECTED_N_GENES || !identical(result$gene, rownames(dge))) {
    stop("Incomplete tested family for ", transition)
  }
  bh <- p.adjust(result$P.Value, method = "BH")
  if (max(abs(bh - result$padj), na.rm = TRUE) > 1e-12) {
    stop("BH rederivation failed for ", transition)
  }
  result[, `:=`(
    transition = transition,
    contrast = contrast,
    method = "limma_voom_qw_C2_true_kleiner_first_biopsy",
    n_cohorts = expected$n_cohorts,
    n_samples = expected$n_samples,
    n_low = expected$n_low,
    n_high = expected$n_high
  )]
  setcolorder(result, c(
    "gene", "logFC", "SE", "t", "P.Value", "padj", "shrunk_logFC",
    "lfsr", "AveExpr", "contrast", "transition", "method", "n_cohorts",
    "n_samples", "n_low", "n_high"
  ))
  main_results[[transition]] <- result

  summary_rows[[transition]] <- data.table(
    transition = transition,
    contrast = contrast,
    n_genes_tested = nrow(result),
    n_deg_05 = sum(result$padj < 0.05),
    n_up = sum(result$padj < 0.05 & result$logFC > 0),
    n_down = sum(result$padj < 0.05 & result$logFC < 0),
    n_low = expected$n_low,
    n_high = expected$n_high,
    n_samples = expected$n_samples,
    n_cohorts = expected$n_cohorts,
    significance_rule = "BH padj < 0.05",
    lfc_gate = "none"
  )

  ordered_samples <- selected[match(colnames(dge), sample_id)]
  sample_rows[[transition]] <- ordered_samples[, .(
    sample_id,
    donor_id,
    donor_id_source,
    dataset,
    source_title,
    biopsy,
    fibrosis_stage,
    transition = transition,
    contrast = contrast,
    arm = ifelse(fibrosis_stage == expected$high_stage, "high", "low"),
    inferred_sex,
    pass_technical = TRUE
  )]
  if (anyDuplicated(sample_rows[[transition]]$donor_id)) {
    stop("A biological donor occurs more than once in ", transition)
  }
  cohort_rows[[transition]] <- sample_rows[[transition]][, .(
    n_low = sum(arm == "low"),
    n_high = sum(arm == "high"),
    n_total = .N,
    both_arms_present = all(c("low", "high") %in% arm)
  ), by = .(transition, contrast, dataset)]
  if (any(!cohort_rows[[transition]]$both_arms_present)) {
    stop("A selected cohort lacks both arms for ", transition)
  }

  design_sample_ids <- rownames(guarded$design)
  model_design <- as.data.table(guarded$design)
  model_design[, sample_id := design_sample_ids]
  model_design[, `:=`(
    transition = transition,
    contrast = contrast,
    arm = as.character(info$fib_group)
  )]
  setcolorder(model_design, c(
    "sample_id", "transition", "contrast", "arm",
    setdiff(names(model_design), c("sample_id", "transition", "contrast", "arm"))
  ))
  model_design_rows[[transition]] <- model_design
  design_rows[[transition]] <- data.table(
    transition = transition,
    contrast = contrast,
    formula = guarded$formula_used,
    coefficient = coefficient,
    dropped_terms = "",
    design_rank = qr(guarded$design)$rank,
    n_design_columns = ncol(guarded$design),
    design_columns = paste(colnames(guarded$design), collapse = ";"),
    n_genes_tested = nrow(result),
    n_samples = ncol(dge),
    n_cohorts = uniqueN(info$dataset),
    cohort_set = paste(sort(unique(as.character(info$dataset))), collapse = ";"),
    expected_cohort_set = paste(sort(expected_cohorts), collapse = ";"),
    n_unique_donors = uniqueN(selected$donor_id)
  )

  permuted <- as.character(info$fib_group)
  for (j in seq_along(expected_cohorts)) {
    cohort <- sort(expected_cohorts)[[j]]
    index <- which(as.character(info$dataset) == cohort)
    set.seed(SHUFFLE_SEED + i * 100L + j)
    shuffled <- sample(permuted[index], length(index), replace = FALSE)
    shuffle_rows[[paste(transition, cohort, sep = "__")]] <- data.table(
      transition = transition,
      contrast = contrast,
      dataset = cohort,
      original_low = sum(permuted[index] == "low"),
      original_high = sum(permuted[index] == "high"),
      shuffled_low = sum(shuffled == "low"),
      shuffled_high = sum(shuffled == "high"),
      n_changed = sum(shuffled != permuted[index]),
      seed = SHUFFLE_SEED + i * 100L + j
    )
    permuted[index] <- shuffled
  }
  if (sum(vapply(
    shuffle_rows[startsWith(names(shuffle_rows), paste0(transition, "__"))],
    function(x) x$n_changed,
    integer(1)
  )) == 0L) stop("Deterministic shuffle changed no labels for ", transition)

  null_info <- info
  null_info$fib_group <- factor(permuted, levels = c("low", "high"))
  null_guarded <- build_design_guarded(
    null_info, c("dataset", "inferred_sex", "fib_group")
  )
  if (length(null_guarded$dropped) != 0L ||
      !identical(null_guarded$formula_used, guarded$formula_used) ||
      qr(null_guarded$design)$rank != ncol(null_guarded$design) ||
      !coefficient %in% colnames(null_guarded$design)) {
    stop("Shuffled-label design contract failed for ", transition)
  }
  null_result <- fit_lvqw(
    dge, null_guarded$design, coefficient, do_ashr = FALSE
  )
  null_bh <- p.adjust(null_result$P.Value, method = "BH")
  if (nrow(null_result) != EXPECTED_N_GENES ||
      max(abs(null_bh - null_result$padj), na.rm = TRUE) > 1e-12) {
    stop("Shuffled-label result contract failed for ", transition)
  }
  null_result <- null_result[, .(gene, logFC, SE, t, P.Value, padj, AveExpr)]
  null_result[, `:=`(
    transition = transition,
    contrast = contrast,
    shuffle_seed_base = SHUFFLE_SEED
  )]
  null_results[[transition]] <- null_result
  null_summary_rows[[transition]] <- data.table(
    transition = transition,
    contrast = contrast,
    n_genes_tested = nrow(null_result),
    n_deg_05 = sum(null_result$padj < 0.05),
    median_logFC = median(null_result$logFC),
    mean_logFC = mean(null_result$logFC),
    median_abs_logFC = median(abs(null_result$logFC)),
    max_abs_logFC = max(abs(null_result$logFC)),
    interpretation = "diagnostic shuffle only; biological positivity is not an acceptance gate"
  )
}

main <- rbindlist(main_results, use.names = TRUE)
null <- rbindlist(null_results, use.names = TRUE)
summary <- rbindlist(summary_rows, use.names = TRUE)
null_summary <- rbindlist(null_summary_rows, use.names = TRUE)
sample_manifest <- rbindlist(sample_rows, use.names = TRUE)
cohort_audit <- rbindlist(cohort_rows, use.names = TRUE)
design_audit <- rbindlist(design_rows, use.names = TRUE)
model_design <- rbindlist(model_design_rows, use.names = TRUE, fill = TRUE)
shuffle_audit <- rbindlist(shuffle_rows, use.names = TRUE)

write_csv(main, file.path(TMP_ROOT, "results/fibrosis_consecutive_lvqw_true_kleiner.csv"))
write_tsv(summary, file.path(TMP_ROOT, "results/transition_summary.tsv"))
write_csv(null, file.path(TMP_ROOT, "audits/null_shuffle_results.csv"))
write_tsv(null_summary, file.path(TMP_ROOT, "audits/null_shuffle_summary.tsv"))
write_tsv(sample_manifest, file.path(TMP_ROOT, "audits/sample_manifest.tsv"))
write_tsv(cohort_audit, file.path(TMP_ROOT, "audits/cohort_arm_audit.tsv"))
write_tsv(design_audit, file.path(TMP_ROOT, "audits/design_audit.tsv"))
write_tsv(model_design, file.path(TMP_ROOT, "audits/model_design.tsv"))
write_tsv(shuffle_audit, file.path(TMP_ROOT, "audits/shuffle_audit.tsv"))
write_tsv(
  gse193066_crosswalk,
  file.path(TMP_ROOT, "audits/gse193066_donor_biopsy_crosswalk.tsv")
)

input_manifest <- file_manifest(
  c(DGE_PATH, META_PATH, QC_PATH, GSE193066_METADATA_PATH),
  c("merged_dge", "meta_matched", "sample_qc", "gse193066_source_metadata")
)
code_manifest <- file_manifest(
  c(PRODUCER, VALIDATOR, WRAPPER, ENGINE),
  c("producer", "validator", "slurm_wrapper", "lvqw_engine")
)
write_tsv(input_manifest, file.path(TMP_ROOT, "manifests/input_manifest.tsv"))
write_tsv(code_manifest, file.path(TMP_ROOT, "manifests/code_manifest.tsv"))

package_names <- c("R", "data.table", "edgeR", "limma", "ashr")
package_versions <- c(
  as.character(getRversion()),
  vapply(package_names[-1L], function(x) as.character(packageVersion(x)), character(1))
)
environment_manifest <- data.table(
  component = package_names,
  version = package_versions,
  hostname = Sys.info()[["nodename"]],
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", "not_slurm"),
  slurm_partition = Sys.getenv("SLURM_JOB_PARTITION", "not_slurm"),
  slurm_cpus_per_task = Sys.getenv("SLURM_CPUS_PER_TASK", "not_slurm"),
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write_tsv(environment_manifest, file.path(TMP_ROOT, "manifests/environment_manifest.tsv"))
write_lines_once(capture.output(sessionInfo()), file.path(TMP_ROOT, "manifests/sessionInfo.txt"))

run_contract <- data.table(
  contract_key = c(
    "candidate_id", "analysis_id", "model", "coefficient", "method",
    "tested_universe", "significance_rule", "lfc_gate", "allowlist",
    "hard_exclusions", "gse193066_timepoint", "gse193066_donor_key",
    "shuffle_seed", "replication_unit", "canonical_write"
  ),
  value = c(
    CANDIDATE_ID,
    ANALYSIS_ID,
    "~ dataset + inferred_sex + fib_group",
    "fib_grouphigh",
    "limma_voom_qw_C2_true_kleiner_first_biopsy",
    as.character(EXPECTED_N_GENES),
    "BH padj < 0.05 per complete contrast family",
    "none",
    paste(TRUE_KLEINER_COHORTS, collapse = ";"),
    paste(HARD_EXCLUDED_COHORTS, collapse = ";"),
    "1st biopsy only",
    "deposited !Sample_title with terminal _1/_2 removed",
    as.character(SHUFFLE_SEED),
    "one first-biopsy biological donor/sample",
    "false"
  )
)
setnames(run_contract, "contract_key", "key")
write_tsv(run_contract, file.path(TMP_ROOT, "manifests/run_contract.tsv"))

producer_status <- data.table(
  status = "PRODUCER_COMPLETE_PENDING_VALIDATION",
  candidate_id = CANDIDATE_ID,
  analysis_id = ANALYSIS_ID,
  n_contrasts = nrow(summary),
  n_result_rows = nrow(main),
  n_null_rows = nrow(null),
  n_unique_primary_donors = uniqueN(meta_pass$donor_id),
  canonical_write = FALSE,
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write_tsv(producer_status, file.path(TMP_ROOT, "manifests/producer_status.tsv"))

artifact_manifest_path <- file.path(TMP_ROOT, "manifests/producer_artifact_manifest.tsv")
artifact_paths <- sort(list.files(TMP_ROOT, recursive = TRUE, full.names = TRUE))
artifact_paths <- artifact_paths[file.info(artifact_paths)$isdir == FALSE]
artifact_roles <- ifelse(
  grepl("/results/", artifact_paths), "result",
  ifelse(grepl("/audits/", artifact_paths), "audit", "manifest")
)
artifact_manifest <- data.table(
  artifact_path = substring(artifact_paths, nchar(TMP_ROOT) + 2L),
  artifact_role = artifact_roles,
  bytes = as.numeric(file.info(artifact_paths)$size),
  sha256 = vapply(artifact_paths, sha256_file, character(1))
)
write_tsv(artifact_manifest, artifact_manifest_path)

if (!file.rename(TMP_ROOT, FINAL_ROOT)) {
  stop("Atomic candidate-root publication failed: ", TMP_ROOT, " -> ", FINAL_ROOT)
}
cat("PRODUCER_COMPLETE_PENDING_VALIDATION:", FINAL_ROOT, "\n")
