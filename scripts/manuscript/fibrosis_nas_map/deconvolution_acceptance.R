#!/usr/bin/env Rscript

# Precontract composition acceptance for BULK-PROGRAM-MAP-v9.
#
# WHY THIS EXISTS
# The v8 composition arm used the deposited MuSiC estimates. Those estimates are
# degenerate on the real bulk substrate: median hepatocyte 0.942, with
# cholangiocytes exactly zero in 97.8% of samples, fibroblasts 87%, macrophages
# 73% and endothelial cells 65%. Because cholangiocyte and fibroblast fractions
# are both identically zero in GSE135251, the 1e-6 pseudocount maps their two
# hepatocyte log-ratios onto the same vector, the design matrix loses rank, and
# GSE135251 (214 of 469 discovery samples) silently left the composition model.
# The project's pseudobulk benchmark does not catch this: MuSiC is well
# calibrated there (hepatocyte estimate 0.626 against a truth of 0.623). The
# benchmark is a cell-number simplex while bulk RNA is mRNA-mass weighted, so it
# does not transfer to the real substrate.
#
# This script therefore selects the composition estimator on real-bulk behaviour
# using criteria fixed below, and seals the decision before any v9 program model
# is fitted.
#
# CRITERIA PROVENANCE, STATED PLAINLY
# Criteria 1, 2 and 4 were evaluated during planning on 2026-08-12 and are not
# blind. The thresholds are set to sit far from every observed value so the
# decision does not turn on where a cut was placed: MuSiC fails criterion 1 at a
# maximum zero fraction of 0.978 against a threshold of 0.20, and BayesPrism
# records 0.000. Criterion 4 uses Kleiner stage, which is an outcome. It
# constrains only the composition estimator, never program selection, but its use
# is disclosed here and must be repeated wherever a composition-adjusted program
# result is reported.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))

options(digits = 17)

# ---------------------------------------------------------------- frozen rules

MAX_ZERO_FRACTION <- 0.20
MIN_COHORT_COVERAGE <- 0.95
HEPATOCYTE_MEDIAN_RANGE <- c(0.60, 0.95)
MIN_DIRECTION_COHORTS <- 3L

REQUIRED_COHORTS <- c(discovery_cohorts, holdout_cohort)
COVARIATE_LINEAGES <- c(
  cholangiocytes = "Cholangiocytes",
  fibroblasts = "Fibroblasts",
  macrophages = "Macrophages",
  tcells = "T cells"
)
DENOMINATOR_LINEAGE <- "Hepatocytes"
DIRECTION_EXPECTED <- c(
  Hepatocytes = -1, Fibroblasts = 1, Cholangiocytes = 1
)

# --------------------------------------------------------------- method inputs

deconvolution_root <- file.path(project_root, "Analysis/Deconvolution/results")

read_per_cohort <- function(suffix) {
  cohorts <- list.dirs(deconvolution_root, recursive = FALSE, full.names = FALSE)
  rows <- list()
  for (cohort in cohorts) {
    path <- file.path(deconvolution_root, cohort, paste0(cohort, suffix))
    if (!file.exists(path)) next
    value <- fread(path, check.names = FALSE)
    if (names(value)[[1L]] %in% c("V1", "")) setnames(value, 1L, "sample_id")
    if (!"sample_id" %in% names(value)) setnames(value, 1L, "sample_id")
    value[, sample_id := as.character(sample_id)]
    value[, dataset := cohort]
    rows[[length(rows) + 1L]] <- value
  }
  if (!length(rows)) return(NULL)
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

read_music_deposited <- function() {
  value <- fread(composition_candidate, check.names = FALSE)
  contract <- fread(composition_testability)
  keep <- contract[testability == "testable", celltype]
  value[, c("sample_id", "dataset", keep), with = FALSE]
}

# --------------------------------------------------------------------- harness

census_samples <- function() {
  dge <- readRDS(input_dge)
  meta <- as.data.table(readRDS(input_meta))
  samples <- colnames(dge)
  dataset <- as.character(dge$samples$dataset)
  if (length(dataset) != length(samples) || anyNA(dataset)) {
    dataset <- meta$dataset[match(samples, meta$sample_id)]
  }
  assert_true(!anyNA(dataset), "Every census sample must carry a dataset")
  rm(dge); gc()
  data.table(sample_id = samples, dataset = dataset)
}

lineage_columns <- function(value) {
  setdiff(names(value)[vapply(value, is.numeric, logical(1))], c("sample_id", "dataset"))
}

evaluate_method <- function(label, proportions, census, phenotype) {
  if (is.null(proportions)) {
    return(list(
      summary = data.table(method = label, eligible = FALSE,
                           failure = "method_output_absent"),
      lineage = data.table(), direction = data.table(), matched = NULL
    ))
  }
  proportions <- unique(proportions, by = c("sample_id", "dataset"))
  matched <- merge(census, proportions, by = c("sample_id", "dataset"), all.x = TRUE)
  matched <- matched[match(census$sample_id, sample_id)]
  lineages <- lineage_columns(matched)
  present <- rowSums(!is.na(matched[, ..lineages])) == length(lineages)

  coverage <- matched[, .(
    n = .N, n_present = sum(present[.I]), coverage = mean(present[.I])
  ), by = dataset]

  observed <- matched[present]
  M <- as.matrix(observed[, ..lineages])
  storage.mode(M) <- "double"

  lineage_audit <- data.table(
    method = label,
    celltype = lineages,
    zero_fraction = round(colMeans(M == 0, na.rm = TRUE), 4),
    median_fraction = round(apply(M, 2L, median, na.rm = TRUE), 5),
    max_fraction = round(apply(M, 2L, max, na.rm = TRUE), 5)
  )

  # criterion 4 uses Kleiner stage and is disclosed at the head of this script
  direction_rows <- list()
  for (lineage in names(DIRECTION_EXPECTED)) {
    if (!lineage %in% lineages) next
    for (cohort in discovery_cohorts) {
      d <- merge(
        observed[dataset == cohort, c("sample_id", lineage), with = FALSE],
        phenotype[, .(sample_id, fibrosis_stage)], by = "sample_id"
      )
      d <- d[is.finite(fibrosis_stage)]
      rho <- if (nrow(d) > 10L && stats::sd(d[[lineage]]) > 0) {
        suppressWarnings(cor(d[[lineage]], d$fibrosis_stage, method = "spearman"))
      } else NA_real_
      direction_rows[[length(direction_rows) + 1L]] <- data.table(
        method = label, celltype = lineage, cohort = cohort, n = nrow(d),
        spearman_fibrosis = round(rho, 4),
        expected_sign = DIRECTION_EXPECTED[[lineage]],
        agrees = is.finite(rho) && sign(rho) == DIRECTION_EXPECTED[[lineage]]
      )
    }
  }
  direction <- if (length(direction_rows)) rbindlist(direction_rows) else data.table()

  rank_rows <- list()
  for (cohort in discovery_cohorts) {
    d <- observed[dataset == cohort]
    full_rank <- NA
    if (nrow(d) > 5L && all(c(DENOMINATOR_LINEAGE, unname(COVARIATE_LINEAGES)) %in% lineages)) {
      denom <- d[[DENOMINATOR_LINEAGE]]
      denom[is.finite(denom) & denom == 0] <- composition_pseudocount
      L <- vapply(unname(COVARIATE_LINEAGES), function(lineage) {
        num <- d[[lineage]]
        num[is.finite(num) & num == 0] <- composition_pseudocount
        log(num / denom)
      }, numeric(nrow(d)))
      full_rank <- qr(cbind(1, L))$rank == (ncol(L) + 1L)
    }
    rank_rows[[length(rank_rows) + 1L]] <- data.table(
      method = label, cohort = cohort, covariate_block_full_rank = full_rank
    )
  }
  rank_audit <- rbindlist(rank_rows)

  covariate_present <- all(
    c(DENOMINATOR_LINEAGE, unname(COVARIATE_LINEAGES)) %in% lineages
  )
  max_zero <- if (covariate_present) {
    max(lineage_audit[celltype %in% unname(COVARIATE_LINEAGES), zero_fraction])
  } else NA_real_
  hep_median <- lineage_audit[celltype == DENOMINATOR_LINEAGE, median_fraction]
  cohort_coverage <- coverage[dataset %in% REQUIRED_COHORTS]

  c1 <- isTRUE(covariate_present) && is.finite(max_zero) && max_zero <= MAX_ZERO_FRACTION
  c2 <- nrow(cohort_coverage) == length(REQUIRED_COHORTS) &&
    all(cohort_coverage$coverage >= MIN_COHORT_COVERAGE)
  c3 <- length(hep_median) == 1L && is.finite(hep_median) &&
    hep_median >= HEPATOCYTE_MEDIAN_RANGE[[1L]] && hep_median <= HEPATOCYTE_MEDIAN_RANGE[[2L]]
  c4 <- nrow(direction) > 0L &&
    all(direction[, .(ok = sum(agrees) >= MIN_DIRECTION_COHORTS), by = celltype]$ok)
  c5 <- nrow(rank_audit) > 0L && all(isTRUE_vector(rank_audit$covariate_block_full_rank))

  failures <- c(
    if (!c1) "degenerate_covariate_lineages",
    if (!c2) "incomplete_required_cohort_coverage",
    if (!c3) "implausible_hepatocyte_fraction",
    if (!c4) "stage_direction_disagreement",
    if (!c5) "rank_deficient_covariate_block"
  )

  list(
    summary = data.table(
      method = label,
      n_lineages = length(lineages),
      n_matched_samples = sum(present),
      max_covariate_zero_fraction = max_zero,
      median_hepatocyte_fraction = if (length(hep_median) == 1L) hep_median else NA_real_,
      min_required_cohort_coverage = if (nrow(cohort_coverage)) min(cohort_coverage$coverage) else NA_real_,
      direction_pass = c4,
      full_rank_all_discovery_cohorts = c5,
      criterion_1_degeneracy = c1,
      criterion_2_coverage = c2,
      criterion_3_plausibility = c3,
      criterion_4_stage_direction = c4,
      criterion_5_rank = c5,
      eligible = c1 && c2 && c3 && c4 && c5,
      failure = if (length(failures)) paste(failures, collapse = ";") else NA_character_
    ),
    lineage = lineage_audit,
    direction = direction,
    rank = rank_audit,
    coverage = cbind(method = label, coverage),
    matched = matched
  )
}

isTRUE_vector <- function(x) vapply(x, isTRUE, logical(1))

# ------------------------------------------------------------------------ main

acceptance_root <- file.path(workstream_root, "precontract", "composition_acceptance")
if (dir.exists(acceptance_root)) fail("Refusing to overwrite acceptance seal: ", acceptance_root)

message("[1/5] Building corrected census")
census <- census_samples()
assert_true(nrow(census) == 1257L, "Corrected census must contain 1,257 QC-passing samples")

phenotype <- as.data.table(readRDS(input_meta))[, .(sample_id, dataset, fibrosis_stage, nas_score)]

message("[2/5] Reading candidate composition estimators")
methods <- list(
  music_deposited = read_music_deposited(),
  bayesprism = read_per_cohort("_bayesprism_proportions.tsv"),
  rectangle = read_per_cohort("_rectangle_proportions.tsv")
)

message("[3/5] Scoring each estimator against the frozen criteria")
evaluated <- lapply(names(methods), function(label) {
  message("      - ", label)
  evaluate_method(label, methods[[label]], census, phenotype)
})
names(evaluated) <- names(methods)

summary_table <- rbindlist(lapply(evaluated, `[[`, "summary"), use.names = TRUE, fill = TRUE)
lineage_table <- rbindlist(lapply(evaluated, `[[`, "lineage"), use.names = TRUE, fill = TRUE)
direction_table <- rbindlist(lapply(evaluated, `[[`, "direction"), use.names = TRUE, fill = TRUE)
rank_table <- rbindlist(lapply(evaluated, function(x) if (is.null(x$rank)) data.table() else x$rank),
                        use.names = TRUE, fill = TRUE)
coverage_table <- rbindlist(lapply(evaluated, function(x) if (is.null(x$coverage)) data.table() else x$coverage),
                            use.names = TRUE, fill = TRUE)

eligible <- summary_table[eligible == TRUE]
if (!nrow(eligible)) {
  fail("No composition estimator satisfies the frozen acceptance criteria; ",
       "composition adjustment cannot be reported for v9")
}
# Decision rule, fixed: among eligible estimators prefer the widest required-cohort
# coverage, then the lowest covariate degeneracy, then alphabetical for determinism.
setorder(eligible, -min_required_cohort_coverage, max_covariate_zero_fraction, method)
accepted_method <- eligible$method[[1L]]
message("[4/5] Accepted estimator: ", accepted_method)

accepted <- evaluated[[accepted_method]]$matched
accepted_lineages <- lineage_columns(accepted)
setcolorder(accepted, c("sample_id", "dataset", accepted_lineages))

tmp <- atomic_dir(acceptance_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

fwrite(accepted, file.path(tmp, "accepted_composition.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA")
write_tsv(summary_table, file.path(tmp, "method_summary.tsv"))
write_tsv(lineage_table, file.path(tmp, "lineage_audit.tsv"))
write_tsv(direction_table, file.path(tmp, "stage_direction_audit.tsv"))
write_tsv(rank_table, file.path(tmp, "covariate_rank_audit.tsv"))
write_tsv(coverage_table, file.path(tmp, "cohort_coverage.tsv"))

criteria <- data.table(
  criterion = c("degeneracy", "coverage", "plausibility", "stage_direction", "rank"),
  rule = c(
    sprintf("max zero fraction across the four covariate lineages <= %.2f", MAX_ZERO_FRACTION),
    sprintf("per-cohort coverage >= %.2f in every discovery cohort and the holdout", MIN_COHORT_COVERAGE),
    sprintf("median hepatocyte fraction within [%.2f, %.2f]", HEPATOCYTE_MEDIAN_RANGE[[1L]],
            HEPATOCYTE_MEDIAN_RANGE[[2L]]),
    sprintf("expected Kleiner-stage sign in at least %d of four discovery cohorts for hepatocytes, fibroblasts and cholangiocytes",
            MIN_DIRECTION_COHORTS),
    "the four hepatocyte log-ratio covariates are full rank in every discovery cohort"
  ),
  uses_outcome = c(FALSE, FALSE, FALSE, TRUE, FALSE)
)
write_tsv(criteria, file.path(tmp, "acceptance_criteria.tsv"))

inputs <- c(input_dge, input_meta, composition_candidate, composition_testability)
write_tsv(data.table(
  input = c("corrected_dge", "corrected_meta", "music_deposited", "composition_testability"),
  path = inputs,
  sha256 = vapply(inputs, sha256_file, character(1))
), file.path(tmp, "input_manifest.tsv"))

writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))

payload <- list(
  state = "COMPOSITION_ESTIMATOR_ACCEPTED",
  release_id = release_id,
  workstream_id = workstream_id,
  accepted_method = accepted_method,
  accepted_table = "accepted_composition.tsv.gz",
  n_census_samples = nrow(census),
  n_lineages = length(accepted_lineages),
  rejected_methods = summary_table[eligible == FALSE, .(method, failure)],
  criteria_are_blind = FALSE,
  criteria_provenance = paste(
    "Criteria 1, 2 and 4 were evaluated during planning on 2026-08-12 and are not blind.",
    "Thresholds were set far from every observed value. Criterion 4 uses Kleiner stage and",
    "constrains only the composition estimator, never program selection."
  ),
  supersedes = "BULK-PROGRAM-MAP-v8 used the deposited MuSiC estimates without a real-bulk acceptance test"
)
jsonlite::write_json(payload, file.path(tmp, "COMPOSITION_ACCEPTED.json"),
                     pretty = TRUE, auto_unbox = TRUE)

artifacts <- setdiff(list.files(tmp), "artifact_manifest.tsv")
write_tsv(data.table(
  artifact = artifacts,
  sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1)),
  size_bytes = file.info(file.path(tmp, artifacts))$size
), file.path(tmp, "artifact_manifest.tsv"))

message("[5/5] Publishing acceptance seal")
publish_dir(tmp, acceptance_root)
Sys.chmod(list.files(acceptance_root, full.names = TRUE), mode = "0440")
cat("COMPOSITION_ACCEPTANCE_COMPLETE\t", accepted_method, "\t", acceptance_root, "\n", sep = "")
