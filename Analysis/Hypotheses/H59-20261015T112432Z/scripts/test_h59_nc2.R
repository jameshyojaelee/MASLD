#!/usr/bin/env Rscript
# Exercise the actual NC2 expression without executing the H59 analysis.
args <- commandArgs(trailingOnly = TRUE)
source_path <- if (length(args)) args[1L] else
  "Analysis/Hypotheses/H59-20261015T112432Z/scripts/01_h59_gates_null.R"

nc2_gls_expression <- function(path) {
  expressions <- parse(path)
  assignments <- Filter(function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
                          identical(e[[2L]], as.name("nc2")), as.list(expressions))
  stopifnot(length(assignments) == 1L)
  expression <- as.list(assignments[[1L]][[3L]])[["gls_rate"]]
  stopifnot(!is.null(expression))
  expression
}

gls_expression <- nc2_gls_expression(source_path)
source_rate <- function(Tn, nfit) eval(gls_expression, list(Tn = Tn, nfit = nfit))
reference_rate <- function(Tn, nfit) {
  vapply(seq_len(ncol(Tn)), function(j)
    mean(2 * pt(abs(Tn[, j]), df = nfit[j] - 3, lower.tail = FALSE) < 0.05,
         na.rm = TRUE), numeric(1))
}
check <- function(Tn, nfit) {
  observed <- source_rate(Tn, nfit)
  expected <- reference_rate(Tn, nfit)
  stopifnot(isTRUE(all.equal(unname(observed), expected, tolerance = 0)))
}

# Unequal degrees of freedom make accidental matrix recycling observable.
different_df <- matrix(c(2.0, 2.1, 2.2, 2.3, 2.4, 2.5,
                         2.5, 2.4, 2.3, 2.2, 2.1, 2.0,
                         1.9, 2.0, 2.1, 2.2, 2.3, 2.4), nrow = 6L)
nfit_fixture <- c(5L, 15L, 55L)
check(different_df, nfit_fixture)
buggy_rate <- colMeans(2 * pt(abs(different_df), nfit_fixture - 3,
                            lower.tail = FALSE) < 0.05)
stopifnot(any(buggy_rate != reference_rate(different_df, nfit_fixture)))

missing <- different_df
missing[cbind(c(1L, 4L, 6L), c(1L, 2L, 3L))] <- c(NA_real_, NaN, NA_real_)
check(missing, nfit_fixture)
check(matrix(c(NA_real_, 2.0, 2.3, NaN), ncol = 1L), 15L)
check(matrix(NA_real_, nrow = 4L, ncol = 1L), 15L)
cat("PASS: actual NC2 source expression; unequal df, missing values, one column, all missing.\n")

# Optional replay reads existing draws but writes only to a new scratch folder.
if (length(args) >= 4L) {
  suppressPackageStartupMessages(library(data.table))
  setDTthreads(1L)
  old_source <- args[2L]
  old_out <- args[3L]
  scratch <- args[4L]
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  Tn <- as.matrix(fread(file.path(old_out, "null_gls_t.tsv.gz")))
  design <- fread(file.path(old_out, "null_design_df.tsv"))
  saved <- fread(file.path(old_out, "NC2_parametric_calibration.tsv"))
  saved_summary <- fread(file.path(old_out, "NC2_summary.tsv"))
  stopifnot(identical(colnames(Tn), design$program), identical(saved$program, design$program))
  nfit <- design$n_fit_gls
  old_rate <- eval(nc2_gls_expression(old_source), list(Tn = Tn, nfit = nfit))
  new_rate <- source_rate(Tn, nfit)
  stopifnot(max(abs(old_rate - saved$gls_rate)) < 1e-14)
  check(Tn, nfit)
  comparison <- data.table(program = design$program, n_fit = nfit,
                           old_gls_rate = old_rate, new_gls_rate = new_rate,
                           difference = new_rate - old_rate)
  corrected <- copy(saved)
  corrected[, gls_rate := new_rate]
  corrected_summary <- copy(saved_summary)
  corrected_summary[, `:=`(
    gls_median = median(new_rate),
    n_gls_above_0.075 = sum(new_rate > 0.075),
    gls_parametric_label_anticonservative = sum(new_rate > 0.075) > 10,
    ductular_gls = new_rate[design$program == "hepatocytes_20"])]
  summary_comparison <- data.table(
    field = names(saved_summary),
    old = vapply(saved_summary, as.character, character(1)),
    new = vapply(corrected_summary, as.character, character(1)))
  summary_comparison[, changed := old != new]
  fwrite(comparison, file.path(scratch, "nc2_program_comparison.tsv"), sep = "\t")
  fwrite(summary_comparison, file.path(scratch, "nc2_summary_comparison.tsv"), sep = "\t")
  fwrite(corrected, file.path(scratch, "NC2_parametric_calibration.corrected.tsv"), sep = "\t")
  fwrite(corrected_summary, file.path(scratch, "NC2_summary.corrected.tsv"), sep = "\t")
  print(summary_comparison)
  cat("PASS: saved rates reproduced for", ncol(Tn), "programs; corrected rates match per-column reference.\n")
  writeLines(capture.output(sessionInfo()), file.path(scratch, "sessionInfo.txt"))
}
