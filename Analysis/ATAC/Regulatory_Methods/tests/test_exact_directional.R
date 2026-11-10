#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Requires a Slurm compute allocation")
arguments <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
source(file.path(dirname(normalizePath(script)), "..", "scripts", "compare_da.R"))
donors <- data.table(donor_id = paste0("D", 1:4), condition = rep(c("NORMAL", "MASH"), 2),
  batch_mm = rep(c("A", "B"), each = 2L), n_cells = 100L, contrast_eligible = TRUE)
design <- restricted_label_design(donors)
choices <- restricted_label_choices(design)
labels <- lapply(seq_len(4L), function(index) restricted_labels(design, choices, index, 4L))
signatures <- vapply(labels, paste, character(1), collapse = ",")
expected <- c("NORMAL,MASH,NORMAL,MASH", "MASH,NORMAL,NORMAL,MASH", "NORMAL,MASH,MASH,NORMAL", "MASH,NORMAL,MASH,NORMAL")
stopifnot(design$assignment_count == 4L, setequal(signatures, expected), !anyDuplicated(signatures),
  sum(signatures == paste(donors$condition, collapse = ",")) == 1L)
set.seed(42)
counts <- matrix(rnbinom(80L * 4L, mu = 20, size = 8), nrow = 4L)
rownames(counts) <- donors$donor_id
colnames(counts) <- paste0("chr1:", 1:80 * 1000L, "-", 1:80 * 1000L + 500L)
discoveries <- setNames(c(1, -1), colnames(counts)[1:2])
actual <- fit_da(counts, donors, "voom_qw")$results
observed <- sum(sign(actual$logFC[match(names(discoveries), actual$peak_coordinate)]) == sign(discoveries))
serial <- directional_da(counts, donors, "voom_qw", discoveries, observed, reps = 1000L, seed = 42L, workers = 1L)
forked <- directional_da(counts, donors, "voom_qw", discoveries, observed, reps = 1000L, seed = 42L, workers = 2L)
independent <- vapply(strsplit(expected, ",", fixed = TRUE), function(labels) {
  shuffled <- copy(donors)
  shuffled$condition <- labels
  value <- fit_da(counts, shuffled, "voom_qw")$results
  sum(sign(value$logFC[match(names(discoveries), value$peak_coordinate)]) == sign(discoveries)) >= observed
}, logical(1))
stopifnot(identical(serial, forked), serial$mode == "exact", serial$n == 4L,
  serial$assignment_count == 4L, serial$p == mean(independent), serial$p >= .25)
monte_carlo <- directional_da(counts, donors, "voom_qw", discoveries, observed, reps = 2L, seed = 42L)
stopifnot(monte_carlo$mode == "monte_carlo", monte_carlo$n == 2L,
  abs(monte_carlo$p * 3 - round(monte_carlo$p * 3)) < 1e-12)
incomplete <- copy(donors)
incomplete$batch_mm[1L] <- NA_character_
untestable <- directional_da(counts, incomplete, "voom_qw", discoveries, observed, reps = 4L, seed = 42L)
stopifnot(is.na(untestable$p), untestable$status == "untestable_incomplete_recorded_batch_metadata")
cat("PASS: exact unique restricted assignments, identity once, exact denominator, Monte Carlo plus-one, workers1/2 equality, incomplete-batch refusal\n")
