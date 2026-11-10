#!/usr/bin/env Rscript
# Run with rnaseq Rscript inside an allocated compute node.
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Tests require a Slurm compute allocation")
arguments <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
source(file.path(dirname(normalizePath(script)), "..", "scripts", "compare_da.R"))
set.seed(42)
donors <- data.table(donor_id = paste0("D", 1:12), condition = rep(c("NORMAL", "MASH"), each = 6),
  n_cells = 100L, contrast_eligible = TRUE)
counts <- matrix(rnbinom(600L * 12L, mu = 25, size = 10), nrow = 12L)
rownames(counts) <- donors$donor_id
colnames(counts) <- paste0("chr1:", seq_len(ncol(counts)) * 1000, "-", seq_len(ncol(counts)) * 1000 + 500)
counts[donors$condition == "MASH", 1:20] <- counts[donors$condition == "MASH", 1:20] * 4L
counts <- as(counts, "dgCMatrix")

# Reordering biological donors cannot change the estimates or hypothesis family.
for (method in c("voom_qw", "edger_ql")) {
  fit <- fit_da(counts, donors, method)$results
  permutation <- c(12:7, 6:1)
  shuffled <- fit_da(counts[permutation, ], donors[permutation], method)$results
  stopifnot(max(abs(fit$logFC - shuffled$logFC)) < 1e-7,
    max(abs(fit$pvalue - shuffled$pvalue)) < 1e-7,
    all(fit$pvalue >= 0 & fit$pvalue <= 1),
    identical(fit$qvalue, p.adjust(fit$pvalue, "BH")), median(fit$logFC[1:20]) > 1)
}

# Reproduce the saved v3 algorithm, particularly the posterior moderated errors.
y <- calcNormFactors(DGEList(counts = t(counts)), method = "TMM")
design <- model.matrix(~ factor(donors$condition, levels = c("NORMAL", "MASH")))
v <- voomWithQualityWeights(y, design, plot = FALSE)
reference_fit <- eBayes(lmFit(v, design), robust = TRUE)
reference <- topTable(reference_fit, coef = 2L, number = Inf, sort.by = "none")
actual <- fit_da(counts, donors, "voom_qw")$results
stopifnot(max(abs(actual$logFC - reference$logFC)) < 1e-10,
  max(abs(actual$standard_error - reference_fit$stdev.unscaled[, 2L] * sqrt(reference_fit$s2.post))) < 1e-10)

expected <- Matrix::colSums(Diagonal(x = 1e6 / Matrix::rowSums(counts)) %*% counts >= 1) >= 4L & Matrix::colSums(counts) >= 20
stopifnot(identical(as.logical(filter_da_peaks(counts, "prevalence")), as.logical(expected)))
stopifnot(inherits(try(make_da_design(copy(donors)[, batch_mm := condition], "batch_mm"), silent = TRUE), "try-error"))

# A completely confounded batch cannot supply an exchangeable permutation null.
blocked <- copy(donors)[, batch_mm := condition]
directional <- directional_da(counts, blocked, "voom_qw", setNames(1, colnames(counts)[1]), 1L, 2L, 42L)
stopifnot(is.na(directional$p), directional$status == "untestable_no_within_batch_label_exchange")
stopifnot(all(is.na(fit_da(counts, donors, "edger_ql")$results$standard_error)))
stopifnot(covariate_availability(donors, "sex")$status == "missing_column",
  covariate_availability(copy(donors)[, sex := NA_character_], "sex")$status == "all_missing",
  covariate_availability(blocked, "batch_mm")$status == "nonestimable")

# Null methods see the same simulated counts and diagnostic feature universe.
null_voom <- null_da(counts, donors, "voom_qw", 2L, 100L, 41L)
null_ql <- null_da(counts, donors, "edger_ql", 2L, 100L, 41L)
stopifnot(identical(null_voom$simulated_count_sum, null_ql$simulated_count_sum),
  identical(null_voom$simulated_count_sum_squares, null_ql$simulated_count_sum_squares))
null_voom_parallel <- null_da(counts, donors, "voom_qw", 2L, 100L, 41L, workers = 2L)
stopifnot(identical(null_voom, null_voom_parallel))
discoveries <- setNames(c(1, -1), colnames(counts)[1:2])
serial_directional <- directional_da(counts, donors, "voom_qw", discoveries, 1L, 2L, 41L)
parallel_directional <- directional_da(counts, donors, "voom_qw", discoveries, 1L, 2L, 41L, workers = 2L)
stopifnot(identical(serial_directional, parallel_directional))
cat("PASS: donor order, effect sign, BH, moderated SE, filtering, rank, exchangeability, covariate availability, paired null draws, serial/parallel identity\n")
