#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
test_dir <- dirname(normalizePath(script_arg))
source(file.path(dirname(test_dir), "lib_molecular_layers.R"))

contract <- ml_read_contract()
stopifnot(
  identical(contract$gene_family_size, 23370L),
  identical(contract$signature_published_size, 145L),
  identical(contract$signature_observed_size, 139L),
  identical(contract$program_family_size, 117L),
  identical(contract$nmf_family_size, 10L),
  identical(contract$tf_family_size, 268L),
  length(contract$window_centers) == 9L,
  identical(contract$window_width, 0.2)
)

p <- c(0.001, 0.02, NA_real_)
adjusted <- ml_complete_bh(p, 117L)
stopifnot(length(adjusted) == 3L, is.na(adjusted[[3L]]), adjusted[[1L]] >= p[[1L]])

meta <- ml_fixed_meta(c(0.2, 0.4), c(0.1, 0.2))
stopifnot(meta$estimable, meta$n_cohorts == 2L, meta$beta > 0.2, meta$beta < 0.4)

x <- c(-1, 0, 1, NA_real_)
stopifnot(all.equal(ml_standardize(x)[1:3], as.numeric(scale(x[1:3]))))
stopifnot(all.equal(ml_percentile(1:3), c(0, 0.5, 1)))

windows <- ml_fixed_windows(seq(0, 1, by = 0.1), contract$window_centers,
                            contract$window_width)
stopifnot(uniqueN(windows$window_id) == 9L, max(windows$upper) == 1)

scripts <- list.files(dirname(test_dir), pattern = "[.](R|py|sh|sbatch)$", full.names = TRUE)
text <- paste(vapply(scripts, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
                     character(1)), collapse = "\n")
stopifnot(
  !grepl("longitudinal progression", text, fixed = TRUE),
  !grepl("patient trajectory", text, fixed = TRUE),
  !grepl("A3_axis.tsv", text, fixed = TRUE)
)

cat("MOLECULAR_LAYER_CONTRACT_TESTS_PASS\n")

extra_tests <- list.files(test_dir, pattern = "^test_.*[.]R$", full.names = TRUE)
for (test_path in sort(extra_tests)) {
  sys.source(test_path, envir = new.env(parent = globalenv()))
}
cat("MOLECULAR_LAYER_ALL_UNIT_TESTS_PASS\n")
