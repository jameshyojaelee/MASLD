#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
stopifnot(length(script_arg) == 1L)
test_dir <- dirname(normalizePath(script_arg))
analysis_dir <- normalizePath(file.path(test_dir, ".."))
ml_dir <- file.path(analysis_dir, "molecular_layers")
project_root <- normalizePath(file.path(analysis_dir, "../../.."))
Sys.setenv(
  HAC_ML_LIBRARY_ONLY = "1",
  HAC_ML_SCRIPT_DIR = ml_dir,
  MASLD_PROJECT_ROOT = project_root
)
source(file.path(ml_dir, "21_tf_continuum.R"))

test_case <- function(name, expression) {
  force(expression)
  cat(sprintf("PASS\t%s\n", name))
}

test_case("Frozen DoRothEA family has 268 TFs and fixed display axes", {
  contract <- ml_read_contract()
  roster_path <- file.path(
    project_root,
    "figures/candidates/pi-figure-redesign-2026-08-17-v13/source_tables/fig4f_all_tf_activity.tsv"
  )
  roster <- tf_fixed_roster(roster_path, as.integer(contract$tf_family_size))
  stopifnot(length(roster) == 268L)
  stopifnot(all(contract$display_tfs %in% roster))
})

test_case("Signed regulon scoring excludes signature targets before L1 coverage", {
  genes <- paste0("G", 1:12)
  set.seed(20260817)
  expression <- matrix(
    rnorm(48), nrow = 12L,
    dimnames = list(genes, paste0("S", 1:4))
  )
  metadata <- data.table(sample_id = colnames(expression), dataset = "C1")
  network <- data.table(
    tf = "TF1", target = genes, mor = rep(c(-1, 1), 6L), confidence = "A"
  )
  scored <- tf_score_regulons(
    expression, metadata, network, roster = "TF1",
    signature = c("G1", "G2"), minimum_fraction = 0.8,
    minimum_targets = 5L
  )
  audit <- scored$coverage
  stopifnot(audit$n_signature_targets == 2L)
  stopifnot(isTRUE(all.equal(audit$retained_l1_fraction, 10 / 12, tolerance = 0)))
  stopifnot(audit$testable)
  stopifnot(identical(rownames(scored$scores$C1), "TF1"))
})

test_case("TF BH correction retains the frozen 268-member family", {
  p <- c(0.001, 0.02, NA_real_)
  observed <- ml_complete_bh(p, 268L)
  expected <- c(p.adjust(p[1:2], method = "BH", n = 268L), NA_real_)
  stopifnot(isTRUE(all.equal(observed, expected, tolerance = 0)))
})

test_case("TF source fixes run_wmean parameters and forbids downloads", {
  source_text <- paste(readLines(file.path(ml_dir, "21_tf_continuum.R")), collapse = "\n")
  stopifnot(grepl("times = permutations, minsize = 5L", source_text, fixed = TRUE))
  stopifnot(grepl("1000L, as.integer(contract$tf_seed)", source_text, fixed = TRUE))
  stopifnot(!grepl("download.file", source_text, fixed = TRUE))
  stopifnot(!grepl("install.packages", source_text, fixed = TRUE))
})

cat("ALL_TF_CONTRACTS_PASS\n")
