#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
stopifnot(length(script_arg) == 1L)
test_dir <- dirname(normalizePath(script_arg))
analysis_dir <- normalizePath(file.path(test_dir, ".."))
ml_dir <- file.path(analysis_dir, "molecular_layers")
Sys.setenv(
  HAC_ML_LIBRARY_ONLY = "1",
  HAC_ML_SCRIPT_DIR = ml_dir,
  MASLD_PROJECT_ROOT = normalizePath(file.path(analysis_dir, "../../.."))
)
source(file.path(ml_dir, "20_pathway_continuum.R"))

test_case <- function(name, expression) {
  force(expression)
  cat(sprintf("PASS\t%s\n", name))
}

test_case("Frozen GMT families retain exact declared cardinality", {
  contract <- ml_read_contract()
  gmt_dir <- ml_resolve(contract$pathway_dir)
  for (collection in names(contract$pathway_collections)) {
    expected <- as.integer(contract$pathway_collections[[collection]])
    observed <- pathway_read_gmt(
      file.path(gmt_dir, paste0(collection, ".gmt")), collection, expected
    )
    stopifnot(length(observed$ids) == expected)
    stopifnot(uniqueN(observed$membership$set_id) == expected)
  }
})

test_case("Signature removal and 80-percent pathway testability use original membership", {
  expression <- matrix(
    seq_len(44), nrow = 11L,
    dimnames = list(LETTERS[1:11], paste0("S", 1:4))
  )
  metadata <- data.table(sample_id = colnames(expression), dataset = "C1")
  membership <- data.table(
    set_id = rep(c("PASS", "FAIL"), each = 10L),
    gene_symbol = c(LETTERS[1:10], LETTERS[2:11])
  )
  scored <- pathway_score_collection(
    expression, metadata, c("PASS", "FAIL"), membership,
    signature_symbols = c("A", "B", "C"), minimum_genes = 7L,
    minimum_fraction = 0.8
  )
  audit <- scored$coverage
  stopifnot(audit[set_id == "PASS", observed_retained_fraction] == 0.7)
  stopifnot(audit[set_id == "FAIL", observed_retained_fraction] == 0.8)
  stopifnot(!audit[set_id == "PASS", testable])
  stopifnot(audit[set_id == "FAIL", testable])
  stopifnot(!any(c("A", "B", "C") %in%
                  membership[!gene_symbol %in% c("A", "B", "C"), gene_symbol]))
})

test_case("Vectorized donor model recovers a positive continuum coefficient", {
  n <- 40L
  axis <- seq(-2, 2, length.out = n)
  stage <- rep(0:3, each = 10L)
  sex <- rep(c("F", "M"), length.out = n)
  outcome <- rbind(
    POSITIVE = 2 * axis + rep(c(-0.1, 0.1), length.out = n),
    NEGATIVE = -axis + rep(c(-0.1, 0.1), length.out = n)
  )
  fit <- pathway_fit_matrix(outcome, axis, stage, sex)
  stopifnot(fit[set_id == "POSITIVE", beta] > 0)
  stopifnot(fit[set_id == "NEGATIVE", beta] < 0)
  stopifnot(all(fit$spline_df == 3L))
})

test_case("Complete-family BH does not shrink to tested subsets", {
  p <- c(0.001, 0.02, NA_real_)
  observed <- ml_complete_bh(p, 50L)
  expected <- c(p.adjust(p[1:2], method = "BH", n = 50L), NA_real_)
  stopifnot(isTRUE(all.equal(observed, expected, tolerance = 0)))
})

cat("ALL_PATHWAY_CONTRACTS_PASS\n")
