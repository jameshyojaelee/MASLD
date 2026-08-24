#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
script_dir <- file.path(root, "scripts/analysis/histology_anchored_continuum/molecular_layers")
Sys.setenv(HAC_ML_LIBRARY_ONLY = "1", HAC_ML_SCRIPT_DIR = script_dir)
source(file.path(script_dir, "20_pathway_continuum.R"))

contract <- ml_read_contract()
for (collection in names(contract$pathway_collections)) {
  gmt <- pathway_read_gmt(
    ml_resolve(file.path(contract$pathway_dir, paste0(collection, ".gmt"))),
    collection, as.integer(contract$pathway_collections[[collection]])
  )
  stopifnot(length(gmt$ids) == as.integer(contract$pathway_collections[[collection]]))
}

set.seed(20260817)
n <- 120L
axis <- rnorm(n)
outcomes <- rbind(
  positive = 0.8 * axis + rnorm(n, sd = 0.5),
  null = rnorm(n)
)
sex <- sample(c("F", "M"), n, replace = TRUE)
fit <- pathway_fit_matrix(
  outcomes, axis, rep(0:3, length.out = n), sex
)
stopifnot(nrow(fit) == 2L, fit[set_id == "positive", beta] > 0)

adjusted <- ml_complete_bh(c(0.001, rep(NA_real_, 49L)), 50L)
stopifnot(abs(adjusted[[1L]] - 0.05) < 1e-12)

source_text <- readLines(file.path(script_dir, "20_pathway_continuum.R"), warn = FALSE)
stopifnot(
  any(grepl("signature_symbols", source_text, fixed = TRUE)),
  any(grepl("voomWithQualityWeights", source_text, fixed = TRUE)),
  !any(grepl("install.packages|BiocManager::install|download.file", source_text))
)
cat("test_pathway_contracts: PASS\n")
