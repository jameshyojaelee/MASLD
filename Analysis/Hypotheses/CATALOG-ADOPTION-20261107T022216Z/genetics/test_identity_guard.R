#!/usr/bin/env Rscript
library(data.table)
here <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[startsWith(commandArgs(FALSE), "--file=")])))
root <- normalizePath(file.path(here, "../../../.."))
source(file.path(root, "GWAS/finemapping/src/coloc_variant_identity.R"))
stopifnot(identical(coloc_variant_identity_collisions(c("17:1", "17:2", "17:1")), c(1L, 3L)),
          length(coloc_variant_identity_collisions(c("17:1", "17:2"))) == 0L)
fixture <- function(duplicate_gwas = FALSE, duplicate_eqtl = FALSE) {
  env <- new.env(parent = globalenv())
  env$chr_num <- 17L
  env$ss_ld <- data.frame(position = if (duplicate_gwas) c(17164840L, 17164840L) else c(17164840L, 17164841L),
    allele1 = c("CTAAA", "C"), allele2 = c("C", "CTAAA"), beta = c(-0.006156, 0.006156),
    standard_error = c(0.002588, 0.002588))
  env$eqtl_snp_ids <- if (duplicate_eqtl) c("17:3", "17:3") else c("17:3", "17:4")
  env$susie_method <- "abf_only"
  env$susie_pp4 <- NA_real_
  env$abf_pp4 <- .8
  env$gene_id <- "ENSG00000133030"
  env$OUT_DIR <- tempfile("identity-fixture-")
  dir.create(env$OUT_DIR)
  env$merged <- data.table(eqtl_pos = 17164840L, gwas_af = .2, eqtl_eaf = .8)
  env$fit_called <- FALSE
  lines <- c('          snp_ids <- paste0(chr_num, ":", ss_ld$position)',
             'fit_called <- TRUE')
  guarded <- coloc_identity_guarded_runner(lines)
  error <- tryCatch({ eval(parse(text = guarded), envir = env); NULL }, error = identity)
  if (duplicate_gwas || duplicate_eqtl) {
    stopifnot(inherits(error, "error"), !env$fit_called, env$abf_pp4 == .8,
      is.na(env$susie_pp4), env$susie_method == COLOC_AMBIGUOUS_VARIANT_METHOD)
    suffix <- if (duplicate_gwas) "_gwas.tsv" else "_eqtl.tsv"
    files <- list.files(env$OUT_DIR, full.names = TRUE)
    stopifnot(any(endsWith(files, suffix)))
    collision <- fread(files[endsWith(files, suffix)])
    stopifnot(nrow(collision) == 2L)
  } else {
    stopifnot(is.null(error), env$fit_called, env$susie_method == "abf_only", env$abf_pp4 == .8)
  }
}
fixture(TRUE, FALSE)
fixture(FALSE, TRUE)
fixture(FALSE, FALSE)
cat("PASS: all collisions retained, no ambiguous fit attempted, ABF preserved, both assay sides guarded, unique-ID path unchanged\n")
