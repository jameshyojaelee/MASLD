#!/usr/bin/env Rscript
# Does the canonical COLOC pipeline's fixed 1e-3 ridge suffice on the shipped
# non-EUR panels at the scale COLOC actually uses?
#
# get_ld_per_locus() reads the SAME block .ld files and takes a principal
# submatrix matched to the cis-region summary stats, then 06_susie_coloc.R adds
# LD_REGULARIZE = 1e-3. A principal submatrix of an indefinite matrix can itself
# be indefinite, so subset size is the thing to test.
suppressPackageStartupMessages({ library(data.table) })
P <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/data/ld_ref"
BLK <- "chr8/9154694.9640787/9154694.9640787"
LD_REGULARIZE <- 1e-3

load_block <- function(root) {
  M <- as.matrix(fread(file.path(P, root, paste0(BLK, ".ld")), header = FALSE))
  M <- (M + t(M)) / 2
  ok <- rowSums(!is.finite(M)) == 0 & colSums(!is.finite(M)) == 0
  M <- M[ok, ok, drop = FALSE]; diag(M) <- 1; M
}
Rs <- load_block("1kg_afr")
Rg <- load_block("1kg_afr_gram")
cat(sprintf("shipped %d x %d | rebuilt %d x %d\n", nrow(Rs), ncol(Rs), nrow(Rg), ncol(Rg)))

pd <- function(M) tryCatch({ invisible(chol(M + LD_REGULARIZE * diag(nrow(M)))); TRUE },
                           error = function(e) FALSE)
set.seed(11)
sizes <- c(200L, 500L, 1000L, 2000L, 3000L)
cat("\ncis-window-sized principal submatrices, 40 random draws each\n")
cat(sprintf("%-8s %-26s %-26s\n", "size", "SHIPPED chol@1e-3 ok", "REBUILT chol@1e-3 ok"))
for (s in sizes) {
  if (s > nrow(Rs)) next
  okS <- okG <- 0L
  for (i in 1:40) {
    idx <- sort(sample(nrow(Rs), s))
    okS <- okS + pd(Rs[idx, idx, drop = FALSE])
    okG <- okG + pd(Rg[idx, idx, drop = FALSE])
  }
  cat(sprintf("%-8d %-26s %-26s\n", s, sprintf("%d/40 (%.0f%%)", okS, 100*okS/40),
              sprintf("%d/40 (%.0f%%)", okG, 100*okG/40)))
}
cat("\nfull block:\n")
cat(sprintf("  SHIPPED chol@1e-3: %s | REBUILT chol@1e-3: %s\n", pd(Rs), pd(Rg)))
