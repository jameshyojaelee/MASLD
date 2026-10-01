#!/usr/bin/env Rscript
# B-COLOC helper: does the marginal sign at a colocalized eQTL lead match the signal's conditional sign?
#
# Input  (arg 1): TSV with columns ensembl, chr, idx2, hit2. idx2 is the SuSiE component index l as
#                 coloc.susie reports it (cs_index value); hit2 is "chr:pos" (hg19).
# Fits   (arg 2): <dir>/chr<chr>/<ensembl>_susie.rds, the eQTL fits COLOC read
#                 (GWAS/finemapping/results/eqtl_susie_polyfun).
# Output (arg 3): TSV, one row per input row: mu[l, hit2] (posterior mean effect of component l at the
#                 lead) and XtXr[hit2] (the fitted marginal, R times the posterior mean of all
#                 components). Both are in the fit's own allele coding, so whether their signs agree does
#                 not depend on allele orientation. Also alpha, PIP, and whether hit2 is the component's
#                 lbf argmax over all fitted variants (coloc takes the argmax over shared variants only).
args <- commandArgs(trailingOnly = TRUE)
inp <- read.delim(args[1], colClasses = "character")
fit_dir <- args[2]
rows <- vector("list", nrow(inp))
for (i in seq_len(nrow(inp))) {
  r <- inp[i, ]
  f <- file.path(fit_dir, paste0("chr", r$chr), paste0(r$ensembl, "_susie.rds"))
  res <- data.frame(ensembl = r$ensembl, idx2 = r$idx2, hit2 = r$hit2, fit_found = file.exists(f),
                    hit2_in_fit = NA, lbf_argmax_is_hit2 = NA, alpha_hit2 = NA_real_, pip_hit2 = NA_real_,
                    mu_hit2 = NA_real_, xtxr_hit2 = NA_real_, conditional_sign_equals_fitted_marginal = NA)
  if (res$fit_found) {
    s <- readRDS(f)
    l <- as.integer(r$idx2)
    j <- match(r$hit2, colnames(s$lbf_variable))
    res$hit2_in_fit <- !is.na(j)
    if (!is.na(j) && l <= nrow(s$mu)) {
      res$lbf_argmax_is_hit2 <- unname(which.max(s$lbf_variable[l, ]) == j)
      res$alpha_hit2 <- s$alpha[l, j]
      res$pip_hit2 <- s$pip[j]
      res$mu_hit2 <- s$mu[l, j]
      res$xtxr_hit2 <- s$XtXr[j]
      res$conditional_sign_equals_fitted_marginal <- sign(s$mu[l, j]) == sign(s$XtXr[j])
    }
  }
  rows[[i]] <- res
}
write.table(do.call(rbind, rows), args[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat("eQTL conditional-sign rows:", nrow(inp), "\n")
