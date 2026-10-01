#!/usr/bin/env Rscript
# Synthetic check of the coloc.susie() shared-posterior eligibility rule
# (review item 2). Two fake SuSiE fits: the GWAS fit covers a subset of the eQTL
# fit's SNPs, as in 06_susie_coloc.R. eQTL signal A sits on a shared SNP; signal
# B sits on an eQTL-only SNP, so almost none of its posterior is shared.

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: smoke_coloc_overlap_eligibility.R <coloc_signal_eligibility.R>")
source(args[[1L]])
stopifnot(as.character(packageVersion("coloc")) == "5.2.3")
set.seed(42)

eqtl_snps <- paste0("1:", 100000L + seq_len(150L))
gwas_snps <- eqtl_snps[1:100]
fake_susie <- function(snps, peaks) {
  lbf <- matrix(0, nrow = length(peaks), ncol = length(snps),
                dimnames = list(NULL, snps))
  for (l in seq_along(peaks)) lbf[l, peaks[[l]]] <- 12
  structure(list(
    lbf_variable = lbf,
    sets = list(cs = stats::setNames(as.list(peaks), paste0("L", seq_along(peaks))),
                cs_index = seq_along(peaks))
  ), class = "susie")
}
s_gwas <- fake_susie(gwas_snps, list(20L))
s_eqtl <- fake_susie(eqtl_snps, list(20L, 130L))   # A shared, B eQTL-only

# 1. Full fits: coloc 5.2.3 keeps A and drops B.
full <- suppressWarnings(coloc.susie(s_gwas, s_eqtl))
stopifnot(nrow(full$summary) == 1L, full$summary$idx2 == 1L,
          is.finite(full$summary$PP.H4.abf))

# 2. Our replicated shares make the same call.
tab <- coloc_susie_pair_table(s_gwas, s_eqtl, full)
stopifnot(
  tab$matches_coloc,
  identical(tab$pairs$eligible, c(TRUE, FALSE)),
  all(abs(tab$pairs$prop_gwas - 1) < 1e-12),
  tab$pairs$prop_eqtl[1] > 0.5,
  tab$pairs$prop_eqtl[2] < 1e-3,
  is.na(tab$pairs$PP.H4[2]),
  coloc_susie_method(tab$kept) == "susie"
)

# 3. The former call (both fits cut to the shared SNPs) passes B by construction,
#    and gives the kept pair the same PP.H4 as the full call.
cut_fit <- function(s, snps) { s$lbf_variable <- s$lbf_variable[, snps, drop = FALSE]; s }
shared <- intersect(gwas_snps, eqtl_snps)
old <- coloc.susie(cut_fit(s_gwas, shared), cut_fit(s_eqtl, shared))
stopifnot(nrow(old$summary) == 2L, all(is.finite(old$summary$PP.H4.abf)),
          isTRUE(all.equal(old$summary[idx2 == 1L]$PP.H4.abf, tab$pairs$PP.H4[1],
                           tolerance = 1e-12)))

# 4. Every pair dropped: coloc returns NA PP.H4 rows; the gene is untestable,
#    not a zero-row best.
s_eqtl_b <- fake_susie(eqtl_snps, list(130L))
none <- suppressWarnings(coloc.susie(s_gwas, s_eqtl_b))
stopifnot(nrow(none$summary) == 1L, all(is.na(none$summary$PP.H4.abf)))
tab_none <- coloc_susie_pair_table(s_gwas, s_eqtl_b, none)
stopifnot(
  tab_none$matches_coloc,
  identical(tab_none$pairs$eligible, FALSE),
  nrow(tab_none$kept) == 0L,
  coloc_susie_method(tab_none$kept) == "susie_untestable_insufficient_shared_posterior"
)

# 5. coloc's last-column-as-null prior: a flat eQTL signal whose LAST SNP is
#    shared keeps the null prior mass on a shared SNP. Our replica must follow
#    coloc here too, even though the uniform-prior share is below 0.5.
flat_snps <- c(paste0("1:", 200000L + seq_len(149L)), gwas_snps[100])
s_gwas_f <- fake_susie(gwas_snps, list(100L))
s_eqtl_f <- s_eqtl_b
s_eqtl_f$lbf_variable <- matrix(0, 1, 150, dimnames = list(NULL, flat_snps))
s_eqtl_f$lbf_variable[1, 150] <- 0.5
flat <- suppressWarnings(coloc.susie(s_gwas_f, s_eqtl_f))
tab_flat <- coloc_susie_pair_table(s_gwas_f, s_eqtl_f, flat)
stopifnot(tab_flat$matches_coloc,
          tab_flat$pairs$prop_eqtl > 0.5,
          tab_flat$pairs$prop_eqtl_uniform < 0.5)

cat("COLOC_OVERLAP_ELIGIBILITY_FIXTURE\tPASS\n")
