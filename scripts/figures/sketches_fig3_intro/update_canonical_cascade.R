#!/usr/bin/env Rscript
# Update method_cascade.csv with the canonical numbers (PolyFun for UKBB EUR,
# 1KG for the rest) and the cross-ancestry physical-locus counts agreed on
# 2026-05-04.
suppressPackageStartupMessages({ library(data.table) })

OUT <- "scripts/figures/sketches_fig3_intro/data"

cascade <- data.table(
  step = c(
    "Lead loci screened",
    "SuSiE converged loci (canonical)",
    "SuSiE-X high-PIP physical loci (cross-ancestry)",
    "meSuSiE shared CS physical loci (cross-ancestry)",
    "Genes — regular (ABF) COLOC PP4 ≥ 0.5",
    "Genes — SuSiE COLOC PP4 ≥ 0.5"),
  count = c(308L, 272L, 91L, 110L, 641L, 394L))
fwrite(cascade, file.path(OUT, "method_cascade.csv"))

cat("Updated method_cascade.csv:\n")
print(cascade)

# Update method_totals (SuSiE-X / meSuSiE physical-locus level)
method_totals <- data.table(
  method = c("SuSiE (single ancestry)",
             "SuSiE-X (cross-ancestry, EUR+EAS)",
             "meSuSiE (multi-ethnic, EUR+EAS)"),
  n_loci_converged    = c(272L, 132L, 122L),
  n_loci_high_pip_or_shared = c(NA_integer_, 91L, 110L),
  n_genes_with_pip05  = c(NA_integer_, 1407L, 1557L),
  n_genes_with_pip09  = c(NA_integer_, 1246L, 1168L))
fwrite(method_totals, file.path(OUT, "method_totals.csv"))

# Update ancestry_summary with canonical counts
ancestry_summary <- data.table(
  ancestry = c("EUR","EAS","AFR","SAS"),
  n_gwas = c(17, 5, 3, 3),
  n_samples_total = c(4872909, 1232652, 19908, 26628),
  n_loci_converged = c(135L, 98L, 30L, 9L),
  n_coloc_genes_pp4_05 = c(364L, 92L, 1L, 5L))
fwrite(ancestry_summary, file.path(OUT, "ancestry_summary.csv"))

cat("\nAll canonical data tables updated.\n")
