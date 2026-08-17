#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: smoke_coloc_susie_export.R <helper.R>")
source(args[[1L]])
set.seed(42)
snps <- paste0("1:", seq_len(100L) + 100000L)
bf1 <- rbind(
  stats::setNames(c(8, rep(0, 99)), snps),
  stats::setNames(c(rep(0, 49), 8, rep(0, 50)), snps)
)
bf2 <- rbind(
  stats::setNames(c(8, rep(0, 99)), snps),
  stats::setNames(c(rep(0, 49), 8, rep(0, 50)), snps)
)
result <- coloc:::coloc.bf_bf(bf1, bf2, trim_by_posterior = FALSE)
stopifnot(nrow(result$summary) == 4L, ncol(result$results) == 5L)
merged <- data.table(
  merge_key = snps,
  eqtl_pos = seq_len(100L) + 100000L,
  gwas_a1 = rep("A", 100L),
  gwas_a2 = rep("C", 100L)
)
output <- tempfile("atac_v3_coloc_export_")
dir.create(output)
export_coloc_susie_result(
  result, merged, "fixture_study", "FIXTURE", "ENSG_FIXTURE", 1L, output
)
signals <- fread(file.path(
  output, "fixture_study", "chr1", "ENSG_FIXTURE", "signal_pairs.tsv"
))
posterior <- fread(file.path(
  output, "fixture_study", "chr1", "ENSG_FIXTURE", "variant_posteriors.tsv.gz"
))
stopifnot(
  nrow(signals) == 4L,
  identical(signals$signal_pair_index, 1:4),
  nrow(posterior) == 400L,
  all(posterior[, .(ok = abs(sum(SNP.PP.H4) - 1) < 1e-8),
    by = signal_pair_index]$ok),
  all(c("hg19_position", "allele1", "allele2") %in% names(posterior))
)
single <- coloc:::coloc.bf_bf(
  bf1[1L, , drop = FALSE], bf2[1L, , drop = FALSE],
  trim_by_posterior = FALSE
)
stopifnot(
  nrow(single$summary) == 1L,
  identical(setdiff(names(single$results), "snp"), "SNP.PP.H4.abf")
)
export_coloc_susie_result(
  single, merged, "fixture_study", "SINGLE", "ENSG_SINGLE", 1L, output
)
ambiguous <- rbind(
  merged,
  merged[1L, .(merge_key, eqtl_pos, gwas_a1 = "G", gwas_a2 = "T")]
)
ambiguous_failed <- inherits(try(
  export_coloc_susie_result(
    single, ambiguous, "fixture_study", "AMBIGUOUS", "ENSG_AMBIGUOUS", 1L, output
  ),
  silent = TRUE
), "try-error")
stopifnot(ambiguous_failed)
cat("COLOC_SUSIE_EXPORT_FIXTURE\tPASS\n")
