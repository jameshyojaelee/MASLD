#!/usr/bin/env Rscript
# Compare the primary-label refits across rosters and count arms.
# Call = canonical BH padj < 0.05 and |logFC| > 0.5 (the refit's canonical family).
suppressPackageStartupMessages(library(data.table))
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 1L) stop("Usage: compare_refit_arms.R output_dir")
out <- argv[[1]]
if (file.exists(out)) stop("Refusing to overwrite output dir")
dir.create(out, recursive = TRUE)
set.seed(20260923)
cand <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates"
fits <- list(
  prior_844       = "bulk-corrected-controls-20260923T231939Z",
  n774_primary    = "bulk-corrected-controls-20261008T162308Z",
  n774_sensitivity = "bulk-corrected-controls-20261008T162309Z",
  n738_primary    = "bulk-corrected-controls-20261008T162310Z",
  n738_sensitivity = "bulk-corrected-controls-20261008T162311Z")
read_fit <- function(root, arm = "primary_source_controls") {
  x <- fread(file.path(cand, root, arm, "deg_results.csv"), select = c("gene", "logFC", "padj", "symbol"))
  x[, call := !is.na(padj) & padj < 0.05 & abs(logFC) > 0.5]
  x[, direction := fifelse(call, fifelse(logFC > 0, "up", "down"), "none")]
  x
}
d <- lapply(fits, read_fit)
cmp <- function(a, b) {
  m <- merge(d[[a]], d[[b]], by = "gene", suffixes = c(".a", ".b"))
  data.table(a = a, b = b, n_common = nrow(m), calls_a_all = sum(d[[a]]$call), calls_b_all = sum(d[[b]]$call),
             called_both = m[, sum(call.a & call.b)], only_a = m[, sum(call.a & !call.b)], only_b = m[, sum(!call.a & call.b)],
             sign_flip_both_called = m[call.a & call.b, sum(sign(logFC.a) != sign(logFC.b))],
             logFC_pearson = m[, cor(logFC.a, logFC.b)], logFC_spearman = m[, cor(logFC.a, logFC.b, method = "spearman")])
}
res <- rbindlist(list(cmp("prior_844", "n774_primary"), cmp("n774_primary", "n774_sensitivity"),
                      cmp("n738_primary", "n738_sensitivity"), cmp("n774_primary", "n738_primary")))
fwrite(res, file.path(out, "refit_comparison.tsv"), sep = "\t")
m <- merge(d$n774_primary, d$n774_sensitivity, by = c("gene", "symbol"), suffixes = c(".primary", ".sensitivity"))
fwrite(m[call.primary != call.sensitivity, .(gene, symbol, logFC.primary, padj.primary, call.primary, logFC.sensitivity, padj.sensitivity, call.sensitivity)],
       file.path(out, "n774_arm_discordant_genes.tsv"), sep = "\t")
g <- merge(d$prior_844, d$n774_primary, by = c("gene", "symbol"), suffixes = c(".prior", ".new"))
fwrite(g[call.prior != call.new, .(gene, symbol, logFC.prior, call.prior, logFC.new, call.new)],
       file.path(out, "n774_vs_prior_844_changed_calls.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(res)
