#!/usr/bin/env Rscript
# _smoke_metafor.R — validate the finalized run_metafor recipe on real data.
# Binary per-study DE for the 5 mega cohorts -> run_metafor -> report the
# tiered DEG counts (REML no-HKSJ vs HKSJ vs FE). One-off; safe to delete.
suppressWarnings(source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT"),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/multimethod_validation/de_validation_helpers.R")))

bp <- bp_param()
d  <- load_mega_data()
mega <- mega_cohorts()
cat("Cohorts:", paste(mega, collapse=", "), "| n =", ncol(d$counts), "\n")

cat("Computing BINARY per-study DE for", length(mega), "cohorts...\n")
ps <- lapply(mega, function(ds) {
  sel <- d$meta$dataset == ds
  r <- run_per_study_voom(ds, d$counts[, sel, drop=FALSE], d$meta[sel], contrast_mode = "binary")
  cat(sprintf("  %s: %d genes, median SE_unmod = %.4f\n", ds, nrow(r), median(r$SE_unmoderated, na.rm=TRUE)))
  r
})

cat("\nRunning finalized run_metafor (REML no-HKSJ + HKSJ + FE)...\n")
t0 <- Sys.time()
res <- run_metafor(ps, K = length(mega), bp)
cat(sprintf("Elapsed %.1f min | genes meta-analyzed: %d\n",
            as.numeric(difftime(Sys.time(), t0, units="mins")), nrow(res)))

f05 <- function(p) res[!is.na(p) & p < 0.05 & abs(logFC) > 0.5, .N]
cat("\n=== TIERED metafor DEG counts (binary estimand) ===\n")
cat(sprintf("  REML no-HKSJ : padj<0.05 = %5d | padj<0.05 & |LFC|>0.5 = %5d\n",
            res[!is.na(padj) & padj<0.05, .N], f05(res$padj)))
cat(sprintf("  REML + HKSJ  : padj<0.05 = %5d | padj<0.05 & |LFC|>0.5 = %5d\n",
            res[!is.na(padj_hksj) & padj_hksj<0.05, .N], f05(res$padj_hksj)))
cat(sprintf("  Fixed-effect : padj<0.05 = %5d | padj<0.05 & |LFC|>0.5 = %5d\n",
            res[!is.na(padj_fe) & padj_fe<0.05, .N], f05(res$padj_fe)))
cat(sprintf("\n  median I2 = %.1f%% | median tau2 = %.4f\n",
            median(res$I2, na.rm=TRUE), median(res$tau2, na.rm=TRUE)))
cat("\nSanity: no-HKSJ DEGs >> HKSJ DEGs (expected at K=5 high-I2)? ",
    res[!is.na(padj) & padj<0.05, .N] > 10 * max(1, res[!is.na(padj_hksj) & padj_hksj<0.05, .N]),
    "\n")
cat("schema cols:", paste(names(res), collapse=", "), "\n")
cat("SMOKE OK\n")
