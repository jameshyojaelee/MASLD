#!/usr/bin/env Rscript
# V6-a Test 3: B6 F2 breakpoint via AUTOMATED change-point detection
# =====================================================================
# Our claim: 36.3% of DEGs are "switch-like" with hardcoded F2 (fib>=2)
# breakpoint. AIC step > linear model.
#
# Competitor: strucchange::breakpoints() finds the OPTIMAL breakpoint
# location per gene via dynamic programming on residual SS.
#
# Test: For each DEG, let strucchange find the best 1-breakpoint on
# expression vs fibrosis_stage. Which F level is most commonly the
# breakpoint? If most genes "prefer" F1 (not F2), our claim fails.
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(segmented)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUT_DIR <- file.path(BASE, "docs/manuscript/verification/method_comparison")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Load data ===\n")
dge <- readRDS(file.path(INT, "results/integration/merged_dge.rds"))
dream <- fread(file.path(INT, "results/integration/dream_results_ashr.csv"))
degs <- dream[padj < 0.05 & abs(logFC) > 0.3]
cat("DEGs:", nrow(degs), "\n")

meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"),
              select = c("sample_id", "fibrosis_stage", "dataset"))
meta <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta[, fibrosis_stage := as.integer(fibrosis_stage)]
shared_samples <- intersect(meta$sample_id, colnames(dge))
meta <- meta[sample_id %in% shared_samples]
shared_genes <- intersect(unique(degs$gene), rownames(dge))
dge_sub <- dge[shared_genes, shared_samples]
logcpm <- edgeR::cpm(dge_sub, log = TRUE, prior.count = 1)
stage_vec <- meta[match(colnames(logcpm), sample_id), fibrosis_stage]

cat("Samples:", length(stage_vec), "  Stage distribution:\n")
print(table(stage_vec))

# Method: For each gene, fit segmented::segmented() to find the optimal
# single breakpoint on expression ~ fibrosis_stage. Returns continuous
# breakpoint psi in [0, 4]; we round to nearest integer F level.

n <- length(stage_vec)
fib_numeric <- as.numeric(stage_vec)

cat("\nRunning segmented::segmented on", nrow(logcpm), "DEGs...\n")
cat("(Optimal single breakpoint per gene)\n")

# Estimate timing on subsample
set.seed(42)
t0 <- Sys.time()
test_idx <- sample(nrow(logcpm), 50)
for (i in test_idx) {
  y <- logcpm[i, ]
  bp <- tryCatch({
    lin_fit <- lm(y ~ fib_numeric)
    seg <- segmented(lin_fit, seg.Z = ~ fib_numeric,
                     psi = 2,
                     control = seg.control(n.boot = 0, it.max = 30))
    seg$psi[1, "Est."]
  }, error = function(e) NA, warning = function(w) NA)
}
dt <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("  Estimate: %.2fs per 50 genes ≈ %.1f min total\n",
            dt, dt * nrow(logcpm) / 50 / 60))

# Full run
cat("\nFull run...\n")
bp_psi <- numeric(nrow(logcpm))
bp_rss_reduction <- numeric(nrow(logcpm))
bp_converged <- logical(nrow(logcpm))

t0 <- Sys.time()
for (i in seq_len(nrow(logcpm))) {
  y <- logcpm[i, ]
  result <- tryCatch({
    lin_fit <- lm(y ~ fib_numeric)
    rss0 <- sum(residuals(lin_fit)^2)
    seg <- suppressWarnings(segmented(lin_fit, seg.Z = ~ fib_numeric,
                     psi = 2,
                     control = seg.control(n.boot = 0, it.max = 30)))
    psi_val <- seg$psi[1, "Est."]
    rss_seg <- sum(residuals(seg)^2)
    list(psi = psi_val, rss_red = (rss0 - rss_seg) / rss0, converged = TRUE)
  }, error = function(e) list(psi = NA, rss_red = NA, converged = FALSE))

  bp_psi[i] <- result$psi
  bp_rss_reduction[i] <- result$rss_red
  bp_converged[i] <- result$converged

  if (i %% 500 == 0) {
    dt <- as.numeric(Sys.time() - t0, units = "mins")
    cat(sprintf("  %d / %d  (%.1f min, %.0f%% converged)\n",
                i, nrow(logcpm), dt, 100*mean(bp_converged[1:i], na.rm=TRUE)))
  }
}
cat(sprintf("Done: %.1f min  Converged: %d / %d\n",
            as.numeric(Sys.time() - t0, units = "mins"),
            sum(bp_converged), length(bp_converged)))

# Classify breakpoint by nearest fibrosis transition
# Breakpoint psi in [0, 4]; assign to nearest of F0/F1, F1/F2, F2/F3, F3/F4
# (boundaries at 0.5, 1.5, 2.5, 3.5)

per_gene <- data.table(
  gene = rownames(logcpm),
  bp_psi = bp_psi,
  rss_reduction = bp_rss_reduction,
  converged = bp_converged
)

classify_psi <- function(psi) {
  if (is.na(psi)) return(NA_character_)
  if (psi < 0.5) return("F0/F1")  # also catches values <= 0
  if (psi < 1.5) return("F0/F1")
  if (psi < 2.5) return("F1/F2")
  if (psi < 3.5) return("F2/F3")
  return("F3/F4")
}
per_gene[, nearest_transition := sapply(bp_psi, classify_psi)]
# Also a finer "round to nearest integer" classification
per_gene[, psi_rounded := round(bp_psi)]

# Summary
cat("\n=== Distribution of nearest fibrosis transition for optimal breakpoint ===\n")
print(table(per_gene$nearest_transition, useNA = "ifany"))

cat("\n=== As % of DEGs ===\n")
tab <- table(per_gene$nearest_transition, useNA = "no")
print(round(100 * tab / sum(tab), 2))

# Restrict to genes with substantial RSS reduction
sig_genes <- per_gene[!is.na(rss_reduction) & rss_reduction > 0.05]
cat(sprintf("\nGenes with >5%% RSS reduction: %d / %d\n",
            nrow(sig_genes), nrow(per_gene)))
cat("Distribution among these:\n")
tab_sig <- table(sig_genes$nearest_transition, useNA = "no")
print(round(100 * tab_sig / sum(tab_sig), 2))

# Strong changepoints (>20% RSS reduction)
strong <- per_gene[!is.na(rss_reduction) & rss_reduction > 0.20]
cat(sprintf("\nGenes with >20%% RSS reduction: %d / %d\n",
            nrow(strong), nrow(per_gene)))
cat("Distribution:\n")
tab_strong <- table(strong$nearest_transition, useNA = "no")
print(round(100 * tab_strong / sum(tab_strong), 2))

# Save
fwrite(per_gene, file.path(OUT_DIR, "test3_B6_auto_breakpoint_per_gene.csv"))

summary_tab <- data.table(
  level = c("all DEGs", ">5% RSS reduction", ">20% RSS reduction"),
  n = c(sum(!is.na(per_gene$bp_index)), nrow(sig_genes), nrow(strong)),
  pct_F0_F1 = c(
    mean(per_gene$nearest_transition == "F0/F1", na.rm = TRUE) * 100,
    mean(sig_genes$nearest_transition == "F0/F1", na.rm = TRUE) * 100,
    mean(strong$nearest_transition == "F0/F1", na.rm = TRUE) * 100
  ),
  pct_F1_F2 = c(
    mean(per_gene$nearest_transition == "F1/F2", na.rm = TRUE) * 100,
    mean(sig_genes$nearest_transition == "F1/F2", na.rm = TRUE) * 100,
    mean(strong$nearest_transition == "F1/F2", na.rm = TRUE) * 100
  ),
  pct_F2_F3 = c(
    mean(per_gene$nearest_transition == "F2/F3", na.rm = TRUE) * 100,
    mean(sig_genes$nearest_transition == "F2/F3", na.rm = TRUE) * 100,
    mean(strong$nearest_transition == "F2/F3", na.rm = TRUE) * 100
  ),
  pct_F3_F4 = c(
    mean(per_gene$nearest_transition == "F3/F4", na.rm = TRUE) * 100,
    mean(sig_genes$nearest_transition == "F3/F4", na.rm = TRUE) * 100,
    mean(strong$nearest_transition == "F3/F4", na.rm = TRUE) * 100
  )
)
print(summary_tab)
fwrite(summary_tab, file.path(OUT_DIR, "test3_B6_auto_breakpoint_summary.csv"))
cat("\nSaved to", OUT_DIR, "\n")
cat("=== Done ===\n")
