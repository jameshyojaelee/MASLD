#!/usr/bin/env Rscript
# 13_nafl_vs_nash_de.R
# ---------------------------------------------------------------------------
# NAFL vs NASH differential expression using diagnosis_harmonized column:
#   A) Per-study limma-voom DE (all cohorts with NAFL+NASH samples)
#   B) Fixed-effects meta-analysis across cohorts
#   C) Dream mega-analysis (NAFL + NASH samples pooled)
#
# Uses NAS-based harmonized diagnosis (Kleiner et al. 2005):
#   NAS < 3 → NAFL; NAS 3-4 → Borderline (grouped with NASH); NAS >= 5 → NASH
# See docs/dataset_labeling_and_harmonization.md for full rationale.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(variancePartition)
  library(BiocParallel)
  library(metafor)
  library(ggplot2)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/disease_signatures")
dir.create(RDIR, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
counts <- readRDS(file.path(INT, "results/integration/merged_counts_raw.rds"))
meta   <- readRDS(file.path(INT, "results/integration/meta_matched.rds"))
qc     <- fread(file.path(INT, "qc/sample_qc_report.csv"))
meta   <- meta[sample_id %in% qc[pass_technical == TRUE, sample_id]]

# PRJNA512027 permanently removed from pipeline 2026-05-15 (L0/S0 library batch
# perfectly confounded with disease severity). Defensive filter retained to
# protect against stale merged_dge.rds inputs.
meta <- meta[dataset != "PRJNA512027"]

# Subset to NAFL + NASH using diagnosis_harmonized
# Borderline (NAS 3-4) is grouped with NASH for this analysis
meta_nn <- meta[diagnosis_harmonized %in% c("NAFL", "NASH", "Borderline")]
meta_nn[, nafl_nash := fifelse(diagnosis_harmonized == "NAFL", "NAFL", "NASH")]
meta_nn[, nafl_nash := factor(nafl_nash, levels = c("NAFL", "NASH"))]
cat("Total NAFL+NASH+Borderline samples:", nrow(meta_nn), "\n")
cat("  Breakdown:\n")
print(meta_nn[, .N, by = .(dataset, diagnosis_harmonized)][order(dataset, diagnosis_harmonized)])
cat("\n  Binary grouping (NAFL vs NASH+Borderline):\n")
print(meta_nn[, .N, by = .(dataset, nafl_nash)][order(dataset, nafl_nash)])

# ============================================================
#  A) Per-study DE (NAFL vs NASH using diagnosis_harmonized)
# ============================================================
cat("\n===== PER-STUDY NAFL vs NASH DE =====\n")

perstudy_results <- list()
datasets <- unique(meta_nn$dataset)

# Skip datasets with < 3 samples in either group
for (ds in datasets) {
  m <- meta_nn[dataset == ds]
  n_nafl <- sum(m$nafl_nash == "NAFL")
  n_nash <- sum(m$nafl_nash == "NASH")
  if (n_nafl < 3 || n_nash < 3) {
    cat(sprintf("\n--- %s --- SKIPPED (NAFL: %d, NASH: %d — need >= 3 each)\n", ds, n_nafl, n_nash))
    next
  }

  cat(sprintf("\n--- %s ---\n", ds))
  cat(sprintf("  NAFL: %d, NASH: %d\n", n_nafl, n_nash))

  # Build DGE
  idx <- colnames(counts) %in% m$sample_id
  dge <- DGEList(counts = counts[, idx])
  dge$samples <- cbind(dge$samples, m[match(colnames(dge), m$sample_id),
    .(nafl_nash, sex, inferred_sex, age)])

  # Normalization
  dge <- calcNormFactors(dge, method = "TMM")
  cat("  Norm: TMM\n")

  # Filter
  keep <- filterByExpr(dge, group = dge$samples$nafl_nash)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("  Genes after filter:", nrow(dge), "\n")

  # Design matrix — dataset-specific covariates
  has_sex <- !all(is.na(m$sex))
  has_inferred_sex <- "inferred_sex" %in% names(m) && !all(is.na(m$inferred_sex))
  has_age <- !all(is.na(m$age))
  sex_col <- if (has_sex) "sex" else if (has_inferred_sex) "inferred_sex" else NULL

  formula_parts <- "~ nafl_nash"
  if (has_age) formula_parts <- paste0(formula_parts, " + age")
  if (!is.null(sex_col)) formula_parts <- paste0(formula_parts, " + ", sex_col)

  cat("  Design:", formula_parts, "\n")
  design <- model.matrix(as.formula(formula_parts), data = dge$samples)

  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  tt <- topTable(fit, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$dataset <- ds
  # Add unmoderated SE for downstream metafor (same logic as Script 02 fix)
  coef_idx <- "nafl_nashNASH"
  stdev_unsc <- if (is.matrix(fit$stdev.unscaled)) fit$stdev.unscaled[, coef_idx] else fit$stdev.unscaled
  tt$SE_unmoderated <- stdev_unsc * fit$sigma
  perstudy_results[[ds]] <- as.data.table(tt)

  sig <- sum(tt$adj.P.Val < 0.05)
  sig_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0)
  sig_down <- sum(tt$adj.P.Val < 0.05 & tt$logFC < 0)
  cat(sprintf("  DEGs (padj<0.05): %d (Up: %d, Down: %d)\n", sig, sig_up, sig_down))
}

# Save per-study results
all_perstudy <- rbindlist(perstudy_results, use.names = TRUE, fill = TRUE)
fwrite(all_perstudy, file.path(RDIR, "nafl_vs_nash_per_study.csv"))
cat("\nSaved: nafl_vs_nash_per_study.csv\n")

# Per-study summary
summary_dt <- all_perstudy[, .(
  total_genes = .N,
  degs_005 = sum(adj.P.Val < 0.05),
  degs_up = sum(adj.P.Val < 0.05 & logFC > 0),
  degs_down = sum(adj.P.Val < 0.05 & logFC < 0),
  mean_lfc = round(mean(logFC), 4)
), by = dataset]
fwrite(summary_dt, file.path(RDIR, "nafl_vs_nash_summary.csv"))
cat("\n=== Per-study summary ===\n")
print(summary_dt)

# ============================================================
#  B) Meta-analysis
# ============================================================
cat("\n===== META-ANALYSIS (NAFL vs NASH) =====\n")

# Get genes present in at least 3 datasets
gene_counts <- all_perstudy[, .N, by = gene]
meta_genes <- gene_counts[N >= 3, gene]
cat("Genes in >= 3 datasets:", length(meta_genes), "\n")

meta_results <- rbindlist(lapply(meta_genes, function(g) {
  d <- all_perstudy[gene == g]
  if (nrow(d) < 2) return(NULL)
  tryCatch({
    # Use UNMODERATED SE (stdev.unscaled * sigma) saved by per-study limma-voom.
    # The old formula SE = |logFC / t| gives moderated SE (eBayes-shrunk),
    # which biases tau^2 downward in metafor's REML estimator.
    # Script 13's per-study DE (section A above) saves topTable output which
    # does not yet carry SE_unmoderated, so we fall back to |logFC/t| for those.
    # NOTE: This is a within-script meta (NAFL-vs-NASH contrast), not the main
    # Disease-vs-Control meta in Script 06, so the per-study results come from
    # section A above, not from Script 02's output files.
    sei_vec <- if ("SE_unmoderated" %in% names(d)) d$SE_unmoderated else abs(d$logFC / d$t)
    fit <- rma(yi = d$logFC, sei = sei_vec,
               method = "REML")
    data.table(
      gene = g,
      meta_logFC = fit$beta[1],
      meta_se = fit$se,
      meta_pval = fit$pval,
      meta_zval = fit$zval,
      n_datasets = nrow(d),
      I2 = fit$I2,
      Q_pval = fit$QEp
    )
  }, error = function(e) NULL)
}))

meta_results[, meta_padj := p.adjust(meta_pval, method = "BH")]
meta_results <- meta_results[order(meta_padj)]

sig_meta <- meta_results[meta_padj < 0.05]
cat(sprintf("Meta-analysis DEGs (padj<0.05): %d\n", nrow(sig_meta)))
cat(sprintf("  Up: %d, Down: %d\n",
  sum(sig_meta$meta_logFC > 0), sum(sig_meta$meta_logFC < 0)))

fwrite(meta_results, file.path(RDIR, "nafl_vs_nash_meta.csv"))
cat("Saved: nafl_vs_nash_meta.csv\n")

# ============================================================
#  C) Dream mega-analysis (NAFL vs NASH)
# ============================================================
cat("\n===== DREAM MEGA-ANALYSIS (NAFL vs NASH) =====\n")

# Pool all NAFL + NASH samples
all_nn <- meta_nn[order(sample_id)]
idx_all <- colnames(counts) %in% all_nn$sample_id
dge_all <- DGEList(counts = counts[, idx_all])
dge_all$samples <- cbind(dge_all$samples,
  all_nn[match(colnames(dge_all), all_nn$sample_id),
    .(nafl_nash, dataset, sex, inferred_sex, age)])

# Dataset is the random-effect grouping (PRJNA512027 L0/S0 subbatch logic
# permanently removed 2026-05-15 with the cohort itself).
dge_all$samples$dataset_subbatch <- as.character(dge_all$samples$dataset)

# Use inferred_sex where annotated sex is missing
dge_all$samples$sex_for_model <- dge_all$samples$sex
na_sex <- is.na(dge_all$samples$sex_for_model) | dge_all$samples$sex_for_model == ""
dge_all$samples$sex_for_model[na_sex] <- as.character(dge_all$samples$inferred_sex[na_sex])
dge_all$samples$sex_for_model <- factor(dge_all$samples$sex_for_model)

dge_all <- calcNormFactors(dge_all, method = "TMM")
keep_all <- filterByExpr(dge_all, group = dge_all$samples$nafl_nash)
dge_all <- dge_all[keep_all, , keep.lib.sizes = FALSE]

cat("Dream samples:", ncol(dge_all), "\n")
cat("Dream genes:", nrow(dge_all), "\n")
cat("Batch levels:", paste(sort(unique(dge_all$samples$dataset_subbatch)), collapse = ", "), "\n")
cat("Sex levels:", paste(levels(dge_all$samples$sex_for_model), collapse = ", "), "\n")

# Dream formula — nafl_nash as fixed, dataset as random
form <- ~ nafl_nash + sex_for_model + (1 | dataset_subbatch)
cat("Formula:", deparse(form), "\n")

n_cores <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")), 32)
cat("Using", n_cores, "CPU cores\n")
BPPARAM <- MulticoreParam(n_cores, progressbar = TRUE)

vobjDream <- voomWithDreamWeights(dge_all, form, dge_all$samples, BPPARAM = BPPARAM)
cat("Running dream...\n")
fitDream <- dream(vobjDream, form, dge_all$samples, BPPARAM = BPPARAM)
# NOTE: Do NOT call eBayes() after dream(). dream() already applies moderated
# t-statistics via Satterthwaite approximation. A second eBayes() would double-shrink
# variance estimates, producing anti-conservative p-values. (Bug fixed 2026-04-02)

dream_tt <- topTable(fitDream, coef = "nafl_nashNASH", number = Inf, sort.by = "none")
dream_tt$gene <- rownames(dream_tt)

sig_dream <- sum(dream_tt$adj.P.Val < 0.05)
sig_dream_up <- sum(dream_tt$adj.P.Val < 0.05 & dream_tt$logFC > 0)
sig_dream_down <- sum(dream_tt$adj.P.Val < 0.05 & dream_tt$logFC < 0)
cat(sprintf("Dream DEGs (padj<0.05): %d (Up: %d, Down: %d)\n",
  sig_dream, sig_dream_up, sig_dream_down))
cat(sprintf("Mean LFC: %.4f\n", mean(dream_tt$logFC)))

fwrite(as.data.table(dream_tt), file.path(RDIR, "nafl_vs_nash_dream.csv"))
cat("Saved: nafl_vs_nash_dream.csv\n")

# ============================================================
#  Consensus DEGs
# ============================================================
cat("\n===== CONSENSUS (NAFL vs NASH) =====\n")

dream_dt <- as.data.table(dream_tt)
meta_dt <- meta_results

# Merge dream + meta
consensus <- merge(
  dream_dt[, .(gene, dream_lfc = logFC, dream_padj = adj.P.Val)],  # C2-OK-sensitivity
  meta_dt[, .(gene, meta_logFC, meta_padj)],
  by = "gene", all = TRUE
)

consensus[, `:=`(
  dream_sig = dream_padj < 0.05,  # C2-OK-sensitivity
  meta_sig = meta_padj < 0.05,
  dream_dir = fifelse(dream_lfc > 0, "up", "down"),  # C2-OK-sensitivity
  meta_dir = fifelse(meta_logFC > 0, "up", "down")
)]

# Tier 1: significant in BOTH and concordant
consensus[, tier := "NS"]
consensus[dream_sig == TRUE & meta_sig == TRUE &  # C2-OK-sensitivity
          !is.na(dream_dir) & !is.na(meta_dir) &  # C2-OK-sensitivity
          dream_dir == meta_dir, tier := "Tier1"]  # C2-OK-sensitivity

# Tier 2: significant in at least one, concordant if both exist
consensus[(dream_sig == TRUE | meta_sig == TRUE) &  # C2-OK-sensitivity
          tier == "NS" &
          (is.na(dream_dir) | is.na(meta_dir) | dream_dir == meta_dir),  # C2-OK-sensitivity
          tier := "Tier2"]

cat("Tier 1 (both sig + concordant):", sum(consensus$tier == "Tier1"), "\n")
cat("Tier 2 (one sig, concordant):", sum(consensus$tier == "Tier2"), "\n")

fwrite(consensus, file.path(RDIR, "nafl_vs_nash_consensus.csv"))
cat("Saved: nafl_vs_nash_consensus.csv\n")

# ============================================================
#  Volcano plot
# ============================================================
cat("\nGenerating volcano plots...\n")

pdf(file.path(RDIR, "nafl_vs_nash_volcanos.pdf"), width = 14, height = 6)

# Dream volcano
p1 <- ggplot(dream_dt, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) > 0.5), size = 0.5, alpha = 0.4) +
  scale_color_manual(values = c("grey70", "#E91E63"), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  labs(title = "NAFL vs NASH — Dream Mega-Analysis",
       subtitle = sprintf("%d samples (%d datasets), %d DEGs (padj<0.05)",
                          ncol(dge_all), length(unique(dge_all$samples$dataset)), sig_dream),
       x = "log2 Fold Change (NASH vs NAFL)", y = "-log10(padj)") +
  theme_minimal(base_size = 11)

# Meta volcano
p2 <- ggplot(meta_dt, aes(x = meta_logFC, y = -log10(meta_padj))) +
  geom_point(aes(color = meta_padj < 0.05 & abs(meta_logFC) > 0.5), size = 0.5, alpha = 0.4) +
  scale_color_manual(values = c("grey70", "#9C27B0"), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  labs(title = "NAFL vs NASH — Fixed-Effects Meta-Analysis",
       subtitle = sprintf("%d genes, %d DEGs (padj<0.05)", nrow(meta_dt), nrow(sig_meta)),
       x = "Meta log2 Fold Change (NASH vs NAFL)", y = "-log10(padj)") +
  theme_minimal(base_size = 11)

library(gridExtra)
grid.arrange(p1, p2, ncol = 2)
dev.off()
cat("Saved: nafl_vs_nash_volcanos.pdf\n")

cat("\n=== Script 13 complete ===\n")
