#!/usr/bin/env Rscript
# 14_fibrosis_progression_de.R
# ---------------------------------------------------------------------------
# Fibrosis progression analysis:
#   A) Pairwise stage contrasts (F0-F1, F1-F2, F2-F3, F3-F4)
#   B) Ordinal DE (fibrosis_stage as continuous covariate)
#   C) Dream mega-analysis with fibrosis as ordinal
#   D) Early vs Late gene classification
# Uses GSE130970 + GSE135251 (both have F0-F4 annotation)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(variancePartition)
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

# Subset to datasets with fibrosis annotation (F0-F4)
meta_fib <- meta[!is.na(fibrosis_stage) & dataset %in% c("GSE130970", "GSE135251", "GSE213621")]
cat("Fibrosis-annotated samples:", nrow(meta_fib), "\n")

cat("  Per dataset:\n")
print(meta_fib[, .N, by = .(dataset, fibrosis_stage)][order(dataset, fibrosis_stage)])

# ============================================================
#  A) Per-study ordinal DE (fibrosis as continuous)
# ============================================================
cat("\n===== ORDINAL FIBROSIS DE =====\n")

ordinal_results <- list()

for (ds in c("GSE130970", "GSE135251", "GSE213621")) {
  cat(sprintf("\n--- %s ---\n", ds))

  m <- meta_fib[dataset == ds]
  
  idx <- colnames(counts) %in% m$sample_id
  dge <- DGEList(counts = counts[, idx])
  dge$samples <- cbind(dge$samples,
    m[match(colnames(dge), m$sample_id), .(fibrosis_stage, condition, sex, inferred_sex, age)])
  dge <- calcNormFactors(dge, method = "TMM")
  
  keep <- filterByExpr(dge, group = factor(dge$samples$fibrosis_stage))
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("  Genes after filter:", nrow(dge), "\n")
  cat("  Samples:", ncol(dge), "\n")
  
  # Model 1: fibrosis only (unadjusted)
  if (ds == "GSE130970") {
    design <- model.matrix(~ fibrosis_stage + age + sex, data = dge$samples)
    cat("  Design: ~ fibrosis_stage + age + sex\n")
  } else if (ds == "GSE213621") {
    design <- model.matrix(~ fibrosis_stage + inferred_sex, data = dge$samples)
    cat("  Design: ~ fibrosis_stage + inferred_sex\n")
  } else {
    design <- model.matrix(~ fibrosis_stage, data = dge$samples)
    cat("  Design: ~ fibrosis_stage\n")
  }

  
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  tt <- topTable(fit, coef = "fibrosis_stage", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$dataset <- ds
  ordinal_results[[ds]] <- as.data.table(tt)
  
  sig <- sum(tt$adj.P.Val < 0.05)
  cat(sprintf("  Ordinal DEGs (padj<0.05): %d\n", sig))
  cat(sprintf("  Mean slope: %.4f (per fibrosis stage)\n", mean(tt$logFC)))
  
  # Model 2: fibrosis + condition covariate (to separate fibrosis from NAFL/NASH)
  cat("  --- With condition covariate ---\n")
  dge$samples$condition_binary <- fifelse(
    dge$samples$condition %in% c("NASH", "NASH_Fibrosis"), "NASH_spectrum", "non_NASH"
  )
  
  if (length(unique(dge$samples$condition_binary)) < 2) {
    cat("  Skipping condition adjustment: only 1 condition level present\n")
  } else {
    if (ds == "GSE130970") {
      design2 <- model.matrix(~ fibrosis_stage + condition_binary + age + sex, data = dge$samples)
    } else if (ds == "GSE213621") {
      design2 <- model.matrix(~ fibrosis_stage + condition_binary + inferred_sex, data = dge$samples)
    } else {
      design2 <- model.matrix(~ fibrosis_stage + condition_binary, data = dge$samples)
    }

    v2 <- voom(dge, design2, plot = FALSE)
    fit2 <- lmFit(v2, design2)
    fit2 <- eBayes(fit2)
    tt2 <- topTable(fit2, coef = "fibrosis_stage", number = Inf, sort.by = "none")
    sig2 <- sum(tt2$adj.P.Val < 0.05)
    cat(sprintf("  Condition-adjusted ordinal DEGs: %d\n", sig2))
  }
}

# Save ordinal results
all_ordinal <- rbindlist(ordinal_results, use.names = TRUE, fill = TRUE)
fwrite(all_ordinal, file.path(RDIR, "fibrosis_ordinal_per_study.csv"))
cat("\nSaved: fibrosis_ordinal_per_study.csv\n")

# Meta-analyze ordinal slopes
cat("\n===== META-ANALYSIS OF FIBROSIS SLOPES =====\n")

gene_counts <- all_ordinal[, .N, by = gene]
common_genes <- gene_counts[N >= 2, gene]
cat("Genes in >= 2 cohorts:", length(common_genes), "\n")


slope_meta <- rbindlist(lapply(common_genes, function(g) {
  d <- all_ordinal[gene == g]
  if (nrow(d) < 2) return(NULL)
  tryCatch({
    se_est <- abs(d$logFC / d$t)
    se_est[se_est == 0 | !is.finite(se_est)] <- 1
    fit <- rma(yi = d$logFC, sei = se_est, method = "FE")
    data.table(
      gene = g,
      slope_meta = fit$beta[1],
      slope_se = fit$se,
      slope_pval = fit$pval,
      I2 = fit$I2
    )
  }, error = function(e) NULL)
}))

slope_meta[, slope_padj := p.adjust(slope_pval, method = "BH")]
slope_meta <- slope_meta[order(slope_padj)]

cat(sprintf("Meta-analyzed fibrosis slope DEGs (padj<0.05): %d\n", sum(slope_meta$slope_padj < 0.05)))
cat(sprintf("  Positive slope (up with fibrosis): %d\n",
  sum(slope_meta$slope_padj < 0.05 & slope_meta$slope_meta > 0)))
cat(sprintf("  Negative slope (down with fibrosis): %d\n",
  sum(slope_meta$slope_padj < 0.05 & slope_meta$slope_meta < 0)))

fwrite(slope_meta, file.path(RDIR, "fibrosis_slopes_meta.csv"))
cat("Saved: fibrosis_slopes_meta.csv\n")

# ============================================================
#  B) Pairwise stage contrasts
# ============================================================
cat("\n===== PAIRWISE FIBROSIS STAGE CONTRASTS =====\n")

pairwise_results <- list()
stage_pairs <- list(c(0, 1), c(1, 2), c(2, 3), c(3, 4))

for (ds in c("GSE130970", "GSE135251", "GSE213621")) {
  m_ds <- meta_fib[dataset == ds]

  
  for (pair in stage_pairs) {
    f_low <- pair[1]
    f_high <- pair[2]
    contrast_name <- sprintf("F%d_vs_F%d", f_high, f_low)
    
    m_pair <- m_ds[fibrosis_stage %in% c(f_low, f_high)]
    if (nrow(m_pair) < 6) {
      cat(sprintf("  %s %s: skipped (n=%d too small)\n", ds, contrast_name, nrow(m_pair)))
      next
    }
    
    m_pair[, fib_group := factor(fifelse(fibrosis_stage == f_high, "high", "low"),
                                  levels = c("low", "high"))]
    
    idx <- colnames(counts) %in% m_pair$sample_id
    dge <- DGEList(counts = counts[, idx])
    dge$samples <- cbind(dge$samples,
      m_pair[match(colnames(dge), m_pair$sample_id), .(fib_group, sex, inferred_sex, age)])
    dge <- calcNormFactors(dge, method = "TMM")
    keep <- filterByExpr(dge, group = dge$samples$fib_group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    if (ds == "GSE130970") {
      design <- model.matrix(~ fib_group + age + sex, data = dge$samples)
    } else if (ds == "GSE213621") {
      design <- model.matrix(~ fib_group + inferred_sex, data = dge$samples)
    } else {
      design <- model.matrix(~ fib_group, data = dge$samples)
    }

    
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    
    tt <- topTable(fit, coef = "fib_grouphigh", number = Inf, sort.by = "none")
    tt$gene <- rownames(tt)
    tt$dataset <- ds
    tt$contrast <- contrast_name
    pairwise_results[[paste(ds, contrast_name)]] <- as.data.table(tt)
    
    sig <- sum(tt$adj.P.Val < 0.05)
    cat(sprintf("  %s %s: n=%d, DEGs=%d\n", ds, contrast_name, nrow(m_pair), sig))
  }
}

all_pairwise <- rbindlist(pairwise_results, use.names = TRUE, fill = TRUE)
fwrite(all_pairwise, file.path(RDIR, "fibrosis_pairwise.csv"))
cat("\nSaved: fibrosis_pairwise.csv\n")

# ============================================================
#  C) Dream mega-analysis
# ============================================================
cat("\n===== DREAM FIBROSIS MEGA-ANALYSIS =====\n")

idx_fib <- colnames(counts) %in% meta_fib$sample_id
dge_fib <- DGEList(counts = counts[, idx_fib])
dge_fib$samples <- cbind(dge_fib$samples,
  meta_fib[match(colnames(dge_fib), meta_fib$sample_id),
    .(fibrosis_stage, condition, dataset, inferred_sex)])
dge_fib <- calcNormFactors(dge_fib, method = "TMM")
keep_fib <- filterByExpr(dge_fib, group = factor(dge_fib$samples$fibrosis_stage))
dge_fib <- dge_fib[keep_fib, , keep.lib.sizes = FALSE]

cat("Dream fibrosis samples:", ncol(dge_fib), "\n")
cat("Dream fibrosis genes:", nrow(dge_fib), "\n")

form_fib <- ~ fibrosis_stage + inferred_sex + (1 | dataset)
n_cores <- min(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")), 32)

BPPARAM <- BiocParallel::SnowParam(n_cores, progressbar = TRUE)

vobjFib <- voomWithDreamWeights(dge_fib, form_fib, dge_fib$samples, BPPARAM = BPPARAM)
cat("Running dream (fibrosis)...\n")
fitFib <- dream(vobjFib, form_fib, dge_fib$samples, BPPARAM = BPPARAM)
# NOTE: Do NOT call eBayes() after dream(). dream() already applies moderated
# t-statistics via Satterthwaite approximation. A second eBayes() would double-shrink
# variance estimates, producing anti-conservative p-values. (Bug fixed 2026-04-02)

dream_fib <- topTable(fitFib, coef = "fibrosis_stage", number = Inf, sort.by = "none")
dream_fib$gene <- rownames(dream_fib)

sig_fib <- sum(dream_fib$adj.P.Val < 0.05)
cat(sprintf("Dream fibrosis DEGs (padj<0.05): %d\n", sig_fib))
cat(sprintf("  Positive slope: %d\n", sum(dream_fib$adj.P.Val < 0.05 & dream_fib$logFC > 0)))
cat(sprintf("  Negative slope: %d\n", sum(dream_fib$adj.P.Val < 0.05 & dream_fib$logFC < 0)))

fwrite(as.data.table(dream_fib), file.path(RDIR, "fibrosis_dream.csv"))
cat("Saved: fibrosis_dream.csv\n")

# ============================================================
#  D) Early vs Late gene classification
# ============================================================
cat("\n===== EARLY vs LATE GENE CLASSIFICATION =====\n")

# For top fibrosis genes, compute mean expression at each stage
sig_genes <- dream_fib$gene[dream_fib$adj.P.Val < 0.05]
if (length(sig_genes) > 0) {
  logcpm <- cpm(dge_fib, log = TRUE, prior.count = 1)
  
  # Compute per-stage mean for sig genes
  stages <- sort(unique(dge_fib$samples$fibrosis_stage))
  stage_means <- do.call(cbind, lapply(stages, function(s) {
    cols <- which(dge_fib$samples$fibrosis_stage == s)
    rowMeans(logcpm[sig_genes, cols, drop = FALSE])
  }))
  colnames(stage_means) <- paste0("F", stages)
  
  # Classify: where is the biggest jump?
  # "Early" = biggest stage-to-stage delta in the first half of the progression axis
  # "Late"  = biggest stage-to-stage delta in the second half
  # The threshold is dynamic based on the number of detected stages, so the
  # classification is valid even when not all 5 fibrosis stages are present.
  n_stages <- length(stages)
  classify_gene <- function(means) {
    deltas <- diff(means)
    max_idx <- which.max(abs(deltas))
    if (max_idx <= floor((n_stages - 1) / 2)) return("Early")
    else return("Late")
  }

  classifications <- apply(stage_means, 1, classify_gene)

  # Build early_late table with dynamic stage columns (no hard-coded F0–F4 assumption)
  early_late <- as.data.table(stage_means)
  setnames(early_late, paste0("F", stages, "_mean"))
  early_late[, gene        := sig_genes]
  early_late[, stage_class := classifications]
  setcolorder(early_late, c("gene", "stage_class"))
  
  cat(sprintf("  Early-response genes: %d\n", sum(early_late$stage_class == "Early")))
  cat(sprintf("  Late-response genes: %d\n", sum(early_late$stage_class == "Late")))
  
  fwrite(early_late, file.path(RDIR, "fibrosis_early_late.csv"))
  cat("Saved: fibrosis_early_late.csv\n")
} else {
  cat("  No significant fibrosis genes to classify.\n")
}

# ============================================================
#  Progression heatmap (top 50 genes)
# ============================================================
cat("\nGenerating fibrosis progression plots...\n")

pdf(file.path(RDIR, "fibrosis_progression.pdf"), width = 12, height = 8)

# Volcano of dream fibrosis slopes
dream_fib_dt <- as.data.table(dream_fib)
p1 <- ggplot(dream_fib_dt, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), size = 0.5, alpha = 0.4) +
  scale_color_manual(values = c("grey70", "#FF5722"), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Fibrosis Progression — Meta-analysis (ordinal slope)",
       subtitle = sprintf("%d samples (GSE130970+GSE135251+GSE213621), %d DEGs", ncol(dge_fib), sig_fib),
       x = "log2FC per fibrosis stage", y = "-log10(padj)") +

  theme_minimal(base_size = 11)
print(p1)

# Top genes heatmap (if we have significant genes)
if (exists("early_late") && nrow(early_late) > 0) {
  top_n <- min(50, nrow(early_late))
  top_genes <- early_late[order(abs(F4_mean - F0_mean), decreasing = TRUE)][1:top_n]
  
  # Melt for ggplot
  hm_data <- melt(top_genes, id.vars = c("gene", "stage_class"),
                   measure.vars = c("F0_mean", "F1_mean", "F2_mean", "F3_mean", "F4_mean"),
                   variable.name = "stage", value.name = "logCPM")
  hm_data[, stage := gsub("_mean", "", stage)]
  hm_data[, gene := factor(gene, levels = rev(top_genes$gene))]
  
  p2 <- ggplot(hm_data, aes(x = stage, y = gene, fill = logCPM)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2196F3", mid = "white", high = "#E91E63",
                          midpoint = median(hm_data$logCPM)) +
    facet_wrap(~ stage_class, scales = "free_y") +
    labs(title = "Fibrosis Progression — Top 50 Genes",
         x = "Fibrosis Stage", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_text(size = 6))
  print(p2)
}

dev.off()
cat("Saved: fibrosis_progression.pdf\n")

cat("\n=== Script 14 complete ===\n")
