#!/usr/bin/env Rscript
# =============================================================================
# Script 195: Critical Slowing Down Analysis at the F2 Tipping Point
# =============================================================================
# Tests dynamical systems hallmarks of a tipping point at the F2→F3 boundary:
#   (a) Increased expression variance (flickering)
#   (b) Jensen-Shannon divergence between adjacent stages
#   (c) Leading eigenvalue dominance (PC1 variance explained)
#   (d) Mean pairwise gene correlation (coordinated response)
#
# Input:
#   - merged_dge.rds (1,444 samples)
#   - dream_results_ashr.csv (DEG annotations)
#   - modeling_metadata.csv (fibrosis_stage)
#
# Output:
#   - figures/supplementary/progression/figS_critical_slowing_down.pdf
#   - RNA-seq/Human/.../results/progression/critical_slowing_down.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

# --- Paths ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")

source(file.path(BASE, "scripts/figures/publication_theme.R"))

FIG_DIR <- file.path(BASE, "figures/supplementary/figS02_progression")
RES_DIR <- file.path(INT, "results/progression")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
cat("Loading merged DGE...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

cat("Loading modeling metadata for fibrosis staging...\n")
modeling_meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"))

cat("Loading dream results...\n")
dream_res <- fread(file.path(RDIR, "dream_results_ashr.csv"))

# --- Merge fibrosis stage ---
sample_meta <- data.table(
  sample_id = colnames(dge),
  dataset   = dge$samples$dataset
)
sample_meta <- merge(sample_meta,
  modeling_meta[, .(sample_id, fibrosis_stage)],
  by = "sample_id", all.x = TRUE)

# Filter to samples with fibrosis staging (F0-F4)
fib_samples <- sample_meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
cat(sprintf("Samples with fibrosis staging: %d\n", nrow(fib_samples)))
cat("Per-stage counts:\n")
print(table(fib_samples$fibrosis_stage))

# --- Prepare expression matrix (logCPM) ---
cat("Computing logCPM...\n")
dge_sub <- dge[, fib_samples$sample_id]
logcpm <- cpm(dge_sub, log = TRUE, prior.count = 1)

# Select top 5,000 variable genes
gene_vars <- apply(logcpm, 1, var)
top5k <- names(sort(gene_vars, decreasing = TRUE))[1:5000]
logcpm_top <- logcpm[top5k, ]
cat(sprintf("Using top %d variable genes\n", length(top5k)))

# Top DEGs for correlation analysis (padj < 0.05, |logFC| > 0.5; migrated 0.3 -> 0.5)
sig_genes <- dream_res[padj < 0.05 & abs(logFC) > 0.5]$gene
sig_in_matrix <- intersect(sig_genes, rownames(logcpm))
cat(sprintf("Significant DEGs in expression matrix: %d\n", length(sig_in_matrix)))

# Use up to 2,000 DEGs for correlation (computational tractability)
if (length(sig_in_matrix) > 2000) {
  # Prioritize by absolute t-statistic
  sig_dt <- dream_res[gene %in% sig_in_matrix][order(-abs(t))]
  sig_in_matrix <- sig_dt$gene[1:2000]
}

# Stage labels
stages <- c("F0", "F1", "F2", "F3", "F4")
fib_samples[, stage_label := paste0("F", fibrosis_stage)]

# ============================================================
# Analysis 1: Per-stage coefficient of variation
# ============================================================
cat("\n--- Analysis 1: Per-stage CV ---\n")

cv_results <- rbindlist(lapply(0:4, function(s) {
  samp <- fib_samples[fibrosis_stage == s]$sample_id
  if (length(samp) < 5) return(NULL)

  mat <- logcpm_top[, samp, drop = FALSE]
  # Per-gene CV across samples within this stage
  gene_means <- rowMeans(mat)
  gene_sds   <- apply(mat, 1, sd)
  # Avoid division by near-zero means: use absolute value
  gene_cv <- gene_sds / pmax(abs(gene_means), 0.1)

  data.table(
    fibrosis_stage = s,
    stage_label    = paste0("F", s),
    mean_cv        = mean(gene_cv),
    median_cv      = median(gene_cv),
    sd_cv          = sd(gene_cv),
    se_cv          = sd(gene_cv) / sqrt(length(gene_cv)),
    n_samples      = length(samp)
  )
}))

cat("Per-stage CV:\n")
print(cv_results)

# ============================================================
# Analysis 2: Jensen-Shannon divergence between adjacent stages
# ============================================================
cat("\n--- Analysis 2: Jensen-Shannon divergence ---\n")

# JSD for two distributions (discretized gene expression)
compute_jsd <- function(x, y, n_bins = 100) {
  # Pool range for consistent binning
  rng <- range(c(x, y))
  breaks <- seq(rng[1], rng[2], length.out = n_bins + 1)

  # Discretize into probability distributions
  p <- tabulate(findInterval(x, breaks, all.inside = TRUE), nbins = n_bins)
  q <- tabulate(findInterval(y, breaks, all.inside = TRUE), nbins = n_bins)

  # Normalize to probabilities
  p <- (p + 1e-10) / sum(p + 1e-10)
  q <- (q + 1e-10) / sum(q + 1e-10)

  # JSD = 0.5 * KL(P||M) + 0.5 * KL(Q||M) where M = (P+Q)/2
  m <- (p + q) / 2
  kl_pm <- sum(p * log2(p / m))
  kl_qm <- sum(q * log2(q / m))

  0.5 * kl_pm + 0.5 * kl_qm
}

# Compute JSD between adjacent stages using all top5k gene expressions
jsd_results <- rbindlist(lapply(0:3, function(s) {
  samp1 <- fib_samples[fibrosis_stage == s]$sample_id
  samp2 <- fib_samples[fibrosis_stage == (s + 1)]$sample_id
  if (length(samp1) < 5 || length(samp2) < 5) return(NULL)

  # Compute JSD per gene, then average
  jsds <- sapply(1:nrow(logcpm_top), function(g) {
    compute_jsd(logcpm_top[g, samp1], logcpm_top[g, samp2])
  })

  data.table(
    transition  = paste0("F", s, "->F", s + 1),
    from_stage  = s,
    to_stage    = s + 1,
    mean_jsd    = mean(jsds),
    median_jsd  = median(jsds),
    sd_jsd      = sd(jsds),
    se_jsd      = sd(jsds) / sqrt(length(jsds)),
    n_genes     = length(jsds)
  )
}))

cat("JSD between adjacent stages:\n")
print(jsd_results)

# ============================================================
# Analysis 3: Leading eigenvalue (PC1 variance explained)
# ============================================================
cat("\n--- Analysis 3: PC1 variance explained ---\n")

pca_results <- rbindlist(lapply(0:4, function(s) {
  samp <- fib_samples[fibrosis_stage == s]$sample_id
  if (length(samp) < 10) return(NULL)

  mat <- logcpm_top[, samp, drop = FALSE]
  # Center genes
  mat_centered <- t(scale(t(mat), center = TRUE, scale = FALSE))

  # PCA on samples (genes as features)
  pca <- prcomp(t(mat_centered), center = FALSE, scale. = FALSE)
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

  data.table(
    fibrosis_stage   = s,
    stage_label      = paste0("F", s),
    pc1_var_pct      = var_explained[1] * 100,
    pc2_var_pct      = var_explained[2] * 100,
    pc1_pc2_ratio    = var_explained[1] / var_explained[2],
    n_samples        = length(samp)
  )
}))

cat("PC1 variance explained per stage:\n")
print(pca_results)

# ============================================================
# Analysis 4: Mean pairwise gene correlation (top DEGs)
# ============================================================
cat("\n--- Analysis 4: Mean pairwise DEG correlation ---\n")

cor_results <- rbindlist(lapply(0:4, function(s) {
  samp <- fib_samples[fibrosis_stage == s]$sample_id
  if (length(samp) < 10) return(NULL)

  mat <- logcpm[sig_in_matrix, samp, drop = FALSE]
  # Pairwise gene correlation
  cor_mat <- cor(t(mat), method = "pearson")
  # Extract upper triangle
  upper <- cor_mat[upper.tri(cor_mat)]

  data.table(
    fibrosis_stage = s,
    stage_label    = paste0("F", s),
    mean_abs_cor   = mean(abs(upper)),
    mean_cor       = mean(upper),
    median_abs_cor = median(abs(upper)),
    sd_abs_cor     = sd(abs(upper)),
    se_abs_cor     = sd(abs(upper)) / sqrt(length(upper)),
    n_pairs        = length(upper),
    n_samples      = length(samp)
  )
}))

cat("Mean pairwise DEG correlation per stage:\n")
print(cor_results)

# ============================================================
# Compile results CSV
# ============================================================
cat("\n--- Compiling results ---\n")

results <- merge(cv_results[, .(fibrosis_stage, stage_label, mean_cv, se_cv, n_samples)],
  pca_results[, .(fibrosis_stage, pc1_var_pct, pc1_pc2_ratio)],
  by = "fibrosis_stage", all = TRUE)
results <- merge(results,
  cor_results[, .(fibrosis_stage, mean_abs_cor, se_abs_cor)],
  by = "fibrosis_stage", all = TRUE)

out_csv <- file.path(RES_DIR, "critical_slowing_down.csv")
fwrite(results, out_csv)
cat(sprintf("Results saved: %s\n", out_csv))

# Also save JSD as separate companion
jsd_csv <- file.path(RES_DIR, "critical_slowing_down_jsd.csv")
fwrite(jsd_results, jsd_csv)
cat(sprintf("JSD results saved: %s\n", jsd_csv))

# ============================================================
# Figure: 4-panel critical slowing down
# ============================================================
cat("\n--- Generating figure ---\n")

# Panel (a): Mean CV per stage
pa <- ggplot(cv_results, aes(x = stage_label, y = mean_cv, group = 1)) +
  geom_line(color = masld_colors$fibrosis, linewidth = 0.6) +
  geom_point(aes(fill = stage_label), shape = 21, size = 2.5, stroke = 0.3) +
  geom_errorbar(aes(ymin = mean_cv - se_cv, ymax = mean_cv + se_cv),
    width = 0.15, linewidth = 0.3) +
  scale_fill_manual(values = fibrosis_stage_colors, guide = "none") +
  labs(x = "Fibrosis stage", y = "Mean CV (top 5K genes)",
    title = "Expression variance (flickering)") +
  theme_masld()

# Panel (b): JSD between adjacent stages
jsd_results[, transition_f := factor(transition,
  levels = c("F0->F1", "F1->F2", "F2->F3", "F3->F4"))]
# Color the maximum JSD bar
jsd_results[, is_max := mean_jsd == max(mean_jsd)]

pb <- ggplot(jsd_results, aes(x = transition_f, y = mean_jsd)) +
  geom_col(aes(fill = is_max), width = 0.6) +
  geom_errorbar(aes(ymin = mean_jsd - se_jsd, ymax = mean_jsd + se_jsd),
    width = 0.15, linewidth = 0.3) +
  scale_fill_manual(values = c("FALSE" = "#90CAF9", "TRUE" = masld_colors$fibrosis),
    guide = "none") +
  labs(x = "Transition", y = "Mean Jensen-Shannon divergence",
    title = "Transcriptomic divergence") +
  theme_masld()

# Panel (c): PC1 variance explained
pc <- ggplot(pca_results, aes(x = stage_label, y = pc1_var_pct, group = 1)) +
  geom_line(color = masld_colors$fibrosis, linewidth = 0.6) +
  geom_point(aes(fill = stage_label), shape = 21, size = 2.5, stroke = 0.3) +
  scale_fill_manual(values = fibrosis_stage_colors, guide = "none") +
  labs(x = "Fibrosis stage", y = "PC1 variance explained (%)",
    title = "Leading eigenvalue dominance") +
  theme_masld()

# Panel (d): Mean pairwise DEG correlation
pd <- ggplot(cor_results, aes(x = stage_label, y = mean_abs_cor, group = 1)) +
  geom_line(color = masld_colors$fibrosis, linewidth = 0.6) +
  geom_point(aes(fill = stage_label), shape = 21, size = 2.5, stroke = 0.3) +
  geom_errorbar(aes(ymin = mean_abs_cor - se_abs_cor, ymax = mean_abs_cor + se_abs_cor),
    width = 0.15, linewidth = 0.3) +
  scale_fill_manual(values = fibrosis_stage_colors, guide = "none") +
  labs(x = "Fibrosis stage", y = "Mean |r| (top DEGs)",
    title = "Gene-gene coordination") +
  theme_masld()

# Combine with patchwork
fig <- (pa | pb) / (pc | pd) +
  plot_annotation(
    tag_levels = "a",
    title = "Critical slowing down at the F2 tipping point",
    theme = theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0,
        family = "Helvetica")
    )
  )

out_pdf <- file.path(FIG_DIR, "figS_critical_slowing_down.pdf")
save_fig(fig, out_pdf, width = fig_full_width, height = 5)
cat(sprintf("Figure saved: %s\n", out_pdf))

cat("\nDone.\n")
