#!/usr/bin/env Rscript
# =============================================================================
# 218_switch_gene_classification.R
# Gap #14: Classify DEGs as switch-like (sigmoid at F2) vs gradual (linear)
#
# Logic:
#   For each of 5,484 DEGs, compute per-fibrosis-stage mean logCPM, then:
#   1. Fit linear model:  expression ~ fibrosis_stage (numeric 0-4)
#   2. Fit step model:    expression ~ I(fibrosis_stage >= 2) (binary F0-F1 vs F2-F4)
#   3. Compare via delta-AIC (AIC_linear - AIC_step); positive = step wins
#   4. Compute switch_ratio = |mean(F2-F4) - mean(F0-F1)| /
#        (|mean(F1)-mean(F0)| + |mean(F4)-mean(F3)| + epsilon)
#   5. Classify: switch_up, switch_down, gradual_up, gradual_down
#
# Outputs:
#   CSV:  RNA-seq/results/stratified_causal/switch_gene_classification.csv
#   Fig:  figures/supplementary/progression/figS_switch_classification.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(edgeR)
  library(patchwork)
  library(fgsea)
  library(msigdbr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUT_CSV <- file.path(BASE, "RNA-seq/results/stratified_causal/switch_gene_classification.csv")
OUT_FIG <- file.path(BASE, "figures/supplementary/figS02_progression/figS_switch_classification.pdf")

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_FIG), recursive = TRUE, showWarnings = FALSE)

cat("=== 218: Switch gene classification ===\n")

# ── 1. Load data ─────────────────────────────────────────────────────────────
cat("Loading merged DGE...\n")
dge <- readRDS(file.path(INT, "results/integration/merged_dge.rds"))
cat("  DGE:", nrow(dge), "genes x", ncol(dge), "samples\n")

# Load dream results for DEG list
dream <- fread(file.path(INT, "results/integration/dream_results_ashr.csv"))
degs <- dream[padj < 0.05 & abs(logFC) > 0.3]
cat("  DEGs:", nrow(degs), "\n")
deg_genes <- unique(degs$gene)

# Load fibrosis staging from modeling metadata
meta <- fread(file.path(INT, "results/staging_classifier/modeling_metadata.csv"),
              select = c("sample_id", "fibrosis_stage", "dataset"))
meta <- meta[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4]
meta[, fibrosis_stage := as.integer(fibrosis_stage)]
cat("  Samples with fibrosis stage:", nrow(meta), "\n")
cat("  Stage distribution:\n")
print(table(meta$fibrosis_stage))

# ── 2. Compute logCPM matrix (subset to DEGs + staged samples) ──────────────
shared_samples <- intersect(meta$sample_id, colnames(dge))
cat("  Shared samples:", length(shared_samples), "\n")
meta <- meta[sample_id %in% shared_samples]

shared_genes <- intersect(deg_genes, rownames(dge))
cat("  DEGs in DGE:", length(shared_genes), "\n")

dge_sub <- dge[shared_genes, shared_samples]
logcpm <- edgeR::cpm(dge_sub, log = TRUE, prior.count = 1)

# ── 3. Per-gene, per-stage mean logCPM ───────────────────────────────────────
cat("Computing per-stage means...\n")
stage_vec <- meta[match(colnames(logcpm), sample_id), fibrosis_stage]

# Pre-compute stage indices
stage_idx <- lapply(0:4, function(s) which(stage_vec == s))
names(stage_idx) <- paste0("F", 0:4)
stage_n <- sapply(stage_idx, length)
cat("  Per-stage N:", paste(names(stage_n), stage_n, sep = "=", collapse = ", "), "\n")

# Stage means matrix: genes x 5 stages
stage_means <- sapply(0:4, function(s) {
  rowMeans(logcpm[, stage_idx[[s + 1]], drop = FALSE])
})
colnames(stage_means) <- paste0("F", 0:4)

# ── 4. Fit linear vs step model per gene ─────────────────────────────────────
cat("Fitting linear vs step models for", nrow(logcpm), "genes...\n")

# Use sample-level data for proper AIC comparison
fib_numeric <- as.numeric(stage_vec)
fib_binary  <- as.integer(stage_vec >= 2)

results_list <- vector("list", nrow(logcpm))
gene_names <- rownames(logcpm)

for (i in seq_len(nrow(logcpm))) {
  y <- logcpm[i, ]

  # Linear model
  fit_lin <- lm(y ~ fib_numeric)
  aic_lin <- AIC(fit_lin)
  r2_lin  <- summary(fit_lin)$r.squared

  # Step model (binary: F0-F1 vs F2-F4)
  fit_step <- lm(y ~ fib_binary)
  aic_step <- AIC(fit_step)
  r2_step  <- summary(fit_step)$r.squared

  # Delta AIC: positive means step model wins
  delta_aic <- aic_lin - aic_step

  # Stage means
  m <- stage_means[i, ]
  mean_early <- mean(m[1:2])  # F0-F1
  mean_late  <- mean(m[3:5])  # F2-F4

  # Switch ratio
  epsilon <- 0.01
  jump <- abs(mean_late - mean_early)
  gradual_denom <- abs(m[2] - m[1]) + abs(m[5] - m[4]) + epsilon
  switch_ratio <- jump / gradual_denom

  # Direction
  direction <- ifelse(mean_late > mean_early, "up", "down")

  results_list[[i]] <- data.table(
    gene = gene_names[i],
    mean_F0 = m[1], mean_F1 = m[2], mean_F2 = m[3],
    mean_F3 = m[4], mean_F4 = m[5],
    mean_early = mean_early, mean_late = mean_late,
    delta_aic = delta_aic,
    r2_linear = r2_lin, r2_step = r2_step,
    switch_ratio = switch_ratio,
    direction = direction
  )

  if (i %% 1000 == 0) cat("  Processed", i, "/", nrow(logcpm), "genes\n")
}

res <- rbindlist(results_list)
cat("  Done. Results:", nrow(res), "genes\n")

# ── 5. Classify genes ────────────────────────────────────────────────────────
# Primary classification by delta_aic (positive = step wins)
# Secondary: switch_ratio for ranking within class
res[, classification := fifelse(
  delta_aic > 0 & direction == "up", "switch_up",
  fifelse(delta_aic > 0 & direction == "down", "switch_down",
  fifelse(delta_aic <= 0 & direction == "up", "gradual_up", "gradual_down"))
)]

# Merge dream statistics
res <- merge(res, degs[, .(gene, logFC, padj, t)], by = "gene", all.x = TRUE)

cat("\nClassification summary:\n")
print(table(res$classification))

cat("\nSwitch ratio summary:\n")
print(summary(res$switch_ratio))

cat("\nDelta AIC summary (positive = step wins):\n")
print(summary(res$delta_aic))

# Save CSV
fwrite(res[order(-switch_ratio)], OUT_CSV)
cat("Saved:", OUT_CSV, "\n")

# ── 6. Figures ───────────────────────────────────────────────────────────────
cat("\nGenerating figures...\n")

# Define class colors
class_colors <- c(
  switch_up    = masld_colors$up,
  switch_down  = masld_colors$down,
  gradual_up   = "#F48FB1",  # soft pink
  gradual_down = "#90CAF9"   # light blue
)

# Panel A: Histogram of switch ratio
p_a <- ggplot(res, aes(x = switch_ratio, fill = classification)) +
  geom_histogram(bins = 60, alpha = 0.8, color = NA) +
  scale_fill_manual(values = class_colors, name = "Class") +
  geom_vline(xintercept = median(res$switch_ratio), linetype = "dashed", color = "grey30") +
  labs(x = "Switch ratio", y = "Number of DEGs",
       title = "Distribution of switch ratio across DEGs") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 2.2,
           label = sprintf("Median = %.2f", median(res$switch_ratio))) +
  theme_masld() +
  theme(legend.position = c(0.8, 0.8))

# Panel B: Example gene trajectories
# Top 3 switch-like + bottom 3 gradual
top_switch <- res[classification %in% c("switch_up", "switch_down")][order(-switch_ratio)][1:3, gene]
top_gradual <- res[classification %in% c("gradual_up", "gradual_down")][order(switch_ratio)][1:3, gene]
example_genes <- c(top_switch, top_gradual)

# Build trajectory data
traj_dt <- rbindlist(lapply(example_genes, function(g) {
  m <- stage_means[g, ]
  cls <- res[gene == g, classification]
  data.table(gene = g, stage = 0:4, logCPM = m,
             class = cls,
             type = ifelse(grepl("switch", cls), "Switch-like", "Gradual"))
}))
traj_dt[, gene_label := paste0(gene, " (", type, ")")]
traj_dt[, gene_label := factor(gene_label, levels = unique(gene_label))]

p_b <- ggplot(traj_dt, aes(x = stage, y = logCPM, color = type, group = gene_label)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.2) +
  facet_wrap(~gene_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Switch-like" = masld_colors$up,
                                 "Gradual" = masld_colors$down)) +
  scale_x_continuous(breaks = 0:4, labels = paste0("F", 0:4)) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(x = "Fibrosis stage", y = "Mean logCPM", title = "Example gene trajectories") +
  theme_masld() +
  theme(legend.position = "none",
        strip.text = element_text(size = 5.5))

# Panel C: Pathway enrichment (switch vs gradual genes)
cat("Running pathway enrichment...\n")

# Rank genes by switch_ratio for fgsea
# Positive = switch-like, negative = gradual
ranked_vec <- res$switch_ratio
names(ranked_vec) <- res$gene

# Load Hallmark + KEGG gene sets
h_sets <- msigdbr(species = "Homo sapiens", collection = "H")
kegg_sets <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  error = function(e) msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
)
all_sets <- rbind(
  as.data.table(h_sets)[, .(gs_name, gene_symbol)],
  as.data.table(kegg_sets)[, .(gs_name, gene_symbol)]
)
pathway_list <- split(all_sets$gene_symbol, all_sets$gs_name)
# Filter to sets with >=15 genes in our data
pathway_list <- pathway_list[sapply(pathway_list, function(x) sum(x %in% res$gene) >= 15)]

fgsea_res <- fgsea(pathways = pathway_list, stats = ranked_vec,
                    minSize = 15, maxSize = 500, nPermSimple = 10000)
fgsea_res <- as.data.table(fgsea_res)
fgsea_sig <- fgsea_res[padj < 0.05][order(NES)]

cat("  Significant pathways (padj < 0.05):", nrow(fgsea_sig), "\n")

# Show top 10 enriched in switch-like (positive NES) and top 10 in gradual (negative NES)
top_switch_pw <- fgsea_sig[NES > 0][order(-NES)][1:min(10, sum(NES > 0))]
top_gradual_pw <- fgsea_sig[NES < 0][order(NES)][1:min(10, sum(NES < 0))]
pw_plot <- rbind(top_switch_pw, top_gradual_pw)

if (nrow(pw_plot) > 0) {
  pw_plot[, pathway_short := gsub("HALLMARK_|KEGG_", "", pathway)]
  pw_plot[, pathway_short := gsub("_", " ", pathway_short)]
  pw_plot[, pathway_short := factor(pathway_short, levels = pathway_short[order(NES)])]
  pw_plot[, enrichment := fifelse(NES > 0, "Switch-like", "Gradual")]

  p_c <- ggplot(pw_plot, aes(x = NES, y = pathway_short, fill = enrichment)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = c("Switch-like" = masld_colors$up,
                                  "Gradual" = masld_colors$down)) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    labs(x = "Normalized Enrichment Score",
         y = NULL,
         title = "Pathway enrichment: switch-like vs gradual DEGs") +
    theme_masld() +
    theme(legend.position = c(0.85, 0.15),
          axis.text.y = element_text(size = 5))
} else {
  p_c <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No significant pathways") +
    theme_masld()
}

# Panel D: Barplot of classification proportions
class_counts <- res[, .N, by = classification]
class_counts[, pct := N / sum(N) * 100]
class_counts[, classification := factor(classification,
  levels = c("switch_up", "switch_down", "gradual_up", "gradual_down"))]

p_d <- ggplot(class_counts, aes(x = classification, y = N, fill = classification)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
            vjust = -0.3, size = 2.2) +
  scale_fill_manual(values = class_colors) +
  labs(x = NULL, y = "Number of DEGs",
       title = "DEG classification by expression dynamics") +
  theme_masld() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  coord_cartesian(ylim = c(0, max(class_counts$N) * 1.15))

# ── 7. Assemble and save ─────────────────────────────────────────────────────
fig <- (p_a | p_d) / (p_b) / (p_c) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave(OUT_FIG, fig, width = 180, height = 240, units = "mm", device = cairo_pdf)
cat("Saved:", OUT_FIG, "\n")

# ── 8. Summary statistics ────────────────────────────────────────────────────
n_switch <- sum(grepl("switch", res$classification))
n_gradual <- sum(grepl("gradual", res$classification))
cat("\n=== SUMMARY ===\n")
cat(sprintf("Total DEGs classified: %d\n", nrow(res)))
cat(sprintf("Switch-like: %d (%.1f%%)\n", n_switch, 100 * n_switch / nrow(res)))
cat(sprintf("Gradual:     %d (%.1f%%)\n", n_gradual, 100 * n_gradual / nrow(res)))
cat(sprintf("Median switch ratio: %.2f\n", median(res$switch_ratio)))
cat(sprintf("Mean delta AIC (step - linear): %.2f\n", mean(res$delta_aic)))

# Top switch genes
cat("\nTop 15 switch-like genes (by switch_ratio):\n")
print(res[classification %in% c("switch_up", "switch_down")][order(-switch_ratio)][1:15,
  .(gene, classification, switch_ratio, delta_aic, logFC)])

cat("\nDone.\n")
