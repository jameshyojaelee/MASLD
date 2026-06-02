#!/usr/bin/env Rscript
# 210b_progression_coloc_plots.R
# Two focused plots for COLOC × fibrosis stage:
#   Plot 1: Per-stage COLOC enrichment (Fisher's exact per F0-F4 OVR contrast)
#   Plot 2: LogFC heatmap across stages for COLOC genes
# DEG THRESHOLD SYSTEM (Two-Tier):
#   Tier 1 (Primary): padj < 0.05, |logFC| > 0.5 — main dream DEGs (Script 05b)
#   Tier 2 (Progression): padj < 0.05, no LFC filter — binary and adjacent contrasts
#   Rationale: Binary contrasts pool multiple stages, diluting per-gene fold changes.
#              Adjacent transitions have lower N (200-400 vs 1,444), making LFC estimates
#              noisier. Cross-contrast comparisons use rank-based enrichment (fgsea)
#              to avoid confounding power with biology.
#   See: figures/supplementary/sensitivity/figS_deg_threshold_landscape.pdf
#
# Usage: Rscript 210b_progression_coloc_plots.R

library(data.table)
library(ggplot2)
library(patchwork)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

source(file.path(BASE, "scripts/figures/publication_theme.R"))

OUT_DIR <- file.path(BASE, "RNA-seq/results/stratified_causal")
FIG_DIR <- file.path(BASE, "figures/supplementary/figS04_coloc")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("210b: Progression × COLOC focused plots\n")
cat("============================================================\n\n")

# ── 1. Load data ─────────────────────────────────────────────────────────────
coloc <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
coloc <- coloc[gene != "" & !is.na(gene)]
cat("COLOC:", nrow(coloc), "genes,", sum(coloc$coloc_best_pp4 > 0.9), "with PP.H4 > 0.9\n")

fib_ovr <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier",
  "one_vs_rest_fibrosis_dream.csv"))
cat("Fibrosis OVR:", nrow(fib_ovr), "rows\n")

# Map Ensembl to symbol
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "human_symbol"))
atlas <- atlas[human_symbol != "" & !is.na(human_symbol)]
atlas[, ensembl_base := sub("\\.[0-9]+$", "", ensembl_id)]
symbol_map <- unique(atlas[, .(ensembl_base, symbol = human_symbol)])

fib_ovr[, ensembl_base := sub("\\.[0-9]+$", "", gene)]
fib_ovr <- merge(fib_ovr, symbol_map, by = "ensembl_base", all.x = TRUE)

# Merge COLOC
coloc_genes <- coloc[coloc_best_pp4 > 0.9, unique(gene)]
cat("COLOC genes (PP.H4 > 0.9):", length(coloc_genes), "\n")

fib_ovr[, is_coloc := symbol %in% coloc_genes]
cat("OVR rows with COLOC gene:", sum(fib_ovr$is_coloc, na.rm = TRUE), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 1: Per-stage COLOC enrichment (Fisher's exact)
# For each fibrosis stage OVR contrast, are COLOC genes over-represented
# among significant DEGs?
# ══════════════════════════════════════════════════════════════════════════════
cat("--- Plot 1: Per-stage enrichment ---\n")

stages <- c("F0_vs_rest", "F1_vs_rest", "F2_vs_rest", "F3_vs_rest", "F4_vs_rest")

# Universe: all genes tested in any OVR contrast
all_genes <- unique(fib_ovr[!is.na(symbol), symbol])
coloc_in_universe <- intersect(coloc_genes, all_genes)
cat("Universe:", length(all_genes), "genes;", length(coloc_in_universe), "COLOC genes\n")

enrich_list <- list()
for (stage in stages) {
  # DEGs at this stage (Tier 2: padj < 0.05, no LFC filter)
  stage_degs <- unique(fib_ovr[stage_label == stage & padj < 0.05 &
                                  !is.na(symbol), symbol])
  # 2x2 table
  a <- length(intersect(coloc_in_universe, stage_degs))
  b <- length(setdiff(coloc_in_universe, stage_degs))
  c <- length(setdiff(stage_degs, coloc_in_universe))
  d <- length(all_genes) - a - b - c

  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))

  enrich_list[[stage]] <- data.table(
    stage       = stage,
    stage_short = gsub("_vs_rest", "", stage),
    n_degs      = length(stage_degs),
    n_coloc_deg = a,
    OR          = ft$estimate,
    ci_lo       = ft$conf.int[1],
    ci_hi       = ft$conf.int[2],
    pvalue      = ft$p.value
  )
  cat(sprintf("  %s: %d DEGs, %d COLOC overlap, OR=%.2f, p=%.3g\n",
              stage, length(stage_degs), a, ft$estimate, ft$p.value))
}
enrich_dt <- rbindlist(enrich_list)
enrich_dt[, padj := p.adjust(pvalue, method = "BH")]
enrich_dt[, sig_label := fifelse(padj < 0.001, "***",
                          fifelse(padj < 0.01, "**",
                          fifelse(padj < 0.05, "*", "ns")))]
enrich_dt[, stage_short := factor(stage_short, levels = c("F0", "F1", "F2", "F3", "F4"))]

# Plot: forest-style OR + 95% CI per stage
plot1 <- ggplot(enrich_dt, aes(x = OR, y = stage_short)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, linewidth = 0.4,
                 color = "grey40") +
  geom_point(aes(size = n_coloc_deg, color = padj < 0.05), shape = 16) +
  geom_text(aes(label = paste0("n=", n_coloc_deg)), hjust = -0.3, vjust = -0.6,
            size = 2, color = "grey30") +
  geom_text(aes(label = sig_label, x = ci_hi + 0.05), hjust = 0, size = 2.5) +
  scale_size_continuous(range = c(2, 6), name = "COLOC genes\nin stage DEGs") +
  scale_color_manual(values = c("TRUE" = "#C2185B", "FALSE" = "#BDBDBD"), guide = "none") +
  labs(x = "Odds ratio (COLOC enrichment among stage DEGs)",
       y = "Fibrosis stage (one-vs-rest)",
       title = "COLOC gene enrichment by fibrosis stage",
       subtitle = "Fisher's exact: COLOC PP.H4 > 0.9 genes among stage DEGs (padj < 0.05)") +
  theme_masld(base_size = 7) +
  theme(plot.subtitle = element_text(size = 5.5, color = "grey40"))

save_fig(plot1, file.path(FIG_DIR, "figS_progression_stage_enrichment.pdf"),
         width = fig_half_width + 0.5, height = 2.8)
cat("Plot 1 saved\n")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 2: LogFC heatmap across stages for COLOC genes
# Show actual expression trajectory across F0-F4 for top COLOC genes
# ══════════════════════════════════════════════════════════════════════════════
cat("\n--- Plot 2: LogFC profiles across stages ---\n")

# Get logFC for each COLOC gene at each stage
coloc_profiles <- fib_ovr[is_coloc == TRUE & !is.na(symbol),
                           .(symbol, stage_label, logFC, padj)]

# Pivot to wide: gene × stage matrix
profiles_wide <- dcast(coloc_profiles, symbol ~ stage_label, value.var = "logFC",
                        fun.aggregate = mean)
padj_wide <- dcast(coloc_profiles, symbol ~ stage_label, value.var = "padj",
                    fun.aggregate = min)

# Merge COLOC PP.H4
profiles_wide <- merge(profiles_wide,
  coloc[, .(gene, coloc_best_pp4)],
  by.x = "symbol", by.y = "gene", all.x = TRUE)

# Identify peak stage for each gene
stage_cols <- c("F0_vs_rest", "F1_vs_rest", "F2_vs_rest", "F3_vs_rest", "F4_vs_rest")
profiles_mat <- as.matrix(profiles_wide[, ..stage_cols])
rownames(profiles_mat) <- profiles_wide$symbol

# Replace NA with 0
profiles_mat[is.na(profiles_mat)] <- 0

# Peak stage = stage with max |logFC|
profiles_wide[, peak_stage := stage_cols[apply(abs(profiles_mat), 1, which.max)]]
profiles_wide[, peak_logFC := apply(abs(profiles_mat), 1, max)]
profiles_wide[, peak_stage_short := gsub("_vs_rest", "", peak_stage)]

cat("COLOC genes with profiles:", nrow(profiles_wide), "\n")
cat("Peak stage distribution:\n")
print(table(profiles_wide$peak_stage_short))

# Select top 30 by PP.H4 for a readable heatmap
top30 <- head(profiles_wide[order(-coloc_best_pp4)], 30)

# Melt for ggplot
top30_long <- melt(top30, id.vars = c("symbol", "coloc_best_pp4", "peak_stage",
                                        "peak_logFC", "peak_stage_short"),
                   measure.vars = stage_cols,
                   variable.name = "stage", value.name = "logFC")
top30_long[, stage_short := gsub("_vs_rest", "", stage)]
top30_long[, stage_short := factor(stage_short, levels = c("F0", "F1", "F2", "F3", "F4"))]

# Add padj for significance marking
padj_long <- melt(padj_wide, id.vars = "symbol", measure.vars = stage_cols,
                   variable.name = "stage", value.name = "padj_val")
top30_long <- merge(top30_long, padj_long, by = c("symbol", "stage"), all.x = TRUE)
top30_long[, is_sig := padj_val < 0.05]

# Order genes: group by peak stage, then by PP.H4 within group
top30[, peak_order := match(peak_stage_short, c("F0", "F1", "F2", "F3", "F4"))]
top30 <- top30[order(peak_order, -coloc_best_pp4)]
gene_order <- top30$symbol
top30_long[, symbol := factor(symbol, levels = gene_order)]

# Add PP.H4 annotation for y-axis
pp4_labels <- top30[, .(symbol, pp4_label = paste0(symbol, "  (", round(coloc_best_pp4, 2), ")"))]
pp4_labels[, pp4_label := factor(pp4_label, levels = pp4_labels$pp4_label)]
top30_long <- merge(top30_long, pp4_labels, by = "symbol")

plot2 <- ggplot(top30_long, aes(x = stage_short, y = pp4_label, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(data = top30_long[is_sig == TRUE],
             aes(x = stage_short, y = pp4_label),
             shape = 8, size = 0.8, color = "black", stroke = 0.3) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C2185B",
                       midpoint = 0, name = "OVR logFC",
                       limits = c(-max(abs(top30_long$logFC), na.rm = TRUE),
                                   max(abs(top30_long$logFC), na.rm = TRUE))) +
  labs(x = "Fibrosis stage (one-vs-rest contrast)",
       y = NULL,
       title = "Expression profiles of top COLOC genes across fibrosis stages",
       subtitle = "PP.H4 > 0.9; * = padj < 0.05; genes grouped by peak stage") +
  theme_masld(base_size = 6) +
  theme(plot.subtitle = element_text(size = 5, color = "grey40"),
        axis.text.y = element_text(size = 5))

save_fig(plot2, file.path(FIG_DIR, "figS_progression_logfc_profiles.pdf"),
         width = fig_half_width + 0.5, height = 5.5)
cat("Plot 2 saved\n")

# ══════════════════════════════════════════════════════════════════════════════
# Combined
# ══════════════════════════════════════════════════════════════════════════════
combined <- plot1 / plot2 +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = "A")

save_fig(combined, file.path(FIG_DIR, "figS_progression_coloc_combined.pdf"),
         width = fig_half_width + 0.5, height = 8)
cat("Combined figure saved\n")

# Save enrichment table
fwrite(enrich_dt, file.path(OUT_DIR, "progression_stage_enrichment_v2.csv"))
cat("\nDone.\n")
