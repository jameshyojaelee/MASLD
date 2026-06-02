#!/usr/bin/env Rscript
##############################################################################
# Supplementary Figure: Interaction-based sex classification stability
#
# Validates the interaction-based sex classification (Script 26 v2) via
# 30-iteration 70% subsampling. Classes: Female_biased, Male_biased,
# Divergent, Concordant. Key result: 378 robust genes (selection_prob >= 0.60)
# with Spearman rho = 0.89 and 100% direction consistency.
#
# Six panels:
#   (a) Selection probability histogram by sex class
#   (b) Full vs subsample interaction logFC correlation
#   (c) Robustness tier stacked bar by sex class
#   (d) Class switching alluvial (full -> modal subsample class)
#   (e) Selection probability vs |interaction logFC|
#   (f) Forest plot: top 30 robust genes with subsampling intervals
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ============================================================================
# Paths
# ============================================================================
SUBSAMP_DIR <- file.path(BASE,
  "RNA-seq/results/audit_sensitivity/sex_interaction_subsampling")
OUT_DIR   <- FIGS_SENS_DIR
PANEL_DIR <- file.path(FIGS_SENS_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load data
# ============================================================================
cat("Loading interaction subsampling results from:", SUBSAMP_DIR, "\n")

per_gene    <- fread(file.path(SUBSAMP_DIR,
                     "sex_interaction_subsampling_per_gene.csv"))
iter_metrics <- fread(file.path(SUBSAMP_DIR,
                     "sex_interaction_subsampling_iter_metrics.csv"))
frac_summary <- fread(file.path(SUBSAMP_DIR,
                     "sex_interaction_subsampling_fraction_summary.csv"))
class_switch <- fread(file.path(SUBSAMP_DIR,
                     "sex_interaction_subsampling_class_switching.csv"))

# Add gene symbols
per_gene <- add_symbols(per_gene, "gene")

cat("Loaded:", nrow(per_gene), "genes,", nrow(iter_metrics), "iterations\n")

# ============================================================================
# Subset to dimorphic genes (full-model classification != Concordant)
# ============================================================================
dimorphic <- per_gene[full_dimorphic == TRUE]
cat("Dimorphic genes:", nrow(dimorphic), "\n")

# ============================================================================
# Color definitions
# ============================================================================
# Use the interaction-based sex class colors from publication_theme.R
class_colors <- c(
  Female_biased = sex_class_colors[["Female_biased"]],
  Male_biased   = sex_class_colors[["Male_biased"]],
  Divergent     = sex_class_colors[["Divergent"]],
  Concordant    = sex_class_colors[["Concordant"]]
)

robust_colors <- c(
  Robust   = "#004D40",
  Stable   = "#00897B",
  Moderate = "#FFA000",
  Fragile  = "#E53935"
)

# ============================================================================
# Panel (a): Selection probability histogram by sex class
# ============================================================================
cat("Panel (a): Selection probability distribution\n")

n_robust <- nrow(dimorphic[selection_prob >= 0.60])
n_robust_by_class <- dimorphic[selection_prob >= 0.60, .N, by = full_class]

p_a <- ggplot(dimorphic, aes(x = selection_prob, fill = full_class)) +
  geom_histogram(binwidth = 0.05, boundary = 0, color = "white",
                 linewidth = 0.2, position = "stack") +
  geom_vline(xintercept = 0.60, linetype = "dashed", color = "grey30",
             linewidth = 0.5) +
  annotate("text", x = 0.62, y = Inf, vjust = 1.5, hjust = 0,
           label = paste0("Robust >= 0.60\n", n_robust, " genes"),
           size = 2, color = "grey30") +
  scale_fill_manual(values = class_colors, name = "Sex class") +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     labels = label_percent(accuracy = 1)) +
  labs(x = "Selection probability (30 iterations)",
       y = "Number of genes",
       title = "Selection probability distribution") +
  theme_masld() +
  theme(legend.position = c(0.75, 0.75),
        legend.background = element_blank(),
        legend.key.size = unit(0.25, "cm"))

# ============================================================================
# Panel (b): Interaction logFC correlation (full vs subsample median)
# ============================================================================
cat("Panel (b): Interaction logFC correlation\n")

# Use all genes for the correlation, color dimorphic ones
per_gene[, plot_class := fifelse(full_dimorphic == TRUE, full_class, "Concordant")]

# Compute Spearman rho on dimorphic genes
rho_val <- cor(dimorphic$full_int_logFC, dimorphic$median_sub_logFC,
               method = "spearman", use = "complete.obs")

# Build scatter — Concordant genes as background, dimorphic on top
bg_dt  <- per_gene[full_dimorphic == FALSE]
fg_dt  <- per_gene[full_dimorphic == TRUE]

p_b <- ggplot() +
  geom_point(data = bg_dt,
             aes(x = full_int_logFC, y = median_sub_logFC),
             color = class_colors[["Concordant"]], alpha = 0.15, size = 0.3) +
  geom_point(data = fg_dt,
             aes(x = full_int_logFC, y = median_sub_logFC, color = full_class),
             alpha = 0.6, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.4) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = paste0("Spearman \u03c1 = ", sprintf("%.2f", rho_val)),
           size = 2.5) +
  scale_color_manual(values = class_colors, name = "Sex class") +
  labs(x = "Full-model interaction logFC",
       y = "Median subsample interaction logFC",
       title = "Interaction logFC concordance") +
  coord_equal() +
  theme_masld() +
  theme(legend.position = c(0.80, 0.20),
        legend.background = element_blank(),
        legend.key.size = unit(0.25, "cm"))

# ============================================================================
# Panel (c): Robustness tier stacked bar by sex class
# ============================================================================
cat("Panel (c): Robustness tier stacked bar\n")

# Order tiers
tier_order <- c("Robust", "Stable", "Moderate", "Fragile")
dimorphic[, robustness := factor(robustness, levels = tier_order)]

tier_counts <- dimorphic[, .N, by = .(full_class, robustness)]
tier_counts[, pct := N / sum(N) * 100, by = full_class]

# Order x-axis
class_order <- c("Female_biased", "Male_biased", "Divergent")
tier_counts[, full_class := factor(full_class, levels = class_order)]

p_c <- ggplot(tier_counts, aes(x = full_class, y = pct, fill = robustness)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(pct >= 5, paste0(round(pct), "%"), "")),
            position = position_stack(vjust = 0.5), size = 2, color = "white") +
  scale_fill_manual(values = robust_colors, name = "Tier",
                    drop = FALSE) +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Percentage of genes",
       title = "Robustness tier by sex class") +
  theme_masld() +
  theme(legend.position = "right")

# ============================================================================
# Panel (d): Class switching alluvial / heatmap
# ============================================================================
cat("Panel (d): Class switching\n")

# Filter to dimorphic genes in full model only (from_full != Concordant)
cs_dimorphic <- class_switch[from_full != "Concordant"]

# Try ggalluvial first; fall back to heatmap
has_alluvial <- requireNamespace("ggalluvial", quietly = TRUE)

if (has_alluvial) {
  library(ggalluvial)

  # Order factors
  cs_dimorphic[, from_full := factor(from_full,
    levels = c("Female_biased", "Male_biased", "Divergent"))]
  cs_dimorphic[, to_subsample_mode := factor(to_subsample_mode,
    levels = c("Female_biased", "Male_biased", "Divergent", "Concordant"))]

  # All colors needed for both axes
  alluvial_colors <- c(class_colors,
    Concordant = sex_class_colors[["Concordant"]])

  p_d <- ggplot(cs_dimorphic[count > 0],
                aes(y = count, axis1 = from_full, axis2 = to_subsample_mode)) +
    geom_alluvium(aes(fill = from_full), width = 1/6, alpha = 0.65,
                  decreasing = FALSE) +
    geom_stratum(width = 1/4, fill = "grey95", color = "grey40",
                 linewidth = 0.3) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 1.8) +
    scale_x_discrete(limits = c("Full model", "Modal subsample"),
                     expand = c(0.15, 0.05)) +
    scale_fill_manual(values = alluvial_colors, guide = "none") +
    labs(y = "Number of genes",
         title = "Class switching (full -> subsample)") +
    theme_masld() +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank())

} else {
  # Fallback: heatmap of transition matrix
  cs_dimorphic[, from_full := factor(from_full,
    levels = c("Female_biased", "Male_biased", "Divergent"))]
  cs_dimorphic[, to_subsample_mode := factor(to_subsample_mode,
    levels = c("Female_biased", "Male_biased", "Divergent", "Concordant"))]

  p_d <- ggplot(cs_dimorphic, aes(x = to_subsample_mode, y = from_full,
                                   fill = count)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = count), size = 2.5) +
    scale_fill_gradient(low = "white", high = "#C2185B", name = "Genes") +
    scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
    scale_y_discrete(labels = function(x) gsub("_", "\n", x)) +
    labs(x = "Modal subsample class", y = "Full-model class",
         title = "Class switching matrix") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ============================================================================
# Panel (e): Selection probability vs |interaction logFC|
# ============================================================================
cat("Panel (e): Selection probability vs effect size\n")

dimorphic[, abs_int_logFC := abs(full_int_logFC)]

# Spearman correlation between |logFC| and selection_prob
rho_eff <- cor(dimorphic$abs_int_logFC, dimorphic$selection_prob,
               method = "spearman", use = "complete.obs")

p_e <- ggplot(dimorphic, aes(x = abs_int_logFC, y = selection_prob,
                              color = full_class)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE,
              color = "grey30", fill = "grey80", linewidth = 0.6,
              alpha = 0.3) +
  geom_hline(yintercept = 0.60, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Spearman \u03c1 = ", sprintf("%.2f", rho_eff)),
           size = 2.5) +
  scale_color_manual(values = class_colors, name = "Sex class") +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(x = "|Interaction logFC| (full model)",
       y = "Selection probability",
       title = "Effect size predicts robustness") +
  theme_masld() +
  theme(legend.position = c(0.75, 0.75),
        legend.background = element_blank(),
        legend.key.size = unit(0.25, "cm"))

# ============================================================================
# Panel (f): Forest plot — top 30 robust genes
# ============================================================================
cat("Panel (f): Forest plot of top robust genes\n")

# Top 30 by selection_prob among dimorphic, breaking ties by |logFC|
top30 <- dimorphic[order(-selection_prob, -abs_int_logFC)][1:30]

# Order by full_int_logFC for visual clarity
top30[, symbol := factor(symbol, levels = rev(top30[order(full_int_logFC)]$symbol))]

p_f <- ggplot(top30, aes(x = full_int_logFC, y = symbol, color = full_class)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey80",
             linewidth = 0.3) +
  geom_errorbar(aes(xmin = subsamp_lower, xmax = subsamp_upper),
                width = 0.3, linewidth = 0.4, orientation = "y") +
  geom_point(size = 1.5) +
  scale_color_manual(values = class_colors, name = "Sex class") +
  labs(x = "Interaction logFC",
       y = NULL,
       title = "Top 30 robust genes") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 5),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"))

# ============================================================================
# Assemble composite figure
# ============================================================================
cat("Assembling composite figure\n")

composite <- (p_a | p_b | p_c) /
             (p_d | p_e | p_f) +
  plot_annotation(
    tag_levels = "a",
    title = "Interaction-based sex classification: subsampling stability",
    subtitle = paste0("3,080 dimorphic genes | 378 robust (selection prob >= 60%) | ",
                      "Spearman \u03c1 = 0.89 | Direction consistency = 100%"),
    theme = theme(
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7, color = "grey40")
    )
  ) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

# ============================================================================
# Save outputs
# ============================================================================
cat("Saving figures\n")

# Combined figure
save_fig_tall(composite,
              file.path(OUT_DIR, "figS_sex_interaction_subsampling.pdf"),
              width = fig_full_width, height = 7.5)

# Individual panels
save_fig(p_a, file.path(PANEL_DIR, "panel_a_selection_prob_hist.pdf"),
         width = fig_half_width, height = 3)
save_fig(p_b, file.path(PANEL_DIR, "panel_b_logfc_correlation.pdf"),
         width = fig_half_width, height = 3.5)
save_fig(p_c, file.path(PANEL_DIR, "panel_c_robustness_tiers.pdf"),
         width = fig_half_width, height = 3)
save_fig(p_d, file.path(PANEL_DIR, "panel_d_class_switching.pdf"),
         width = fig_half_width, height = 3)
save_fig(p_e, file.path(PANEL_DIR, "panel_e_effect_size_robustness.pdf"),
         width = fig_half_width, height = 3)
save_fig(p_f, file.path(PANEL_DIR, "panel_f_forest_top30.pdf"),
         width = fig_half_width, height = 4)

cat("Done. Composite saved to:", file.path(OUT_DIR,
    "figS_sex_interaction_subsampling.pdf"), "\n")
