#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure 3: Deconvolution Methods Comparison
# Comprehensive MuSiC + BayesPrism deconvolution panels
#
# Outputs individual PDFs to figures/supplementary/figS03_deconvolution/
#
#  1. music_all_celltypes.pdf     — MuSiC 16 cell types, Control vs Disease
#  2. bp_all_celltypes.pdf        — BayesPrism 16 cell types, Control vs Disease
#  3. method_comparison_scatter.pdf — MuSiC vs BayesPrism per-sample scatter
#  4. method_comparison_bland_altman.pdf — Bland-Altman agreement plots
#  5. stacked_composition.pdf     — Stacked bars: per-dataset, per-condition
#  6. per_dataset_boxplots.pdf    — Per-dataset hepatocyte + macrophage fractions
#  7. attribution_sankey.pdf      — Attribution class concordance (MuSiC vs BP)
#  8. hep_fraction_vs_disease.pdf — Hepatocyte fraction vs disease severity
# ==========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_DIR <- FIGS03_DIR

DECONV_DIR <- file.path(BASE, "Analysis/Deconvolution/results")
META_PATH  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                        "metadata/unified_metadata.csv")

# 9 cohorts in the cohort presentation; PRJNA512027 (Gerhard 2018) excluded
# for L0/S0 library-prep batch confound with diagnosis.
datasets_10 <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE167523",
                  "GSE174478", "GSE193066", "GSE213621", "GSE240729")

disease_fill <- c(Control = masld_colors$control, Disease = masld_colors$masld)

# ============================================================================
# Load both methods for the 9 cohorts in the cohort presentation
# ============================================================================
load_method <- function(method_suffix) {
  rbindlist(lapply(datasets_10, function(ds) {
    f <- file.path(DECONV_DIR, ds, paste0(ds, method_suffix))
    if (!file.exists(f)) return(NULL)
    raw <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    mat <- as.data.table(raw)
    mat[, sample_id := rownames(raw)]
    mat[, dataset := ds]
    mat
  }), fill = TRUE)
}

cat("Loading BayesPrism proportions...\n")
bp_all <- load_method("_bayesprism_proportions.tsv")
cat("  BayesPrism:", nrow(bp_all), "samples\n")

cat("Loading MuSiC proportions (weighted)...\n")
music_all <- load_method("_music_prop_weighted.tsv")
cat("  MuSiC:", nrow(music_all), "samples\n")

# Merge metadata
meta <- fread(META_PATH, select = c("sample_id", "group_binary", "dataset", "condition"))

bp_all    <- merge(bp_all, meta, by = c("sample_id", "dataset"), all.x = TRUE)
music_all <- merge(music_all, meta, by = c("sample_id", "dataset"), all.x = TRUE)
bp_all    <- bp_all[!is.na(group_binary)]
music_all <- music_all[!is.na(group_binary)]

# Cell type columns
meta_cols <- c("sample_id", "dataset", "group_binary", "condition")
ct_cols   <- setdiff(intersect(names(bp_all), names(music_all)), meta_cols)
cat("Cell types:", length(ct_cols), "\n")
cat(" ", paste(ct_cols, collapse = ", "), "\n\n")

# Order cell types by mean BayesPrism abundance
ct_means <- sapply(ct_cols, function(x) mean(bp_all[[x]], na.rm = TRUE))
ct_order_desc <- names(sort(ct_means, decreasing = TRUE))
ct_order_asc  <- rev(ct_order_desc)

# ============================================================================
# 1. MuSiC — all 16 cell types, Control vs Disease
# ============================================================================
cat("--- Panel 1: MuSiC all cell types ---\n")

music_long <- melt(music_all, id.vars = meta_cols, measure.vars = ct_cols,
                   variable.name = "cell_type", value.name = "proportion")
music_long[, cell_type := factor(cell_type, levels = ct_order_desc)]
music_long[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]

p1 <- ggplot(music_long, aes(x = group_binary, y = proportion, fill = group_binary)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.65, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 0.05, alpha = 0.1, color = "gray30") +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = disease_fill, guide = "none") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Estimated proportion",
       title = "MuSiC cell-type fractions across all 16 cell types") +
  theme_masld() +
  theme(strip.text = element_text(size = 6, face = "bold"),
        plot.title = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 5))

save_fig_tall(p1, file.path(OUT_DIR, "music_all_celltypes.pdf"),
              width = fig_full_width, height = 9)

# ============================================================================
# 2. BayesPrism — all 16 cell types, Control vs Disease
# ============================================================================
cat("--- Panel 2: BayesPrism all cell types ---\n")

bp_long <- melt(bp_all, id.vars = meta_cols, measure.vars = ct_cols,
                variable.name = "cell_type", value.name = "proportion")
bp_long[, cell_type := factor(cell_type, levels = ct_order_desc)]
bp_long[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]

p2 <- ggplot(bp_long, aes(x = group_binary, y = proportion, fill = group_binary)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.65, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 0.05, alpha = 0.1, color = "gray30") +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = disease_fill, guide = "none") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Estimated proportion",
       title = "BayesPrism cell-type fractions across all 16 cell types") +
  theme_masld() +
  theme(strip.text = element_text(size = 6, face = "bold"),
        plot.title = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 5))

save_fig_tall(p2, file.path(OUT_DIR, "bp_all_celltypes.pdf"),
              width = fig_full_width, height = 9)

# ============================================================================
# 3. Method comparison — MuSiC vs BayesPrism scatter (per cell type)
# ============================================================================
cat("--- Panel 3: Method comparison scatter ---\n")

# Merge both methods by sample_id
bp_slim <- melt(bp_all, id.vars = "sample_id", measure.vars = ct_cols,
                variable.name = "cell_type", value.name = "BayesPrism")
music_slim <- melt(music_all, id.vars = "sample_id", measure.vars = ct_cols,
                   variable.name = "cell_type", value.name = "MuSiC")
comparison <- merge(bp_slim, music_slim, by = c("sample_id", "cell_type"))

# Compute per-cell-type correlations
cor_dt <- comparison[, .(
  pearson = cor(BayesPrism, MuSiC, use = "complete.obs"),
  spearman = cor(BayesPrism, MuSiC, use = "complete.obs", method = "spearman"),
  n = .N
), by = cell_type]
cor_dt[, label := sprintf("r=%.2f\nrho=%.2f", pearson, spearman)]

comparison[, cell_type := factor(cell_type, levels = ct_order_desc)]
cor_dt[, cell_type := factor(cell_type, levels = ct_order_desc)]

p3 <- ggplot(comparison, aes(x = MuSiC, y = BayesPrism)) +
  geom_point(size = 0.1, alpha = 0.08, color = "gray30", shape = 16) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
  geom_text(data = cor_dt, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.3, size = 1.8, color = masld_colors$hep_intrinsic) +
  facet_wrap(~ cell_type, scales = "free", ncol = 4) +
  scale_x_continuous(labels = percent_format()) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "MuSiC proportion", y = "BayesPrism proportion",
       title = "Method concordance: MuSiC vs BayesPrism per cell type") +
  theme_masld() +
  theme(strip.text = element_text(size = 6, face = "bold"),
        plot.title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 4.5))

save_fig_tall(p3, file.path(OUT_DIR, "method_comparison_scatter.pdf"),
              width = fig_full_width, height = 9)

# ============================================================================
# 4. Bland-Altman agreement (key cell types)
# ============================================================================
cat("--- Panel 4: Bland-Altman ---\n")

key_cts <- c("Hepatocytes", "Macrophages", "Fibroblasts", "Endothelial cells",
             "T cells", "Cholangiocytes")
ba_dt <- comparison[cell_type %in% key_cts]
ba_dt[, mean_prop := (BayesPrism + MuSiC) / 2]
ba_dt[, diff_prop := BayesPrism - MuSiC]
ba_dt[, cell_type := factor(cell_type, levels = key_cts)]

# Per-facet limits
ba_stats <- ba_dt[, .(mean_diff = mean(diff_prop, na.rm = TRUE),
                       sd_diff = sd(diff_prop, na.rm = TRUE)),
                  by = cell_type]
ba_stats[, upper := mean_diff + 1.96 * sd_diff]
ba_stats[, lower := mean_diff - 1.96 * sd_diff]

p4 <- ggplot(ba_dt, aes(x = mean_prop, y = diff_prop)) +
  geom_point(size = 0.1, alpha = 0.1, color = "gray30", shape = 16) +
  geom_hline(data = ba_stats, aes(yintercept = mean_diff),
             linetype = "solid", color = masld_colors$hep_intrinsic, linewidth = 0.3) +
  geom_hline(data = ba_stats, aes(yintercept = upper),
             linetype = "dashed", color = "red", linewidth = 0.25) +
  geom_hline(data = ba_stats, aes(yintercept = lower),
             linetype = "dashed", color = "red", linewidth = 0.25) +
  facet_wrap(~ cell_type, scales = "free", ncol = 3) +
  scale_x_continuous(labels = percent_format()) +
  labs(x = "Mean proportion (MuSiC + BP) / 2",
       y = "Difference (BP - MuSiC)",
       title = "Bland-Altman agreement: BayesPrism vs MuSiC") +
  theme_masld() +
  theme(strip.text = element_text(size = 6.5, face = "bold"),
        plot.title = element_text(size = 8, face = "bold"))

save_fig(p4, file.path(OUT_DIR, "method_comparison_bland_altman.pdf"),
         width = fig_full_width, height = 5)

# ============================================================================
# 5. Stacked composition — per-dataset and per-condition
# ============================================================================
cat("--- Panel 5: Stacked composition ---\n")

# Build fill palette from ct_palette, fill missing
fill_vals <- ct_palette[ct_cols]
missing <- is.na(fill_vals)
if (any(missing)) {
  gray_seq <- colorRampPalette(c("#78909C", "#CFD8DC"))(sum(missing))
  fill_vals[missing] <- gray_seq
}
names(fill_vals) <- ct_cols

# 5a: Per-condition (BayesPrism)
bp_comp <- melt(bp_all, id.vars = meta_cols, measure.vars = ct_cols,
                variable.name = "cell_type", value.name = "proportion")
bp_comp_mean <- bp_comp[, .(mean_prop = mean(proportion, na.rm = TRUE)),
                        by = .(group_binary, cell_type)]
bp_comp_mean[, cell_type := factor(cell_type, levels = ct_order_asc)]
bp_comp_mean[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]

p5a <- ggplot(bp_comp_mean, aes(x = group_binary, y = mean_prop, fill = cell_type)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = fill_vals, name = "Cell type") +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Mean proportion", title = "BayesPrism: by condition") +
  theme_masld() +
  theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"), plot.title = element_text(size = 7, face = "bold")) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

# 5b: Same for MuSiC
music_comp <- melt(music_all, id.vars = meta_cols, measure.vars = ct_cols,
                   variable.name = "cell_type", value.name = "proportion")
music_comp_mean <- music_comp[, .(mean_prop = mean(proportion, na.rm = TRUE)),
                              by = .(group_binary, cell_type)]
music_comp_mean[, cell_type := factor(cell_type, levels = ct_order_asc)]
music_comp_mean[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]

p5b <- ggplot(music_comp_mean, aes(x = group_binary, y = mean_prop, fill = cell_type)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = fill_vals, name = "Cell type") +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Mean proportion", title = "MuSiC: by condition") +
  theme_masld() +
  theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"), plot.title = element_text(size = 7, face = "bold")) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

# 5c: Per-dataset mean (BayesPrism, both conditions)
bp_ds_comp <- bp_comp[, .(mean_prop = mean(proportion, na.rm = TRUE)),
                      by = .(dataset, group_binary, cell_type)]
bp_ds_comp[, cell_type := factor(cell_type, levels = ct_order_asc)]
bp_ds_comp[, ds_cond := paste0(dataset, "\n", group_binary)]
# Order: grouped by dataset, Control then Disease
ds_levels <- unlist(lapply(datasets_10, function(d) paste0(d, "\n", c("Control", "Disease"))))
bp_ds_comp[, ds_cond := factor(ds_cond, levels = ds_levels)]

p5c <- ggplot(bp_ds_comp, aes(x = ds_cond, y = mean_prop, fill = cell_type)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = fill_vals, name = "Cell type") +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Mean proportion",
       title = "BayesPrism composition per dataset and condition") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 3.5, angle = 90, hjust = 1, vjust = 0.5),
        legend.text = element_text(size = 5), legend.title = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"), plot.title = element_text(size = 7, face = "bold")) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

# Assemble stacked composition figure
p5_top <- p5a + p5b + plot_layout(guides = "collect") &
  theme(legend.position = "right")
p5_full <- p5_top / p5c + plot_layout(heights = c(1, 1.3)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig_tall(p5_full, file.path(OUT_DIR, "stacked_composition.pdf"),
              width = fig_full_width, height = 8)

# ============================================================================
# 6. Per-dataset boxplots — Hepatocyte + Macrophage, both methods
# ============================================================================
cat("--- Panel 6: Per-dataset boxplots ---\n")

# BayesPrism hepatocytes
bp_hep <- bp_all[, .(sample_id, dataset, group_binary, proportion = Hepatocytes, method = "BayesPrism")]
music_hep <- music_all[, .(sample_id, dataset, group_binary, proportion = Hepatocytes, method = "MuSiC")]
both_hep <- rbind(bp_hep, music_hep)
both_hep[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]
ds_order <- both_hep[method == "BayesPrism", .(med = median(proportion, na.rm = TRUE)), by = dataset][order(-med)]$dataset
both_hep[, dataset := factor(dataset, levels = ds_order)]

p6a <- ggplot(both_hep, aes(x = dataset, y = proportion, fill = group_binary)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.7, alpha = 0.8,
               position = position_dodge(0.8)) +
  facet_wrap(~ method, ncol = 1) +
  scale_fill_manual(values = disease_fill, name = "Group") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Hepatocyte proportion",
       title = "Hepatocyte fraction per dataset") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 8, face = "bold"),
        legend.position = "bottom")

# Macrophages
bp_mac <- bp_all[, .(sample_id, dataset, group_binary, proportion = Macrophages, method = "BayesPrism")]
music_mac <- music_all[, .(sample_id, dataset, group_binary, proportion = Macrophages, method = "MuSiC")]
both_mac <- rbind(bp_mac, music_mac)
both_mac[, group_binary := factor(group_binary, levels = c("Control", "Disease"))]
both_mac[, dataset := factor(dataset, levels = ds_order)]

p6b <- ggplot(both_mac, aes(x = dataset, y = proportion, fill = group_binary)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.7, alpha = 0.8,
               position = position_dodge(0.8)) +
  facet_wrap(~ method, ncol = 1) +
  scale_fill_manual(values = disease_fill, name = "Group") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Macrophage proportion",
       title = "Macrophage fraction per dataset") +
  theme_masld() +
  theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        strip.text = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 8, face = "bold"),
        legend.position = "bottom")

p6_full <- p6a + p6b + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"),
        legend.position = "bottom")

save_fig_tall(p6_full, file.path(OUT_DIR, "per_dataset_boxplots.pdf"),
              width = fig_full_width, height = 8)

# ============================================================================
# 7. Attribution class concordance — MuSiC vs BayesPrism
# ============================================================================
cat("--- Panel 7: Attribution concordance ---\n")

attr_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                       "results/deconvolution/bayesprism/attribution_cross_method_comparison.csv")
if (file.exists(attr_path)) {
  attr_dt <- fread(attr_path)
  # Standardize class names
  if ("bp_class" %in% names(attr_dt) && "music_class" %in% names(attr_dt)) {
    # Clean NS → Not_significant for consistency
    attr_dt[bp_class == "NS", bp_class := "Not_significant"]
    attr_dt[music_class == "NS", music_class := "Not_significant"]

    # Confusion matrix as heatmap
    conf <- attr_dt[, .N, by = .(music_class, bp_class)]

    # Order classes
    class_order <- c("Hepatocyte_intrinsic", "Composition_driven", "Unmasked", "Not_significant")
    conf[, music_class := factor(music_class, levels = class_order)]
    conf[, bp_class := factor(bp_class, levels = class_order)]
    conf <- conf[!is.na(music_class) & !is.na(bp_class)]

    # Percentage within MuSiC class
    conf[, total_music := sum(N), by = music_class]
    conf[, pct := N / total_music * 100]

    p7 <- ggplot(conf, aes(x = bp_class, y = music_class, fill = pct)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = paste0(N, "\n(", round(pct, 1), "%)")),
                size = 2.2, color = "black") +
      scale_fill_gradient(low = "white", high = masld_colors$hep_intrinsic,
                          name = "% of\nMuSiC class") +
      scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
      scale_y_discrete(labels = function(x) gsub("_", "\n", x)) +
      labs(x = "BayesPrism attribution class",
           y = "MuSiC attribution class",
           title = "Attribution class concordance: MuSiC vs BayesPrism") +
      theme_masld() +
      theme(plot.title = element_text(size = 8, face = "bold"),
            axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5),
            axis.text.y = element_text(size = 6),
            panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.3))

    save_fig(p7, file.path(OUT_DIR, "attribution_concordance.pdf"),
             width = fig_half_width + 1, height = 4.5)
  }
} else {
  cat("  attribution comparison file not found, skipping\n")
}

# ============================================================================
# 8. Hepatocyte fraction vs disease severity
# ============================================================================
cat("--- Panel 8: Hepatocyte fraction vs disease severity ---\n")

# Load full metadata for condition (includes NAFL, NASH, fibrosis stages etc.)
meta_full <- fread(META_PATH)
# Focus on datasets with graded disease (condition column)
bp_sev <- merge(bp_all[, .(sample_id, dataset, Hepatocytes, Macrophages)],
                meta_full[, .(sample_id, condition, group_binary)],
                by = "sample_id", all.x = TRUE)
bp_sev <- bp_sev[!is.na(condition)]

# Simplify conditions into ordered severity
bp_sev[, severity := fcase(
  condition %in% c("Control", "Healthy", "Normal"), "Control",
  condition %in% c("NAFL", "Steatosis", "MASL"), "MASL",
  condition %in% c("NASH", "MASH", "Steatohepatitis"), "MASH",
  condition %in% c("NASH_Fibrosis", "Advanced_Fibrosis", "Fibrosis",
                   "Fibrosis_F3", "Fibrosis_F4", "Cirrhosis"), "Advanced",
  default = NA_character_
)]
bp_sev <- bp_sev[!is.na(severity)]
bp_sev[, severity := factor(severity, levels = c("Control", "MASL", "MASH", "Advanced"))]

sev_colors <- c(Control = masld_colors$control, MASL = masld_colors$masl,
                MASH = masld_colors$mash, Advanced = masld_colors$fibrosis)

p8a <- ggplot(bp_sev, aes(x = severity, y = Hepatocytes, fill = severity)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.65, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 0.1, alpha = 0.15, color = "gray30") +
  scale_fill_manual(values = sev_colors, guide = "none") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Disease stage", y = "Hepatocyte proportion",
       title = "Hepatocyte fraction by disease severity") +
  theme_masld() +
  theme(plot.title = element_text(size = 7, face = "bold"))

p8b <- ggplot(bp_sev, aes(x = severity, y = Macrophages, fill = severity)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.65, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 0.1, alpha = 0.15, color = "gray30") +
  scale_fill_manual(values = sev_colors, guide = "none") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Disease stage", y = "Macrophage proportion",
       title = "Macrophage fraction by disease severity") +
  theme_masld() +
  theme(plot.title = element_text(size = 7, face = "bold"))

p8_full <- p8a + p8b +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(p8_full, file.path(OUT_DIR, "hep_fraction_vs_disease.pdf"),
         width = fig_full_width, height = 4)

cat("\n=== All deconvolution supplementary figures saved to:", OUT_DIR, "===\n")
