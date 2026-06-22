# Supplementary Figure: ComBat-seq Sensitivity Analysis
# Panels: (a) LFC scatter primary vs ComBat, (b) DEG overlap, (c) Direction concordance

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dir.create(file.path(FIGS01_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(FIGS_SENS_DIR, "figS_combat_sensitivity.pdf")

dream_primary <- load_dream_results()
dream_combat <- load_combat_dream()
concordance <- load_combat_concordance()

# ---- Panel (a): LFC scatter primary vs ComBat (rasterized) ----
if (!is.null(dream_primary) && !is.null(dream_combat)) {
  merged <- merge(
    dream_primary[, .(gene, lfc_primary = bulk_logFC, padj_primary = bulk_padj, symbol)],
    dream_combat[, .(gene, lfc_combat = logFC, padj_combat = padj)],
    by = "gene"
  )
  cor_val <- cor(merged$lfc_primary, merged$lfc_combat, use = "complete.obs")

  p_a <- ggplot(merged, aes(x = lfc_primary, y = lfc_combat)) +
    rasterize_layer(geom_point(size = 0.2, alpha = 0.3, color = masld_colors$ns, shape = 16)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, color = "gray50") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = masld_colors$up) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
             label = paste0("r = ", round(cor_val, 3)), size = 2.5, color = masld_colors$up) +
    labs(x = "Primary meta-analysis logFC", y = "ComBat-seq meta-analysis logFC") +
    theme_masld()
} else {
  p_a <- placeholder("Panel a: dream or ComBat data not found")
}

# ---- Panel (b): DEG overlap Venn-style bar ----
if (!is.null(concordance)) {
  # concordance is key-value: metric, value
  conc_wide <- dcast(concordance, . ~ metric, value.var = "value")
  sig_primary <- as.numeric(conc_wide$sig_primary)
  sig_combat  <- as.numeric(conc_wide$sig_combat)
  sig_overlap <- as.numeric(conc_wide$sig_overlap)
  jaccard     <- as.numeric(conc_wide$jaccard)

  overlap_dt <- data.table(
    Category = c("Primary only", "Overlap", "ComBat only"),
    Count = c(sig_primary - sig_overlap, sig_overlap, sig_combat - sig_overlap)
  )
  overlap_dt[, Category := factor(Category, levels = c("Primary only", "Overlap", "ComBat only"))]

  p_b <- ggplot(overlap_dt, aes(x = Category, y = Count, fill = Category)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = comma(Count)), vjust = -0.3, size = 2.2) +
    scale_fill_manual(values = c("Primary only" = masld_colors$down,
                                  "Overlap" = masld_colors$up,
                                  "ComBat only" = "#F48FB1"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    annotate("text", x = 2, y = max(overlap_dt$Count) * 0.9,
             label = paste0("Jaccard = ", round(jaccard, 3)), size = 2.5) +
    labs(x = NULL, y = "DEG count") +
    theme_masld()
} else if (!is.null(dream_primary) && !is.null(dream_combat)) {
  # Compute from data
  sig_p <- dream_primary[bulk_padj < 0.1]$gene
  sig_c <- dream_combat[padj < 0.1]$gene
  overlap_genes <- intersect(sig_p, sig_c)
  overlap_dt <- data.table(
    Category = c("Primary only", "Overlap", "ComBat only"),
    Count = c(length(setdiff(sig_p, sig_c)), length(overlap_genes), length(setdiff(sig_c, sig_p)))
  )
  overlap_dt[, Category := factor(Category, levels = Category)]
  jac <- length(overlap_genes) / length(union(sig_p, sig_c))

  p_b <- ggplot(overlap_dt, aes(x = Category, y = Count, fill = Category)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = comma(Count)), vjust = -0.3, size = 2.2) +
    scale_fill_manual(values = c("Primary only" = masld_colors$down,
                                  "Overlap" = masld_colors$up,
                                  "ComBat only" = "#F48FB1"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    annotate("text", x = 2, y = max(overlap_dt$Count) * 0.9,
             label = paste0("Jaccard = ", round(jac, 3)), size = 2.5) +
    labs(x = NULL, y = "DEG count") +
    theme_masld()
} else {
  p_b <- placeholder("Panel b: concordance data not found")
}

# ---- Panel (c): Direction concordance summary ----
if (!is.null(concordance)) {
  conc_wide <- dcast(concordance, . ~ metric, value.var = "value")
  dir_conc <- as.numeric(conc_wide$direction_concordance)
  lfc_r    <- as.numeric(conc_wide$lfc_pearson_r)
  jac      <- as.numeric(conc_wide$jaccard)

  summary_dt <- data.table(
    Metric = c("LFC Pearson r", "Direction concordance", "Jaccard index"),
    Value = c(lfc_r, dir_conc, jac)
  )
  summary_dt[, Metric := factor(Metric, levels = rev(Metric))]

  p_c <- ggplot(summary_dt, aes(x = Value, y = Metric)) +
    geom_col(fill = masld_colors$up, width = 0.5) +
    geom_text(aes(label = sprintf("%.3f", Value)), hjust = -0.1, size = 2.5) +
    scale_x_continuous(limits = c(0, 1.15), expand = c(0, 0)) +
    labs(x = "Value", y = NULL, title = "Concordance metrics") +
    theme_masld()
} else {
  p_c <- placeholder("Panel c: concordance metrics not found")
}

# ---- Assemble ----
figS1 <- (p_a | p_b | p_c) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(figS1, OUT, height = 3.5)
message("FigS1 saved to ", OUT)
