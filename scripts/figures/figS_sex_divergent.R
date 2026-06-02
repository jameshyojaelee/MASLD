# Supplementary Figure: Sex-Divergent Sensitivity Analysis
# Panels: (a) sex class counts, (b) |LFC diff| histogram, (c) top 20 genes

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

dir.create(file.path(FIGS03_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(FIGS03_DIR, "panels", "figS_sex_divergent.pdf")

sex_class <- load_sex_classification()
sex_perm  <- load_sex_divergent_perm()

# ---- Panel (a): Sex class counts ----
if (!is.null(sex_class) && "sex_class" %in% names(sex_class)) {
  class_counts <- sex_class[!sex_class %in% c("Not_significant", ""), .N, by = sex_class]
  # Script 26 v2 uses interaction-based names: Female_biased, Male_biased, Concordant, Divergent
  class_counts[, sex_class := factor(sex_class,
    levels = c("Female_biased", "Male_biased", "Concordant", "Divergent"))]
  class_counts <- class_counts[!is.na(sex_class)]

  p_a <- ggplot(class_counts, aes(x = sex_class, y = N, fill = sex_class)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = comma(N)), vjust = -0.3, size = 2.2) +
    scale_fill_manual(values = sex_class_colors, guide = "none") +
    scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Number of DEGs", title = "Stringent classification") +
    theme_masld()
} else {
  p_a <- placeholder("Panel a: sex classification not found")
}

# ---- Panel (b): |LFC diff| histogram ----
if (!is.null(sex_perm) && "abs_lfc_diff" %in% names(sex_perm)) {
  lfc_diff_col <- "abs_lfc_diff"
} else if (!is.null(sex_perm) && "lfc_diff_abs" %in% names(sex_perm)) {
  lfc_diff_col <- "lfc_diff_abs"
} else {
  lfc_diff_col <- NULL
}

if (!is.null(sex_perm) && !is.null(lfc_diff_col)) {
  p_b <- ggplot(sex_perm, aes(x = get(lfc_diff_col))) +
    geom_histogram(bins = 40, fill = masld_colors$up, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "gray40") +
    geom_vline(xintercept = 1.0, linetype = "dashed", linewidth = 0.3, color = masld_colors$up) +
    annotate("text", x = 0.55, y = Inf, vjust = 1.5, hjust = 0,
             label = "Permissive (0.5)", size = 2, color = "gray40") +
    annotate("text", x = 1.05, y = Inf, vjust = 1.5, hjust = 0,
             label = "Stringent (1.0)", size = 2, color = masld_colors$up) +
    labs(x = "|Male LFC - Female LFC|", y = "Genes",
         title = paste0("Permissive divergent (n=", nrow(sex_perm), ")")) +
    theme_masld()
} else {
  p_b <- placeholder("Panel b: LFC diff data not found")
}

# ---- Panel (c): Top 20 sex-divergent genes ----
if (!is.null(sex_perm) && !is.null(lfc_diff_col)) {
  top20 <- sex_perm[order(-get(lfc_diff_col))][1:min(20, .N)]
  top20[, symbol := factor(symbol, levels = rev(symbol))]

  p_c <- ggplot(top20) +
    geom_segment(aes(x = meta_logFC_F, xend = meta_logFC_M, y = symbol, yend = symbol),
                 color = "gray70", linewidth = 0.4) +
    geom_point(aes(x = meta_logFC_F, y = symbol), color = masld_colors$female, size = 1.5, shape = 16) +
    geom_point(aes(x = meta_logFC_M, y = symbol), color = masld_colors$male, size = 1.5, shape = 16) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "gray50") +
    labs(x = "Meta-analysis logFC", y = NULL,
         title = "Top 20 sex-divergent genes") +
    annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1,
             label = "Female", color = masld_colors$female, size = 2) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1,
             label = "Male", color = masld_colors$male, size = 2) +
    theme_masld() +
    theme(axis.text.y = element_text(face = "italic", size = 6))
} else {
  p_c <- placeholder("Panel c: sex-divergent genes not found")
}

# ---- Assemble ----
figS3 <- (p_a | p_b | p_c) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(figS3, OUT, height = 4)
message("FigS3 saved to ", OUT)
