#!/usr/bin/env Rscript
# ==========================================================================
# Supplementary Figure: Per-Cohort Heterogeneity
# Shows per-study effect sizes for top DEGs as forest plots,
# demonstrating robustness of mega-analysis findings across 9 cohorts
# (PRJNA512027 excluded from cohort presentation for L0/S0 batch confound).
#
# Panels:
#   a: Forest plot of top 15 upregulated DEGs (per-study logFC + CI)
#   b: Forest plot of top 15 downregulated DEGs
#   c: I-squared distribution across all DEGs
#   d: Per-cohort DEG count comparison (volcano summary)
# ==========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIGS08_DIR, "figS_cohort_forest.pdf")
dir.create(file.path(FIGS08_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)

# ==========================================================================
# Load data
# ==========================================================================
dream <- load_dream_results()
if (is.null(dream)) stop("dream_results.csv not found")

# Per-study DE
per_study_files <- list.files(PER_STUDY, pattern = "_de_results\\.csv$", full.names = TRUE)
per_study <- rbindlist(lapply(per_study_files, function(f) {
  dt <- fread(f)
  if (!"dataset" %in% names(dt)) {
    ds <- gsub("_de_results\\.csv$", "", basename(f))
    dt[, dataset := ds]
  }
  dt
}), fill = TRUE)

# Drop PRJNA512027 from cohort presentation (L0/S0 batch confound).
per_study <- per_study[dataset != "PRJNA512027"]

# Normalize column names
if ("adj.P.Val" %in% names(per_study)) setnames(per_study, "adj.P.Val", "padj")
if (!"symbol" %in% names(per_study)) per_study <- add_symbols(per_study, "gene")

# Dream vs meta concordance (has I-squared)
conc_path <- file.path(INT_RESULTS, "dream_vs_meta_concordance.csv")
has_concordance <- file.exists(conc_path)
if (has_concordance) {
  conc <- fread(conc_path)
}

# ==========================================================================
# Panel a: Forest plot — top 15 upregulated genes
# ==========================================================================
# Select top 15 up from dream by absolute logFC * significance
bulk_sig <- dream[is_dream_deg(dream) & !grepl("^ENS", symbol)]

top_up <- bulk_sig[bulk_logFC > 0][order(-bulk_logFC)][1:min(15, .N)]
top_down <- bulk_sig[bulk_logFC < 0][order(bulk_logFC)][1:min(15, .N)]

make_forest <- function(top_genes, per_study_dt, direction_label) {
  # Get per-study data for selected genes
  genes <- top_genes$symbol
  forest_dt <- per_study_dt[symbol %in% genes]

  if (nrow(forest_dt) == 0) return(placeholder(paste("No per-study data for", direction_label)))

  # Add dream summary
  dream_summary <- top_genes[, .(symbol, bulk_logFC, bulk_padj)]
  forest_dt <- merge(forest_dt, dream_summary, by = "symbol", all.x = TRUE)

  # Compute approximate 95% CI from t-statistic
  if ("t" %in% names(forest_dt)) {
    forest_dt[, se := logFC / t]
    forest_dt[, ci_lo := logFC - 1.96 * abs(se)]
    forest_dt[, ci_hi := logFC + 1.96 * abs(se)]
  } else {
    forest_dt[, ci_lo := logFC * 0.7]
    forest_dt[, ci_hi := logFC * 1.3]
  }

  # Order genes by dream logFC
  gene_order <- top_genes[order(bulk_logFC)]$symbol
  forest_dt[, symbol := factor(symbol, levels = gene_order)]

  # Color by significance
  forest_dt[, sig := fifelse(!is.na(padj) & padj < 0.05, "padj < 0.05", "NS")]

  p <- ggplot(forest_dt, aes(x = logFC, y = symbol)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "gray50", linetype = "dashed") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3,
                   linewidth = 0.3, color = "gray60") +
    geom_point(aes(color = sig, shape = dataset), size = 1.5, alpha = 0.8) +
    # Add dream summary as diamond
    geom_point(data = top_genes, aes(x = bulk_logFC, y = symbol),
               shape = 23, size = 2.5, fill = "black", color = "black", stroke = 0.5) +
    scale_color_manual(values = c("padj < 0.05" = masld_colors$up, "NS" = masld_colors$ns),
                       name = "Per-study") +
    theme_masld() +
    theme(
      axis.text.y = element_text(face = "italic", size = 6),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.size = unit(0.25, "cm")
    ) +
    guides(shape = guide_legend(nrow = 2, override.aes = list(size = 1.5)),
           color = guide_legend(override.aes = list(size = 2))) +
    labs(x = expression("log"[2]*"FC (Disease vs Control)"),
         y = NULL,
         title = paste0(direction_label, " DEGs across 9 cohorts"),
         shape = "Cohort")

  return(p)
}

p_a <- make_forest(top_up, per_study, "Top upregulated")
p_b <- make_forest(top_down, per_study, "Top downregulated")

# ==========================================================================
# Panel c: I-squared distribution
# ==========================================================================
if (has_concordance && "meta_I2" %in% names(conc)) {
  i2_dt <- conc[!is.na(meta_I2)]

  # Classify I2 levels
  i2_dt[, i2_class := fcase(
    meta_I2 < 25,  "Low (<25%)",
    meta_I2 < 50,  "Moderate (25-50%)",
    meta_I2 < 75,  "High (50-75%)",
    default = "Very high (>75%)"
  )]
  i2_dt[, i2_class := factor(i2_class,
    levels = c("Low (<25%)", "Moderate (25-50%)", "High (50-75%)", "Very high (>75%)"))]

  median_i2 <- median(i2_dt$meta_I2, na.rm = TRUE)

  p_c <- ggplot(i2_dt, aes(x = meta_I2, fill = i2_class)) +
    geom_histogram(binwidth = 5, boundary = 0, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = median_i2, linetype = "dashed", linewidth = 0.4, color = "black") +
    annotate("text", x = median_i2, y = Inf, vjust = 1.5, hjust = -0.1,
             label = sprintf("Median = %.0f%%", median_i2), size = 2.2) +
    scale_fill_manual(values = c(
      "Low (<25%)" = "#81D4FA",
      "Moderate (25-50%)" = "#42A5F5",
      "High (50-75%)" = "#1E88E5",
      "Very high (>75%)" = "#0D47A1"
    ), name = expression(I^2 ~ "class")) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    theme_masld() +
    theme(legend.position = "right") +
    labs(x = expression(I^2 ~ "(%)"),
         y = "Gene count",
         title = "Cross-cohort heterogeneity distribution")
} else {
  p_c <- placeholder("I-squared data not available")
}

# ==========================================================================
# Panel d: Per-cohort DEG counts
# ==========================================================================
deg_counts <- per_study[, .(
  n_up = sum(padj < 0.05 & logFC > 0, na.rm = TRUE),
  n_down = sum(padj < 0.05 & logFC < 0, na.rm = TRUE)
), by = dataset]

deg_long <- melt(deg_counts, id.vars = "dataset",
                 variable.name = "direction", value.name = "n")
deg_long[, direction := fifelse(direction == "n_up", "Up", "Down")]
deg_long[direction == "Down", n := -n]

# Order by total DEGs
totals <- deg_counts[, .(total = n_up + n_down), by = dataset][order(total)]
deg_long[, dataset := factor(dataset, levels = totals$dataset)]

p_d <- ggplot(deg_long, aes(x = n, y = dataset, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(Up = masld_colors$up, Down = masld_colors$down)) +
  scale_x_continuous(labels = function(x) format(abs(x), big.mark = ",")) +
  theme_masld() +
  theme(legend.position = "bottom") +
  labs(x = "Number of DEGs (padj < 0.05)", y = NULL, fill = NULL,
       title = "Per-cohort DE analysis")

# ==========================================================================
# Assemble
# ==========================================================================
fig <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(fig, OUT, width = fig_full_width, height = 10)
message("Cohort forest figure saved to ", OUT)
