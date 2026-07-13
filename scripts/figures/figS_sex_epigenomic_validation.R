#!/usr/bin/env Rscript
# figS_sex_epigenomic_validation.R — Supplementary figure panels for
# sex-stratified epigenomic validation (scATAC-seq)
#
# 4 panels:
#   A: XIST vs Y-score scatter (sex inference validation)
#   B: Promoter enrichment bar plot (Fisher OR + CI per sex class)
#   C: chromVAR TF heatmap (top sex-differential TFs)
#   D: RNA vs ATAC sex LFC scatter (concordance)
#
# Reads CSVs from Analysis/ATAC/Human_Multiome/results/sex_validation/
# Outputs to figures/supplementary/figS05_epigenomic_spatial/

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# I/O paths
RESULTS_DIR <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/sex_validation")
FIG_DIR     <- FIGS05_DIR
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Panel A: Sex inference scatter
# ---------------------------------------------------------------------------
plot_sex_inference <- function() {
  df <- read.csv(file.path(RESULTS_DIR, "donor_sex_inference.csv"),
                 stringsAsFactors = FALSE)

  p <- ggplot(df, aes(x = y_score, y = xist_score,
                       color = inferred_sex, shape = condition)) +
    geom_point(size = 3, stroke = 0.5) +
    ggrepel::geom_text_repel(aes(label = donor_id), size = 2, max.overlaps = 20,
                              show.legend = FALSE) +
    scale_color_manual(values = c(Female = masld_colors$female,
                                   Male   = masld_colors$male),
                        name = "Inferred sex") +
    scale_shape_manual(values = c(NORMAL = 16, MASL = 17, MASH = 15),
                        name = "Condition") +
    labs(x = "Mean Y-linked gene activity",
         y = "Mean XIST gene activity",
         title = "Sex inference from scATAC gene activity") +
    theme_masld() +
    theme(legend.position = "right")

  out <- file.path(FIG_DIR, "panel_sex_validation_A.pdf")
  ggsave(out, p, width = fig_half_width, height = 3, device = cairo_pdf)
  message("Saved: ", out)
}

# ---------------------------------------------------------------------------
# Panel B: Promoter enrichment bar plot
# ---------------------------------------------------------------------------
plot_promoter_enrichment <- function() {
  df <- read.csv(file.path(RESULTS_DIR, "promoter_enrichment_summary.csv"),
                 stringsAsFactors = FALSE)

  if (nrow(df) == 0) {
    message("No enrichment results, skipping panel B")
    return(invisible(NULL))
  }

  # Use disease-only results (properly controlled)
  if ("subset" %in% colnames(df)) {
    df <- df[df$subset == "disease_only", ]
  }

  df$sex_class <- factor(df$sex_class,
                          levels = c("Female_biased", "Male_biased", "Divergent"))

  # Cap CI for plotting
  df$ci_upper_capped <- pmin(df$ci_upper, 10)

  # Significance labels
  df$sig_label <- ifelse(df$fisher_p < 0.001, "***",
                   ifelse(df$fisher_p < 0.01, "**",
                    ifelse(df$fisher_p < 0.05, "*", "n.s.")))

  p <- ggplot(df, aes(x = sex_class, y = odds_ratio, fill = sex_class)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper_capped),
                  width = 0.2, linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_text(aes(label = sig_label, y = ci_upper_capped + 0.1),
              size = 2.5, vjust = 0) +
    geom_text(aes(label = sprintf("n=%d", n_genes), y = 0.05),
              size = 2, color = "white", vjust = 0) +
    scale_fill_manual(values = sex_class_colors, guide = "none") +
    scale_x_discrete(labels = c("Female_biased" = "Female-\nbiased",
                                 "Male_biased"   = "Male-\nbiased",
                                 "Divergent"     = "Divergent")) +
    labs(x = NULL, y = "Odds ratio (vs Concordant)",
         title = "Sex-preferential promoter accessibility\n(CPM-normalized, disease only)") +
    theme_masld() +
    theme(axis.text.x = element_text(size = 6))

  out <- file.path(FIG_DIR, "panel_sex_validation_B.pdf")
  ggsave(out, p, width = fig_half_width * 0.7, height = 3, device = cairo_pdf)
  message("Saved: ", out)
}

# ---------------------------------------------------------------------------
# Panel C: chromVAR TF heatmap
# ---------------------------------------------------------------------------
plot_chromvar_sex <- function() {
  df <- read.csv(file.path(RESULTS_DIR, "chromvar_sex_differential.csv"),
                 stringsAsFactors = FALSE)

  if (nrow(df) == 0) {
    message("No chromVAR results, skipping panel C")
    return(invisible(NULL))
  }

  # Select top TFs: prioritize disease regulon TFs, then by p-value
  df$in_disease_regulon <- as.logical(df$in_disease_regulon)
  df <- df %>% arrange(wilcoxon_p)

  # Take all disease regulon TFs + top non-regulon TFs
  disease_tfs <- df %>% filter(in_disease_regulon) %>% head(24)
  other_tfs   <- df %>% filter(!in_disease_regulon) %>% head(max(0, 20 - nrow(disease_tfs)))
  top <- bind_rows(disease_tfs, other_tfs) %>%
    arrange(wilcoxon_p) %>%
    head(25)

  if (nrow(top) == 0) {
    message("No TFs to plot, skipping panel C")
    return(invisible(NULL))
  }

  # Prepare for heatmap-style dot plot
  top$tf_label <- paste0(
    top$tf_name,
    ifelse(top$in_disease_regulon, " *", "")
  )
  top$tf_label <- factor(top$tf_label, levels = rev(top$tf_label))
  top$neg_log10_p <- -log10(top$padj + 1e-300)

  p <- ggplot(top, aes(x = delta_dev, y = tf_label)) +
    geom_point(aes(size = neg_log10_p,
                    color = ifelse(delta_dev > 0, "Female > Male", "Male > Female")),
               stroke = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
    scale_color_manual(values = c("Female > Male" = masld_colors$female,
                                   "Male > Female" = masld_colors$male),
                        name = "Direction") +
    scale_size_continuous(name = expression(-log[10]~padj),
                           range = c(1, 4)) +
    labs(x = expression(Delta~"deviation (F - M)"),
         y = NULL,
         title = "Sex-differential TF motif activity\n(hepatocytes, * = disease regulon)") +
    theme_masld() +
    theme(axis.text.y = element_text(size = 5.5),
          legend.position = "right")

  out <- file.path(FIG_DIR, "panel_sex_validation_C.pdf")
  ggsave(out, p, width = fig_half_width, height = 4.5, device = cairo_pdf)
  message("Saved: ", out)
}

# ---------------------------------------------------------------------------
# Panel D: RNA vs ATAC sex LFC concordance
# ---------------------------------------------------------------------------
plot_concordance <- function() {
  df <- read.csv(file.path(RESULTS_DIR, "gene_activity_sex_concordance.csv"),
                 stringsAsFactors = FALSE)

  if (nrow(df) == 0) {
    message("No concordance results, skipping panel D")
    return(invisible(NULL))
  }

  # Compute global correlation for annotation
  cor_test <- cor.test(df$rna_sex_lfc, df$atac_sex_lfc, method = "spearman")
  rho_label <- sprintf("rho = %.3f\np = %.2e\nn = %d",
                        cor_test$estimate, cor_test$p.value, nrow(df))

  # Downsample concordant for visual clarity
  set.seed(42)
  concordant <- df %>% filter(sex_class == "Concordant")
  others     <- df %>% filter(sex_class != "Concordant")
  if (nrow(concordant) > 2000) {
    concordant <- concordant %>% slice_sample(n = 2000)
  }
  plot_df <- bind_rows(concordant, others)

  # Order so sex-dimorphic genes plot on top
  plot_df$sex_class <- factor(plot_df$sex_class,
                               levels = c("Concordant", "Divergent",
                                          "Male_biased", "Female_biased"))

  p <- ggplot(plot_df, aes(x = rna_sex_lfc, y = atac_sex_lfc, color = sex_class)) +
    rasterize_layer(
      geom_point(alpha = 0.3, size = 0.5, stroke = 0)
    ) +
    geom_smooth(data = df, inherit.aes = FALSE,
                aes(x = rna_sex_lfc, y = atac_sex_lfc),
                method = "lm", se = TRUE, color = "grey30",
                linewidth = 0.4, linetype = "dashed") +
    scale_color_manual(values = sex_class_colors, name = "Sex class") +
    annotate("text", x = Inf, y = Inf, label = rho_label,
             hjust = 1.1, vjust = 1.3, size = 2.5, fontface = "italic") +
    labs(x = expression(RNA~sex~LFC~(logFC[F] - logFC[M])),
         y = expression(ATAC~sex~LFC~(log[2]~(F/M))),
         title = "Transcriptomic vs epigenomic\nsex effect concordance") +
    theme_masld() +
    theme(legend.position = "right")

  out <- file.path(FIG_DIR, "panel_sex_validation_D.pdf")
  ggsave(out, p, width = fig_half_width, height = 3.5, device = cairo_pdf)
  message("Saved: ", out)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
message("=== figS_sex_epigenomic_validation.R ===")
message("Results dir: ", RESULTS_DIR)
message("Figure dir:  ", FIG_DIR)

if (!dir.exists(RESULTS_DIR)) {
  stop("Results directory not found: ", RESULTS_DIR,
       "\nRun 13_sex_epigenomic_validation.py first.")
}

plot_sex_inference()
plot_promoter_enrichment()
plot_chromvar_sex()
plot_concordance()

message("=== All panels complete ===")
