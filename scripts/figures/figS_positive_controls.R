# Supplementary Figure: Positive Controls Validation
# Panels: (a) Expression-driven controls, (b) Genetic-risk controls

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

dir.create(file.path(FIGS01_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(FIGS01_DIR, "panels", "figS_positive_controls.pdf")

pos_ctrl <- load_positive_controls()

if (!is.null(pos_ctrl)) {
  lfc_col  <- intersect(c("bulk_logFC", "logFC"), names(pos_ctrl))[1]
  padj_col <- intersect(c("bulk_padj", "padj"), names(pos_ctrl))[1]
  sym_col  <- intersect(c("symbol", "gene"), names(pos_ctrl))[1]

  # Determine control type split
  type_col <- intersect(c("control_type", "type"), names(pos_ctrl))[1]
  if (!is.na(type_col)) {
    expr_ctrl <- pos_ctrl[grepl("express", get(type_col), ignore.case = TRUE)]
    gene_ctrl <- pos_ctrl[grepl("genetic|risk", get(type_col), ignore.case = TRUE)]
  } else {
    # Fallback: split by known genetic-risk genes
    genetic_genes <- c("PNPLA3", "TM6SF2", "GCKR", "MBOAT7", "HSD17B13")
    expr_ctrl <- pos_ctrl[!get(sym_col) %in% genetic_genes]
    gene_ctrl <- pos_ctrl[get(sym_col) %in% genetic_genes]
  }

  # ---- Panel (a): Expression-driven dot plot ----
  if (nrow(expr_ctrl) > 0 && !is.na(lfc_col) && !is.na(padj_col)) {
    expr_ctrl[, neg_log10p := -log10(get(padj_col))]
    expr_ctrl[, sig := get(padj_col) < 0.1]
    expr_ctrl[, gene_label := get(sym_col)]
    expr_ctrl <- expr_ctrl[order(-abs(get(lfc_col)))]
    expr_ctrl[, gene_label := factor(gene_label, levels = rev(gene_label))]

    p_a <- ggplot(expr_ctrl, aes(x = get(lfc_col), y = gene_label, color = sig, size = neg_log10p)) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "gray50") +
      geom_point(shape = 16) +
      scale_color_manual(values = c(`TRUE` = masld_colors$up, `FALSE` = masld_colors$ns),
                         labels = c("NS", "padj < 0.1"), name = NULL) +
      scale_size_continuous(range = c(0.5, 3), name = expression(-log[10]~padj)) +
      labs(x = "Meta-analysis logFC", y = NULL, title = "Expression-driven controls") +
      theme_masld() +
      theme(axis.text.y = element_text(face = "italic", size = 6))
  } else {
    p_a <- placeholder("Panel a: expression controls insufficient")
  }

  # ---- Panel (b): Genetic-risk controls ----
  if (nrow(gene_ctrl) > 0 && !is.na(lfc_col) && !is.na(padj_col)) {
    gene_ctrl[, neg_log10p := -log10(get(padj_col))]
    gene_ctrl[, sig := get(padj_col) < 0.1]
    gene_ctrl[, gene_label := get(sym_col)]
    gene_ctrl <- gene_ctrl[order(-abs(get(lfc_col)))]
    gene_ctrl[, gene_label := factor(gene_label, levels = rev(gene_label))]

    # Add rationale annotation if available
    rationale_col <- intersect(c("de_expected_rationale", "rationale"), names(gene_ctrl))[1]

    p_b <- ggplot(gene_ctrl, aes(x = get(lfc_col), y = gene_label, color = sig)) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "gray50") +
      geom_point(size = 2.5, shape = 16) +
      scale_color_manual(values = c(`TRUE` = masld_colors$up, `FALSE` = masld_colors$ns),
                         labels = c("NS", "padj < 0.1"), name = NULL) +
      labs(x = "Meta-analysis logFC", y = NULL,
           title = "Genetic-risk controls",
           subtitle = "DE not expected (protein-altering variants)") +
      theme_masld() +
      theme(axis.text.y = element_text(face = "italic", size = 6),
            plot.subtitle = element_text(size = 5, color = "gray40"))
  } else {
    p_b <- placeholder("Panel b: genetic-risk controls not found")
  }
} else {
  p_a <- placeholder("Panel a: positive controls not found")
  p_b <- placeholder("Panel b: positive controls not found")
}

# ---- Assemble ----
figS2 <- (p_a | p_b) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig(figS2, OUT, height = 5)
message("FigS2 saved to ", OUT)
