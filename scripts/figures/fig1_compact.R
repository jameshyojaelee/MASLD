##############################################################################
# Figure 1 (Compact): Multi-Cohort MASLD Transcriptomic Atlas
# 6 panels: (a) study design [Illustrator], (b) sunburst wheel,
#   (c) metadata overview, (d) UpSet (per-study vs integrated DEGs),
#   (e) LOO-CV, (f) AUROC heatmap
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
  library(png)
  library(grid)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIG1_DIR, "fig1_compact.pdf")

# Helper: embed an external PDF as a raster grob wrapped for patchwork
embed_pdf <- function(pdf_path, dpi = 600) {
  tmp_png <- tempfile(fileext = ".png")
  cmd <- sprintf("convert -density %d '%s' -flatten '%s'", dpi, pdf_path, tmp_png)
  system(cmd)
  img <- png::readPNG(tmp_png)
  unlink(tmp_png)
  grob <- grid::rasterGrob(img, interpolate = TRUE)
  wrap_elements(full = grob)
}

# Pre-initialize all 6 panels
p_a <- placeholder("(a) Study design\n(Illustrator)")
p_b <- placeholder("(b) Sunburst wheel")
p_c <- placeholder("(c) Metadata overview")
p_d <- placeholder("(d) UpSet plot")
p_e <- placeholder("(e) LOO-CV robustness")
p_f <- placeholder("(f) AUROC heatmap")

# Study name mapping (shared) — PRJNA512027 (Gerhard 2018) excluded from
# cohort presentation due to L0/S0 library-prep / diagnosis confound.
STUDY_NAMES <- c(
  GSE213621 = "Chen", GSE135251 = "Govaere", GSE130970 = "Hoang",
  GSE162694 = "Bril", GSE174478 = "Kawamura", GSE193066 = "Hoshida",
  GSE240729 = "Verschuren", GSE126848 = "Suppli",
  GSE167523 = "Pantano"
)

# ==========================================================================
# (b) Sunburst wheel (embed pre-rendered PDF)
# ==========================================================================
sunburst_pdf <- file.path(FIG1_DIR, "fig1b_wheel_sunburst.pdf")
if (file.exists(sunburst_pdf)) {
  p_b <- embed_pdf(sunburst_pdf)
}

# ==========================================================================
# (c) Metadata overview (embed pre-rendered PDF)
# ==========================================================================
meta_pdf <- file.path(FIG1_DIR, "bulkrna_metadata_panel.pdf")
if (file.exists(meta_pdf)) {
  p_c <- embed_pdf(meta_pdf)
}

# ==========================================================================
# (d) UpSet — per-study DEGs vs dream integrated DEGs
#     Only cohorts with healthy controls (Disease vs Control contrast)
# ==========================================================================
dream_masld <- load_dream_results()
per_study <- load_per_study_de()

if (!is.null(per_study) && !is.null(dream_masld)) {
  # 5 studies with Disease vs Control contrast (healthy controls)
  upset_studies <- c("GSE126848", "GSE130970", "GSE135251",
                     "GSE162694", "GSE213621")

  upset_names <- c(
    GSE126848 = "Suppli",
    GSE130970 = "Hoang",
    GSE135251 = "Govaere",
    GSE162694 = "Bril",
    GSE213621 = "Chen"
  )

  # Consistent thresholds
  padj_thr <- 0.1
  lfc_thr  <- 0.5   # 2026-05-28 (P0-D): canonical Tier-1 cutoff (STAR -s 2; PI decision, within CV stability plateau)

  # Build DEG lists per study
  deg_lists <- list()
  for (ds in upset_studies) {
    ds_data <- per_study[dataset == ds]
    pcol <- intersect(c("padj", "adj.P.Val"), names(ds_data))[1]
    lcol <- intersect(c("logFC", "study_logFC"), names(ds_data))[1]
    if (is.na(pcol) || is.na(lcol)) next
    degs <- ds_data[get(pcol) < padj_thr & abs(get(lcol)) > lfc_thr, gene]
    degs <- sub("\\..*", "", degs)
    deg_lists[[upset_names[ds]]] <- unique(degs)
  }

  # Dream integrated DEGs
  dream_degs <- dream_masld[dream_padj < padj_thr & abs(dream_logFC) > lfc_thr, gene]
  dream_degs <- sub("\\..*", "", dream_degs)
  deg_lists[["Integrated"]] <- unique(dream_degs)

  set_names <- names(deg_lists)

  # Binary membership matrix
  all_genes_upset <- unique(unlist(deg_lists))
  bin_dt <- data.table(gene = all_genes_upset)
  for (s in set_names) {
    bin_dt[, (s) := as.integer(gene %in% deg_lists[[s]])]
  }

  # Compute intersections
  bin_dt[, int_key := do.call(paste0, .SD), .SDcols = set_names]
  int_counts <- bin_dt[, .N, by = int_key][order(-N)]
  top_n <- min(25, nrow(int_counts))
  top_ints <- int_counts[1:top_n]
  top_ints[, rank := .I]

  # Decode which sets are active per intersection
  for (i in seq_along(set_names)) {
    top_ints[, (set_names[i]) := as.integer(substr(int_key, i, i))]
  }

  # Set sizes for axis labels
  set_sizes <- vapply(set_names, function(s) sum(bin_dt[[s]]), integer(1))
  set_labels_vec <- paste0(set_names, " (", format(set_sizes, big.mark = ","), ")")

  # Melt for dot matrix
  dot_long <- melt(top_ints, id.vars = c("int_key", "N", "rank"),
                   measure.vars = set_names,
                   variable.name = "set", value.name = "active")
  dot_long[, set := factor(set, levels = set_names)]

  # Segments connecting active dots within each intersection
  seg_data <- dot_long[active == 1, .(
    ymin = min(as.numeric(set)),
    ymax = max(as.numeric(set))
  ), by = rank]
  seg_data <- seg_data[ymin != ymax]

  # Highlight Integrated-only bars
  top_ints[, has_dream := as.logical(get("Integrated"))]
  n_active <- rowSums(as.matrix(top_ints[, ..set_names]))
  top_ints[, bar_fill := fifelse(has_dream & n_active == 1,
                                  masld_colors$up,
                                  fifelse(has_dream, masld_colors$mash, masld_colors$down))]

  # ---- Bar chart (intersection sizes) ----
  p_d_bars <- ggplot(top_ints, aes(x = rank, y = N)) +
    geom_col(fill = top_ints$bar_fill, width = 0.7) +
    geom_text(aes(label = ifelse(N >= max(N) * 0.03,
                                  format(N, big.mark = ","), "")),
              vjust = -0.3, size = 1.5) +
    scale_x_continuous(limits = c(0.4, top_n + 0.6), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                       labels = comma) +
    labs(y = "Intersection\nsize",
         # 2026-05-28 (P0-D): title derived from the actual thresholds so it can't drift from the code
         title = sprintf("Per-study vs integrated DEGs (padj<%g, |LFC|>%g)", padj_thr, lfc_thr)) +
    theme_masld() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          plot.margin = margin(5, 5, 0, 5))

  # ---- Dot matrix (set membership indicators) ----
  p_d_dots <- ggplot() +
    geom_segment(data = seg_data,
                 aes(x = rank, xend = rank, y = ymin, yend = ymax),
                 color = "gray30", linewidth = 0.4) +
    geom_point(data = dot_long,
               aes(x = rank, y = as.numeric(set),
                   fill = factor(active)),
               shape = 21, size = 2.5, stroke = 0.3, color = "gray40") +
    scale_fill_manual(values = c("0" = "#E8E8E8", "1" = "gray20"),
                      guide = "none") +
    scale_x_continuous(limits = c(0.4, top_n + 0.6), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(set_names),
                       labels = set_labels_vec,
                       limits = c(0.5, length(set_names) + 0.5)) +
    labs(x = NULL, y = NULL) +
    theme_masld() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(0, 5, 5, 5))

  # Combine into single panel
  p_d <- p_d_bars / p_d_dots +
    plot_layout(heights = c(2, 2))
}

# ==========================================================================
# (e) LOO-CV stability (compact horizontal bar)
# ==========================================================================
LOO_DIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/loo_cv")
loo_summary_f <- file.path(LOO_DIR, "loo_cv_summary.csv")
if (file.exists(loo_summary_f)) {
  loo <- fread(loo_summary_f)
  loo[, author := STUDY_NAMES[held_out]]

  # Load sample counts for sorting
  qc_file <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")
  if (file.exists(qc_file)) {
    qc <- fread(qc_file)
    cohort_n <- qc[pass_technical == TRUE, .(n_samples = .N), by = dataset]
    loo <- merge(loo, cohort_n, by.x = "held_out", by.y = "dataset", all.x = TRUE)
  } else {
    loo[, n_samples := 0L]
  }
  loo[is.na(n_samples), n_samples := 0L]
  loo[, author := factor(author, levels = author[order(n_samples)])]

  p_e <- ggplot(loo, aes(x = author, y = pct_full_recovered)) +
    geom_col(fill = masld_colors$down, width = 0.65) +
    geom_hline(yintercept = 100, linetype = "dashed", linewidth = 0.3, color = "gray50") +
    geom_text(aes(label = sprintf("%.0f%%", pct_full_recovered)),
              hjust = -0.15, size = 1.6, color = "gray30") +
    labs(x = NULL, y = "Recovery (%)",
         title = sprintf("LOO-CV (mean %.0f%%)", mean(loo$pct_full_recovered))) +
    coord_flip(ylim = c(0, 115)) +
    theme_masld() +
    theme(axis.text.y = element_text(size = 5),
          plot.margin = margin(2, 2, 2, 2))
}

# ==========================================================================
# (f) AUROC heatmap (cross-study prediction)
# ==========================================================================
auroc_matrix_f <- file.path(LOO_DIR, "auroc_matrix.csv")
if (file.exists(auroc_matrix_f)) {
  auroc_mat <- fread(auroc_matrix_f)
  auroc_mat[, train_author := STUDY_NAMES[train]]
  auroc_mat[, test_author := STUDY_NAMES[test]]

  # Add integrated signature row if available
  auroc_int_f <- file.path(LOO_DIR, "auroc_integrated.csv")
  if (file.exists(auroc_int_f)) {
    auroc_int <- fread(auroc_int_f)
    auroc_int[, test_author := STUDY_NAMES[held_out]]
    auroc_int[, train_author := "Integrated"]
    auroc_int[, train := "Integrated"]
    auroc_int[, test := held_out]
    int_rows <- auroc_int[, .(train, test, auroc, n_test, n_disease, n_control,
                               train_author, test_author)]
    auroc_all <- rbind(auroc_mat, int_rows, fill = TRUE)
  } else {
    auroc_all <- auroc_mat
  }

  # Fibrosis-contrast cohorts (no healthy controls) — train rows only
  FIBROSIS_ONLY <- c("Kawamura", "Hoshida", "Verschuren")

  # Factor levels: single cohorts alphabetical, Integrated at bottom
  train_levels <- c(sort(unique(auroc_mat$train_author)), "Integrated")
  test_levels <- sort(unique(auroc_all$test_author))
  auroc_all[, train_author := factor(train_author, levels = rev(train_levels))]
  auroc_all[, test_author := factor(test_author, levels = test_levels)]

  # Dagger annotation for fibrosis-contrast training signatures
  train_labels_display <- rev(train_levels)
  train_labels_styled <- ifelse(
    train_labels_display %in% FIBROSIS_ONLY,
    paste0(train_labels_display, "\u2020"),
    train_labels_display
  )

  # Gray out self-predictions (diagonal)
  auroc_all[, fill_val := fifelse(train == test, NA_real_, auroc)]

  p_f <- ggplot(auroc_all, aes(x = test_author, y = train_author, fill = fill_val)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(train == test, "-", sprintf("%.2f", auroc))),
              size = 1.5, color = "black") +
    scale_fill_gradient(low = "white", high = masld_colors$up,
                        limits = c(0.5, 1), oob = squish,
                        na.value = "gray90", name = "AUROC") +
    scale_y_discrete(labels = setNames(train_labels_styled, rev(train_levels))) +
    labs(x = "Test cohort", y = "Training signature",
         title = "Cross-study prediction") +
    theme_masld() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          axis.text.y = element_text(size = 5),
          panel.grid = element_blank())

  # Bold rectangle around integrated row
  if ("Integrated" %in% train_levels) {
    int_y <- which(rev(train_levels) == "Integrated")
    p_f <- p_f +
      annotate("rect", xmin = 0.5, xmax = length(test_levels) + 0.5,
               ymin = int_y - 0.5, ymax = int_y + 0.5,
               fill = NA, color = masld_colors$up, linewidth = 0.8)
  }
}

# ==========================================================================
# Compose
# Row 1: (a) study design | (b) sunburst
# Row 2: (c) metadata overview (full width)
# Row 3: (d) UpSet (full width)
# Row 4: (e) LOO-CV (narrow) | (f) AUROC heatmap (wide)
# ==========================================================================
fig1 <- (p_a | p_b) /
         p_c /
         wrap_elements(full = p_d) /
         (p_e + p_f + plot_layout(widths = c(1, 2)))

fig1 <- fig1 +
  plot_layout(heights = c(1, 1, 1.2, 1.2)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

save_fig_tall(fig1, OUT, width = fig_full_width, height = 13)
cat("Saved:", OUT, "\n")
