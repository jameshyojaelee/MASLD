##############################################################################
# Figure 1: Atlas Overview (4 panels, remaining panels TBD)
#   Row 1: (a) Study design schematic (~80%),  (b) Convergence wheel (~20%)
#   Row 2: (c) Data & metadata overview (full width)
#   Row 3: (d) UpSet: per-study vs integrated DEGs (full width)
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
  library(edgeR)
  library(png)
  library(grid)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT <- file.path(FIG1_DIR, "fig1_atlas_overview.pdf")

# Helper: embed an external PDF as a raster grob wrapped for patchwork
# Uses ImageMagick convert (system) → PNG → grid::rasterGrob
embed_pdf <- function(pdf_path, dpi = 600) {
  tmp_png <- tempfile(fileext = ".png")
  cmd <- sprintf("convert -density %d '%s' -flatten '%s'", dpi, pdf_path, tmp_png)
  system(cmd)
  img <- png::readPNG(tmp_png)
  unlink(tmp_png)
  grob <- grid::rasterGrob(img, interpolate = TRUE)
  wrap_elements(full = grob)
}

# Pre-initialize all panels with placeholders
p_a <- placeholder("Panel a: Study design schematic\n(Insert via Illustrator)")
p_b <- placeholder("Panel b:\nConvergence\nwheel")
p_c <- placeholder("Panel c: Data & metadata overview")
p_d <- placeholder("Panel d: UpSet plot")

# ==========================================================================
# Panel (b): Convergence wheel (embed from pre-rendered PDF)
# ==========================================================================
wheel_pdf <- file.path(FIG1_DIR, "fig1b_wheel_notitle.pdf")
if (file.exists(wheel_pdf)) {
  message("Embedding convergence wheel (no title) from ", wheel_pdf)
  p_b <- embed_pdf(wheel_pdf)
} else {
  message("WARNING: convergence wheel PDF not found at ", wheel_pdf)
}

# ==========================================================================
# Panel (c): Metadata matrix (embed from pre-rendered PDF)
# ==========================================================================
meta_pdf <- file.path(FIG1_DIR, "bulkrna_metadata_panel.pdf")
if (file.exists(meta_pdf)) {
  message("Embedding metadata matrix from ", meta_pdf)
  p_c <- embed_pdf(meta_pdf)
} else {
  message("WARNING: metadata matrix PDF not found at ", meta_pdf)
}

# ==========================================================================
# Panel (d): UpSet — per-study DEGs vs dream integrated DEGs
# ==========================================================================
dream_masld <- load_dream_results()
per_study <- load_per_study_de()

if (!is.null(per_study) && !is.null(dream_masld)) {
  # Only use studies with Disease vs Control contrast (those with controls)
  # PRJNA512027 excluded: L0/S0 library-prep batch confounded with disease severity
  upset_studies <- c("GSE126848", "GSE130970", "GSE135251",
                     "GSE162694", "GSE213621")

  # Short display names (first author)
  upset_names <- c(
    GSE126848   = "Suppli",
    GSE130970   = "Hoang",
    GSE135251   = "Govaere",
    GSE162694   = "Bril",
    GSE213621   = "Chen"
  )

  # Consistent thresholds: padj < 0.1, |LFC| > 0.58
  padj_thr <- 0.1
  lfc_thr  <- 0.58

  # Build DEG lists per study
  deg_lists <- list()
  for (ds in upset_studies) {
    ds_data <- per_study[dataset == ds]
    pcol <- intersect(c("padj", "adj.P.Val"), names(ds_data))[1]
    lcol <- intersect(c("logFC", "study_logFC"), names(ds_data))[1]
    if (is.na(pcol) || is.na(lcol)) next
    degs <- ds_data[get(pcol) < padj_thr & abs(get(lcol)) > lfc_thr, gene]
    degs <- sub("\\..*", "", degs)   # strip Ensembl version
    deg_lists[[upset_names[ds]]] <- unique(degs)
  }

  # Dream DEGs
  dream_masld[, sig_cat := fifelse(dream_padj < 0.1 & dream_logFC > 0.5, "Up",
                     fifelse(dream_padj < 0.1 & dream_logFC < -0.5, "Down", "NS"))]
  dream_degs <- dream_masld[dream_padj < padj_thr & abs(dream_logFC) > lfc_thr, gene]
  dream_degs <- sub("\\..*", "", dream_degs)
  deg_lists[["Integrated"]] <- unique(dream_degs)

  # Set names in display order (studies first, Integrated last = top of y-axis)
  set_names <- c(names(deg_lists))

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
         title = "Per-study vs integrated DEGs (padj<0.1, |LFC|>0.5)") +
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

  # Combine into single panel (shorter bars relative to dots)
  p_d <- p_d_bars / p_d_dots +
    plot_layout(heights = c(2, 2))
}

# ==========================================================================
# Assemble figure
#   Row 1: A (~80%) + B (~20%)   — schematic + convergence wheel
#   Row 2: C (full width)        — metadata matrix
#   Row 3: D (full width)        — UpSet
# ==========================================================================
design <- "
AAABBB
CCCCCC
DDDDDD
"

fig1 <- p_a + p_b + p_c + wrap_elements(full = p_d) +
  plot_layout(design = design, heights = c(1.3, 0.8, 0.7)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_fig_tall(fig1, OUT, height = 10)
message("Fig 1 (atlas overview) saved to ", OUT)
