#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Fig 3 (RNA-seq) panel: Cell-type bulk attribution
#
# Left:  Stacked bar — MuSiC category breakdown of bulk DEGs
# Right: Horizontal bar — cell-type attribution for intrinsic DEGs
#
# Data: C2 recount attribution matrix (canonical limma-voom-qw Tier-1).
# Output: figures/main/fig3_RNAseq/panels/celltype_bulk_attribution.pdf
# ---------------------------------------------------------------------------

# --- Shared theme + data loaders -----------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

# --- Load data -----------------------------------------------------------
# C2 recount attribution matrix (Script 80_*_C2.R): bulk side = canonical
# limma-voom-qw Tier-1, MuSiC + per-cell-type pseudobulk concordance all
# recomputed on C2. Supersedes the dream-era top-level matrix. The matrix
# carries its own canonical-C2 bulk_lfc/bulk_padj, so no external join needed.
attr_dt <- fread(file.path(BASE,
  "RNA-seq/results/celltype_attribution/c2_recount/celltype_attribution_matrix.csv"))

# Bulk DEGs: C2 Tier-1 (padj < 0.05 & |logFC| > 0.5) using the matrix's own
# canonical bulk columns — same definition used to assign attribution_class.
merged <- attr_dt[bulk_padj < 0.05 & abs(bulk_lfc) > 0.5,
                  .(ensg_base, music_category, attribution_class)]

# =========================================================================
# LEFT PANEL — MuSiC category stacked bar
# =========================================================================

# Collapse Unmasked into Hepatocyte_intrinsic (1 gene; minor)
merged[music_category == "Unmasked", music_category := "Hepatocyte_intrinsic"]

# Summarize
cat_summary <- merged[, .N, by = music_category]
cat_summary[, pct := 100 * N / sum(N)]
cat_summary[, label := sprintf("%s\n(n = %s, %.1f%%)",
  c(Hepatocyte_intrinsic = "Hepatocyte-\nintrinsic",
    Composition_driven   = "Composition-\ndriven",
    Not_significant      = "NS in scRNA")[music_category],
  formatC(N, big.mark = ","), pct)]

# Order for stacking (bottom to top: Hep-intrinsic, Composition, NS).
# Under C2 Tier-1 the DEG set partitions into Hepatocyte_intrinsic +
# Composition_driven only (Not_significant applies to non-DEG genes), so drop
# any absent level to keep the legend honest.
stack_levels <- intersect(
  c("Not_significant", "Composition_driven", "Hepatocyte_intrinsic"),
  unique(cat_summary$music_category))
cat_summary[, music_category := factor(music_category, levels = stack_levels)]
# Order rows to match position_stack, which builds the bar from the HIGHEST
# factor level at the base. Computing cumsum in this same (descending) order
# keeps the text-label midpoints aligned with their colored segments.
cat_summary <- cat_summary[order(-as.integer(music_category))]

# Cumulative positions for label placement
cat_summary[, cumN := cumsum(N)]
cat_summary[, midpoint := cumN - N / 2]

# Colors per spec
cat_colors <- c(
  Hepatocyte_intrinsic = "#0D47A1",
  Composition_driven   = "#F57F17",
  Not_significant      = "#9E9E9E"
)

# Display labels (cleaner, without linebreak for legend)
cat_labels <- c(
  Hepatocyte_intrinsic = "Hepatocyte-intrinsic",
  Composition_driven   = "Composition-driven",
  Not_significant      = "NS in scRNA"
)

total_n <- sum(cat_summary$N)

p_left <- ggplot(cat_summary, aes(x = 1, y = N, fill = music_category)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(y = midpoint,
                label = sprintf("n=%s\n%.1f%%",
                                formatC(N, big.mark = ","), pct)),
            size = PUB_GEOM_TEXT, color = "white", fontface = "plain",
            lineheight = 0.85) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = cat_colors,
                    breaks = rev(levels(cat_summary$music_category)),
                    labels = cat_labels[rev(levels(cat_summary$music_category))]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL,
       y = paste0("Bulk DEGs (n = ", formatC(total_n, big.mark = ","), ")"),
       fill = NULL) +
  theme_masld() + theme_pub() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.margin = margin(t = 0),
    legend.box.margin = margin(t = -2),
    plot.margin = margin(3, 6, 3, 3)
  ) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 1,
                             override.aes = list(color = NA)))

# =========================================================================
# RIGHT PANEL — cell-type attribution for intrinsic DEGs
# =========================================================================

# Keep only genes with explicit cell-type attribution
ct_genes <- merged[grepl("_intrinsic_", attribution_class)]

# Parse cell type from attribution_class
ct_genes[, cell_type := sub("_intrinsic_(likely|strong)$", "",
                            attribution_class)]
# Clean underscores to match ct_palette names
ct_name_map <- c(
  "Hepatocytes"        = "Hepatocytes",
  "Cholangiocytes"     = "Cholangiocytes",
  "Endothelial_cells"  = "Endothelial cells",
  "Fibroblasts"        = "Fibroblasts",
  "Macrophages"        = "Macrophages",
  "T_cells"            = "T cells",
  "B_cells"            = "B cells",
  "Plasma_cells"       = "Plasma cells",
  "Resident_NK"        = "Resident NK",
  "Circulating_NK_NKT" = "Circulating NK/NKT",
  "Mono+mono_derived_cells" = "Mono+mono derived cells"
)
ct_genes[, cell_type_clean := ifelse(cell_type %in% names(ct_name_map),
                                     ct_name_map[cell_type],
                                     gsub("_", " ", cell_type))]

# Summarize per cell type
ct_summary <- ct_genes[, .N, by = cell_type_clean]
ct_summary <- ct_summary[order(-N)]
ct_summary[, cell_type_clean := factor(cell_type_clean,
                                       levels = rev(ct_summary$cell_type_clean))]

# Map colors from ct_palette; fall back to gray
ct_summary[, bar_color := ifelse(cell_type_clean %in% names(ct_palette),
                                 ct_palette[as.character(cell_type_clean)],
                                 "#9E9E9E")]
bar_colors <- setNames(ct_summary$bar_color,
                       as.character(ct_summary$cell_type_clean))

p_right <- ggplot(ct_summary, aes(x = cell_type_clean, y = N,
                                  fill = cell_type_clean)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = N), hjust = -0.2, size = PUB_GEOM_TEXT) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL,
       y = paste0("Cell-type-attributed DEGs (n = ",
                  sum(ct_summary$N), ")"),
       title = NULL) +
  theme_masld() + theme_pub() +
  theme(
    plot.margin = margin(3, 3, 3, 3)
  )

# =========================================================================
# Composite
# =========================================================================
p_composite <- p_left + p_right +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = list(c("", "")))

# Output
out_dir <- file.path(FIG2_DIR, "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "celltype_bulk_attribution.pdf")

# ggsave(out_file, p_composite,
#        width = 180 / 25.4, height = 60 / 25.4,
#        device = cairo_pdf)
# 
# cat("Saved:", out_file, "\n")
# cat("File size:", file.size(out_file), "bytes\n")
