#!/usr/bin/env Rscript
# Fig 1 — scATAC UMAP triad (cell_type | F_stage_augmented | condition)
#
# Inputs:
#   Analysis/ATAC/Human_Multiome/results/snapatac2/umap_export.tsv.gz
#     produced by Analysis/ATAC/Human_Multiome/scripts/30_export_umap_triad_data.py
#
# Output:
#   figures/main/fig1_atlas_overview/panels/fig1_scatac_umap_triad.pdf
#
# Liang-style: bold black axis titles outside plot, compact legend below.
# 30% random subsample (~27K cells) at set.seed(42) for vector PDF size.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

UMAP_TSV <- file.path(
  BASE, "Analysis/ATAC/Human_Multiome/results/snapatac2/umap_export.tsv.gz"
)
OUT_PDF <- file.path(FIGS05_DIR, "figS05_scatac_umap_triad.pdf")
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load + subsample (30%, seed 42)
# ---------------------------------------------------------------------------
stopifnot(file.exists(UMAP_TSV))
dat <- fread(UMAP_TSV)
n_full <- nrow(dat)
set.seed(42)
idx <- sample.int(n_full, size = round(0.30 * n_full))
dat_sub <- dat[idx]
n_sub <- nrow(dat_sub)

message(sprintf("[load] %d cells, sub-sampling to %d (30%%)", n_full, n_sub))

# ---------------------------------------------------------------------------
# Palettes — label-transferred cell types (matching snRNA atlas):
#   Hepatocytes, Fibroblasts, Macrophages, Plasma cells, Cholangiocytes,
#   Endothelial cells, T cells, Circulating NK/NKT, Resident NK, Low_confidence
# ---------------------------------------------------------------------------
ct_colors <- c(
  "Hepatocytes"   = "#C9265E",  # magenta — hep (canonical positive)
  "Fibroblasts"   = "#FF9800",  # orange — fib lineage
  "Macrophages"   = "#9C27B0",  # purple — macs
  "Endothelial"   = "#1565C0",  # blue — endo (shortened)
  "Cholangiocytes"= "#4CAF50",  # green — chol
  "Lymphocytes"   = "#9E9E9E",  # gray — merged T + NK + NKT + Plasma
  "Low conf."     = "#E0E0E0"
)
# Map full label-transferred names to display labels (shorter for legend)
ct_rename <- c(
  "Hepatocytes"        = "Hepatocytes",
  "Fibroblasts"        = "Fibroblasts",
  "Macrophages"        = "Macrophages",
  "Endothelial cells"  = "Endothelial",
  "Cholangiocytes"     = "Cholangiocytes",
  "Plasma cells"       = "Lymphocytes",
  "T cells"            = "Lymphocytes",
  "Resident NK"        = "Lymphocytes",
  "Circulating NK/NKT" = "Lymphocytes",
  "Low_confidence"     = "Low conf."
)
dat_sub[, cell_type := ct_rename[as.character(cell_type)]]
ct_levels <- names(ct_colors)
dat_sub[, cell_type := factor(cell_type, levels = ct_levels)]

# F_stage_augmented: 5 distinct colors (gray -> blue -> amber -> orange -> magenta)
# Sequential but more perceptually distinct than a single gradient
fstage_levels <- c("0", "1", "2", "3", "4")
fstage_colors <- c(
  "0" = "#BDBDBD",  # gray (baseline)
  "1" = "#1565C0",  # blue
  "2" = "#FFB300",  # amber
  "3" = "#FF6F00",  # deep orange
  "4" = "#C9265E"   # magenta (advanced disease)
)
dat_sub[, F_stage_augmented := as.character(F_stage_augmented)]
dat_sub[, F_stage_augmented := factor(F_stage_augmented, levels = fstage_levels)]

# condition: Normal=gray, MASL=blue, MASH=magenta
cond_colors <- c(
  NORMAL = "#9E9E9E",
  MASL   = "#1565C0",
  MASH   = "#C9265E"
)
dat_sub[, condition := factor(condition, levels = names(cond_colors))]

# ---------------------------------------------------------------------------
# Shared UMAP scaffold: strip axes, add small bottom-left arrow labels
# ---------------------------------------------------------------------------
xrng <- range(dat_sub$x, na.rm = TRUE)
yrng <- range(dat_sub$y, na.rm = TRUE)
xpad <- 0.04 * diff(xrng)
ypad <- 0.04 * diff(yrng)
arrow_len_x <- 0.18 * diff(xrng)
arrow_len_y <- 0.18 * diff(yrng)
ax0_x <- xrng[1] - xpad
ax0_y <- yrng[1] - ypad

axis_arrow_layers <- list(
  annotate("segment",
           x = ax0_x, xend = ax0_x + arrow_len_x,
           y = ax0_y, yend = ax0_y,
           arrow = arrow(length = unit(0.06, "cm"), type = "closed"),
           linewidth = 0.3, color = "black"),
  annotate("segment",
           x = ax0_x, xend = ax0_x,
           y = ax0_y, yend = ax0_y + arrow_len_y,
           arrow = arrow(length = unit(0.06, "cm"), type = "closed"),
           linewidth = 0.3, color = "black"),
  annotate("text",
           x = ax0_x + arrow_len_x / 2, y = ax0_y - 0.025 * diff(yrng),
           label = "UMAP1", size = GEOM_TEXT_6PT, fontface = "plain", hjust = 0.5,
           vjust = 1, color = "black"),
  annotate("text",
           x = ax0_x - 0.025 * diff(xrng), y = ax0_y + arrow_len_y / 2,
           label = "UMAP2", size = GEOM_TEXT_6PT, fontface = "plain", angle = 90,
           hjust = 0.5, vjust = 0, color = "black")
)

theme_umap <- theme_void(base_size = 7, base_family = "Helvetica") +
  theme(
    plot.title       = element_text(size = 6, face = "plain", hjust = 0,
                                    margin = margin(b = 1)),
    legend.position  = "bottom",
    legend.title     = element_text(size = 6, face = "plain"),
    legend.text      = element_text(size = 6),
    legend.key.size  = unit(0.18, "cm"),
    legend.box.spacing = unit(1, "pt"),
    legend.margin    = margin(0, 0, 0, 0),
    legend.justification = "center",
    legend.box.just = "center",
    plot.margin      = margin(2, 2, 2, 2)
  )

make_umap <- function(df, color_col, color_scale, tag, legend_nrow = 2) {
  ggplot(df, aes(x = x, y = y, color = .data[[color_col]])) +
    rasterize_layer(
      geom_point(size = 0.1, alpha = 0.4, stroke = 0, shape = 16),
      dpi = 350
    ) +
    color_scale +
    axis_arrow_layers +
    coord_fixed(clip = "off") +
    labs(tag = tag) +
    guides(color = guide_legend(
      title = NULL,
      nrow = legend_nrow,
      override.aes = list(size = 1.4, alpha = 1)
    )) +
    theme_umap +
    theme(plot.tag = element_text(size = 9, face = "plain"),
          plot.tag.position = c(0.02, 0.98))
}

# ---------------------------------------------------------------------------
# Build the three panels
# ---------------------------------------------------------------------------
# Drop empty levels in this draw so the legend reflects only present cell types.
ct_present <- intersect(ct_levels, unique(as.character(dat_sub$cell_type)))
p_ct <- make_umap(
  dat_sub,
  color_col   = "cell_type",
  color_scale = scale_color_manual(values = ct_colors[ct_present],
                                   drop = TRUE,
                                   na.value = "#E0E0E0"),
  tag         = "a",
  legend_nrow = 3
)

p_fs <- make_umap(
  dat_sub[!is.na(F_stage_augmented)],
  color_col   = "F_stage_augmented",
  color_scale = scale_color_manual(values = fstage_colors,
                                   labels = paste0("F", fstage_levels),
                                   drop = TRUE,
                                   na.value = "#E0E0E0"),
  tag         = "b",
  legend_nrow = 1
)

cond_present <- intersect(names(cond_colors),
                          unique(as.character(dat_sub$condition)))
cond_labels <- c(NORMAL = "Normal", MASL = "MASL", MASH = "MASH")
p_cond <- make_umap(
  dat_sub[!is.na(condition)],
  color_col   = "condition",
  color_scale = scale_color_manual(values = cond_colors[cond_present],
                                   labels = cond_labels[cond_present],
                                   drop = TRUE,
                                   na.value = "#E0E0E0"),
  tag         = "c",
  legend_nrow = 1
)

# ---------------------------------------------------------------------------
# Compose with patchwork (~7.2in x 2.5in)
# ---------------------------------------------------------------------------
combo <- (p_ct | p_fs | p_cond) +
  plot_layout(ncol = 3, guides = "keep") &
  theme(plot.tag = element_text(size = 7, face = "plain"))

message(sprintf("[save] %s", OUT_PDF))
ggsave(
  filename = OUT_PDF,
  plot     = combo,
  width    = 7.2,
  height   = 2.5,
  device   = cairo_pdf
)
message("[done]")
