#!/usr/bin/env Rscript

# ⛔ DEPRECATED 2026-07-08 → moved to _legacy. SUPERSEDED by scripts/figures/progression_cascade.R,
# which now builds the combined NAS + fibrosis cascade (Fig 3E) directly. This script is redundant
# (it also wrote fig3e_progression_cascade.pdf, causing a double-write) — do NOT run.
#
# Figure 3E — compact NAS activity | fibrosis stage cascade.
# Reuses the two standalone cascade builders and composes their four matched
# tracks with shared row labels and one collected legend block.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# Load the component builders without writing their standalone files. Keeping
# each script in a separate environment prevents their p1/p2/p3a/p3b objects
# and data tables from overwriting one another.
load_cascade <- function(script_name) {
  old_write <- Sys.getenv("CASCADE_WRITE_OUTPUT", unset = NA_character_)
  on.exit({
    if (is.na(old_write)) {
      Sys.unsetenv("CASCADE_WRITE_OUTPUT")
    } else {
      Sys.setenv(CASCADE_WRITE_OUTPUT = old_write)
    }
  }, add = TRUE)
  Sys.setenv(CASCADE_WRITE_OUTPUT = "false")
  env <- new.env(parent = globalenv())
  sys.source(file.path(BASE, "scripts/figures", script_name), envir = env)
  env
}

nas <- load_cascade("nas_activity_cascade.R")
fib <- load_cascade("progression_cascade.R")

# Shared legend definitions. Both columns use the same visible labels and
# ordering so patchwork can collect each duplicated guide exactly once.
direction_labels <- c("Upregulated", "Downregulated")
nmf_breaks <- c(
  "Inflammatory-EMT", "Fibrotic-ECM", "Hepatocyte-Metabolic",
  "Unresolved", "Noncoding", "Skeletal-muscle"
)
nmf_colors <- c(
  "Inflammatory-EMT" = cat_palette[3],
  "Fibrotic-ECM" = cat_palette[6],
  "Hepatocyte-Metabolic" = "#BDBDBD",
  "Unresolved" = "#D6D6D6",
  "Noncoding" = "#9E9E9E",
  "Skeletal-muscle" = "#7D7D7D"
)
celltype_breaks <- c(
  "T cells", "Macrophages", "Fibroblasts", "Endothelial cells"
)
celltype_colors <- c(
  "T cells" = ct_palette[["T cells"]],
  "Macrophages" = ct_palette[["Macrophages"]],
  "Fibroblasts" = ct_palette[["Fibroblasts"]],
  "Endothelial cells" = ct_palette[["Endothelial cells"]]
)

# A single horizontal title on the NAS plot labels each complete two-column
# row. The fibrosis plot retains its numeric ticks but has no duplicate title.
shared_y_theme <- theme(
  axis.title.y = element_text(
    size = 6, face = "plain", angle = 0, hjust = 1, vjust = 0.5,
    margin = margin(r = 4)
  ),
  plot.margin = margin(1, 2, 1, 1)
)
right_y_theme <- theme(
  axis.title.y = element_blank(),
  plot.margin = margin(1, 2, 1, 1)
)
legend_theme <- theme(
  legend.position = "right",
  legend.justification = "top",
  legend.box = "vertical",
  legend.box.spacing = unit(0.5, "mm"),
  legend.spacing.y = unit(0.5, "mm"),
  legend.margin = margin(0, 0, 0, 1),
  legend.title = element_text(size = 6, face = "plain"),
  legend.text = element_text(size = 6, face = "plain"),
  legend.key.height = unit(2.2, "mm"),
  legend.key.width = unit(2.2, "mm")
)

drop_text_layers <- function(plot) {
  plot$layers <- Filter(
    function(layer) !inherits(layer$geom, "GeomText"),
    plot$layers
  )
  plot
}

nas_p1 <- nas$p1 +
  scale_fill_manual(
    values = nas$deg_colors,
    breaks = levels(nas$deg_long$direction),
    labels = direction_labels,
    name = NULL
  ) +
  labs(x = NULL, y = "DEGs") +
  shared_y_theme

fib_p1 <- fib$p1 +
  scale_fill_manual(
    values = fib$DEG_COLORS,
    breaks = levels(fib$deg_long$direction),
    labels = direction_labels,
    name = NULL
  ) +
  labs(x = NULL, y = NULL) +
  guides(fill = "none") +
  right_y_theme

nas_p2 <- nas$p2 +
  scale_fill_manual(values = nmf_colors, breaks = nmf_breaks,
                    name = "NMF program") +
  labs(x = NULL, y = "NMF %") +
  shared_y_theme

fib_p2 <- fib$p2 +
  scale_fill_manual(values = nmf_colors, breaks = nmf_breaks,
                    name = "NMF program") +
  labs(x = NULL, y = NULL) +
  right_y_theme

nas_p3a <- drop_text_layers(nas$p3a) +
  labs(x = NULL, y = "Hepatocytes %") +
  shared_y_theme

fib_p3a <- drop_text_layers(fib$p3a) +
  labs(x = NULL, y = NULL) +
  right_y_theme

nas_p3b <- nas$p3b +
  scale_color_manual(values = celltype_colors, breaks = celltype_breaks,
                     name = "Cell type") +
  scale_x_discrete(
    drop = FALSE,
    limits = nas$GROUPS,
    labels = sub("^NAS", "", nas$GROUPS),
    expand = nas$X_EXPAND
  ) +
  labs(x = "NAS activity group", y = "Non-hepatocytes %") +
  shared_y_theme

fib_p3b <- fib$p3b +
  scale_color_manual(values = celltype_colors, breaks = celltype_breaks,
                     name = "Cell type") +
  labs(x = "Fibrosis stage", y = NULL) +
  right_y_theme

# Row-major placement keeps matched NAS/fibrosis tracks adjacent. The 180 x
# 96 mm canvas matches the compact aspect of the supplied layout while all
# text remains exactly 6 pt.
combined <-
  nas_p1 + fib_p1 +
  nas_p2 + fib_p2 +
  nas_p3a + fib_p3a +
  nas_p3b + fib_p3b +
  plot_layout(
    ncol = 2,
    widths = c(1, 1),
    heights = c(1.0, 1.0, 0.68, 0.82),
    guides = "collect"
  ) &
  legend_theme

out <- file.path(PANEL_DIR, "fig3e_progression_cascade.pdf")
save_fig(combined, out, width = 180 / 25.4, height = 96 / 25.4)
cat("Wrote combined NAS/fibrosis cascade:", out, "\n")
