#!/usr/bin/env Rscript
# Render the endpoint PCA from fixed coordinates. Endpoint definitions remain
# in source metadata; the reader-facing control legend is intentionally terse.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
grDevices::pdf.options(useDingbats = FALSE)

render_endpoint_pca <- function(coords, pve, out_path) {
  base <- Sys.getenv("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
  if (!exists("theme_masld") || !exists("masld_colors") || !exists("red_gradient")) {
    source(file.path(base, "scripts/figures/publication_theme.R"))
  }
  coords <- as.data.table(coords)
  coords[, endpoint := factor(endpoint,
    levels = c("Strict control", "Advanced disease"))]
  if (!"endpoint_code" %in% names(coords)) {
    coords[, endpoint_code := fifelse(endpoint == "Strict control",
                                      "StrictControl", "Advanced")]
  }
  fibrosis_ramp <- colorRampPalette(c(red_gradient[[1L]], masld_colors$mash))(5)
  density <- ggplot(coords, aes(PC1, color = endpoint, fill = endpoint)) +
    geom_density(alpha = 0.25, linewidth = 0.4) +
    scale_color_manual(values = c(`Strict control` = "#9E9E9E",
                                  `Advanced disease` = masld_colors$mash), guide = "none") +
    scale_fill_manual(values = c(`Strict control` = "#9E9E9E",
                                 `Advanced disease` = masld_colors$mash), guide = "none") +
    labs(x = NULL, y = NULL) + theme_masld(base_size = 6) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          plot.margin = margin(1, 1, 0, 1))
  scatter <- ggplot() +
    geom_point(data = coords[endpoint_code == "StrictControl"],
               aes(PC1, PC2, shape = "Control"), fill = "white",
               color = "#9E9E9E", size = 0.85, stroke = 0.3, alpha = 0.6) +
    geom_point(data = coords[endpoint_code == "Advanced"],
               aes(PC1, PC2, color = fibrosis_stage), size = 0.85, alpha = 0.6) +
    scale_color_gradientn(colors = fibrosis_ramp, limits = c(0, 4), name = "Fibrosis") +
    scale_shape_manual(name = NULL, values = c(Control = 1)) +
    guides(color = guide_colorbar(order = 1, barwidth = 0.4, barheight = 2.4),
           shape = guide_legend(order = 2, override.aes = list(size = 1.6))) +
    labs(x = sprintf("PC1 (%.1f%%)", pve[[1L]]),
         y = sprintf("PC2 (%.1f%%)", pve[[2L]])) +
    theme_masld(base_size = 6) +
    theme(legend.position = "right", axis.text = element_text(size = 6, face = "plain"),
          legend.title = element_text(size = 6, face = "plain"),
          legend.text = element_text(size = 6, face = "plain"),
          legend.key.size = unit(0.22, "cm"))
  panel <- (density / scatter) + plot_layout(heights = c(1, 4))
  ggsave(out_path, panel, width = 2.42, height = 2.23, device = cairo_pdf)
  invisible(out_path)
}

if (sys.nframe() == 0L) {
  candidate_root <- normalizePath(Sys.getenv("FIGURE_CANDIDATE_ROOT", ""), mustWork = TRUE)
  coords <- fread(file.path(candidate_root, "source_tables",
                            "fig3d_endpoint_pca_coordinates.tsv"))
  pve <- c(as.numeric(Sys.getenv("FIG3D_PC1_PVE", "7.6")),
           as.numeric(Sys.getenv("FIG3D_PC2_PVE", "3.9")))
  render_endpoint_pca(coords, pve,
    file.path(candidate_root, "figure3", "panels", "fig3d_pca_fibrosis_gradient.pdf"))
}
