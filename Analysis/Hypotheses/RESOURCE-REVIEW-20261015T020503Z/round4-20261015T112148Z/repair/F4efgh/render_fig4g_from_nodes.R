# Copy of figures/candidates/continuum-774-panels-20261011T152106Z/fig4_compact_layout/render_G_774.R
# (review round 3 repair, F4efgh). Only changes: the node table is passed as argument 1
# (fig4g_nodes_<roster>.tsv.gz from build_fig4g_disease_only_nodes.py) and the output
# directory and file name as arguments 2 and 3. Geometry, labels, colours, limits rule,
# size (220 x 176 pt) and seed are unchanged.
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrepel)})
grDevices::pdf.options(useDingbats=FALSE)
theme_systems <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(size = 6, face = "plain"),
      axis.text = element_text(size = 6, color = "black", face = "plain"),
      axis.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      legend.title = element_text(size = 6, face = "plain"),
      strip.text = element_text(size = 6, face = "plain"),
      strip.background = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_text(size = 6, face = "plain", hjust = 0),
      plot.caption = element_text(size = 6, face = "plain", hjust = 0),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.key.size = grid::unit(2.5, "mm"),
      panel.spacing = grid::unit(0.9, "mm"),
      plot.margin = margin(1.5, 1.5, 1.5, 1.5)
    )
}

system_short_label <- function(system_id) sprintf('S%02d', as.integer(sub('system_', '', system_id)))
compact_numeric_labels <- function(values) formatC(values, format='f', digits=1)
contract <- list(seed=20260817)
nodes <- fread(args[1])
communities <- fread("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/continuum-774-panels-20261011T152106Z/molecular_systems_atlas/geometry/community_registry.tsv")
system_nodes <- nodes[component_role=='main_system_component']
stopifnot(nrow(system_nodes)==5439L, nrow(communities)==43L)
map_labels <- copy(communities)[order(-n_features)][1:15]
map_labels[, map_label := system_short_label(community_id)]
map_label_layer <- function() {
  geom_label_repel(
    data = map_labels,
    aes(centroid_x, centroid_y, label = map_label),
    inherit.aes = FALSE, size = 6 / ggplot2::.pt, color = "black",
    fill = scales::alpha("white", 0.82), label.size = NA,
    label.padding = grid::unit(0.35, "mm"),
    seed = contract$seed, box.padding = 0.11, point.padding = 0.06,
    min.segment.length = 0, segment.size = 0.14, max.overlaps = Inf,
    max.time = 10, max.iter = 100000
  )
}

effect_limits <- function(values) {
  bound <- as.numeric(quantile(
    abs(values[is.finite(values)]), 0.98, names = FALSE, type = 8
  ))
  c(-bound, bound)
}

render_effect_map <- function(effect_column, support_column, legend_title,
                              subtitle, filename) {
  limits <- effect_limits(system_nodes[[effect_column]])
  n_supported <- sum(system_nodes[[support_column]] %in% TRUE)
  plot <- ggplot(system_nodes, aes(umap_1, umap_2)) +
    geom_point(
      data = system_nodes[feature_type == "pathway"],
      aes(color = .data[[effect_column]]), size = 0.43, alpha = 0.78, stroke = 0
    ) +
    geom_point(
      data = system_nodes[feature_type == "hotspot_program"],
      aes(color = .data[[effect_column]]), shape = 18, size = 1.05, alpha = 0.95
    ) +
    map_label_layer() +
    scale_color_gradient2(
      low = "#1565C0", mid = "#F7F7F7", high = "#C9265E", midpoint = 0,
      limits = limits, oob = scales::squish,
      breaks = scales::breaks_pretty(n = 3), labels = compact_numeric_labels,
      name = legend_title
    ) +
    coord_equal() +
    labs(
      x = "Membership UMAP 1 (arbitrary)",
      y = "Membership UMAP 2 (arbitrary)",
      subtitle = subtitle
    ) +
    theme_systems() +
    theme(legend.position = "bottom") +
    guides(color = guide_colorbar(
      title.position = "top", title.hjust = 0.5, display="rectangles",
      barwidth = grid::unit(30, "mm"), barheight = grid::unit(2, "mm")
    ))
  ggsave(file.path(panel_dir, filename), plot, width=220/72, height=176/72, device=cairo_pdf)
}

panel_dir <- args[2]
render_effect_map('continuum_beta', 'continuum_supported',
  'Stage-adjusted continuum beta
(score SD / continuum SD)',
  'Stage-adjusted continuum effects
on the same fixed atlas',
  args[3])
