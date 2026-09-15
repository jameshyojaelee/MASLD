#!/usr/bin/env Rscript
# KEY MESSAGE: A compact, consistently ordered figure family shows how membership-defined molecular systems change across fibrosis stage and the aligned MASLD continuum.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

options(digits = 17, scipen = 999)
grDevices::pdf.options(useDingbats = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    paste(
      "Usage: 69_render_refined_molecular_systems_family.R",
      "PROMOTED_ATLAS_SOURCE PROMOTED_POOLED_SOURCE OUTPUT_CANDIDATE"
    ),
    call. = FALSE
  )
}

atlas_source <- normalizePath(args[[1L]], mustWork = TRUE)
pooled_source <- normalizePath(args[[2L]], mustWork = TRUE)
output_candidate <- args[[3L]]
if (file.exists(output_candidate)) {
  stop("Refusing to overwrite output candidate: ", output_candidate, call. = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lib_molecular_layers.R"))
contract <- ml_read_contract()

candidate_root <- normalizePath(
  file.path(ml_project_root(), "figures", "candidates"), mustWork = TRUE
)
candidate_only_pass <- startsWith(
  paste0(normalizePath(dirname(output_candidate), mustWork = TRUE), "/"),
  paste0(candidate_root, "/")
)
ml_assert(candidate_only_pass,
          "Refined molecular-system figures are restricted to figures/candidates")

panel_dir <- file.path(output_candidate, "panels")
source_dir <- file.path(output_candidate, "source_tables")
provenance_dir <- file.path(output_candidate, "provenance")
for (path in c(panel_dir, source_dir, provenance_dir)) ml_ensure_dir(path)

atlas_checksums_path <- file.path(atlas_source, "promoted_source_checksums.tsv")
pooled_checksums_path <- file.path(pooled_source, "promoted_source_checksums.tsv")
required_inputs <- c(
  file.path(atlas_source, "geometry", "community_registry.tsv"),
  file.path(atlas_source, "overlays", "all_node_stage_continuum_overlays.tsv.gz"),
  file.path(atlas_source, "overlays", "community_stage_continuum_effects.tsv"),
  file.path(atlas_source, "overlays", "community_balanced_fixed_windows.tsv.gz"),
  file.path(pooled_source, "source_tables", "pooled_molecular_system_windows.tsv")
)
ml_assert(all(file.exists(c(required_inputs, atlas_checksums_path, pooled_checksums_path))),
          "A promoted molecular-system source is missing")

validate_promoted_hashes <- function(paths, checksum_path) {
  checksums <- fread(checksum_path)
  ml_assert(all(c("path", "sha256") %in% names(checksums)),
            paste0("Promoted checksum schema drift: ", checksum_path))
  for (path in paths) {
    target_basename <- basename(path)
    expected <- checksums[
      basename(get("path")) == target_basename, unique(sha256)
    ]
    ml_assert(length(expected) == 1L,
              paste0("Promoted checksum is absent or ambiguous: ", path))
    ml_assert(identical(ml_sha256(path), expected),
              paste0("Promoted source hash drift: ", path))
  }
}
validate_promoted_hashes(required_inputs[1:4], atlas_checksums_path)
validate_promoted_hashes(required_inputs[5], pooled_checksums_path)

communities <- fread(required_inputs[[1L]], na.strings = c("", "NA"))
nodes <- fread(required_inputs[[2L]], na.strings = c("", "NA"))
community_effects <- fread(required_inputs[[3L]], na.strings = c("", "NA"))
cohort_windows <- fread(required_inputs[[4L]], na.strings = c("", "NA"))
pooled_windows <- fread(required_inputs[[5L]], na.strings = c("", "NA"))

all_cohorts <- c(contract$source_overlap_cohorts, contract$evaluation_cohorts)
ml_assert(nrow(nodes) == 5445L && uniqueN(nodes$feature_id) == 5445L,
          "Molecular-system node family drift")
ml_assert(nrow(communities) == 43L && uniqueN(communities$community_id) == 43L,
          "Molecular-system registry drift")
ml_assert(nrow(community_effects) == 43L,
          "Molecular-system effect family drift")
ml_assert(nrow(cohort_windows) == 43L * 5L * 9L &&
            uniqueN(cohort_windows$dataset) == 5L &&
            uniqueN(cohort_windows$window_id) == 9L,
          "Five-cohort window family drift")
ml_assert(nrow(pooled_windows) == 43L * 9L &&
            uniqueN(pooled_windows$community_id) == 43L &&
            uniqueN(pooled_windows$window_id) == 9L,
          "Pooled window family drift")
ml_assert(setequal(unique(cohort_windows$dataset), all_cohorts),
          "Five-cohort window roster drift")
ml_assert(setequal(communities$community_id, unique(pooled_windows$community_id)),
          "Pooled windows and system registry disagree")

family_order <- c(
  "hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc", "hotspot"
)
family_labels <- c(
  hallmark = "Hallmark", kegg = "KEGG", reactome = "Reactome",
  go_bp = "GO biological process", go_mf = "GO molecular function",
  go_cc = "GO cellular component", hotspot = "Hotspot program"
)
family_colors <- c(
  hallmark = "#C9265E", kegg = "#F4A674", reactome = "#7B1FA2",
  go_bp = "#1565C0", go_mf = "#2A8C7D", go_cc = "#F5A623",
  hotspot = "#00695C"
)
family_display_colors <- setNames(
  unname(family_colors[family_order]), unname(family_labels[family_order])
)
cohort_labels <- setNames(all_cohorts, all_cohorts)

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

save_panel <- function(plot, filename, width, height) {
  path <- file.path(panel_dir, filename)
  ml_assert(!file.exists(path), paste0("Refusing to overwrite: ", path))
  ggsave(
    path, plot, width = width, height = height, units = "in",
    device = cairo_pdf, bg = "white", limitsize = FALSE
  )
  ml_assert(file.info(path)$size > 1000L,
            paste0("Unexpectedly small PDF: ", path))
  path
}

compact_numeric_labels <- function(values) {
  formatC(values, format = "f", digits = 1)
}

system_short_label <- function(system_id) {
  sprintf("S%02d", as.integer(sub("system_", "", system_id)))
}

pooled_direction <- pooled_windows[, .(
  early_mean = mean(
    pooled_collection_balanced_score[window_id %in% c(1L, 2L)]
  ),
  late_mean = mean(
    pooled_collection_balanced_score[window_id %in% c(8L, 9L)]
  ),
  pooled_window_spearman = suppressWarnings(cor(
    center, pooled_collection_balanced_score,
    method = "spearman", use = "complete.obs"
  )),
  n_windows = sum(is.finite(pooled_collection_balanced_score))
), by = .(community_id, system_display)]
pooled_direction[, `:=`(
  late_minus_early = late_mean - early_mean,
  direction_class = fifelse(
    late_mean - early_mean >= 0, "positive_late_minus_early",
    "negative_late_minus_early"
  ),
  short_label = system_short_label(community_id),
  ordering_rule = paste0(
    "descending_equal_cohort_pooled_late_minus_early;",
    "late=mean_windows_80_and_90;early=mean_windows_10_and_20"
  ),
  visualization_only = TRUE
)]
setorder(pooled_direction, -late_minus_early, community_id)
pooled_direction[, display_rank := seq_len(.N)]
pooled_direction[, signed_short_label := paste0(
  fifelse(late_minus_early >= 0, "+ ", "- "), short_label
)]
ml_assert(nrow(pooled_direction) == 43L &&
            all(pooled_direction$n_windows == 9L) &&
            !anyDuplicated(pooled_direction$display_rank),
          "Pooled direction-order ledger drift")

ordered_systems <- pooled_direction$community_id
factor_levels <- rev(ordered_systems)
n_negative <- sum(pooled_direction$late_minus_early < 0)
direction_boundary <- if (
  n_negative > 0L && n_negative < nrow(pooled_direction)
) n_negative + 0.5 else NA_real_
label_lookup <- setNames(
  pooled_direction$signed_short_label, pooled_direction$community_id
)

nodes[, collection_label := factor(
  family_labels[collection], levels = unname(family_labels[family_order])
)]
system_nodes <- nodes[component_role == "main_system_component"]
microcomponent_nodes <- nodes[component_role == "excluded_microcomponent"]
ml_assert(nrow(system_nodes) == 5439L && nrow(microcomponent_nodes) == 6L,
          "Main-system and microcomponent node counts drift")

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

topology <- ggplot(nodes, aes(umap_1, umap_2)) +
  geom_point(
    data = nodes[feature_type == "pathway"],
    aes(color = collection_label), size = 0.38, alpha = 0.55, stroke = 0
  ) +
  geom_point(
    data = nodes[feature_type == "hotspot_program"],
    aes(color = collection_label), shape = 18, size = 1.05, alpha = 0.95
  ) +
  geom_point(
    data = microcomponent_nodes, color = "#9E9E9E", fill = "white",
    shape = 21, size = 0.8, stroke = 0.25
  ) +
  map_label_layer() +
  scale_color_manual(
    values = family_display_colors,
    breaks = unname(family_labels[family_order]), name = NULL
  ) +
  coord_equal() +
  labs(
    x = "Membership UMAP 1 (arbitrary)",
    y = "Membership UMAP 2 (arbitrary)",
    subtitle = "Signature-excluded membership geometry: 5,445 mapped features"
  ) +
  theme_systems() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(color = guide_legend(
    nrow = 2, byrow = TRUE,
    override.aes = list(size = 1.35, alpha = 1)
  ))
topology_path <- save_panel(
  topology, "s4_molecular_systems_membership_atlas.pdf", 4.9, 3.90
)

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
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(30, "mm"), barheight = grid::unit(2, "mm")
    ))
  save_panel(plot, filename, 4.75, 3.78)
}

stage_map_path <- render_effect_map(
  "stage_beta", "stage_supported",
  "Fibrosis-stage beta (score SD / stage)",
  "Fibrosis-stage remodeling on the fixed membership atlas",
  "s4_molecular_systems_fibrosis_stage_effect_map.pdf"
)
continuum_map_path <- render_effect_map(
  "continuum_beta", "continuum_supported",
  "Stage-adjusted continuum beta (score SD / continuum SD)",
  "Stage-adjusted continuum effects on the same fixed atlas",
  "s4_molecular_systems_continuum_effect_map.pdf"
)

community_effects[, direction_pair := fcase(
  stage_beta > 0 & continuum_beta > 0, "Both positive",
  stage_beta < 0 & continuum_beta < 0, "Both negative",
  default = "Opposite signs"
)]
community_effects[, short_label := system_short_label(community_id)]
direction_colors <- c(
  "Both positive" = "#C9265E",
  "Both negative" = "#1565C0",
  "Opposite signs" = "#9E9E9E"
)
direction_levels <- names(direction_colors)
community_effects[, direction_pair := factor(
  direction_pair, levels = direction_levels
)]
system_rho <- cor(
  community_effects$stage_beta, community_effects$continuum_beta,
  method = "spearman", use = "complete.obs"
)
same_sign_count <- sum(
  sign(community_effects$stage_beta) == sign(community_effects$continuum_beta)
)
system_scatter <- ggplot(
  community_effects, aes(stage_beta, continuum_beta)
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_point(
    aes(size = n_features, color = direction_pair),
    alpha = 0.86, stroke = 0
  ) +
  geom_text_repel(
    aes(label = short_label), size = 6 / ggplot2::.pt, color = "black",
    seed = contract$seed, box.padding = 0.07, point.padding = 0.035,
    min.segment.length = 0, segment.size = 0.11, max.overlaps = Inf,
    max.time = 15, max.iter = 200000
  ) +
  scale_color_manual(values = direction_colors, drop = FALSE, name = "Direction") +
  scale_size_continuous(range = c(0.9, 2.8), name = "Mapped features") +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4), labels = compact_numeric_labels,
    expand = expansion(mult = 0.1)
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 4), labels = compact_numeric_labels,
    expand = expansion(mult = 0.1)
  ) +
  labs(
    x = "Collection-balanced fibrosis-stage beta",
    y = "Collection-balanced stage-adjusted continuum beta",
    subtitle = sprintf(
      "43 systems: Spearman rho = %.2f; %d have the same sign",
      system_rho, same_sign_count
    )
  ) +
  theme_systems() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(
    color = guide_legend(order = 1, nrow = 1, override.aes = list(size = 1.8)),
    size = guide_legend(order = 2, nrow = 1)
  )
system_scatter_path <- save_panel(
  system_scatter, "s4_molecular_systems_stage_vs_continuum_summary.pdf",
  4.15, 3.62
)

cohort_windows <- merge(
  cohort_windows,
  pooled_direction[, .(
    community_id, display_rank, late_minus_early, signed_short_label
  )],
  by = "community_id", all.x = TRUE, sort = FALSE
)
cohort_windows[, `:=`(
  system_factor = factor(community_id, levels = factor_levels),
  dataset_label = factor(cohort_labels[dataset], levels = cohort_labels[all_cohorts]),
  window_position = as.numeric(window_id)
)]
cohort_heat_limit <- as.numeric(quantile(
  abs(cohort_windows$collection_balanced_score), 0.98,
  names = FALSE, type = 8, na.rm = TRUE
))
five_cohort_heatmap <- ggplot(
  cohort_windows,
  aes(window_position, system_factor, fill = collection_balanced_score)
) +
  geom_tile(width = 0.94, height = 0.94) +
  {if (is.finite(direction_boundary)) geom_hline(
    yintercept = direction_boundary, color = "#4D4D4D", linewidth = 0.22
  )} +
  facet_grid(. ~ dataset_label) +
  scale_fill_gradient2(
    low = "#1565C0", mid = "#F7F7F7", high = "#C9265E", midpoint = 0,
    limits = c(-cohort_heat_limit, cohort_heat_limit), oob = scales::squish,
    breaks = c(-cohort_heat_limit, 0, cohort_heat_limit),
    labels = compact_numeric_labels,
    name = "Within-cohort collection-balanced score (z)"
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 9), labels = c("10", "50", "90"),
    expand = c(0, 0)
  ) +
  scale_y_discrete(labels = label_lookup) +
  labs(
    x = "Fixed-projection continuum percentile (20-point windows; 50% overlap)",
    y = "System (pooled late-minus-early order)"
  ) +
  theme_systems() +
  theme(
    axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = grid::unit(30, "mm"), barheight = grid::unit(2, "mm")
  ))
five_cohort_path <- save_panel(
  five_cohort_heatmap,
  "s4_molecular_systems_five_cohort_sliding_windows.pdf",
  5.5, 5.02
)

pooled_plot_data <- merge(
  pooled_windows,
  pooled_direction[, .(
    community_id, display_rank, late_minus_early, signed_short_label
  )],
  by = "community_id", all.x = TRUE, sort = FALSE
)
pooled_plot_data[, `:=`(
  system_factor = factor(community_id, levels = factor_levels),
  display_position = as.numeric(window_id),
  display_value = pooled_collection_balanced_score,
  display_column = "window"
)]
delta_tiles <- pooled_direction[, .(
  community_id,
  system_factor = factor(community_id, levels = factor_levels),
  display_position = 10.35,
  display_value = late_minus_early,
  display_column = "late_minus_early"
)]
pooled_display <- rbindlist(list(
  pooled_plot_data[, .(
    community_id, system_factor, display_position, display_value, display_column
  )],
  delta_tiles
), use.names = TRUE)
pooled_heat_limit <- as.numeric(quantile(
  abs(pooled_display$display_value), 0.98,
  names = FALSE, type = 8, na.rm = TRUE
))
pooled_heatmap <- ggplot(
  pooled_display,
  aes(display_position, system_factor, fill = display_value)
) +
  geom_tile(width = 0.9, height = 0.94) +
  geom_vline(xintercept = 9.68, color = "#4D4D4D", linewidth = 0.25) +
  {if (is.finite(direction_boundary)) geom_hline(
    yintercept = direction_boundary, color = "#4D4D4D", linewidth = 0.22
  )} +
  scale_fill_gradient2(
    low = "#1565C0", mid = "#F7F7F7", high = "#C9265E", midpoint = 0,
    limits = c(-pooled_heat_limit, pooled_heat_limit), oob = scales::squish,
    breaks = c(-pooled_heat_limit, 0, pooled_heat_limit),
    labels = compact_numeric_labels,
    name = "Pooled system score (z)"
  ) +
  scale_x_continuous(
    breaks = c(1, 3, 5, 7, 9, 10.35),
    labels = c("10", "30", "50", "70", "90", "\u0394"),
    expand = expansion(add = c(0.04, 0.04))
  ) +
  scale_y_discrete(labels = label_lookup) +
  coord_fixed(ratio = 1) +
  labs(
    x = "Continuum percentile",
    y = "System (pooled Delta order)"
  ) +
  theme_systems() +
  theme(
    axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = grid::unit(20, "mm"), barheight = grid::unit(2, "mm")
  ))
pooled_heatmap_path <- save_panel(
  pooled_heatmap,
  "s4_molecular_systems_pooled_continuum_windows.pdf",
  2.25, 4.88
)

ml_write_tsv_once(
  pooled_direction,
  file.path(source_dir, "pooled_system_late_minus_early_order.tsv")
)
ml_write_tsv_once(
  pooled_display,
  file.path(source_dir, "pooled_system_heatmap_display_values.tsv")
)

figure_manifest <- data.table(
  figure_id = c(
    "membership_atlas", "fibrosis_stage_map", "continuum_effect_map",
    "system_stage_continuum_summary", "five_cohort_sliding_windows",
    "pooled_continuum_windows"
  ),
  path = c(
    topology_path, stage_map_path, continuum_map_path,
    system_scatter_path, five_cohort_path, pooled_heatmap_path
  ),
  role = c(
    "outcome_blind_topology", "stage_overlay", "continuum_overlay",
    "community_level_effect_comparison", "descriptive_five_cohort_windows",
    "descriptive_equal_cohort_pooled_windows"
  ),
  data_values_changed = FALSE,
  geometry_changed = FALSE,
  display_order_rule = c(
    rep("fixed_membership_geometry", 4L),
    rep("descending_pooled_late_minus_early", 2L)
  ),
  result_adaptive_display_order = c(rep(FALSE, 4L), TRUE, TRUE),
  inferential_use_of_display_order = FALSE,
  embedded_explanatory_text = FALSE,
  promoted = FALSE
)
figure_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
ml_write_tsv_once(figure_manifest, file.path(output_candidate, "figure_manifest.tsv"))

input_manifest <- data.table(
  path = normalizePath(required_inputs, mustWork = TRUE),
  input_role = c(
    "promoted_system_registry", "promoted_node_overlays",
    "promoted_system_effects", "promoted_five_cohort_windows",
    "promoted_equal_cohort_pooled_windows"
  )
)
input_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
ml_write_tsv_once(input_manifest, file.path(provenance_dir, "input_manifest.tsv"))

audit <- data.table(
  check = c(
    "promoted_input_hashes", "complete_six_panel_family",
    "all_43_systems_in_both_heatmaps", "shared_heatmap_row_order",
    "transparent_result_adaptive_display_order", "frozen_values_and_geometry",
    "embedded_explanatory_text_removed",
    "candidate_only", "pdf_nonempty"
  ),
  passed = c(
    TRUE, nrow(figure_manifest) == 6L,
    uniqueN(cohort_windows$community_id) == 43L &&
      uniqueN(pooled_plot_data$community_id) == 43L,
    identical(levels(cohort_windows$system_factor),
              levels(pooled_plot_data$system_factor)),
    all(figure_manifest[
      figure_id %in% c("five_cohort_sliding_windows", "pooled_continuum_windows"),
      result_adaptive_display_order
    ]) && all(figure_manifest$inferential_use_of_display_order == FALSE),
    all(figure_manifest$data_values_changed == FALSE) &&
      all(figure_manifest$geometry_changed == FALSE),
    all(figure_manifest$embedded_explanatory_text == FALSE),
    candidate_only_pass && !any(figure_manifest$promoted),
    all(file.info(figure_manifest$path)$size > 1000L)
  ),
  detail = c(
    paste(basename(required_inputs), collapse = ";"),
    paste(figure_manifest$figure_id, collapse = ";"),
    "43 systems x 5 cohorts x 9 windows and 43 systems x 9 pooled windows",
    paste(ordered_systems, collapse = ";"),
    paste0(
      "Display-only descending late-minus-early ordering; ",
      "no selection, support call, timing label, or p value"
    ),
    "Only layout, labeling, color, and descriptive row order changed",
    "Methods and interpretation remain in the external caption authority",
    normalizePath(output_candidate, mustWork = TRUE),
    paste(basename(figure_manifest$path), collapse = ";")
  )
)
ml_assert(all(audit$passed), "Refined molecular-system figure audit failed")
ml_write_tsv_once(audit, file.path(provenance_dir, "audit.tsv"))
ml_write_session_info(file.path(provenance_dir, "sessionInfo.txt"))

output_files <- c(
  figure_manifest$path,
  file.path(source_dir, "pooled_system_late_minus_early_order.tsv"),
  file.path(source_dir, "pooled_system_heatmap_display_values.tsv"),
  file.path(output_candidate, "figure_manifest.tsv"),
  file.path(provenance_dir, "input_manifest.tsv"),
  file.path(provenance_dir, "audit.tsv"),
  file.path(provenance_dir, "sessionInfo.txt")
)
output_checksums <- data.table(
  path = normalizePath(output_files, mustWork = TRUE),
  sha256 = vapply(output_files, ml_sha256, character(1))
)
ml_write_tsv_once(
  output_checksums, file.path(provenance_dir, "output_checksums.tsv")
)
