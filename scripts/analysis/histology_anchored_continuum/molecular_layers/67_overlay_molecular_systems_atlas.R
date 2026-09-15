#!/usr/bin/env Rscript
# KEY MESSAGE: The same outcome-blind molecular systems reveal where fibrosis-stage remodeling, stage-adjusted continuum effects, and sliding-window trajectories converge across MASLD biopsies.

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
      "Usage: 67_overlay_molecular_systems_atlas.R",
      "ANALYSIS_CANDIDATE COMPARISON_CANDIDATE OUTPUT_CANDIDATE"
    ),
    call. = FALSE
  )
}
analysis_candidate <- normalizePath(args[[1L]], mustWork = TRUE)
comparison_candidate <- normalizePath(args[[2L]], mustWork = TRUE)
output_candidate <- normalizePath(args[[3L]], mustWork = TRUE)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lib_molecular_layers.R"))
contract <- ml_read_contract()

candidate_root <- normalizePath(
  file.path(ml_project_root(), "figures", "candidates"), mustWork = TRUE
)
candidate_only_pass <- startsWith(
  paste0(output_candidate, "/"), paste0(candidate_root, "/")
)
ml_assert(candidate_only_pass,
          "Molecular-system rendering is restricted to figures/candidates")

geometry_ok <- file.path(output_candidate, "geometry", "GEOMETRY_FROZEN.ok")
ml_assert(file.exists(geometry_ok), "Outcome-blind geometry is not frozen")
panel_dir <- file.path(output_candidate, "panels")
overlay_dir <- file.path(output_candidate, "overlays")
provenance_dir <- file.path(output_candidate, "provenance")
ml_assert(!dir.exists(panel_dir) && !dir.exists(overlay_dir),
          "Refusing to overwrite molecular-system overlays")
ml_ensure_dir(panel_dir)
ml_ensure_dir(overlay_dir)

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
cohort_labels <- c(
  GSE126848 = "GSE126848", GSE130970 = "GSE130970",
  GSE135251 = "GSE135251", GSE162694 = "GSE162694",
  GSE213621 = "GSE213621"
)
all_cohorts <- c(contract$source_overlap_cohorts, contract$evaluation_cohorts)

theme_atlas <- function() {
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
      legend.key.size = grid::unit(2.8, "mm"),
      panel.spacing = grid::unit(1.2, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

save_panel <- function(plot, filename, width, height) {
  path <- file.path(panel_dir, filename)
  ml_assert(!file.exists(path), paste0("Refusing to overwrite: ", path))
  ggsave(path, plot, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white", limitsize = FALSE)
  ml_assert(file.info(path)$size > 1000L, paste0("Unexpectedly small PDF: ", path))
  path
}

compact_numeric_labels <- function(values) {
  formatC(values, format = "f", digits = 1)
}

coordinates_path <- file.path(output_candidate, "geometry", "latent_coordinates.tsv.gz")
community_path <- file.path(output_candidate, "geometry", "community_registry.tsv")
geometry_checksum_path <- file.path(
  output_candidate, "provenance", "geometry_output_checksums.tsv"
)
geometry_checksums <- fread(geometry_checksum_path)
expected_geometry_hash <- geometry_checksums[
  basename(path) == basename(coordinates_path), unique(sha256)
]
ml_assert(length(expected_geometry_hash) == 1L,
          "Geometry checksum manifest lacks one coordinate checksum")
geometry_hash_at_overlay_entry <- ml_sha256(coordinates_path)
ml_assert(identical(geometry_hash_at_overlay_entry, expected_geometry_hash),
          "Geometry differs from the serial pre-outcome checksum")
coordinates <- fread(coordinates_path, na.strings = c("", "NA"))
communities <- fread(community_path, na.strings = c("", "NA"))
complete_ledger <- fread(file.path(
  output_candidate, "source_tables", "complete_feature_testability_ledger.tsv.gz"
))
ml_assert(nrow(coordinates) == 5445L && uniqueN(coordinates$feature_id) == 5445L,
          "Frozen molecular-system coordinate count drift")
ml_assert(sum(coordinates$component_role == "main_system_component") == 5439L &&
            sum(coordinates$component_role == "excluded_microcomponent") == 6L,
          "Main-component and microcomponent node counts drifted")
ml_assert(nrow(communities) == 43L,
          "Expected 43 main-component molecular systems")
ml_assert(all(family_order %in% coordinates$collection),
          "Molecular-system collection coverage drift")

pathway_effect_path <- file.path(
  comparison_candidate, "source_tables", "pathway_stage_vs_continuum.tsv.gz"
)
hotspot_effect_path <- file.path(
  comparison_candidate, "source_tables", "hotspot_stage_vs_continuum.tsv"
)
ml_assert(all(file.exists(c(pathway_effect_path, hotspot_effect_path))),
          "Missing frozen stage-versus-continuum coordinates")
pathway_effects <- fread(pathway_effect_path, na.strings = c("", "NA"))
setnames(pathway_effects, "set_id", "feature_id")
pathway_effects[, continuum_supported_normalized := continuum_supported %in% TRUE]
hotspot_effects <- fread(hotspot_effect_path, na.strings = c("", "NA"))
setnames(hotspot_effects, "program_uid", "feature_id")
hotspot_effects[, `:=`(
  collection = "hotspot",
  continuum_supported_normalized = continuum_associated %in% TRUE
)]
effects <- rbindlist(list(
  pathway_effects[, .(
    feature_id, collection, stage_beta = beta, stage_q = q_value,
    stage_supported = stage_supported %in% TRUE,
    continuum_beta, continuum_q,
    continuum_supported = continuum_supported_normalized,
    support_class = as.character(support_class)
  )],
  hotspot_effects[, .(
    feature_id, collection, stage_beta = beta, stage_q = q_value,
    stage_supported = stage_supported %in% TRUE,
    continuum_beta, continuum_q,
    continuum_supported = continuum_supported_normalized,
    support_class = as.character(support_class)
  )]
), use.names = TRUE, fill = TRUE)
ml_assert(!anyDuplicated(effects$feature_id), "Effect feature IDs are duplicated")
node_results <- merge(coordinates, effects, by = c("feature_id", "collection"),
                      all.x = TRUE, sort = FALSE)
ml_assert(nrow(node_results) == nrow(coordinates), "Effect overlay changed node count")
ml_assert(all(is.finite(node_results$stage_beta) &
                is.finite(node_results$continuum_beta)),
          "A geometry-eligible node lacks a stage or continuum estimate")
node_results[, collection_label := factor(
  family_labels[collection], levels = unname(family_labels[family_order])
)]
system_node_results <- node_results[component_role == "main_system_component"]
microcomponent_results <- node_results[component_role == "excluded_microcomponent"]

canonical_coordinate_path <- file.path(
  overlay_dir, "canonical_atlas_coordinates.tsv"
)
canonical_coordinates <- coordinates[, .(
  feature_id, umap_1, umap_2, component_id, component_role,
  microcomponent_id, community_id
)]
ml_write_tsv_once(canonical_coordinates, canonical_coordinate_path)
canonical_coordinate_hash <- ml_sha256(canonical_coordinate_path)

stage_map_source <- system_node_results[, .(
  feature_id, collection, feature_type, community_id, umap_1, umap_2,
  effect_name = "fibrosis_stage_beta", effect_beta = stage_beta,
  complete_family_supported = stage_supported
)]
continuum_map_source <- system_node_results[, .(
  feature_id, collection, feature_type, community_id, umap_1, umap_2,
  effect_name = "stage_adjusted_fixed_projection_beta",
  effect_beta = continuum_beta,
  complete_family_supported = continuum_supported
)]
stage_map_source_path <- file.path(overlay_dir, "stage_map_node_source.tsv.gz")
continuum_map_source_path <- file.path(overlay_dir, "continuum_map_node_source.tsv.gz")
ml_write_tsv_once(stage_map_source, stage_map_source_path)
ml_write_tsv_once(continuum_map_source, continuum_map_source_path)

community_labels <- communities[order(-n_features)]
maximum_map_labels <- 15L
community_labels[, label_on_map := seq_len(.N) <= maximum_map_labels]
map_labels <- community_labels[label_on_map == TRUE]
map_labels[, map_label := sub("system_", "S", community_id)]

topology <- ggplot(node_results, aes(umap_1, umap_2)) +
  geom_point(
    data = node_results[feature_type == "pathway"],
    aes(color = collection_label), size = 0.36, alpha = 0.48, stroke = 0
  ) +
  geom_point(
    data = node_results[feature_type == "hotspot_program"],
    aes(color = collection_label), shape = 18, size = 1.15, alpha = 0.95
  ) +
  geom_point(
    data = microcomponent_results, color = "#9E9E9E", fill = "white",
    shape = 21, size = 0.8, stroke = 0.25
  ) +
  geom_text_repel(
    data = map_labels,
    aes(centroid_x, centroid_y, label = map_label),
    inherit.aes = FALSE, size = 6 / ggplot2::.pt, color = "black",
    seed = contract$seed, box.padding = 0.12, point.padding = 0.08,
    min.segment.length = 0, segment.size = 0.15, max.overlaps = Inf,
    max.time = 10, max.iter = 100000
  ) +
  scale_color_manual(
    values = family_display_colors,
    breaks = unname(family_labels[family_order]), name = NULL
  ) +
  coord_equal() +
  labs(
    x = "Membership UMAP 1 (arbitrary)",
    y = "Membership UMAP 2 (arbitrary)",
    subtitle = sprintf(
      "%s pathways and %s Hotspot programs; %s main-component systems",
      format(sum(node_results$feature_type == "pathway"), big.mark = ","),
      format(sum(node_results$feature_type == "hotspot_program"), big.mark = ","),
      nrow(communities)
    ),
    caption = paste0(
      "Coordinates and systems use signature-excluded gene membership only. Six grey nodes",
      "\nin three disconnected microcomponents remain explicit but are not summarized as systems."
    )
  ) +
  theme_atlas() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(color = guide_legend(
    nrow = 2, byrow = TRUE,
    override.aes = list(size = 1.4, alpha = 1)
  ))
topology_path <- save_panel(
  topology, "s4_molecular_systems_membership_atlas.pdf", 5.5, 4.8
)

effect_limits <- function(values) {
  bound <- as.numeric(quantile(abs(values[is.finite(values)]), 0.98,
                               names = FALSE, type = 8))
  c(-bound, bound)
}

render_effect_map <- function(effect_column, support_column, legend_title,
                              subtitle, filename) {
  values <- system_node_results[[effect_column]]
  limits <- effect_limits(values)
  n_supported <- sum(system_node_results[[support_column]] %in% TRUE)
  plot <- ggplot(system_node_results, aes(umap_1, umap_2)) +
    geom_point(
      data = system_node_results[feature_type == "pathway"],
      aes(color = .data[[effect_column]]), size = 0.42, alpha = 0.72, stroke = 0
    ) +
    geom_point(
      data = system_node_results[feature_type == "hotspot_program"],
      aes(color = .data[[effect_column]]), shape = 18, size = 1.05, alpha = 0.9
    ) +
    geom_text_repel(
      data = map_labels,
      aes(centroid_x, centroid_y, label = map_label), inherit.aes = FALSE,
      size = 6 / ggplot2::.pt, color = "black", seed = contract$seed,
      box.padding = 0.12, point.padding = 0.08, min.segment.length = 0,
      segment.size = 0.15, max.overlaps = Inf, max.time = 10,
      max.iter = 100000
    ) +
    scale_color_gradient2(
      low = "#1565C0", mid = "#E5E5E5", high = "#C9265E", midpoint = 0,
      limits = limits, oob = scales::squish,
      breaks = scales::breaks_pretty(n = 3), labels = compact_numeric_labels,
      name = legend_title
    ) +
    coord_equal() +
    labs(
      x = "Membership UMAP 1 (arbitrary)",
      y = "Membership UMAP 2 (arbitrary)",
      subtitle = subtitle,
      caption = paste0(
        format(n_supported, big.mark = ","),
        " displayed nodes pass the declared complete-family support gate; support is tabulated",
        "\nbut not encoded by opacity. Both overlays use one checksummed coordinate table."
      )
    ) +
    theme_atlas() +
    theme(legend.position = "bottom") +
    guides(color = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(32, "mm"), barheight = grid::unit(2, "mm")
    ))
  save_panel(plot, filename, 5.2, 4.55)
}

stage_map_path <- render_effect_map(
  "stage_beta", "stage_supported",
  "Fibrosis-stage trend beta\n(score SD / stage increment)",
  "Fibrosis-stage remodeling over the fixed membership atlas",
  "s4_molecular_systems_fibrosis_stage_effect_map.pdf"
)
continuum_map_path <- render_effect_map(
  "continuum_beta", "continuum_supported",
  "Stage-adjusted continuum beta\n(score SD / continuum SD)",
  "Stage-adjusted fixed-projection effects over the same fixed atlas",
  "s4_molecular_systems_continuum_effect_map.pdf"
)

collection_balanced_summary <- function(data, value_column, output_name) {
  family_medians <- data[!is.na(community_id), .(
    family_median = median(get(value_column), na.rm = TRUE),
    n_features = sum(is.finite(get(value_column)))
  ), by = .(community_id, collection)]
  family_medians[!is.finite(family_median), family_median := NA_real_]
  systems <- family_medians[is.finite(family_median), .(
    community_effect = mean(family_median),
    n_collections = .N,
    n_features = sum(n_features)
  ), by = community_id]
  setnames(systems, "community_effect", output_name)
  list(family = family_medians, system = systems)
}

feature_balanced_summary <- function(data, value_column, output_name) {
  systems <- data[!is.na(community_id) & is.finite(get(value_column)), .(
    community_effect = median(get(value_column)),
    n_features = .N
  ), by = community_id]
  setnames(systems, "community_effect", output_name)
  systems
}

stage_systems <- collection_balanced_summary(node_results, "stage_beta", "stage_beta")
continuum_systems <- collection_balanced_summary(
  node_results, "continuum_beta", "continuum_beta"
)
stage_systems_feature <- feature_balanced_summary(
  node_results, "stage_beta", "stage_beta_feature_balanced"
)
continuum_systems_feature <- feature_balanced_summary(
  node_results, "continuum_beta", "continuum_beta_feature_balanced"
)
community_effects <- Reduce(
  function(x, y) merge(x, y, by = "community_id", all = TRUE),
  list(
    communities,
    stage_systems$system[, .(community_id, stage_beta)],
    continuum_systems$system[, .(community_id, continuum_beta)],
    stage_systems_feature,
    continuum_systems_feature
  )
)
community_effects[, same_direction := sign(stage_beta) == sign(continuum_beta)]
community_effects[, same_direction_feature_balanced :=
                    sign(stage_beta_feature_balanced) ==
                    sign(continuum_beta_feature_balanced)]
community_effects[, figure_label := sub("system_", "S", community_id)]
aggregation_effect_sensitivity <- rbindlist(list(
  community_effects[, .(
    community_id, aggregation_method = "collection_balanced",
    stage_beta, continuum_beta, same_direction
  )],
  community_effects[, .(
    community_id, aggregation_method = "feature_balanced_median",
    stage_beta = stage_beta_feature_balanced,
    continuum_beta = continuum_beta_feature_balanced,
    same_direction = same_direction_feature_balanced
  )]
))
aggregation_effect_summary <- aggregation_effect_sensitivity[, .(
  n_systems = .N,
  stage_continuum_spearman = cor(
    stage_beta, continuum_beta, method = "spearman", use = "complete.obs"
  ),
  same_direction_fraction = mean(same_direction)
), by = aggregation_method]
system_scatter <- ggplot(community_effects, aes(stage_beta, continuum_beta)) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_point(aes(size = n_features, color = same_direction), alpha = 0.8) +
  geom_text_repel(
    aes(label = figure_label), size = 6 / ggplot2::.pt,
    seed = contract$seed, box.padding = 0.08, point.padding = 0.04,
    min.segment.length = 0, segment.size = 0.12, max.overlaps = Inf,
    max.time = 15, max.iter = 200000
  ) +
  scale_color_manual(values = c(`TRUE` = "#00695C", `FALSE` = "#9E9E9E"),
                     name = "Same direction") +
  scale_size_continuous(range = c(1, 3), name = "Mapped features") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 4),
                     labels = compact_numeric_labels,
                     expand = expansion(mult = 0.1)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4),
                     labels = compact_numeric_labels,
                     expand = expansion(mult = 0.1)) +
  labs(
    x = "Collection-balanced fibrosis-stage beta",
    y = "Collection-balanced stage-adjusted continuum beta",
    subtitle = sprintf(
      "Descriptive alignment across outcome-blind systems; Spearman rho = %.2f",
      cor(community_effects$stage_beta, community_effects$continuum_beta,
          method = "spearman", use = "complete.obs")
    ),
    caption = paste0(
      "Each point is one membership-defined system. This descriptive comparison contains",
      "\nnested, correlated gene sets; feature-balanced sensitivity is reported separately."
    )
  ) +
  theme_atlas() +
  theme(legend.position = "bottom")
system_scatter_path <- save_panel(
  system_scatter, "s4_molecular_systems_stage_vs_continuum_summary.pdf", 4.8, 4.4
)

message("Computing fixed sliding windows for all outcome-blind atlas nodes")
axis_path <- file.path(ml_source_root(contract), "projection", "participant_scores.tsv")
axes <- fread(axis_path, na.strings = c("", "NA"))[
  axis_id == "fixed_projection" & dataset %in% all_cohorts,
  .(sample_id, dataset, axis_percentile)
]
ml_assert(!anyDuplicated(axes[, .(sample_id, dataset)]),
          "Fixed-projection participant rows are duplicated")
ml_assert(setequal(unique(axes$dataset), all_cohorts),
          "Fixed projection does not cover all five descriptive cohorts")

score_rows <- list()
for (collection in setdiff(family_order, "hotspot")) {
  collection_id <- collection
  score_path <- file.path(
    analysis_candidate, "pathway_tf", "pathways", collection_id,
    "donor_scores.tsv.gz"
  )
  scores <- fread(score_path, na.strings = c("", "NA"))[
    set_id %in% coordinates[collection == collection_id, feature_id] &
      dataset %in% all_cohorts
  ]
  scores[, `:=`(
    feature_id = set_id, collection = collection_id,
    score_z = as.numeric(pathway_score)
  )]
  score_rows[[collection_id]] <- scores[, .(
    feature_id, collection, sample_id, dataset, score_z
  )]
}
hotspot_score_path <- file.path(
  analysis_candidate, "programs", "hotspot", "hotspot_participant_scores.tsv.gz"
)
hotspot_scores <- fread(hotspot_score_path, na.strings = c("", "NA"))[
  program_uid %in% coordinates[collection == "hotspot", feature_id] &
    dataset %in% all_cohorts
]
hotspot_scores[, `:=`(
  feature_id = program_uid, collection = "hotspot", score_z = as.numeric(outcome_z)
)]
score_rows[["hotspot"]] <- hotspot_scores[, .(
  feature_id, collection, sample_id, dataset, score_z
)]
scores <- rbindlist(score_rows, use.names = TRUE, fill = TRUE)
scores <- merge(scores, axes, by = c("sample_id", "dataset"), all = FALSE,
                sort = FALSE)
ml_assert(all(is.finite(scores$axis_percentile)),
          "Sliding-window score rows lack fixed-projection percentiles")

centers <- as.numeric(contract$window_centers)
window_width <- as.numeric(contract$window_width)
window_rows <- vector("list", length(centers))
for (index in seq_along(centers)) {
  lower <- centers[[index]] - window_width / 2
  upper <- centers[[index]] + window_width / 2
  window_data <- if (index == length(centers)) {
    scores[axis_percentile >= lower & axis_percentile <= upper]
  } else {
    scores[axis_percentile >= lower & axis_percentile < upper]
  }
  window_rows[[index]] <- window_data[, .(
    mean_score = mean(score_z, na.rm = TRUE),
    se_score = sd(score_z, na.rm = TRUE) / sqrt(sum(is.finite(score_z))),
    n_participants = sum(is.finite(score_z))
  ), by = .(feature_id, collection, dataset)][, `:=`(
    window_id = index, center = centers[[index]], lower = lower, upper = upper,
    visualization_only = TRUE
  )]
}
feature_windows_observed <- rbindlist(window_rows, use.names = TRUE, fill = TRUE)
window_grid <- CJ(
  feature_id = coordinates$feature_id,
  dataset = all_cohorts,
  window_id = seq_along(centers),
  unique = TRUE
)
window_grid[, `:=`(
  center = centers[window_id],
  lower = centers[window_id] - window_width / 2,
  upper = centers[window_id] + window_width / 2
)]
feature_windows <- merge(
  window_grid,
  coordinates[, .(feature_id, collection, community_id, system_display)],
  by = "feature_id", all.x = TRUE, sort = FALSE
)
feature_windows <- merge(
  feature_windows,
  feature_windows_observed[, .(
    feature_id, dataset, window_id, mean_score, se_score, n_participants,
    visualization_only
  )],
  by = c("feature_id", "dataset", "window_id"), all.x = TRUE, sort = FALSE
)
feature_windows[is.na(visualization_only), visualization_only := TRUE]
ml_assert(nrow(feature_windows) == 5445L * length(all_cohorts) * length(centers),
          "Complete fixed-window feature grid drift")

window_family <- feature_windows[!is.na(community_id) & is.finite(mean_score), .(
  family_median_score = median(mean_score),
  n_features = .N,
  n_participants_min = min(n_participants),
  n_participants_max = max(n_participants)
), by = .(
  community_id, system_display, collection, dataset, window_id, center, lower, upper
)]
community_windows <- window_family[, .(
  collection_balanced_score = mean(family_median_score),
  n_collections = .N,
  n_features = sum(n_features),
  n_participants_min = min(n_participants_min),
  n_participants_max = max(n_participants_max)
), by = .(
  community_id, system_display, dataset, window_id, center, lower, upper
)]
feature_balanced_windows <- feature_windows[
  !is.na(community_id) & is.finite(mean_score), .(
    feature_balanced_score = median(mean_score),
    n_features_feature_balanced = .N
  ), by = .(
    community_id, system_display, dataset, window_id, center, lower, upper
  )
]
community_windows <- merge(
  community_windows, feature_balanced_windows,
  by = c(
    "community_id", "system_display", "dataset", "window_id",
    "center", "lower", "upper"
  ),
  all = TRUE, sort = FALSE
)
community_windows[, `:=`(
  dataset_label = factor(cohort_labels[dataset], levels = unname(cohort_labels[all_cohorts])),
  center_percent = center * 100,
  visualization_only = TRUE,
  collection_balancing = "median_within_collection_then_equal_mean_across_collections"
)]
community_windows_long <- melt(
  community_windows,
  id.vars = c(
    "community_id", "system_display", "dataset", "window_id", "center",
    "lower", "upper", "dataset_label", "center_percent", "visualization_only"
  ),
  measure.vars = c("collection_balanced_score", "feature_balanced_score"),
  variable.name = "aggregation_method", value.name = "system_score"
)
community_windows_long[, aggregation_method := fifelse(
  aggregation_method == "collection_balanced_score",
  "collection_balanced", "feature_balanced_median"
)]
community_window_trends <- community_windows_long[, .(
  window_spearman = suppressWarnings(cor(
    center, system_score, method = "spearman", use = "complete.obs"
  )),
  window_dynamic_range = max(system_score, na.rm = TRUE) -
    min(system_score, na.rm = TRUE),
  n_windows = sum(is.finite(system_score))
), by = .(community_id, system_display, dataset, aggregation_method)]
community_trend_concordance <- community_window_trends[, {
  evaluation <- window_spearman[dataset %in% contract$evaluation_cohorts]
  source_overlap <- window_spearman[dataset %in% contract$source_overlap_cohorts]
  all_values <- window_spearman[dataset %in% all_cohorts]
  same_direction <- function(values) {
    length(values) > 1L && all(is.finite(values)) &&
      (all(values > 0) || all(values < 0))
  }
  .(
    evaluation_same_direction = same_direction(evaluation),
    source_overlap_same_direction = same_direction(source_overlap),
    all_five_same_direction = same_direction(all_values),
    evaluation_median_rho = median(evaluation, na.rm = TRUE),
    source_overlap_median_rho = median(source_overlap, na.rm = TRUE),
    all_five_median_rho = median(all_values, na.rm = TRUE)
  )
}, by = .(community_id, system_display, aggregation_method)]

aggregation_window_summary <- community_trend_concordance[, .(
  n_systems = .N,
  evaluation_direction_concordant = sum(evaluation_same_direction),
  all_five_direction_concordant = sum(all_five_same_direction),
  median_evaluation_rho = median(evaluation_median_rho, na.rm = TRUE),
  median_all_five_rho = median(all_five_median_rho, na.rm = TRUE)
), by = aggregation_method]

collection_window_check <- feature_windows[
  !is.na(community_id) & is.finite(mean_score), .(
    independently_recomputed_collection_median = median(mean_score)
  ), by = .(community_id, collection, dataset, window_id)
][, .(
  independently_recomputed_score = mean(independently_recomputed_collection_median)
), by = .(community_id, dataset, window_id)]
collection_window_check <- merge(
  collection_window_check,
  community_windows[, .(
    community_id, dataset, window_id, collection_balanced_score
  )],
  by = c("community_id", "dataset", "window_id"), all = TRUE
)
collection_balance_max_difference <- max(abs(
  collection_window_check$independently_recomputed_score -
    collection_window_check$collection_balanced_score
), na.rm = TRUE)

atlas_summary <- data.table(
  metric = c(
    "complete_feature_family", "mapped_atlas_nodes", "mapped_pathways",
    "mapped_hotspot_programs", "main_component_nodes",
    "excluded_microcomponent_nodes", "membership_defined_systems",
    "node_stage_continuum_spearman", "node_same_direction_fraction",
    "collection_balanced_system_stage_continuum_spearman",
    "collection_balanced_system_same_direction_fraction",
    "feature_balanced_system_stage_continuum_spearman",
    "feature_balanced_system_same_direction_fraction",
    "systems_with_evaluation_window_direction_concordance",
    "systems_with_all_five_window_direction_concordance"
  ),
  value = as.character(c(
    nrow(complete_ledger), nrow(node_results),
    sum(node_results$feature_type == "pathway"),
    sum(node_results$feature_type == "hotspot_program"),
    nrow(system_node_results), nrow(microcomponent_results), nrow(communities),
    cor(system_node_results$stage_beta, system_node_results$continuum_beta,
        method = "spearman", use = "complete.obs"),
    mean(sign(system_node_results$stage_beta) ==
           sign(system_node_results$continuum_beta)),
    cor(community_effects$stage_beta, community_effects$continuum_beta,
        method = "spearman", use = "complete.obs"),
    mean(community_effects$same_direction),
    cor(community_effects$stage_beta_feature_balanced,
        community_effects$continuum_beta_feature_balanced,
        method = "spearman", use = "complete.obs"),
    mean(community_effects$same_direction_feature_balanced),
    sum(community_trend_concordance[
      aggregation_method == "collection_balanced", evaluation_same_direction
    ]),
    sum(community_trend_concordance[
      aggregation_method == "collection_balanced", all_five_same_direction
    ])
  )),
  interpretation_scope = c(
    rep("outcome_blind_geometry", 7L),
    rep("descriptive_correlated_feature_alignment_no_inferential_p_value", 6L),
    rep("descriptive_overlapping_windows_no_inferential_p_value", 2L)
  )
)
system_order <- communities[
  order(as.integer(sub("system_", "", community_id))), system_display
]
community_windows[, system_factor := factor(system_display, levels = rev(system_order))]
heat_limit <- as.numeric(quantile(
  abs(community_windows$collection_balanced_score), 0.98,
  names = FALSE, type = 8, na.rm = TRUE
))
heatmap_height <- max(4.2, 1.2 + 0.135 * nrow(communities))
window_heatmap <- ggplot(
  community_windows,
  aes(center_percent, system_factor, fill = collection_balanced_score)
) +
  geom_tile(width = 10, height = 0.92) +
  facet_grid(. ~ dataset_label) +
  scale_fill_gradient2(
    low = "#1565C0", mid = "white", high = "#C9265E", midpoint = 0,
    limits = c(-heat_limit, heat_limit), oob = scales::squish,
    breaks = scales::breaks_pretty(n = 3), labels = compact_numeric_labels,
    name = "Collection-balanced\nprogram score (z)"
  ) +
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90), expand = c(0, 0)) +
  labs(
    x = "Fixed-projection continuum percentile (20-percentile windows; 50% overlap)",
    y = "Membership-defined molecular system",
    caption = paste0(
      "Descriptive five-cohort sliding windows. Within each system, feature medians are",
      "\ncalculated per collection and collections receive equal weight; rows follow fixed size rank."
    )
  ) +
  theme_atlas() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = grid::unit(30, "mm"), barheight = grid::unit(2, "mm")
  ))
window_heatmap_path <- save_panel(
  window_heatmap, "s4_molecular_systems_five_cohort_sliding_windows.pdf",
  5.5, heatmap_height
)

rendered_label_manifest <- rbindlist(list(
  CJ(
    figure_id = c("membership_atlas", "fibrosis_stage_map", "continuum_effect_map"),
    community_id = map_labels$community_id,
    unique = TRUE
  )[, `:=`(
    rendered_label = sub("system_", "S", community_id),
    selection_rule = "15_largest_systems_fixed_before_outcome_overlay"
  )],
  data.table(
    figure_id = "system_stage_continuum_summary",
    community_id = communities$community_id,
    rendered_label = sub("system_", "S", communities$community_id),
    selection_rule = "all_main_component_systems"
  ),
  data.table(
    figure_id = "five_cohort_sliding_windows",
    community_id = communities$community_id,
    rendered_label = communities$system_display,
    selection_rule = "all_main_component_systems_as_y_axis_labels"
  )
), use.names = TRUE, fill = TRUE)
rendered_label_manifest[, intended_and_rendered := TRUE]

ml_write_tsv_once(node_results,
                  file.path(overlay_dir, "all_node_stage_continuum_overlays.tsv.gz"))
ml_write_tsv_once(community_effects,
                  file.path(overlay_dir, "community_stage_continuum_effects.tsv"))
ml_write_tsv_once(rbindlist(list(
  stage_systems$family[, effect_type := "fibrosis_stage"],
  continuum_systems$family[, effect_type := "stage_adjusted_continuum"]
), use.names = TRUE, fill = TRUE),
file.path(overlay_dir, "community_collection_effect_components.tsv"))
ml_write_tsv_once(aggregation_effect_sensitivity,
                  file.path(overlay_dir, "community_effect_aggregation_sensitivity.tsv"))
ml_write_tsv_once(aggregation_effect_summary,
                  file.path(overlay_dir, "effect_aggregation_method_summary.tsv"))
ml_write_tsv_once(feature_windows,
                  file.path(overlay_dir, "all_node_five_cohort_fixed_windows.tsv.gz"))
ml_write_tsv_once(window_family,
                  file.path(overlay_dir, "community_collection_fixed_windows.tsv.gz"))
ml_write_tsv_once(community_windows,
                  file.path(overlay_dir, "community_balanced_fixed_windows.tsv.gz"))
ml_write_tsv_once(community_windows_long,
                  file.path(overlay_dir, "community_fixed_windows_by_aggregation.tsv.gz"))
ml_write_tsv_once(community_window_trends,
                  file.path(overlay_dir, "community_five_cohort_window_trends.tsv"))
ml_write_tsv_once(community_trend_concordance,
                  file.path(overlay_dir, "community_window_direction_concordance.tsv"))
ml_write_tsv_once(aggregation_window_summary,
                  file.path(overlay_dir, "window_aggregation_method_summary.tsv"))
ml_write_tsv_once(collection_window_check,
                  file.path(overlay_dir, "collection_balance_independent_check.tsv.gz"))
ml_write_tsv_once(atlas_summary,
                  file.path(overlay_dir, "atlas_summary.tsv"))
ml_write_tsv_once(rendered_label_manifest,
                  file.path(overlay_dir, "rendered_label_manifest.tsv"))

figure_manifest <- data.table(
  figure_id = c(
    "membership_atlas", "fibrosis_stage_map", "continuum_effect_map",
    "system_stage_continuum_summary", "five_cohort_sliding_windows"
  ),
  path = c(
    topology_path, stage_map_path, continuum_map_path,
    system_scatter_path, window_heatmap_path
  ),
  role = c(
    "outcome_blind_topology", "stage_overlay", "continuum_overlay",
    "community_level_effect_comparison", "descriptive_five_cohort_windows"
  ),
  geometry_outcome_blind = TRUE,
  main_inference_cohorts = c(NA_integer_, 2L, 2L, 2L, NA_integer_),
  descriptive_cohorts = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, 5L),
  result_adaptive_geometry = FALSE,
  promoted = FALSE,
  coordinate_source = c(
    canonical_coordinate_path, canonical_coordinate_path,
    canonical_coordinate_path, NA_character_, NA_character_
  ),
  rendered_label_count = c(
    nrow(map_labels), nrow(map_labels), nrow(map_labels),
    nrow(communities), nrow(communities)
  )
)
figure_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
ml_write_tsv_once(figure_manifest, file.path(output_candidate, "figure_manifest.tsv"))

input_manifest <- data.table(
  path = normalizePath(c(
    coordinates_path, community_path, pathway_effect_path, hotspot_effect_path,
    axis_path, hotspot_score_path,
    file.path(
      analysis_candidate, "pathway_tf", "pathways", family_order[family_order != "hotspot"],
      "donor_scores.tsv.gz"
    )
  ), mustWork = TRUE),
  input_role = c(
    "frozen_outcome_blind_coordinates", "frozen_outcome_blind_communities",
    "pathway_stage_continuum_overlay", "hotspot_stage_continuum_overlay",
    "fixed_projection_percentiles", "hotspot_participant_scores",
    rep("pathway_donor_scores", 6L)
  )
)
input_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
ml_write_tsv_once(input_manifest, file.path(provenance_dir, "overlay_input_manifest.tsv"))

geometry_hash_at_overlay_exit <- ml_sha256(coordinates_path)
same_coordinate_source_pass <-
  uniqueN(figure_manifest[
    figure_id %in% c("fibrosis_stage_map", "continuum_effect_map"),
    coordinate_source
  ]) == 1L &&
  identical(
    figure_manifest[figure_id == "fibrosis_stage_map", coordinate_source],
    figure_manifest[figure_id == "continuum_effect_map", coordinate_source]
  ) &&
  identical(ml_sha256(canonical_coordinate_path), canonical_coordinate_hash)
complete_geometry_pass <-
  nrow(node_results) == sum(complete_ledger$geometry_eligible) &&
  nrow(system_node_results) + nrow(microcomponent_results) == nrow(node_results)
fixed_windows_pass <-
  identical(sort(unique(feature_windows$window_id)), 1:9) &&
  identical(sort(unique(feature_windows$center)), as.numeric(contract$window_centers))
five_cohort_pass <- setequal(unique(community_windows$dataset), all_cohorts)
collection_balance_pass <-
  is.finite(collection_balance_max_difference) &&
  collection_balance_max_difference < 1e-12
labels_pass <-
  all(rendered_label_manifest$intended_and_rendered) &&
  rendered_label_manifest[, uniqueN(community_id), by = figure_id][
    figure_id %in% c(
      "membership_atlas", "fibrosis_stage_map", "continuum_effect_map"
    ), all(V1 == maximum_map_labels)
  ]
candidate_manifest_pass <- candidate_only_pass && !any(figure_manifest$promoted)
overlay_audit <- data.table(
  check = c(
    "geometry_unchanged_by_outcome_overlay", "same_coordinates_for_stage_and_continuum",
    "complete_geometry_node_family", "fixed_nine_windows",
    "five_descriptive_cohorts", "collection_balanced_system_summaries",
    "rendered_labels_match_manifest", "candidate_only"
  ),
  passed = c(
    identical(expected_geometry_hash, geometry_hash_at_overlay_entry) &&
      identical(expected_geometry_hash, geometry_hash_at_overlay_exit),
    same_coordinate_source_pass, complete_geometry_pass, fixed_windows_pass,
    five_cohort_pass, collection_balance_pass, labels_pass,
    candidate_manifest_pass
  ),
  detail = c(
    paste(expected_geometry_hash, geometry_hash_at_overlay_exit, sep = " -> "),
    paste0("Canonical coordinate SHA-256: ", canonical_coordinate_hash),
    paste0(
      nrow(node_results), " mapped nodes; ", nrow(system_node_results),
      " in 43 systems and 6 in three explicit microcomponents"
    ),
    "20-percentile-wide windows centered at 10 through 90 with 50% overlap",
    paste(all_cohorts, collapse = ";"),
    paste0("Independent recomputation maximum absolute difference: ",
           format(collection_balance_max_difference, scientific = TRUE)),
    paste0(nrow(rendered_label_manifest), " figure-label rows rendered without check_overlap"),
    paste0("Candidate root asserted: ", output_candidate, "; promoted=FALSE")
  )
)
ml_assert(all(overlay_audit$passed), "Molecular-system overlay audit failed")
ml_write_tsv_once(overlay_audit, file.path(provenance_dir, "overlay_audit.tsv"))
ml_write_session_info(file.path(provenance_dir, "overlay_sessionInfo.txt"))

output_files <- c(
  figure_manifest$path,
  file.path(overlay_dir, "all_node_stage_continuum_overlays.tsv.gz"),
  file.path(overlay_dir, "community_stage_continuum_effects.tsv"),
  file.path(overlay_dir, "community_effect_aggregation_sensitivity.tsv"),
  file.path(overlay_dir, "effect_aggregation_method_summary.tsv"),
  file.path(overlay_dir, "all_node_five_cohort_fixed_windows.tsv.gz"),
  file.path(overlay_dir, "community_balanced_fixed_windows.tsv.gz"),
  file.path(overlay_dir, "community_fixed_windows_by_aggregation.tsv.gz"),
  file.path(overlay_dir, "community_window_direction_concordance.tsv"),
  file.path(overlay_dir, "window_aggregation_method_summary.tsv"),
  file.path(overlay_dir, "rendered_label_manifest.tsv"),
  canonical_coordinate_path, stage_map_source_path, continuum_map_source_path,
  file.path(overlay_dir, "atlas_summary.tsv")
)
output_checksums <- data.table(
  path = normalizePath(output_files, mustWork = TRUE),
  sha256 = vapply(output_files, ml_sha256, character(1))
)
ml_write_tsv_once(output_checksums,
                  file.path(provenance_dir, "atlas_output_checksums.tsv"))
