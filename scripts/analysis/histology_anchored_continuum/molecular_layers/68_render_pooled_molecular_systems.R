#!/usr/bin/env Rscript
# KEY MESSAGE: Equal weighting across cohorts reveals the shared molecular-system trajectories along the aligned MASLD continuum without letting the largest cohort dominate.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(digits = 17, scipen = 999)
grDevices::pdf.options(useDingbats = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: 68_render_pooled_molecular_systems.R PROMOTED_ATLAS_SOURCE OUTPUT_CANDIDATE",
    call. = FALSE
  )
}

atlas_source <- normalizePath(args[[1L]], mustWork = TRUE)
output_candidate <- args[[2L]]
if (file.exists(output_candidate)) {
  stop("Refusing to overwrite output candidate: ", output_candidate, call. = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lib_molecular_layers.R"))
contract <- ml_read_contract()

panel_dir <- file.path(output_candidate, "panels")
source_dir <- file.path(output_candidate, "source_tables")
provenance_dir <- file.path(output_candidate, "provenance")
for (path in c(panel_dir, source_dir, provenance_dir)) ml_ensure_dir(path)

window_path <- file.path(
  atlas_source, "overlays", "community_balanced_fixed_windows.tsv.gz"
)
community_path <- file.path(atlas_source, "geometry", "community_registry.tsv")
promoted_checksums_path <- file.path(atlas_source, "promoted_source_checksums.tsv")
required_inputs <- c(window_path, community_path, promoted_checksums_path)
ml_assert(all(file.exists(required_inputs)), "A promoted atlas input is missing")

promoted_checksums <- fread(promoted_checksums_path)
ml_assert(all(c("path", "sha256") %in% names(promoted_checksums)),
          "Promoted source checksum schema drift")
for (path in c(window_path, community_path)) {
  target_basename <- basename(path)
  expected <- promoted_checksums[
    basename(get("path")) == target_basename, unique(sha256)
  ]
  ml_assert(length(expected) == 1L,
            paste0("Promoted checksum is absent or ambiguous: ", path))
  ml_assert(identical(ml_sha256(path), expected),
            paste0("Promoted source hash drift: ", path))
}

windows <- fread(window_path, na.strings = c("", "NA"))
communities <- fread(community_path, na.strings = c("", "NA"))
all_cohorts <- c(contract$source_overlap_cohorts, contract$evaluation_cohorts)
expected_rows <- 43L * length(all_cohorts) * length(contract$window_centers)
required_columns <- c(
  "community_id", "system_display", "dataset", "window_id", "center",
  "lower", "upper", "collection_balanced_score", "feature_balanced_score",
  "n_collections", "n_features", "n_participants_min", "n_participants_max"
)
ml_assert(all(required_columns %in% names(windows)),
          "Community-window source schema drift")
ml_assert(nrow(windows) == expected_rows,
          "Community-window source family size drift")
ml_assert(!anyDuplicated(windows[, .(community_id, dataset, window_id)]),
          "Community-window rows are duplicated")
ml_assert(setequal(unique(windows$dataset), all_cohorts),
          "Community-window cohort roster drift")
ml_assert(setequal(unique(windows$community_id), communities$community_id),
          "Community registry and window systems disagree")
ml_assert(
  identical(sort(unique(windows$window_id)), 1:9) &&
    identical(sort(unique(windows$center)), as.numeric(contract$window_centers)),
  "Fixed-window contract drift"
)
ml_assert(all(is.finite(windows$collection_balanced_score)) &&
            all(is.finite(windows$feature_balanced_score)),
          "A system window score is non-finite")
ml_assert(all(windows$n_participants_min == windows$n_participants_max),
          "Feature-level participant counts differ within a cohort window")

pooled <- windows[, .(
  pooled_collection_balanced_score = mean(collection_balanced_score),
  pooled_collection_balanced_median = median(collection_balanced_score),
  pooled_feature_balanced_score = mean(feature_balanced_score),
  cohort_sd = sd(collection_balanced_score),
  cohort_min = min(collection_balanced_score),
  cohort_max = max(collection_balanced_score),
  source_overlap_mean = mean(
    collection_balanced_score[dataset %in% contract$source_overlap_cohorts]
  ),
  evaluation_mean = mean(
    collection_balanced_score[dataset %in% contract$evaluation_cohorts]
  ),
  n_cohorts = uniqueN(dataset),
  n_source_overlap_cohorts = sum(dataset %in% contract$source_overlap_cohorts),
  n_evaluation_cohorts = sum(dataset %in% contract$evaluation_cohorts),
  n_window_participants = sum(n_participants_min),
  cohort_roster = paste(sort(dataset), collapse = ";")
), by = .(
  community_id, system_display, window_id, center, lower, upper
)]
setorder(pooled, community_id, window_id)
pooled[, `:=`(
  center_percent = center * 100,
  pooling_estimand = paste0(
    "mean_of_five_cohort_specific_collection_balanced_scores;",
    "equal_cohort_weight;within_cohort_percentiles_and_feature_standardization"
  ),
  visualization_only = TRUE,
  overlapping_windows = TRUE
)]

ml_assert(nrow(pooled) == 43L * 9L, "Pooled system-window family size drift")
ml_assert(all(pooled$n_cohorts == 5L) &&
            all(pooled$n_source_overlap_cohorts == 3L) &&
            all(pooled$n_evaluation_cohorts == 2L),
          "Pooled windows do not give every cohort equal representation")

independent_pool <- windows[, .(
  independently_recomputed = sum(collection_balanced_score / .N)
), by = .(community_id, window_id)]
pool_check <- merge(
  pooled[, .(community_id, window_id, pooled_collection_balanced_score)],
  independent_pool,
  by = c("community_id", "window_id"), all = TRUE
)
pool_check[, absolute_difference := abs(
  pooled_collection_balanced_score - independently_recomputed
)]
pooling_max_difference <- max(pool_check$absolute_difference)
ml_assert(is.finite(pooling_max_difference) && pooling_max_difference < 1e-12,
          "Independent equal-cohort pooling check failed")

system_order <- communities[
  order(as.integer(sub("system_", "", community_id))), system_display
]
pooled[, system_factor := factor(system_display, levels = rev(system_order))]
heat_limit <- as.numeric(quantile(
  abs(pooled$pooled_collection_balanced_score), 0.98,
  names = FALSE, type = 8, na.rm = TRUE
))
heatmap_height <- max(4.2, 1.2 + 0.135 * nrow(communities))

theme_atlas <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(size = 6, face = "plain"),
      axis.text = element_text(size = 6, color = "black", face = "plain"),
      axis.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      legend.title = element_text(size = 6, face = "plain"),
      plot.caption = element_text(size = 6, face = "plain", hjust = 0),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.key.size = grid::unit(2.8, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

compact_numeric_labels <- function(values) {
  formatC(values, format = "f", digits = 1)
}

pooled_heatmap <- ggplot(
  pooled,
  aes(center_percent, system_factor, fill = pooled_collection_balanced_score)
) +
  geom_tile(width = 10, height = 0.92) +
  scale_fill_gradient2(
    low = "#1565C0", mid = "white", high = "#C9265E", midpoint = 0,
    limits = c(-heat_limit, heat_limit), oob = scales::squish,
    breaks = c(-heat_limit, 0, heat_limit), labels = compact_numeric_labels,
    name = "Equal-cohort pooled\nsystem score (z)"
  ) +
  scale_x_continuous(
    breaks = c(10, 30, 50, 70, 90),
    expand = c(0, 0)
  ) +
  labs(
    x = paste0(
      "Aligned within-cohort fixed-projection percentile ",
      "(20-percentile windows; 50% overlap)"
    ),
    y = "Membership-defined molecular system",
    caption = paste0(
      "Equal cohort weight; percentiles and feature z scores remain within cohort.",
      "\nDescriptive overlapping windows; fixed size-rank row order."
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

panel_path <- file.path(
  panel_dir, "s4_molecular_systems_pooled_continuum_windows.pdf"
)
ggsave(
  panel_path, pooled_heatmap,
  width = 4.4, height = heatmap_height, units = "in",
  device = cairo_pdf, bg = "white", limitsize = FALSE
)
ml_assert(file.info(panel_path)$size > 10000L,
          "Pooled molecular-system PDF is unexpectedly small")

ml_write_tsv_once(pooled,
                  file.path(source_dir, "pooled_molecular_system_windows.tsv"))
ml_write_tsv_once(pool_check,
                  file.path(source_dir, "equal_cohort_pooling_check.tsv"))

input_manifest <- data.table(
  path = normalizePath(required_inputs, mustWork = TRUE),
  input_role = c(
    "promoted_five_cohort_system_windows",
    "promoted_outcome_blind_system_registry",
    "promoted_source_checksum_registry"
  )
)
input_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
ml_write_tsv_once(input_manifest,
                  file.path(provenance_dir, "input_manifest.tsv"))

audit <- data.table(
  check = c(
    "promoted_input_hashes", "complete_system_window_family",
    "fixed_nine_windows", "equal_weight_five_cohort_pooling",
    "fixed_system_order", "descriptive_only", "pdf_nonempty"
  ),
  passed = c(
    TRUE,
    nrow(windows) == 43L * 5L * 9L && nrow(pooled) == 43L * 9L,
    identical(sort(unique(pooled$window_id)), 1:9),
    all(pooled$n_cohorts == 5L) && pooling_max_difference < 1e-12,
    identical(levels(pooled$system_factor), rev(system_order)),
    all(pooled$visualization_only) && all(pooled$overlapping_windows),
    file.info(panel_path)$size > 10000L
  ),
  detail = c(
    paste(basename(required_inputs), collapse = ";"),
    "43 systems x 5 cohorts x 9 windows pooled to 43 x 9",
    "20-percentile-wide windows centered at 10 through 90 with 50% overlap",
    paste0("Maximum independent recomputation difference: ",
           format(pooling_max_difference, scientific = TRUE)),
    "Rows retain fixed system_01 through system_43 size-rank order",
    "No p values, support calls, timing labels, or result-adaptive ordering",
    normalizePath(panel_path, mustWork = TRUE)
  )
)
ml_assert(all(audit$passed), "Pooled molecular-system audit failed")
ml_write_tsv_once(audit, file.path(provenance_dir, "audit.tsv"))

figure_manifest <- data.table(
  figure_id = "pooled_continuum_windows",
  path = normalizePath(panel_path, mustWork = TRUE),
  role = "descriptive_equal_cohort_pooled_system_windows",
  n_systems = 43L,
  n_windows = 9L,
  n_cohorts = 5L,
  pooling_weight = "equal_cohort",
  continuum_alignment = "within_cohort_percentile",
  result_adaptive_ordering = FALSE,
  promoted = FALSE,
  sha256 = ml_sha256(panel_path)
)
ml_write_tsv_once(figure_manifest,
                  file.path(output_candidate, "figure_manifest.tsv"))
ml_write_session_info(file.path(provenance_dir, "sessionInfo.txt"))

output_files <- c(
  panel_path,
  file.path(source_dir, "pooled_molecular_system_windows.tsv"),
  file.path(source_dir, "equal_cohort_pooling_check.tsv"),
  file.path(provenance_dir, "input_manifest.tsv"),
  file.path(provenance_dir, "audit.tsv"),
  file.path(provenance_dir, "sessionInfo.txt"),
  file.path(output_candidate, "figure_manifest.tsv")
)
output_checksums <- data.table(
  path = normalizePath(output_files, mustWork = TRUE),
  sha256 = vapply(output_files, ml_sha256, character(1))
)
ml_write_tsv_once(output_checksums,
                  file.path(provenance_dir, "output_checksums.tsv"))
