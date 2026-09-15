#!/usr/bin/env Rscript
# KEY MESSAGE: A single strict donor-level marker identifies the molecular systems that remain continuum-associated after recorded fibrosis stage and sex adjustment.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(digits = 17, scipen = 999)
grDevices::pdf.options(useDingbats = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    paste(
      "Usage: 71_render_system_inference_strip.R",
      "ACTIVE_REFINED_SOURCE SYSTEM_INFERENCE_CANDIDATE OUTPUT_CANDIDATE"
    ),
    call. = FALSE
  )
}

refined_source <- normalizePath(args[[1L]], mustWork = TRUE)
inference_candidate <- normalizePath(args[[2L]], mustWork = TRUE)
output_candidate <- args[[3L]]
if (file.exists(output_candidate)) {
  stop("Refusing to overwrite output candidate: ", output_candidate, call. = FALSE)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lib_molecular_layers.R"))

candidate_root <- normalizePath(
  file.path(ml_project_root(), "figures", "candidates"), mustWork = TRUE
)
candidate_only_pass <- startsWith(
  paste0(normalizePath(dirname(output_candidate), mustWork = TRUE), "/"),
  paste0(candidate_root, "/")
)
ml_assert(candidate_only_pass,
          "System-inference figures are restricted to figures/candidates")

panel_dir <- file.path(output_candidate, "panels")
source_dir <- file.path(output_candidate, "source_tables")
provenance_dir <- file.path(output_candidate, "provenance")
for (path in c(panel_dir, source_dir, provenance_dir)) ml_ensure_dir(path)

order_path <- file.path(
  refined_source, "source_tables", "pooled_system_late_minus_early_order.tsv"
)
display_path <- file.path(
  refined_source, "source_tables", "pooled_system_heatmap_display_values.tsv"
)
membership_path <- file.path(
  inference_candidate, "tables", "system_inference_membership.tsv"
)
refined_checksums_path <- file.path(
  refined_source, "provenance", "output_checksums.tsv"
)
inference_checksums_path <- file.path(
  inference_candidate, "provenance", "output_checksums.tsv"
)
required_inputs <- c(order_path, display_path, membership_path)
ml_assert(all(file.exists(c(
  required_inputs, refined_checksums_path, inference_checksums_path
))), "A system-inference figure input is missing")

validate_output_hash <- function(path, registry_path) {
  registry <- fread(registry_path)
  ml_assert(all(c("path", "sha256") %in% names(registry)),
            paste0("Checksum registry schema drift: ", registry_path))
  target_basename <- basename(path)
  expected <- registry[
    basename(get("path")) == target_basename, unique(sha256)
  ]
  ml_assert(length(expected) == 1L,
            paste0("Input hash is absent or ambiguous: ", path))
  ml_assert(identical(ml_sha256(path), expected),
            paste0("Input hash drift: ", path))
}
validate_output_hash(order_path, refined_checksums_path)
validate_output_hash(display_path, refined_checksums_path)
validate_output_hash(membership_path, inference_checksums_path)

input_manifest <- data.table(
  path = normalizePath(required_inputs, mustWork = TRUE),
  input_role = c(
    "active_pooled_display_order", "active_pooled_heatmap_values",
    "donor_level_system_inference"
  )
)
input_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
input_manifest[, hash_verified_before_render := TRUE]
ml_write_tsv_once(input_manifest, file.path(provenance_dir, "input_manifest.tsv"))

pooled_order <- fread(order_path, na.strings = c("", "NA"))
pooled_display <- fread(display_path, na.strings = c("", "NA"))
membership <- fread(membership_path, na.strings = c("", "NA"))
ml_assert(nrow(pooled_order) == 43L &&
            identical(pooled_order$display_rank, seq_len(43L)) &&
            uniqueN(pooled_order$community_id) == 43L,
          "Active pooled system order drift")
ml_assert(nrow(pooled_display) == 43L * 10L &&
            uniqueN(pooled_display$community_id) == 43L &&
            uniqueN(pooled_display$display_position) == 10L,
          "Active pooled heatmap family drift")
ml_assert(nrow(membership) == 43L &&
            uniqueN(membership$community_id) == 43L,
          "Donor-level system inference family drift")
ml_assert(setequal(pooled_order$community_id, membership$community_id),
          "Display and inference system rosters disagree")

ordered_systems <- pooled_order$community_id
factor_levels <- rev(ordered_systems)
label_lookup <- setNames(
  pooled_order$signed_short_label, pooled_order$community_id
)
pooled_display[, system_factor := factor(community_id, levels = factor_levels)]
plot_display <- pooled_display[display_position <= 9L]
ml_assert(nrow(plot_display) == 43L * 9L &&
            uniqueN(plot_display$display_position) == 9L,
          "Displayed fixed-window family drift")

primary_supported <- sum(membership$continuum_associated_system)
maxT_supported <- sum(membership$permutation_maxT_supported_all_four)
ml_assert(primary_supported == 33L && maxT_supported == 18L,
          "System-support census drift")

inference_strip <- rbindlist(list(
  membership[, .(
    community_id,
    inference_column = "BH",
    display_position = 11.65,
    supported = continuum_associated_system,
    displayed_metric = maximum_meta_q,
    metric_name = "maximum co-primary fixed-effect meta-analysis BH q",
    criterion = paste0(
      "both co-primary meta-analysis BH q<0.05; ",
      "same direction in both evaluation cohorts for both axes"
    ),
    effect_direction = fifelse(
      mean_meta_beta > 0, "positive", fifelse(mean_meta_beta < 0, "negative", "zero")
    ),
    effect_size = mean_meta_beta,
    outline_shape = 1L,
    supported_shape = 16L,
    rendered_in_panel = FALSE,
    render_symbol = NA_character_
  )],
  membership[, .(
    community_id,
    inference_column = "maxT",
    display_position = 13.00,
    supported = permutation_maxT_supported_all_four,
    displayed_metric = maximum_maxT_fwer_p,
    metric_name = "maximum within-stage/sex maxT family-wise P across four tests",
    criterion = paste0(
      "maxT family-wise P<0.05 in both evaluation cohorts for both axes; ",
      "same residual direction in all four tests"
    ),
    effect_direction = fifelse(
      mean_meta_beta > 0, "positive", fifelse(mean_meta_beta < 0, "negative", "zero")
    ),
    effect_size = mean_meta_beta,
    outline_shape = 5L,
    supported_shape = 18L,
    rendered_in_panel = permutation_maxT_supported_all_four,
    render_symbol = fifelse(permutation_maxT_supported_all_four, "*", NA_character_)
  )]
), use.names = TRUE)
inference_strip[, system_factor := factor(community_id, levels = factor_levels)]
setorder(inference_strip, inference_column, community_id)
ml_assert(nrow(inference_strip) == 43L * 2L &&
            inference_strip[inference_column == "BH", sum(supported)] == 33L &&
            inference_strip[inference_column == "maxT", sum(supported)] == 18L &&
            inference_strip[inference_column == "BH", sum(rendered_in_panel)] == 0L &&
            inference_strip[inference_column == "maxT", sum(rendered_in_panel)] == 18L,
          "System-inference strip family drift")

strict_markers <- inference_strip[
  inference_column == "maxT" & rendered_in_panel
]

theme_systems <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(size = 6, face = "plain"),
      axis.text = element_text(size = 6, color = "black", face = "plain"),
      axis.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      legend.title = element_text(size = 6, face = "plain"),
      plot.caption = element_text(size = 6, face = "plain", hjust = 0),
      plot.caption.position = "plot",
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      plot.margin = margin(1.5, 1.5, 1.5, 1.5)
    )
}

compact_numeric_labels <- function(values) {
  formatC(values, format = "f", digits = 1)
}

pooled_heat_limit <- as.numeric(quantile(
  abs(pooled_display$display_value), 0.98,
  names = FALSE, type = 8, na.rm = TRUE
))
pooled_heatmap <- ggplot(
  plot_display,
  aes(display_position, system_factor, fill = display_value)
) +
  geom_tile(width = 0.9, height = 0.94) +
  geom_text(
    data = strict_markers,
    aes(x = 9.80, y = system_factor, label = render_symbol),
    inherit.aes = FALSE, family = "Helvetica", fontface = "plain",
    size = 6 / ggplot2::.pt, color = "#4D4D4D", vjust = 0.35
  ) +
  scale_fill_gradient2(
    low = "#1565C0", mid = "#F7F7F7", high = "#C9265E", midpoint = 0,
    limits = c(-pooled_heat_limit, pooled_heat_limit), oob = scales::squish,
    breaks = c(-pooled_heat_limit, 0, pooled_heat_limit),
    labels = compact_numeric_labels,
    name = "Pooled system score (z)"
  ) +
  scale_x_continuous(
    breaks = c(1, 3, 5, 7, 9),
    labels = c("10", "30", "50", "70", "90"),
    expand = expansion(add = c(0.04, 0.04))
  ) +
  scale_y_discrete(labels = label_lookup) +
  coord_fixed(ratio = 1) +
  labs(
    x = "Continuum percentile",
    y = "System (late-early order)"
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

panel_path <- file.path(
  panel_dir, "s4_molecular_systems_pooled_continuum_windows_with_inference.pdf"
)
ggsave(
  panel_path, pooled_heatmap, width = 2.15, height = 4.92, units = "in",
  device = cairo_pdf, bg = "white", limitsize = FALSE
)
ml_assert(file.info(panel_path)$size > 1000L,
          "Unexpectedly small system-inference PDF")

display_output_path <- file.path(
  source_dir, "pooled_system_heatmap_display_values.tsv"
)
order_output_path <- file.path(
  source_dir, "pooled_system_late_minus_early_order.tsv"
)
membership_output_path <- file.path(
  source_dir, "system_inference_membership.tsv"
)
ml_assert(file.copy(display_path, display_output_path, overwrite = FALSE),
          "Failed to preserve the active heatmap display source byte-for-byte")
ml_assert(file.copy(order_path, order_output_path, overwrite = FALSE),
          "Failed to preserve the active row-order source byte-for-byte")
ml_assert(file.copy(membership_path, membership_output_path, overwrite = FALSE),
          "Failed to preserve the donor-level inference registry byte-for-byte")
ml_write_tsv_once(
  inference_strip,
  file.path(source_dir, "system_inference_strip.tsv")
)

figure_manifest <- data.table(
  figure_id = "pooled_continuum_windows_with_inference",
  path = normalizePath(panel_path, mustWork = TRUE),
  role = "descriptive_equal_cohort_windows_with_strict_donor_support_asterisk",
  window_values_changed = FALSE,
  display_order_changed = FALSE,
  delta_tile_displayed = FALSE,
  delta_order_retained = TRUE,
  inferential_model_changed = FALSE,
  windows_inferentially_tested = FALSE,
  embedded_explanatory_text = FALSE,
  promoted = FALSE
)
figure_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
ml_write_tsv_once(figure_manifest, file.path(output_candidate, "figure_manifest.tsv"))

audit <- data.table(
  check = c(
    "input_hashes_verified", "complete_43_system_family",
    "active_window_values_byte_identical", "active_display_order_byte_identical",
    "system_inference_registry_byte_identical",
    "primary_bh_census", "strict_maxT_census",
    "strict_marker_census", "primary_bh_not_rendered",
    "delta_tile_removed", "delta_order_retained",
    "embedded_explanatory_text_removed",
    "windows_remain_descriptive", "candidate_only", "pdf_nonempty"
  ),
  passed = c(
    all(input_manifest$hash_verified_before_render),
    uniqueN(inference_strip$community_id) == 43L,
    identical(ml_sha256(display_path), ml_sha256(display_output_path)),
    identical(ml_sha256(order_path), ml_sha256(order_output_path)),
    identical(ml_sha256(membership_path), ml_sha256(membership_output_path)),
    primary_supported == 33L,
    maxT_supported == 18L,
    nrow(strict_markers) == 18L,
    inference_strip[inference_column == "BH", sum(rendered_in_panel)] == 0L,
    !figure_manifest$delta_tile_displayed,
    figure_manifest$delta_order_retained,
    !figure_manifest$embedded_explanatory_text,
    !figure_manifest$windows_inferentially_tested,
    candidate_only_pass && !figure_manifest$promoted,
    file.info(panel_path)$size > 1000L
  ),
  detail = c(
    paste(basename(required_inputs), collapse = ";"),
    "43 systems x 2 inferential criteria",
    basename(display_path), basename(order_path),
    basename(membership_path),
    "Both co-primary fixed-effect meta q<0.05 and four-way direction agreement",
    "Within-stage/sex maxT FWER P<0.05 in all four axis-by-cohort tests",
    "One neutral asterisk per strict maxT-supported system",
    "Primary BH membership is retained in source tables but omitted from the panel",
    "Late-minus-early contrast is not rendered",
    "Late-minus-early contrast remains the display-only row-order statistic",
    "Methods and significance definitions remain in the external caption authority",
    "Nine overlapping windows are display only",
    normalizePath(output_candidate, mustWork = TRUE),
    basename(panel_path)
  )
)
ml_assert(all(audit$passed), "System-inference figure audit failed")
ml_write_tsv_once(audit, file.path(provenance_dir, "audit.tsv"))
ml_write_session_info(file.path(provenance_dir, "sessionInfo.txt"))

output_files <- c(
  panel_path,
  display_output_path,
  order_output_path,
  membership_output_path,
  file.path(source_dir, "system_inference_strip.tsv"),
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
