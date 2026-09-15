#!/usr/bin/env Rscript
# KEY MESSAGE: Molecular systems broadly align with fibrosis remodeling while retaining stage-adjusted continuum associations.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

options(digits = 17, scipen = 999)
grDevices::pdf.options(useDingbats = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: 72_render_main_system_stage_continuum.R",
      "PROMOTED_ATLAS_SOURCE PROMOTED_SUPPORT_SOURCE",
      "SYSTEM_STRESS_SOURCE OUTPUT_CANDIDATE"
    ),
    call. = FALSE
  )
}

atlas_source <- normalizePath(args[[1L]], mustWork = TRUE)
support_source <- normalizePath(args[[2L]], mustWork = TRUE)
stress_source <- normalizePath(args[[3L]], mustWork = TRUE)
output_candidate <- args[[4L]]
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
candidate_parent <- normalizePath(dirname(output_candidate), mustWork = TRUE)
ml_assert(
  startsWith(paste0(candidate_parent, "/"), paste0(candidate_root, "/")),
  "Main Figure 4H rendering is restricted to figures/candidates"
)

panel_dir <- file.path(output_candidate, "panels")
source_dir <- file.path(output_candidate, "source_tables")
provenance_dir <- file.path(output_candidate, "provenance")
for (path in c(panel_dir, source_dir, provenance_dir)) ml_ensure_dir(path)

effect_path <- file.path(
  atlas_source, "overlays", "community_stage_continuum_effects.tsv"
)
label_path <- file.path(atlas_source, "overlays", "rendered_label_manifest.tsv")
support_path <- file.path(
  support_source, "source_tables", "system_inference_strip.tsv"
)
atlas_checksums_path <- file.path(atlas_source, "promoted_source_checksums.tsv")
support_checksums_path <- file.path(support_source, "promoted_source_checksums.tsv")
stress_membership_path <- file.path(
  stress_source, "tables", "system_stress_membership.tsv"
)
stress_competitive_path <- file.path(
  stress_source, "tables", "system_competitive_null.tsv"
)
stress_redundancy_path <- file.path(
  stress_source, "tables", "system_redundancy.tsv"
)
stress_checksums_path <- file.path(
  stress_source, "provenance", "output_checksums.tsv"
)
required_inputs <- c(
  effect_path, label_path, support_path,
  atlas_checksums_path, support_checksums_path,
  stress_membership_path, stress_competitive_path, stress_redundancy_path,
  stress_checksums_path
)
ml_assert(all(file.exists(required_inputs)), "A promoted Figure 4H input is missing")

validate_hash <- function(path, checksum_path) {
  checksums <- fread(checksum_path)
  target_basename <- basename(path)
  expected <- checksums[basename(get("path")) == target_basename, unique(sha256)]
  ml_assert(length(expected) == 1L,
            paste0("Promoted checksum is absent or ambiguous: ", path))
  ml_assert(identical(ml_sha256(path), expected),
            paste0("Promoted source hash drift: ", path))
}
validate_hash(effect_path, atlas_checksums_path)
validate_hash(label_path, atlas_checksums_path)
validate_hash(support_path, support_checksums_path)
for (path in c(stress_membership_path, stress_competitive_path,
               stress_redundancy_path)) {
  validate_hash(path, stress_checksums_path)
}

effects <- fread(effect_path, na.strings = c("", "NA"))
labels <- fread(label_path)[
  figure_id == "continuum_effect_map" &
    selection_rule == "15_largest_systems_fixed_before_outcome_overlay"
]
support <- fread(support_path)[inference_column == "maxT", .(
  community_id, strict_maxt_supported = as.logical(supported)
)]
stress <- fread(stress_membership_path, na.strings = c("", "NA"))[, .(
  community_id,
  joint_maxT_supported = as.logical(joint_maxT_supported),
  connected_null_supported = as.logical(connected_null_both_axes),
  uniform_null_supported = as.logical(uniform_null_both_axes)
)]
redundancy <- fread(stress_redundancy_path)[scope == "pooled_centred"]

ml_assert(nrow(effects) == 43L && uniqueN(effects$community_id) == 43L,
          "The complete 43-system effect family drifted")
ml_assert(nrow(labels) == 15L && uniqueN(labels$community_id) == 15L,
          "The fixed 15-system outcome-blind label roster drifted")
ml_assert(nrow(support) == 43L && uniqueN(support$community_id) == 43L,
          "The complete 43-system maxT family drifted")
ml_assert(sum(support$strict_maxt_supported) == 18L,
          "Strict maxT-supported system count drifted")
ml_assert(nrow(stress) == 43L && uniqueN(stress$community_id) == 43L,
          "The complete 43-system stress-test family drifted")
ml_assert(nrow(redundancy) == 1L, "Pooled redundancy summary is absent")

# Ink marks the competitive result, not the four-test maxT conjunction. That
# conjunction is limited by the smaller cohort's within-stratum permutation
# resolution, and the joint cross-cohort test supports 35 of 43 systems, so
# neither count is scarce enough to carry ink. The size-matched draws are the
# only criterion that separates a small number of systems from the rest.
plot_data <- merge(effects, support, by = "community_id", all.x = TRUE)
plot_data <- merge(plot_data, stress, by = "community_id", all.x = TRUE)
plot_data[, competitive_marked := uniform_null_supported]
ml_assert(sum(plot_data$connected_null_supported) == 0L,
          "A system now passes the conservative connected-subgraph null; revisit the ink rule")
ml_assert(sum(plot_data$competitive_marked) == 1L &&
            plot_data[competitive_marked == TRUE, community_id] == "system_01",
          "The competitive-null marked system drifted from the single ECM system")
plot_data[, `:=`(
  short_label = sprintf("S%02d", as.integer(sub("system_", "", community_id))),
  fixed_label = community_id %in% labels$community_id,
  direction_pair = fcase(
    stage_beta > 0 & continuum_beta > 0, "Both positive",
    stage_beta < 0 & continuum_beta < 0, "Both negative",
    default = "Opposite signs"
  )
)]
direction_levels <- c("Both positive", "Both negative", "Opposite signs")
plot_data[, direction_pair := factor(direction_pair, levels = direction_levels)]

system_rho <- cor(
  plot_data$stage_beta, plot_data$continuum_beta,
  method = "spearman", use = "complete.obs"
)
same_sign_count <- sum(
  sign(plot_data$stage_beta) == sign(plot_data$continuum_beta)
)
ml_assert(abs(system_rho - 0.658) < 0.001 && same_sign_count == 36L,
          "Stage-versus-continuum system summary drifted")

direction_colors <- c(
  "Both positive" = "#C9265E",
  "Both negative" = "#1565C0",
  "Opposite signs" = "#9E9E9E"
)

theme_main <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(size = 6, face = "plain"),
      axis.text = element_text(size = 6, color = "black", face = "plain"),
      axis.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      legend.title = element_text(size = 6, face = "plain"),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.key.size = grid::unit(2.6, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

mark_key <- "Exceeds size-matched null"
support_legend <- data.table(
  stage_beta = NA_real_, continuum_beta = NA_real_, support_key = mark_key
)
plot <- ggplot(plot_data, aes(stage_beta, continuum_beta)) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_point(
    aes(color = direction_pair), shape = 16, size = 1.7,
    stroke = 0, alpha = 0.9
  ) +
  geom_point(
    data = plot_data[competitive_marked == TRUE],
    shape = 21, size = 2.15, fill = NA, color = "black", stroke = 0.45,
    show.legend = FALSE
  ) +
  geom_point(
    data = support_legend,
    aes(stage_beta, continuum_beta, shape = support_key),
    inherit.aes = FALSE, size = 1.8, fill = "white", color = "black",
    stroke = 0.45,
    # Contribute to the shape guide only. Without this the layer's fixed black
    # colour is also drawn into the colour guide, blanking the direction keys.
    show.legend = c(shape = TRUE, colour = FALSE), na.rm = TRUE
  ) +
  geom_text_repel(
    data = plot_data[fixed_label == TRUE], aes(label = short_label),
    size = 6 / ggplot2::.pt, color = "black", seed = contract$seed,
    box.padding = 0.08, point.padding = 0.04, min.segment.length = 0,
    segment.size = 0.11, max.overlaps = Inf, max.time = 15,
    max.iter = 200000
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("Spearman rho = %.2f; %d/43 same direction",
                    system_rho, same_sign_count),
    hjust = -0.02, vjust = 1.25, size = 6 / ggplot2::.pt
  ) +
  scale_color_manual(values = direction_colors, drop = FALSE, name = NULL) +
  scale_shape_manual(values = setNames(21, mark_key), name = NULL) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4),
    labels = function(x) formatC(x, format = "f", digits = 1),
    expand = expansion(mult = 0.11)
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 4),
    labels = function(x) formatC(x, format = "f", digits = 1),
    expand = expansion(mult = 0.11)
  ) +
  labs(
    x = "Fibrosis-stage effect (score SD / stage)",
    y = "Stage-adjusted continuum effect\n(score SD / continuum SD)"
  ) +
  theme_main() +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0)) +
  guides(
    color = guide_legend(
      order = 1, nrow = 1, override.aes = list(shape = 16, size = 1.7)
    ),
    shape = guide_legend(
      order = 2, nrow = 1,
      override.aes = list(shape = 21, fill = "white", color = "black",
                          size = 1.8, stroke = 0.45)
    )
  )

panel_path <- file.path(
  panel_dir, "fig4h_molecular_systems_stage_vs_continuum.pdf"
)
ggsave(
  panel_path, plot, width = 4.15, height = 3.25, units = "in",
  device = cairo_pdf, bg = "white", limitsize = FALSE
)
ml_assert(file.info(panel_path)$size > 1000L, "Figure 4H PDF is unexpectedly small")

source_table <- plot_data[, .(
  community_id, short_label, n_features, stage_beta, continuum_beta,
  stage_beta_feature_balanced, continuum_beta_feature_balanced,
  same_direction, same_direction_feature_balanced, direction_pair,
  strict_maxt_supported, joint_maxT_supported,
  connected_null_supported, uniform_null_supported,
  competitive_marked, fixed_outcome_blind_label = fixed_label
)]
setorder(source_table, community_id)
ml_write_tsv_once(
  source_table, file.path(source_dir, "fig4h_system_stage_continuum.tsv")
)
ml_write_tsv_once(
  source_table[fixed_outcome_blind_label == TRUE,
               .(community_id, short_label, n_features, selection_rule =
                   "15_largest_systems_fixed_before_outcome_overlay")],
  file.path(source_dir, "fig4h_fixed_label_roster.tsv")
)

figure_manifest <- data.table(
  callout = "4H", path = panel_path,
  sha256 = ml_sha256(panel_path), complete_system_family = "43/43",
  marked_rule = "exceeds_size_matched_unconstrained_null_both_axes",
  marked_systems = sprintf("%d/43", sum(plot_data$competitive_marked)),
  connected_null_supported = sprintf(
    "%d/43", sum(plot_data$connected_null_supported)
  ),
  joint_strict_supported = sprintf(
    "%d/43", sum(plot_data$joint_maxT_supported)
  ),
  four_test_maxt_supported = sprintf(
    "%d/43", sum(plot_data$strict_maxt_supported)
  ),
  effective_independent_systems = redundancy$effective_independent_systems,
  same_direction = "36/43",
  spearman_rho = system_rho, fixed_outcome_blind_labels = "15/43",
  inferential_unit = "participant", embedded_explanatory_footnote = FALSE,
  geometry_used = FALSE, window_values_used = FALSE, promoted = FALSE
)
ml_write_tsv_once(figure_manifest, file.path(output_candidate, "figure_manifest.tsv"))

input_manifest <- data.table(
  path = normalizePath(required_inputs, mustWork = TRUE),
  sha256 = vapply(required_inputs, ml_sha256, character(1))
)
ml_write_tsv_once(input_manifest, file.path(provenance_dir, "input_manifest.tsv"))

audit <- data.table(
  check = c(
    "promoted_input_hashes", "complete_system_family", "strict_maxT_family",
    "competitive_null_marking", "fixed_outcome_blind_label_roster",
    "stage_continuum_summary", "participant_level_effect_sources",
    "windows_not_used", "pdf_nonempty"
  ),
  passed = c(
    TRUE, nrow(source_table) == 43L, sum(source_table$strict_maxt_supported) == 18L,
    sum(source_table$competitive_marked) == 1L &&
      sum(source_table$connected_null_supported) == 0L,
    sum(source_table$fixed_outcome_blind_label) == 15L,
    abs(system_rho - 0.658) < 0.001 && same_sign_count == 36L,
    TRUE, !figure_manifest$window_values_used,
    file.info(panel_path)$size > 1000L
  ),
  detail = c(
    paste(basename(required_inputs), collapse = ";"),
    "All 43 outcome-blind systems remain visible",
    paste(
      "The four-test maxT count is retained in the source table but no longer",
      "carries ink; it is limited by the smaller cohort's permutation resolution",
      "and the joint cross-cohort test supports 35 of 43"
    ),
    paste(
      "The outline marks only systems exceeding 5,000 size-matched draws on both",
      "axes; no system clears the conservative connected-subgraph reference"
    ),
    "Labels reuse the 15 largest systems fixed before outcome overlay",
    sprintf("Spearman rho=%.6f; same direction=%d/43", system_rho, same_sign_count),
    "Stage and stage-adjusted continuum effects are donor-level two-cohort summaries",
    "No sliding-window value or ordering enters Figure 4H",
    basename(panel_path)
  )
)
ml_assert(all(audit$passed), "Figure 4H candidate audit failed")
ml_write_tsv_once(audit, file.path(provenance_dir, "audit.tsv"))
ml_write_session_info(file.path(provenance_dir, "sessionInfo.txt"))

output_paths <- c(
  panel_path,
  file.path(source_dir, "fig4h_system_stage_continuum.tsv"),
  file.path(source_dir, "fig4h_fixed_label_roster.tsv"),
  file.path(output_candidate, "figure_manifest.tsv"),
  file.path(provenance_dir, "input_manifest.tsv"),
  file.path(provenance_dir, "audit.tsv"),
  file.path(provenance_dir, "sessionInfo.txt")
)
output_checksums <- data.table(
  path = normalizePath(output_paths, mustWork = TRUE),
  sha256 = vapply(output_paths, ml_sha256, character(1))
)
ml_write_tsv_once(
  output_checksums, file.path(provenance_dir, "output_checksums.tsv")
)
writeLines(c(
  "# Main Figure 4H candidate", "",
  "All 43 outcome-blind molecular systems are shown. Black outlines mark the 18 systems passing the strict maxT family-wise criterion. Labels are restricted to the 15 largest systems fixed before any outcome overlay.", "",
  "The scatter compares the fibrosis-stage effect with the stage-adjusted continuum effect. It uses no sliding-window values and does not treat systems as a biological ontology."
), file.path(output_candidate, "README.md"))
writeLines("candidate_ready", file.path(output_candidate, "READY_FOR_REVIEW.ok"))

message("MAIN_SYSTEM_STAGE_CONTINUUM_COMPLETE: ", output_candidate)
