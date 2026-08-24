#!/usr/bin/env Rscript

# KEY MESSAGE: The externally anchored molecular continuum adds donor-level
# resolution within recorded fibrosis stages without replacing histology.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "lib_continuum.R"))
source(file.path(project_root(), "scripts/figures/publication_theme.R"))

pre <- read_prespec()
source_root <- read_env_path("HAC_SOURCE_CANDIDATE")
review_root_raw <- Sys.getenv("HAC_FIGURE_REVIEW_ROOT", unset = "")
assert_true(nzchar(review_root_raw), "HAC_FIGURE_REVIEW_ROOT is unset")
review_root <- normalizePath(review_root_raw, mustWork = FALSE)
ensure_new_dir(review_root)

panel_dir <- file.path(review_root, "panels")
source_dir <- file.path(review_root, "source_tables")
provenance_dir <- file.path(review_root, "provenance")
ensure_new_dir(panel_dir)
ensure_new_dir(source_dir)
ensure_new_dir(provenance_dir)

required <- c(
  file.path(source_root, "unsupervised", "participant_scores.tsv"),
  file.path(source_root, "projection", "participant_scores.tsv"),
  file.path(source_root, "programs", "fixed_program_windows.tsv.gz"),
  file.path(source_root, "programs", "program_fibrosis_models.tsv"),
  file.path(source_root, "programs", "program_meta_analysis.tsv"),
  file.path(source_root, "paired", "donor_deltas.tsv"),
  file.path(source_root, "paired", "paired_tests.tsv"),
  file.path(source_root, "figures", "figure_s3_controls", "all_scorer_agreement.pdf"),
  file.path(source_root, "figures", "figure_s3_controls", "all_117_programs.pdf")
)
assert_true(all(file.exists(required)), paste0(
  "Review input missing: ", paste(required[!file.exists(required)], collapse = ", ")
))

current_fig3d <- file.path(
  project_root(), "figures", "main", "fig3_bulk_transcriptomics", "panels",
  "fig3d_pca_fibrosis_gradient.pdf"
)
current_fig3d_source <- file.path(
  project_root(), "figures", "main", "fig3_bulk_transcriptomics", "source_tables",
  "fig3d_endpoint_pca_coordinates.tsv"
)
assert_true(file.exists(current_fig3d) && file.exists(current_fig3d_source),
            "Current Figure 3D reference or source table is missing")

axis_labels <- c(
  signature_pc1 = "145-gene cohort PC1",
  fixed_projection = "Fixed discovery projection",
  full_transcriptome_pc1 = "Full-transcriptome PC1",
  rf_fibrosis = "Discovery-trained RF",
  consensus_rank = "Consensus rank"
)
axis_roles <- c(
  signature_pc1 = "Co-primary",
  fixed_projection = "Co-primary",
  full_transcriptome_pc1 = "Unsupervised control",
  rf_fibrosis = "Histology-trained control",
  consensus_rank = "Presentation score"
)
axis_colors <- c(
  signature_pc1 = "#C9265E",
  fixed_projection = "#1565C0",
  full_transcriptome_pc1 = "#9E9E9E",
  rf_fibrosis = "#7B1FA2",
  consensus_rank = "#00695C"
)
cohort_colors <- c(GSE162694 = palette1[[1L]], GSE213621 = palette1[[6L]])
stage_colors <- c(
  "0" = "#9E9E9E", "1" = palette2[[12L]], "2" = palette2[[10L]],
  "3" = palette2[[7L]], "4" = "#C9265E"
)
independent <- as.character(pre$cohorts$independent_primary)
co_primary <- as.character(pre$axes$co_primary)
all_axes <- names(axis_labels)

theme_review <- function() {
  theme_masld(base_size = 6) +
    theme(
      text = element_text(size = 6, face = "plain"),
      plot.title = element_blank(),
      strip.text = element_text(size = 6, face = "plain"),
      legend.position = "top",
      legend.justification = "center",
      legend.margin = margin(0, 0, 1, 0),
      panel.spacing = grid::unit(0.8, "lines")
    )
}

save_panel <- function(plot, file_name, width, height) {
  path <- file.path(panel_dir, file_name)
  assert_true(!file.exists(path), paste0("Refusing to overwrite panel: ", path))
  grDevices::cairo_pdf(path, width = width, height = height)
  print(plot)
  grDevices::dev.off()
  assert_true(file.info(path)$size > 1000L, paste0("Panel is unexpectedly small: ", path))
  path
}

copy_panel <- function(source, file_name) {
  destination <- file.path(panel_dir, file_name)
  assert_true(!file.exists(destination), paste0("Refusing to overwrite panel: ", destination))
  assert_true(file.copy(source, destination, overwrite = FALSE, copy.mode = TRUE),
              paste0("Could not copy panel: ", source))
  destination
}

manifest_path <- read_env_path(
  "HAC_MANIFEST_PATH",
  file.path(project_root(), "figures", "candidates", "pi-figure-redesign-2026-08-13-v3",
            "analysis", "stage_extensions", "five_cohort_sample_manifest.tsv")
)
required <- c(required, manifest_path)
unsupervised_scores <- fread(
  file.path(source_root, "unsupervised", "participant_scores.tsv"),
  na.strings = c("", "NA")
)
projection_scores <- fread(
  file.path(source_root, "projection", "participant_scores.tsv"),
  na.strings = c("", "NA")
)
scores <- unique(rbindlist(list(
  unsupervised_scores[, .(sample_id, dataset, axis_id, axis_raw, axis_percentile)],
  projection_scores[, .(sample_id, dataset, axis_id, axis_raw, axis_percentile)]
), fill = TRUE), by = c("sample_id", "dataset", "axis_id"))
metadata <- fread(manifest_path, na.strings = c("", "NA"))[
  , .(sample_id, dataset, fibrosis_stage)
]
scores <- merge(scores, metadata, by = c("sample_id", "dataset"), all.x = TRUE)
scores <- scores[
  dataset %in% independent & axis_id %in% all_axes &
    is.finite(axis_percentile) & !is.na(fibrosis_stage)
]
scores[, `:=`(
  fibrosis_factor = factor(as.integer(fibrosis_stage), levels = 0:4),
  axis_label = factor(axis_labels[axis_id], levels = axis_labels[all_axes]),
  axis_role = axis_roles[axis_id],
  dataset = factor(dataset, levels = independent)
)]
score_census <- scores[, .(n_donors = uniqueN(sample_id)), by = dataset]
assert_true(
  identical(as.integer(score_census[match(independent, dataset), n_donors]), c(109L, 361L)),
  "Fibrosis-complete independent score roster no longer contains 109 and 361 participants"
)

# Main-sized addition candidate: both transported scores retain wide donor-level
# variation within every recorded stage while their medians rise with fibrosis.
calibration <- scores[axis_id %in% co_primary]
calibration_medians <- calibration[, .(
  median_percentile = median(axis_percentile),
  n_donors = .N
), by = .(dataset, axis_id, axis_label, fibrosis_factor)]
write_tsv_once(calibration, file.path(source_dir, "continuum_by_stage.tsv"))
write_tsv_once(calibration_medians, file.path(source_dir, "continuum_by_stage_medians.tsv"))

set.seed(pre$seeds$master)
p_calibration <- ggplot(
  calibration,
  aes(x = fibrosis_factor, y = axis_percentile, color = fibrosis_factor)
) +
  geom_point(
    position = position_jitter(width = 0.13, height = 0, seed = pre$seeds$master),
    size = 0.72, alpha = 0.52, stroke = 0
  ) +
  geom_point(
    data = calibration_medians,
    aes(y = median_percentile), shape = 95, size = 4.0, color = "black"
  ) +
  facet_grid(axis_label ~ dataset) +
  scale_color_manual(values = stage_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "50", "100"),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(x = "Recorded fibrosis stage", y = "Within-cohort continuum percentile") +
  theme_review() +
  theme(legend.position = "none")
panel_calibration <- save_panel(
  p_calibration, "fig3_candidate_continuum_by_stage.pdf", 5.5, 3.65
)

# Supplement-sized comparison: the co-primary scores, controls, and presentation
# score are shown under identical scales and donor rosters.
score_medians <- scores[, .(
  median_percentile = median(axis_percentile),
  n_donors = .N
), by = .(dataset, axis_id, axis_label, fibrosis_factor)]
write_tsv_once(scores, file.path(source_dir, "all_scores_by_stage.tsv"))
write_tsv_once(score_medians, file.path(source_dir, "all_scores_by_stage_medians.tsv"))

p_all_stage <- ggplot(
  scores,
  aes(x = fibrosis_factor, y = axis_percentile, color = fibrosis_factor)
) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = pre$seeds$master),
    size = 0.48, alpha = 0.35, stroke = 0
  ) +
  geom_point(
    data = score_medians,
    aes(y = median_percentile), shape = 95, size = 3.3, color = "black"
  ) +
  facet_grid(axis_label ~ dataset) +
  scale_color_manual(values = stage_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "50", "100"),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(x = "Recorded fibrosis stage", y = "Within-cohort score percentile") +
  theme_review() +
  theme(legend.position = "none")
panel_all_stage <- save_panel(
  p_all_stage, "figs3_candidate_all_scores_by_stage.pdf", 5.5, 7.4
)

# Recompute equal-footing donor-bootstrap intervals for all five axes. The
# co-primary rows reproduce the gate file; controls gain intervals for display.
bootstrap_rows <- list()
bootstrap_i <- 0L
for (axis in all_axes) {
  for (cohort in independent) {
    d <- scores[axis_id == axis & as.character(dataset) == cohort]
    assert_true(nrow(d) >= 100L, paste0("Too few donors for ", axis, " in ", cohort))
    set.seed(pre$seeds$bootstrap + match(axis, all_axes) * 100L + match(cohort, independent))
    boot <- replicate(pre$resampling$bootstrap_replicates, {
      idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
      suppressWarnings(cor(d$axis_percentile[idx], d$fibrosis_stage[idx], method = "spearman"))
    })
    bootstrap_i <- bootstrap_i + 1L
    bootstrap_rows[[bootstrap_i]] <- data.table(
      axis_id = axis,
      dataset = cohort,
      n_donors = nrow(d),
      spearman_rho = suppressWarnings(cor(
        d$axis_percentile, d$fibrosis_stage, method = "spearman"
      )),
      ci_low = quantile(boot[is.finite(boot)], 0.025, names = FALSE, type = 8),
      ci_high = quantile(boot[is.finite(boot)], 0.975, names = FALSE, type = 8),
      bootstrap_replicates = pre$resampling$bootstrap_replicates
    )
  }
}
anchors <- rbindlist(bootstrap_rows)
anchors[, `:=`(
  axis_label = factor(axis_labels[axis_id], levels = rev(axis_labels[all_axes])),
  axis_role = axis_roles[axis_id],
  dataset = factor(dataset, levels = independent)
)]
write_tsv_once(anchors, file.path(source_dir, "all_score_fibrosis_anchors.tsv"))

p_anchors <- ggplot(
  anchors,
  aes(x = spearman_rho, y = axis_label, color = axis_id, shape = dataset)
) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = axis_label),
    position = position_dodge(width = 0.45), linewidth = 0.45
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 1.8, stroke = 0.35) +
  scale_color_manual(values = axis_colors, guide = "none") +
  scale_shape_manual(values = c(GSE162694 = 16, GSE213621 = 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75)) +
  coord_cartesian(xlim = c(-0.15, 0.75)) +
  labs(x = "Spearman correlation with recorded fibrosis stage", y = NULL,
       shape = "Cohort") +
  theme_review()
panel_anchors <- save_panel(
  p_anchors, "figs3_candidate_score_fibrosis_anchor_forest.pdf", 4.5, 2.55
)

# Fixed visualization windows show timing only; no window enters inference.
windows <- fread(file.path(source_root, "programs", "fixed_program_windows.tsv.gz"))
focal_ids <- unname(unlist(pre$programs$focal_programs, use.names = FALSE))
trajectory <- windows[
  dataset %in% independent & axis_id == "consensus_rank" & program_uid %in% focal_ids
]
assert_true(uniqueN(trajectory$window_center) == 9L,
            "Focal trajectories no longer contain the nine fixed windows")
trajectory[, `:=`(
  ci_low = mean_score - 1.96 * se_score,
  ci_high = mean_score + 1.96 * se_score,
  program_label = factor(module_name, levels = unique(module_name)),
  dataset = factor(dataset, levels = independent)
)]
write_tsv_once(trajectory, file.path(source_dir, "focal_program_trajectories.tsv"))

p_trajectory <- ggplot(
  trajectory,
  aes(x = 100 * window_center, y = mean_score, color = dataset, fill = dataset)
) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "#9E9E9E") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.14,
              linewidth = 0, color = NA) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 0.8, stroke = 0) +
  facet_wrap(~program_label, scales = "free_y", nrow = 1) +
  scale_color_manual(values = cohort_colors) +
  scale_fill_manual(values = cohort_colors) +
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90), limits = c(0, 100)) +
  labs(x = "Consensus continuum percentile", y = "Program score (cohort z)",
       color = "Cohort", fill = "Cohort") +
  theme_review()
panel_trajectory <- save_panel(
  p_trajectory, "fig3_candidate_focal_program_trajectories.pdf", 5.5, 2.55
)

# This is the inferential focal-program display: coefficients are conditional
# on recorded fibrosis stage and inferred sex in each independent cohort.
models <- fread(file.path(source_root, "programs", "program_fibrosis_models.tsv"))[
  program_uid %in% focal_ids & axis_id %in% co_primary & dataset %in% independent
]
meta <- fread(file.path(source_root, "programs", "program_meta_analysis.tsv"))[
  program_uid %in% focal_ids & axis_id %in% co_primary
]
models_display <- models[, .(
  program_uid, module_name, axis_id, source = dataset, beta, ci_low, ci_high,
  q_value = NA_real_, direction_concordant = NA
)]
meta_display <- meta[, .(
  program_uid, module_name, axis_id, source = "Fixed-effect meta-analysis",
  beta, ci_low, ci_high, q_value, direction_concordant
)]
focal_forest <- rbindlist(list(models_display, meta_display), fill = TRUE)
focal_forest[, `:=`(
  program_label = factor(module_name, levels = unique(module_name)),
  axis_label = factor(axis_labels[axis_id], levels = axis_labels[co_primary]),
  source = factor(source, levels = c(independent, "Fixed-effect meta-analysis"))
)]
write_tsv_once(focal_forest, file.path(source_dir, "focal_program_adjusted_forest.tsv"))

p_focal_forest <- ggplot(
  focal_forest,
  aes(x = beta, y = source, color = source, shape = source)
) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = source), linewidth = 0.5) +
  geom_point(size = 1.8, stroke = 0.35) +
  facet_grid(program_label ~ axis_label) +
  scale_color_manual(values = c(
    GSE162694 = cohort_colors[["GSE162694"]],
    GSE213621 = cohort_colors[["GSE213621"]],
    "Fixed-effect meta-analysis" = "black"
  ), guide = "none") +
  scale_shape_manual(values = c(
    GSE162694 = 16, GSE213621 = 1, "Fixed-effect meta-analysis" = 18
  ), guide = "none") +
  labs(x = "Continuum coefficient adjusted for fibrosis stage and sex", y = NULL) +
  theme_review() +
  theme(legend.position = "none")
panel_focal_forest <- save_panel(
  p_focal_forest, "fig3_candidate_focal_program_adjusted_forest.pdf", 5.5, 3.15
)

# Paired replication retains stable-fibrosis participants; points are biological
# donors and black ticks are median continuum changes at each fibrosis change.
deltas <- fread(file.path(source_root, "paired", "donor_deltas.tsv"))[
  axis_id %in% co_primary & is.finite(delta_continuum) & is.finite(delta_fibrosis)
]
paired_tests <- fread(file.path(source_root, "paired", "paired_tests.tsv"))[
  axis_id %in% co_primary & endpoint == "delta_fibrosis"
]
deltas[, axis_label := factor(axis_labels[axis_id], levels = axis_labels[co_primary])]
delta_medians <- deltas[, .(
  median_delta_continuum = median(delta_continuum),
  n_donors = .N
), by = .(axis_id, axis_label, delta_fibrosis)]
paired_annotations <- paired_tests[, .(
  axis_id,
  axis_label = factor(axis_labels[axis_id], levels = axis_labels[co_primary]),
  label = sprintf("rho = %.2f; Holm P = %.3f", spearman_rho, holm_p_value)
)]
write_tsv_once(deltas, file.path(source_dir, "paired_delta_continuum.tsv"))
write_tsv_once(delta_medians, file.path(source_dir, "paired_delta_continuum_medians.tsv"))
write_tsv_once(paired_annotations, file.path(source_dir, "paired_delta_annotations.tsv"))

y_limits <- range(deltas$delta_continuum, finite = TRUE)
y_pad <- diff(y_limits) * 0.12
p_paired <- ggplot(
  deltas,
  aes(x = delta_fibrosis, y = delta_continuum)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#9E9E9E") +
  geom_point(
    aes(color = factor(sign(delta_fibrosis), levels = c(-1, 0, 1))),
    position = position_jitter(width = 0.09, height = 0, seed = pre$seeds$master),
    size = 0.8, alpha = 0.6, stroke = 0
  ) +
  geom_point(
    data = delta_medians,
    aes(y = median_delta_continuum), shape = 95, size = 4.0, color = "black"
  ) +
  geom_text(
    data = paired_annotations,
    aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -0.02, vjust = 1.2, size = PUB_GEOM_TEXT
  ) +
  facet_wrap(~axis_label, nrow = 1) +
  scale_color_manual(
    values = c("-1" = "#1565C0", "0" = "#9E9E9E", "1" = "#C9265E"),
    labels = c("Fibrosis decreased", "Fibrosis stable", "Fibrosis increased"),
    name = NULL
  ) +
  scale_x_continuous(breaks = sort(unique(deltas$delta_fibrosis))) +
  coord_cartesian(ylim = c(y_limits[[1L]] - y_pad, y_limits[[2L]] + y_pad)) +
  labs(x = "Within-participant change in recorded fibrosis stage",
       y = "Within-participant change in continuum z-score") +
  theme_review()
panel_paired <- save_panel(
  p_paired, "fig3_candidate_paired_delta_replication.pdf", 5.5, 2.75
)

# Validated controls and the current 3D panel are copied byte-for-byte so the
# review suite compares new additions against the exact existing artwork.
panel_agreement <- copy_panel(
  file.path(source_root, "figures", "figure_s3_controls", "all_scorer_agreement.pdf"),
  "figs3_reference_all_scorer_agreement.pdf"
)
panel_programs <- copy_panel(
  file.path(source_root, "figures", "figure_s3_controls", "all_117_programs.pdf"),
  "figs3_reference_all_117_programs.pdf"
)
panel_current <- copy_panel(current_fig3d, "reference_current_fig3d_endpoint_pca.pdf")
assert_true(file.copy(
  current_fig3d_source,
  file.path(source_dir, "reference_current_fig3d_endpoint_pca.tsv"),
  overwrite = FALSE, copy.mode = TRUE
), "Could not copy the current Figure 3D source table")

review <- data.table(
  panel_id = c(
    "continuum_by_stage", "focal_program_adjusted_forest", "paired_delta_replication",
    "focal_program_trajectories", "score_fibrosis_anchor_forest",
    "all_scores_by_stage", "all_scorer_agreement", "all_117_programs",
    "current_fig3d_endpoint_pca"
  ),
  role = c(
    "compact main-figure addition", "compact main-figure addition",
    "main-or-supplement replication", "main-or-supplement presentation",
    "supplementary scorer control", "supplementary scorer control",
    "supplementary method control", "supplementary complete family",
    "current-panel reference"
  ),
  replacement_potential = c(
    "high for current 3D only if space is required", "low", "medium", "low",
    "low", "low", "none", "none", "already current"
  ),
  decision_question = c(
    "Does the continuum add donor resolution within fibrosis stages?",
    "Do focal programs track the continuum after recorded-stage adjustment?",
    "Do within-participant molecular changes follow fibrosis changes?",
    "Where along the continuum do the focal programs rise?",
    "Do co-primary scores outperform the unsupervised control without hiding the RF benchmark?",
    "How do all five scores distribute across recorded stages?",
    "Are scorer relationships stable across cohorts?",
    "What happens across the complete frozen 117-program family?",
    "What does the current endpoint PCA already establish?"
  )
)
write_tsv_once(review, file.path(review_root, "panel_review.tsv"))

manifest <- data.table(
  panel_id = review$panel_id,
  role = review$role,
  path = c(
    panel_calibration, panel_focal_forest, panel_paired, panel_trajectory,
    panel_anchors, panel_all_stage, panel_agreement, panel_programs, panel_current
  ),
  source_kind = c(rep("new_render", 6L), rep("validated_byte_copy", 3L)),
  current_figure_modified = FALSE,
  automatic_promotion = FALSE,
  useDingbats = FALSE
)
manifest[, `:=`(
  bytes = as.numeric(file.size(path)),
  sha256 = vapply(path, sha256_file, character(1L))
)]
write_tsv_once(manifest, file.path(review_root, "figure_manifest.tsv"))

input_manifest <- data.table(
  input = c(required, current_fig3d, current_fig3d_source),
  sha256 = vapply(c(required, current_fig3d, current_fig3d_source),
                  sha256_file, character(1L))
)
write_tsv_once(input_manifest, file.path(provenance_dir, "input_manifest.tsv"))
write_session_info(file.path(provenance_dir, "sessionInfo.txt"))

message("FIGURE3_REVIEW_SUITE_COMPLETE: ", review_root)
