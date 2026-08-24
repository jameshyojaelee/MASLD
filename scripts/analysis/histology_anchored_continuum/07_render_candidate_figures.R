#!/usr/bin/env Rscript

# Key message: the molecular continuum is calibrated to, and complements,
# recorded fibrosis while retaining within-stage donor-level variation.

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
out <- out_root()
figures_root <- file.path(out, "figures")
ensure_new_dir(figures_root)

decision_path <- file.path(out, "decision", "promotion_decision.tsv")
manifest_path <- read_env_path(
  "HAC_MANIFEST_PATH",
  file.path(project_root(), "figures/candidates/pi-figure-redesign-2026-08-13-v3",
            "analysis/stage_extensions/five_cohort_sample_manifest.tsv")
)
required <- c(
  decision_path,
  file.path(out, "unsupervised", "participant_scores.tsv"),
  file.path(out, "projection", "participant_scores.tsv"),
  file.path(out, "programs", "fixed_program_windows.tsv.gz"),
  file.path(out, "programs", "all_axis_score_agreement.tsv"),
  file.path(out, "programs", "program_meta_analysis.tsv"),
  manifest_path
)
assert_true(all(file.exists(required)),
            paste0("Figure input missing: ", paste(required[!file.exists(required)], collapse = ", ")))

decision <- fread(decision_path)
assert_true(nrow(decision) == 1L &&
              all(c("overall_pass", "figure_route") %in% names(decision)),
            "Promotion decision schema drift")
overall_pass <- as.logical(decision$overall_pass[[1L]])
route <- as.character(decision$figure_route[[1L]])
expected_route <- if (overall_pass) "figure3_proposal" else "figure_s3_only"
assert_true(identical(route, expected_route), "Figure route disagrees with the gate decision")

main_dir <- file.path(figures_root, route)
control_dir <- file.path(figures_root, "figure_s3_controls")
ensure_new_dir(main_dir)
ensure_new_dir(control_dir)

save_cairo_once <- function(plot, path, width, height) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite figure: ", path))
  grDevices::cairo_pdf(path, width = width, height = height)
  print(plot)
  grDevices::dev.off()
  assert_true(file.info(path)$size > 1000L, paste0("Figure is unexpectedly small: ", path))
  invisible(path)
}

axis_labels <- c(
  signature_pc1 = "145-gene cohort PC1",
  fixed_projection = "Fixed discovery projection",
  full_transcriptome_pc1 = "Full-transcriptome PC1",
  rf_fibrosis = "Discovery-trained RF",
  consensus_rank = "Consensus rank"
)
independent <- as.character(pre$cohorts$independent_primary)
co_primary <- as.character(pre$axes$co_primary)

unsupervised <- fread(file.path(out, "unsupervised", "participant_scores.tsv"))
projection <- fread(file.path(out, "projection", "participant_scores.tsv"))
scores <- unique(rbindlist(list(unsupervised, projection), fill = TRUE)[
  , .(sample_id, dataset, axis_id, axis_raw, axis_percentile)
], by = c("sample_id", "axis_id"))
metadata <- fread(manifest_path, na.strings = c("", "NA"))[
  , .(sample_id, dataset, fibrosis_stage)
]

# Candidate Fig. 3D replacement panel, or the same panel routed to Fig. S3 when
# any prespecified gate fails. Points are biological donors; black ticks are
# stage medians. F0 remains the invariant control gray.
calibration <- merge(
  scores[dataset %in% independent & axis_id %in% co_primary],
  metadata,
  by = c("sample_id", "dataset"), all.x = TRUE
)
calibration <- calibration[!is.na(fibrosis_stage) & is.finite(axis_percentile)]
calibration[, `:=`(
  fibrosis_factor = factor(fibrosis_stage, levels = 0:4),
  axis_label = factor(axis_labels[axis_id], levels = axis_labels[co_primary]),
  dataset = factor(dataset, levels = independent)
)]
calibration_medians <- calibration[, .(
  median_percentile = median(axis_percentile),
  n_donors = .N
), by = .(dataset, axis_label, fibrosis_factor)]
stage_colors <- c(
  "0" = "#9E9E9E", "1" = palette1[[12L]], "2" = palette1[[10L]],
  "3" = palette1[[4L]], "4" = "#C9265E"
)
set.seed(pre$seeds$master)
p_calibration <- ggplot(
  calibration,
  aes(x = fibrosis_factor, y = axis_percentile, color = fibrosis_factor)
) +
  geom_point(
    position = position_jitter(width = 0.13, height = 0, seed = pre$seeds$master),
    size = 0.75, alpha = 0.55, stroke = 0
  ) +
  geom_point(
    data = calibration_medians,
    aes(y = median_percentile), inherit.aes = TRUE,
    shape = 95, size = 4.2, color = "black"
  ) +
  facet_grid(axis_label ~ dataset) +
  scale_color_manual(values = stage_colors, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", "50", "100"), expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = "Recorded fibrosis stage", y = "Within-cohort continuum percentile", color = NULL) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(legend.position = "none")
calibration_path <- file.path(main_dir, "donor_continuum_calibration.pdf")
save_cairo_once(p_calibration, calibration_path, 5.5, 3.7)

# Prespecified 20-percentile windows are a display only. No window boundary or
# focal program is selected from differential-expression counts.
windows <- fread(file.path(out, "programs", "fixed_program_windows.tsv.gz"))
focal_ids <- unname(unlist(pre$programs$focal_programs, use.names = FALSE))
trajectory <- windows[
  dataset %in% independent & axis_id == "consensus_rank" & program_uid %in% focal_ids
]
assert_true(uniqueN(trajectory$window_center) == 9L,
            "Program trajectory does not contain the nine fixed windows")
trajectory[, `:=`(
  ci_low = mean_score - 1.96 * se_score,
  ci_high = mean_score + 1.96 * se_score,
  program_label = factor(module_name, levels = unique(module_name)),
  dataset = factor(dataset, levels = independent)
)]
cohort_colors <- setNames(c(palette1[[1L]], palette1[[6L]]), independent)
p_trajectory <- ggplot(
  trajectory,
  aes(x = 100 * window_center, y = mean_score, color = dataset, fill = dataset)
) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "#9E9E9E") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15,
              linewidth = 0, color = NA) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 0.8, stroke = 0) +
  facet_wrap(~program_label, scales = "free_y", nrow = 1) +
  scale_color_manual(values = cohort_colors) +
  scale_fill_manual(values = cohort_colors) +
  scale_x_continuous(breaks = c(10, 30, 50, 70, 90), limits = c(0, 100)) +
  labs(x = "Consensus continuum percentile", y = "Program score (cohort z)",
       color = "Cohort", fill = "Cohort") +
  theme_masld(base_size = 6) + theme_pub() +
  theme(legend.position = "top")
trajectory_path <- file.path(main_dir, "focal_program_trajectories.pdf")
save_cairo_once(p_trajectory, trajectory_path, 5.5, 2.6)

# All scorer comparisons, including the full-transcriptome control and the
# histology-anchored RF benchmark, remain supplementary regardless of gates.
agreement <- fread(file.path(out, "programs", "all_axis_score_agreement.tsv"))
axes <- names(axis_labels)
reverse_agreement <- agreement[, .(
  dataset, evaluation_role, axis_x = axis_y, axis_y = axis_x,
  n_donors, spearman_rho, pearson_r
)]
diagonal <- CJ(dataset = unique(agreement$dataset), axis_x = axes)[
  , `:=`(axis_y = axis_x, evaluation_role = ifelse(
    dataset %in% independent, "independent_validation", "source_overlap_replication"
  ), n_donors = NA_integer_, spearman_rho = 1, pearson_r = 1)]
agreement_square <- rbindlist(list(agreement, reverse_agreement, diagonal), fill = TRUE)
agreement_square[, `:=`(
  axis_x_label = factor(axis_labels[axis_x], levels = rev(axis_labels[axes])),
  axis_y_label = factor(axis_labels[axis_y], levels = axis_labels[axes])
)]
p_agreement <- ggplot(agreement_square,
                      aes(x = axis_y_label, y = axis_x_label, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = PUB_GEOM_TEXT) +
  facet_wrap(~dataset, nrow = 1) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C9265E",
                       midpoint = 0, limits = c(-1, 1), breaks = c(-1, 0, 1),
                       name = "Spearman rho",
                       guide = guide_colorbar(
                         direction = "horizontal", title.position = "left",
                         title.hjust = 0.5, barwidth = grid::unit(2.4, "cm"),
                         barheight = grid::unit(0.16, "cm")
                       )) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", legend.justification = "center")
agreement_path <- file.path(control_dir, "all_scorer_agreement.pdf")
save_cairo_once(p_agreement, agreement_path, 8.0, 3.5)

# Every frozen program remains visible. Missing tiles denote explicit
# untestability; black rings identify meta-analysis BH q < 0.05 (n = 117).
meta <- fread(file.path(out, "programs", "program_meta_analysis.tsv"))
assert_true(uniqueN(meta$program_uid) == pre$programs$family_size,
            "Program meta-analysis does not contain all 117 frozen programs")
meta[, label := paste0(module_name, " | ", cell_type, " | ",
                        substr(program_uid, nchar(program_uid) - 5L, nchar(program_uid)))]
program_order <- meta[axis_id %in% co_primary, .(
  order_value = {
    finite_beta <- beta[is.finite(beta)]
    if (length(finite_beta)) max(abs(finite_beta)) else NA_real_
  }
), by = .(program_uid, label)][order(order_value, na.last = TRUE), label]
meta[, `:=`(
  label = factor(label, levels = program_order),
  axis_label = factor(axis_labels[axis_id], levels = axis_labels[axes]),
  significant = is.finite(q_value) & q_value < pre$promotion_gates$program_meta_bh_maximum
)]
fill_limit <- max(abs(meta$beta[is.finite(meta$beta)]), na.rm = TRUE)
p_programs <- ggplot(meta, aes(x = axis_label, y = label, fill = beta)) +
  geom_tile(color = "white", linewidth = 0.15) +
  geom_point(data = meta[significant == TRUE], shape = 21, size = 1.5,
             stroke = 0.3, fill = NA, color = "black") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C9265E",
                       midpoint = 0, limits = c(-fill_limit, fill_limit),
                       breaks = c(-fill_limit, 0, fill_limit),
                       labels = sprintf("%.1f", c(-fill_limit, 0, fill_limit)),
                       na.value = "#E6E6E6", name = "Meta beta",
                       guide = guide_colorbar(
                         direction = "horizontal", title.position = "top",
                         title.hjust = 0.5, barwidth = grid::unit(3.0, "cm"),
                         barheight = grid::unit(0.16, "cm")
                       )) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top", legend.justification = "center")
programs_path <- file.path(control_dir, "all_117_programs.pdf")
save_cairo_once(p_programs, programs_path, 7.0, 18.0)

figure_manifest <- data.table(
  figure_id = c("donor_continuum_calibration", "focal_program_trajectories",
                "all_scorer_agreement", "all_117_programs"),
  role = c(if (overall_pass) "figure3_proposal" else "figure_s3_only",
           if (overall_pass) "figure3_proposal" else "figure_s3_only",
           "figure_s3_control", "figure_s3_control"),
  path = c(calibration_path, trajectory_path, agreement_path, programs_path),
  promotion_gate_passed = overall_pass,
  automatic_promotion = FALSE,
  useDingbats = FALSE
)
figure_manifest[, `:=`(
  bytes = as.numeric(file.size(path)),
  sha256 = vapply(path, sha256_file, character(1L))
)]
write_tsv_once(figure_manifest, file.path(figures_root, "figure_manifest.tsv"))
write_session_info(file.path(figures_root, "sessionInfo.txt"))
message("CANDIDATE_FIGURES_COMPLETE: ", figures_root)
