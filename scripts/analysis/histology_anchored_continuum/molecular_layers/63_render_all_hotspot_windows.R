#!/usr/bin/env Rscript
# KEY MESSAGE: The complete frozen Hotspot registry shows coordinated, cohort-replicated molecular remodeling across the externally anchored continuum.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

options(digits = 17, scipen = 999)
grDevices::pdf.options(useDingbats = FALSE)

fail <- function(...) stop(paste0(...), call. = FALSE)
assert_true <- function(value, message) {
  if (!isTRUE(value)) fail(message)
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
assert_true(length(args) == 2L,
            "Usage: 63_render_all_hotspot_windows.R ANALYSIS_CANDIDATE OUTPUT_CANDIDATE")
analysis_candidate <- normalizePath(args[[1L]], mustWork = TRUE)
output_candidate <- args[[2L]]
assert_true(!file.exists(output_candidate),
            paste0("Refusing to overwrite output candidate: ", output_candidate))

project_root <- Sys.getenv("MASLD_PROJECT_ROOT", unset = "")
assert_true(nzchar(project_root) && dir.exists(project_root),
            "MASLD_PROJECT_ROOT is unset or invalid")

contract_path <- file.path(
  project_root,
  "scripts/analysis/histology_anchored_continuum/molecular_layers/00_contract.json"
)
assert_true(requireNamespace("jsonlite", quietly = TRUE), "jsonlite is unavailable")
contract <- jsonlite::fromJSON(contract_path, simplifyVector = TRUE)

program_root <- file.path(analysis_candidate, "programs", "hotspot")
input_paths <- c(
  windows = file.path(program_root, "hotspot_fixed_windows.tsv.gz"),
  membership = file.path(program_root, "hotspot_continuum_membership.tsv"),
  testability = file.path(program_root, "hotspot_testability.tsv"),
  meta = file.path(program_root, "hotspot_meta_analysis.tsv"),
  registry = file.path(project_root, contract$program_registry),
  contract = contract_path
)
assert_true(all(file.exists(input_paths)), paste0(
  "Missing inputs: ", paste(names(input_paths)[!file.exists(input_paths)], collapse = ", ")
))

dir.create(file.path(output_candidate, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_candidate, "source_tables"), showWarnings = FALSE)
dir.create(file.path(output_candidate, "provenance"), showWarnings = FALSE)

sha256 <- function(path) {
  value <- system2("sha256sum", path, stdout = TRUE)
  assert_true(length(value) == 1L, paste0("sha256sum failed: ", path))
  sub("[[:space:]].*$", "", value)
}

write_tsv_once <- function(x, path) {
  assert_true(!file.exists(path), paste0("Refusing to overwrite: ", path))
  fwrite(x, path, sep = "\t", na = "NA", quote = FALSE)
  invisible(path)
}

windows <- fread(input_paths[["windows"]])
membership <- fread(input_paths[["membership"]])
testability <- fread(input_paths[["testability"]])
meta <- fread(input_paths[["meta"]])
registry <- fread(input_paths[["registry"]])

expected_programs <- as.integer(contract$program_family_size)
expected_cohorts <- as.character(contract$evaluation_cohorts)
expected_axes <- as.character(contract$co_primary_axes)
expected_centers <- as.numeric(contract$window_centers)

assert_true(nrow(registry) == expected_programs && uniqueN(registry$program_uid) == expected_programs,
            "Frozen registry is not the complete 117-program family")
assert_true(nrow(membership) == expected_programs && uniqueN(membership$program_uid) == expected_programs,
            "Continuum membership is not the complete 117-program family")
assert_true(setequal(unique(windows$dataset), expected_cohorts),
            "Window cohorts differ from the two evaluation cohorts")
assert_true(setequal(unique(windows$axis_id), "fixed_projection"),
            "Window table is not restricted to fixed projection")
assert_true(setequal(sort(unique(windows$window_center)), expected_centers),
            "Window centers differ from the nine frozen centers")
assert_true(all(abs(windows$window_upper - windows$window_lower - contract$window_width) < 1e-12),
            "Window width differs from the frozen 20-percentile width")
assert_true(!anyDuplicated(windows[, .(program_uid, dataset, window_id)]),
            "Duplicate program/cohort/window rows")
assert_true(nrow(windows) == expected_programs * length(expected_cohorts) * length(expected_centers),
            "Window table does not contain 117 programs x 2 cohorts x 9 windows")
assert_true(all(windows$visualization_only), "A fixed-window row is not marked visualization-only")
assert_true(nrow(testability) == expected_programs * length(expected_cohorts),
            "Testability table is incomplete")
assert_true(nrow(meta) == expected_programs * length(expected_axes),
            "Meta-analysis table is incomplete")
assert_true(setequal(unique(meta$axis_id), expected_axes), "Co-primary meta-analysis axes drifted")

meta_check <- meta[, .(
  both_axes_meta_q_check = .N == length(expected_axes) &&
    all(is.finite(fixed_q_value) & fixed_q_value < 0.05),
  all_four_direction_check = .N == length(expected_axes) && all(direction_concordant)
), by = program_uid]
meta_check[, continuum_associated_check :=
             both_axes_meta_q_check & all_four_direction_check]
membership_check <- merge(
  membership[, .(program_uid, continuum_associated)], meta_check,
  by = "program_uid", all = TRUE
)
assert_true(all(membership_check$continuum_associated ==
                  membership_check$continuum_associated_check),
            "Continuum membership does not rederive from the two-axis gate")

registry <- registry[, .(
  program_uid, registry_order = .I, cell_type, module,
  registry_module_name = module_name
)]
membership <- merge(membership, registry, by = c("program_uid", "cell_type", "module"),
                    all.x = TRUE, sort = FALSE)
assert_true(!anyNA(membership$registry_order), "A continuum program is absent from the frozen registry")
assert_true(all(membership$module_name == membership$registry_module_name),
            "Program names differ between continuum results and frozen registry")

testability_program <- testability[, .(
  testable_both_cohorts = .N == length(expected_cohorts) && all(testable),
  testable_cohort_count = sum(testable),
  minimum_retained_l1_fraction = min(retained_l1_fraction),
  minimum_observed_genes = min(n_observed_genes)
), by = program_uid]

effect_wide <- dcast(
  meta,
  program_uid ~ axis_id,
  value.var = c("fixed_beta", "fixed_se", "fixed_ci_low", "fixed_ci_high",
                "fixed_q_value", "direction_concordant")
)
programs <- Reduce(
  function(x, y) merge(x, y, by = "program_uid", all = TRUE, sort = FALSE),
  list(membership, testability_program, effect_wide)
)
assert_true(nrow(programs) == expected_programs && !anyNA(programs$registry_order),
            "Program annotation join changed the frozen family")

cell_labels <- c(
  cholangiocytes = "Cholangiocytes",
  fibroblasts = "Fibroblasts",
  hepatocytes = "Hepatocytes",
  macrophages = "Macrophages",
  tcells = "T cells"
)
assert_true(all(programs$cell_type %in% names(cell_labels)), "Unexpected cell-type label")
programs[, program_label := sprintf(
  "%s %02d · %s", cell_labels[cell_type], as.integer(module), module_name
)]
programs[, fixed_beta := fixed_beta_fixed_projection]
programs[, status := fifelse(
  !testable_both_cohorts, "Untestable",
  fifelse(
    continuum_associated & fixed_beta > 0, "Supported ↑",
    fifelse(continuum_associated & fixed_beta < 0,
            "Supported ↓", "Not supported")
  )
)]
programs[, status := factor(
  status,
  levels = c("Supported ↑", "Supported ↓", "Not supported", "Untestable")
)]

windows_for_plot <- copy(windows)
windows_for_plot[, c("module_name", "robust_display") := NULL]
plot_data <- merge(
  windows_for_plot,
  programs[, .(
    program_uid, registry_order, program_label, cell_type, module, module_name,
    continuum_associated, testable_both_cohorts, status,
    fixed_beta_signature_pc1, fixed_ci_low_signature_pc1,
    fixed_ci_high_signature_pc1, fixed_q_value_signature_pc1,
    direction_concordant_signature_pc1, fixed_beta_fixed_projection,
    fixed_ci_low_fixed_projection, fixed_ci_high_fixed_projection,
    fixed_q_value_fixed_projection, direction_concordant_fixed_projection
  )],
  by = "program_uid", all.x = TRUE, sort = FALSE
)
plot_data <- merge(
  plot_data,
  testability[, .(program_uid, dataset, cohort_testable = testable)],
  by = c("program_uid", "dataset"), all.x = TRUE, sort = FALSE
)
assert_true(nrow(plot_data) == nrow(windows) && !anyNA(plot_data$registry_order),
            "Plot-data joins changed the complete window table")
plot_data[, plotted_mean_score := fifelse(cohort_testable, mean_score, NA_real_)]

window_counts <- plot_data[cohort_testable == TRUE, .(
  n_participants = max(n_participants)
), by = .(dataset, window_id, window_center)]
assert_true(nrow(window_counts) == length(expected_cohorts) * length(expected_centers),
            "Window participant-count header is incomplete")
window_counts[, window_label := sprintf("%d%%\nn=%d", round(100 * window_center), n_participants)]
plot_data <- merge(
  plot_data, window_counts[, .(dataset, window_id, window_label)],
  by = c("dataset", "window_id"), all.x = TRUE, sort = FALSE
)
plot_data[, window_label := factor(
  window_label,
  levels = unique(window_counts[order(window_id), window_label])
)]

source_columns <- c(
  "program_uid", "registry_order", "program_label", "cell_type", "module", "module_name",
  "dataset", "window_id", "window_center", "window_lower", "window_upper",
  "n_participants", "mean_score", "se_score", "cohort_testable",
  "continuum_associated", "status", "fixed_beta_signature_pc1",
  "fixed_ci_low_signature_pc1", "fixed_ci_high_signature_pc1",
  "fixed_q_value_signature_pc1", "direction_concordant_signature_pc1",
  "fixed_beta_fixed_projection", "fixed_ci_low_fixed_projection",
  "fixed_ci_high_fixed_projection", "fixed_q_value_fixed_projection",
  "direction_concordant_fixed_projection", "visualization_only"
)
write_tsv_once(
  plot_data[, ..source_columns],
  file.path(output_candidate, "source_tables", "all117_hotspot_continuum_windows.tsv.gz")
)
write_tsv_once(
  programs[order(registry_order)],
  file.path(output_candidate, "source_tables", "hotspot_program_display_registry.tsv")
)

font_family <- "Helvetica"
text_size <- 6
theme_hotspot <- function() {
  theme_classic(base_size = text_size, base_family = font_family) %+replace%
    theme(
      text = element_text(size = text_size, face = "plain"),
      axis.text = element_text(size = text_size, color = "black", face = "plain"),
      axis.title = element_text(size = text_size, face = "plain"),
      legend.text = element_text(size = text_size, face = "plain"),
      legend.title = element_text(size = text_size, face = "plain"),
      strip.text = element_text(size = text_size, face = "plain"),
      strip.background = element_blank(),
      plot.title = element_blank(),
      axis.line = element_line(linewidth = 0.25, color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      legend.key.size = grid::unit(2.5, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

status_colors <- c(
  "Supported ↑" = "#C9265E",
  "Supported ↓" = "#1565C0",
  "Not supported" = "#9E9E9E",
  "Untestable" = "#E0E0E0"
)
axis_colors <- c("Cohort PC1" = "#9E9E9E", "Fixed projection" = "#00695C")

render_package <- function(selected_ids, filename, height) {
  selected_programs <- programs[program_uid %in% selected_ids][order(registry_order)]
  selected_windows <- plot_data[program_uid %in% selected_ids]
  labels <- selected_programs$program_label
  selected_windows[, program_factor := factor(program_label, levels = rev(labels))]
  selected_programs[, program_factor := factor(program_label, levels = rev(labels))]

  finite_values <- selected_windows[is.finite(plotted_mean_score), plotted_mean_score]
  assert_true(length(finite_values) > 0L, "No finite program scores to plot")
  clip <- as.numeric(quantile(abs(finite_values), 0.98, names = FALSE, type = 8))
  assert_true(is.finite(clip) && clip > 0, "Invalid heatmap clipping limit")

  heatmap <- ggplot(
    selected_windows,
    aes(window_label, program_factor, fill = plotted_mean_score)
  ) +
    geom_tile(color = "white", linewidth = 0.08) +
    facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = "#1565C0", mid = "white", high = "#C9265E", midpoint = 0,
      limits = c(-clip, clip), oob = scales::squish,
      na.value = "#E0E0E0", name = "Mean program score\n(cohort z)"
    ) +
    labs(
      x = "Fixed-projection window center (overlapping 20-percentile windows)",
      y = "Frozen Hotspot program"
    ) +
    theme_hotspot() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.spacing.x = grid::unit(1.5, "mm"),
      legend.position = "bottom"
    )

  effect <- rbindlist(list(
    selected_programs[, .(
      program_factor, axis = "Cohort PC1", beta = fixed_beta_signature_pc1,
      ci_low = fixed_ci_low_signature_pc1, ci_high = fixed_ci_high_signature_pc1
    )],
    selected_programs[, .(
      program_factor, axis = "Fixed projection", beta = fixed_beta_fixed_projection,
      ci_low = fixed_ci_low_fixed_projection, ci_high = fixed_ci_high_fixed_projection
    )]
  ))
  effect[, axis := factor(axis, levels = c("Cohort PC1", "Fixed projection"))]
  forest <- ggplot(effect, aes(beta, program_factor, color = axis, shape = axis)) +
    geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.25) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.25,
                  position = position_dodge(width = 0.45), orientation = "y",
                  na.rm = TRUE) +
    geom_point(size = 0.8, position = position_dodge(width = 0.45), na.rm = TRUE) +
    facet_grid(. ~ "Adjusted meta-effect") +
    scale_color_manual(values = axis_colors, name = NULL) +
    scale_shape_manual(values = c("Cohort PC1" = 16, "Fixed projection" = 17), name = NULL) +
    labs(x = "Meta-analysis beta\n(program SD / continuum SD)", y = NULL) +
    theme_hotspot() +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      legend.position = "bottom", legend.box = "horizontal"
    )

  status <- ggplot(selected_programs, aes("Gate", program_factor, fill = status)) +
    geom_tile(color = "white", linewidth = 0.08) +
    facet_grid(. ~ "Gate") +
    scale_fill_manual(values = status_colors, drop = FALSE, name = "Continuum gate") +
    labs(x = NULL, y = NULL) +
    theme_hotspot() +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
      legend.position = "bottom"
    )

  assembled <- heatmap + forest + status +
    plot_layout(widths = c(5.0, 1.65, 0.38), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.box.just = "left")
  path <- file.path(output_candidate, "panels", filename)
  assert_true(!file.exists(path), paste0("Refusing to overwrite: ", path))
  ggsave(path, assembled, width = 7.0, height = height, units = "in",
         device = cairo_pdf, bg = "white", limitsize = FALSE)
  data.table(
    artifact = sub("\\.pdf$", "", filename), path = normalizePath(path),
    n_programs = length(selected_ids), n_cohorts = length(expected_cohorts),
    n_windows = length(expected_centers), color_clip_abs_z = clip,
    source_rows = nrow(selected_windows)
  )
}

all_manifest <- render_package(
  programs[order(registry_order), program_uid],
  "s4_all117_hotspot_continuum_windows_labeled.pdf",
  13.0
)
supported_manifest <- render_package(
  programs[continuum_associated == TRUE & testable_both_cohorts == TRUE][
    order(registry_order), program_uid
  ],
  "s4_continuum_associated_hotspot_windows_labeled.pdf",
  10.5
)
figure_manifest <- rbindlist(list(all_manifest, supported_manifest))
figure_manifest[, sha256 := vapply(path, sha256, character(1))]
write_tsv_once(figure_manifest, file.path(output_candidate, "figure_manifest.tsv"))

input_manifest <- data.table(
  input = names(input_paths), path = normalizePath(input_paths),
  sha256 = vapply(input_paths, sha256, character(1))
)
write_tsv_once(input_manifest, file.path(output_candidate, "provenance", "input_manifest.tsv"))

summary <- data.table(
  metric = c(
    "frozen_programs", "continuum_associated_programs", "testable_programs",
    "untestable_programs", "evaluation_cohorts", "windows_per_cohort",
    "window_width_percentile", "adjacent_window_overlap", "window_role"
  ),
  value = c(
    expected_programs,
    sum(programs$continuum_associated),
    sum(programs$testable_both_cohorts),
    sum(!programs$testable_both_cohorts),
    length(expected_cohorts),
    length(expected_centers),
    contract$window_width,
    0.5,
    "descriptive_visualization_only"
  )
)
write_tsv_once(summary, file.path(output_candidate, "summary.tsv"))

readme <- c(
  "# Complete Hotspot continuum-window display candidate",
  "",
  "The complete panel labels all 117 frozen Hotspot programs and shows their mean cohort-standardized scores across the nine prespecified fixed-projection windows in each evaluation cohort. The second panel retains only programs that pass the frozen continuum-associated criterion.",
  "",
  "Window columns are overlapping descriptive summaries: each window spans 20 percentile points, adjacent windows overlap by 50%, and participants may therefore contribute to neighboring columns. Window means do not define significance.",
  "",
  "The support rail uses the prespecified complete-family criteria: BH q<0.05 for both co-primary axes and concordant direction in both cohorts for both axes. The effect forest reports fixed-effect inverse-variance meta-analysis coefficients from `program_z ~ continuum_z + factor(fibrosis_stage) + inferred_sex`.",
  "",
  "This is a timestamped candidate for Figure S4. It does not alter the frozen program registry, program weights, Figure 3, or the synchronized manuscript release."
)
writeLines(readme, file.path(output_candidate, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(output_candidate, "provenance", "sessionInfo.txt"))

message("Rendered complete and continuum-associated Hotspot window panels: ", output_candidate)
