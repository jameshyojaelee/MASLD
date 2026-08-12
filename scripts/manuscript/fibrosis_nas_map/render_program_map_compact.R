#!/usr/bin/env Rscript

# KEY MESSAGE: Outcome-locked programs show broader conditional fibrosis support
# than NAS support, with no testable program supported on both axes.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

fail <- function(...) stop(paste0(...), call. = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  fail("Usage: render_program_map_compact.R WORKSTREAM_ROOT OUTPUT_ROOT")
}

workstream_root <- normalizePath(args[[1L]], mustWork = TRUE)
output_root <- normalizePath(args[[2L]], mustWork = FALSE)
project_root <- normalizePath(file.path(dirname(workstream_root), "../../../../../.."),
                              mustWork = TRUE)
source(file.path(project_root, "scripts/figures/publication_theme.R"))
grDevices::pdf.options(useDingbats = FALSE)

input_path <- file.path(workstream_root, "discovery", "program_map.tsv")
validation_path <- file.path(workstream_root, "discovery", "validation_summary.tsv")
paired_path <- file.path(workstream_root, "holdout", "paired_validation_summary.tsv")
if (!file.exists(input_path)) fail("Missing sealed source table: ", input_path)
if (!file.exists(validation_path)) fail("Missing sealed validation table: ", validation_path)
if (!file.exists(paired_path)) fail("Missing sealed paired-validation table: ", paired_path)
if (dir.exists(output_root) || file.exists(output_root)) {
  fail("Refusing to overwrite compact figure release: ", output_root)
}

tmp_root <- paste0(output_root, ".tmp.", Sys.getpid())
if (!dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)) {
  fail("Could not create temporary output directory: ", tmp_root)
}
published <- FALSE
on.exit(if (!published && dir.exists(tmp_root)) unlink(tmp_root, recursive = TRUE), add = TRUE)

effect <- fread(input_path)
validation <- fread(validation_path)
paired <- fread(paired_path)
required <- c(
  "feature_id", "cell_type", "module", "module_name", "robust_display",
  "weight_coverage", "testable", "beta_meta_fibrosis", "beta_meta_nas",
  "se_meta_fibrosis", "se_meta_nas", "ci_lower_fibrosis", "ci_lower_nas",
  "ci_upper_fibrosis", "ci_upper_nas", "p_value_meta_fibrosis",
  "p_value_meta_nas", "q_value_fibrosis", "q_value_nas", "evidence_state",
  "robust_fibrosis", "robust_nas", "linearity_reversal_fibrosis",
  "linearity_reversal_nas"
)
missing <- setdiff(required, names(effect))
if (length(missing)) fail("Program map is missing columns: ", paste(missing, collapse = ", "))
if (nrow(effect) != 117L) fail("Expected the complete 117-program registry; found ", nrow(effect))
if (anyDuplicated(effect$feature_id)) fail("Duplicate feature_id in program map")
if (sum(effect$robust_display, na.rm = TRUE) != 2L) {
  fail("Expected exactly two frozen robust_display programs")
}

loco_global <- validation[
  split == "leave_one_cohort_out_global" & view == "hotspot_program"
]
if (nrow(loco_global) != 2L || !setequal(loco_global$axis, c("fibrosis", "nas"))) {
  fail("Expected one global leave-one-cohort-out result per conditional axis")
}
if (any(loco_global$positive_folds != 4L) || any(loco_global$p_holm >= 0.05)) {
  fail("The impact panel requires both locked axes to pass all four held-out cohorts")
}
discovery_n <- unique(loco_global$biological_n)
if (length(discovery_n) != 1L || discovery_n != 469L) {
  fail("Expected 469 discovery biological samples; found ", paste(discovery_n, collapse = ", "))
}
paired_primary <- paired[
  split == "paired_biopsy_changed" & view == "hotspot_program" &
    score_definition == "weighted_mean_z" & axis == "joint_fibrosis_nas"
]
if (nrow(paired_primary) != 1L || paired_primary$n_participants != 46L) {
  fail("Expected the locked primary paired validation for 46 changed participants")
}

effect[, display_group := fcase(
  !testable, "Untestable",
  evidence_state == "fibrosis_supported_nas_unsupported", "Fibrosis",
  evidence_state == "nas_supported_fibrosis_unsupported", "NAS",
  evidence_state == "both_axes_same_direction", "Both axes",
  evidence_state == "both_axes_opposing", "Opposing axes",
  evidence_state == "unsupported", "Unsupported",
  default = NA_character_
)]
if (anyNA(effect$display_group)) {
  fail("Unmapped evidence state(s): ",
       paste(unique(effect[is.na(display_group), evidence_state]), collapse = ", "))
}

map_data <- effect[
  testable == TRUE & is.finite(beta_meta_fibrosis) & is.finite(beta_meta_nas)
]
if (nrow(map_data) != sum(effect$testable)) {
  fail("Every testable program must have finite fibrosis and NAS coefficients")
}

group_order <- c("Both axes", "Opposing axes", "Fibrosis", "NAS", "Unsupported")
observed_groups <- group_order[group_order %in% map_data$display_group]
counts <- map_data[, .N, by = display_group]
legend_group_names <- c(
  "Both axes" = "Both axes",
  "Opposing axes" = "Opposing axes",
  Fibrosis = "Fibrosis only",
  NAS = "NAS only",
  Unsupported = "Unsupported"
)
legend_labels <- setNames(
  sprintf("%s (%d)", unname(legend_group_names[observed_groups]),
          counts$N[match(observed_groups, counts$display_group)]),
  observed_groups
)

group_colors <- c(
  "Both axes" = "#00695C",
  "Opposing axes" = "#7B1FA2",
  Fibrosis = "#C9265E",
  NAS = "#1565C0",
  Unsupported = "#BDBDBD"
)
group_shapes <- c(
  "Both axes" = 15,
  "Opposing axes" = 18,
  Fibrosis = 16,
  NAS = 17,
  Unsupported = 1
)

map_data[, display_group := factor(display_group, levels = observed_groups)]
map_data[, display_rank := fifelse(display_group == "Unsupported", 1L, 2L)]
setorder(map_data, display_rank, feature_id)
map_data[, point_size := fifelse(display_group == "Unsupported", 1.25, 1.8)]
map_data[, point_alpha := fifelse(display_group == "Unsupported", 0.58, 0.95)]

named <- copy(map_data[robust_display == TRUE])
if (nrow(named) != 2L || any(!named$display_group %in% c("Fibrosis", "NAS", "Both axes", "Opposing axes"))) {
  fail("Both frozen display programs must be testable and supported")
}
named[, display_label := sub("\\s*\\([^()]++\\)\\s*$", "", module_name, perl = TRUE)]
named[, label_x := beta_meta_fibrosis + fifelse(grepl("Ductular injury", module_name), -0.010, 0.014)]
named[, label_y := beta_meta_nas + fifelse(grepl("Ductular injury", module_name), 0.021, -0.026)]
named[, label_hjust := fifelse(grepl("Ductular injury", module_name), 1.0, 0.0)]
named[, segment_xend := label_x + fifelse(label_hjust == 1.0, 0.006, -0.006)]

x_range <- range(map_data$beta_meta_fibrosis, finite = TRUE)
y_range <- range(map_data$beta_meta_nas, finite = TRUE)
x_limits <- c(x_range[[1L]] - 0.04 * diff(x_range),
              x_range[[2L]] + 0.12 * diff(x_range))
y_limits <- y_range + c(-1, 1) * 0.08 * diff(y_range)

p <- ggplot(map_data, aes(beta_meta_fibrosis, beta_meta_nas)) +
  geom_hline(yintercept = 0, color = "grey72", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey72", linewidth = 0.25) +
  geom_point(
    aes(color = display_group, shape = display_group,
        size = point_size, alpha = point_alpha),
    stroke = 0.45
  ) +
  geom_segment(
    data = named,
    aes(xend = segment_xend, yend = label_y),
    inherit.aes = TRUE, color = "grey45", linewidth = 0.25,
    show.legend = FALSE
  ) +
  geom_text(
    data = named,
    aes(x = label_x, y = label_y, label = display_label, hjust = label_hjust),
    inherit.aes = FALSE, size = 6 / ggplot2::.pt, color = "black",
    vjust = 0.5, show.legend = FALSE
  ) +
  scale_color_manual(
    values = group_colors, breaks = observed_groups, labels = legend_labels,
    name = NULL, drop = TRUE
  ) +
  scale_shape_manual(
    values = group_shapes, breaks = observed_groups, labels = legend_labels,
    name = NULL, drop = TRUE
  ) +
  scale_size_identity(guide = "none") +
  scale_alpha_identity(guide = "none") +
  coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "on") +
  labs(
    title = "Frozen programs reveal broader fibrosis-associated remodeling",
    subtitle = "No testable program supported on both; both axes reproduced in 4/4 held-out cohorts",
    x = "Fibrosis association\n(adjusted for NAS)",
    y = "NAS association\n(adjusted for fibrosis)"
  ) +
  theme_masld(base_size = 6) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.margin = margin(0, 0, 1, 0),
    legend.spacing.x = unit(3, "pt"),
    legend.key.width = unit(8, "pt"),
    legend.key.height = unit(7, "pt"),
    panel.grid = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(margin = margin(b = 1)),
    plot.subtitle = element_text(margin = margin(b = 1)),
    plot.margin = margin(2, 3, 2, 2)
  ) +
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(
        shape = unname(group_shapes[observed_groups]), size = 1.8, alpha = 1
      )
    ),
    shape = "none"
  )

map_pdf_path <- file.path(tmp_root, "fig3_program_effect_map_impact.pdf")
ggsave(map_pdf_path, p, width = 3.6, height = 3.0, device = cairo_pdf)

method_steps <- data.table(
  x = 1:4,
  number = as.character(1:4),
  label = c(
    "117 frozen scRNA programs\nfixed genes + L1 weights",
    "L1-weighted bulk score\nwithin-cohort gene z-scores",
    "Mutually adjusted effects\nfibrosis | NAS; NAS | fibrosis",
    "Random-effects meta-analysis\n4 cohorts; n = 469"
  )
)
p_method <- ggplot(method_steps, aes(x, 0)) +
  geom_segment(
    data = data.table(x = 1:3, xend = 2:4, y = 0.18, yend = 0.18),
    aes(x = x + 0.12, xend = xend - 0.12, y = y, yend = yend),
    inherit.aes = FALSE, color = "grey55", linewidth = 0.3,
    arrow = arrow(length = unit(3, "pt"), type = "closed")
  ) +
  geom_point(y = 0.18, shape = 21, size = 4.0, stroke = 0.4,
             fill = "white", color = "black") +
  geom_text(aes(y = 0.18, label = number), size = 6 / ggplot2::.pt) +
  geom_text(aes(y = -0.04, label = label), size = 6 / ggplot2::.pt,
            vjust = 1, lineheight = 0.95) +
  annotate(
    "text", x = 0.55, y = 0.42, hjust = 0,
    label = "Outcome-locked program scoring: no histology-based feature selection or score fitting",
    size = 6 / ggplot2::.pt
  ) +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(-0.35, 0.48), clip = "off") +
  theme_void(base_size = 6, base_family = "Helvetica") +
  theme(
    text = element_text(size = 6, family = "Helvetica", face = "plain"),
    plot.margin = margin(2, 2, 2, 2)
  )
method_pdf_path <- file.path(tmp_root, "fig3_program_score_method.pdf")
ggsave(method_pdf_path, p_method, width = 5.2, height = 1.15, device = cairo_pdf)

effect[, display_label := fifelse(
  robust_display,
  sub("\\s*\\([^()]++\\)\\s*$", "", module_name, perl = TRUE),
  NA_character_
)]
effect[, legend_label := fifelse(
  display_group %in% names(legend_labels),
  unname(legend_labels[display_group]),
  display_group
)]
source_columns <- c(
  "feature_id", "cell_type", "module", "module_name", "display_label",
  "robust_display", "weight_coverage", "testable", "beta_meta_fibrosis",
  "beta_meta_nas", "se_meta_fibrosis", "se_meta_nas", "ci_lower_fibrosis",
  "ci_lower_nas", "ci_upper_fibrosis", "ci_upper_nas",
  "p_value_meta_fibrosis", "p_value_meta_nas", "q_value_fibrosis",
  "q_value_nas", "evidence_state", "display_group", "legend_label",
  "robust_fibrosis", "robust_nas", "linearity_reversal_fibrosis",
  "linearity_reversal_nas"
)
fwrite(effect[, ..source_columns], file.path(tmp_root, "source_table.tsv"), sep = "\t")

summary_table <- data.table(
  metric = c(
    "frozen_programs", "testable_programs", "untestable_programs",
    "fibrosis_supported_only", "nas_supported_only", "both_axes_supported",
    "unsupported", "discovery_cohorts", "discovery_biological_n",
    "fibrosis_positive_loco_folds", "nas_positive_loco_folds",
    "fibrosis_loco_holm_p", "nas_loco_holm_p", "paired_participants",
    "paired_median_cosine", "paired_ci_lower",
    "paired_ci_upper", "paired_empirical_p"
  ),
  value = c(
    nrow(effect), sum(effect$testable), sum(!effect$testable),
    sum(effect$display_group == "Fibrosis"), sum(effect$display_group == "NAS"),
    sum(effect$display_group %in% c("Both axes", "Opposing axes")),
    sum(effect$display_group == "Unsupported"), 4L, discovery_n,
    loco_global[axis == "fibrosis", positive_folds],
    loco_global[axis == "nas", positive_folds],
    loco_global[axis == "fibrosis", p_holm],
    loco_global[axis == "nas", p_holm],
    paired_primary$n_participants, paired_primary$estimate,
    paired_primary$ci_lower, paired_primary$ci_upper, paired_primary$empirical_p
  ),
  source = c(
    rep(input_path, 7), rep(validation_path, 6), rep(paired_path, 5)
  )
)
fwrite(summary_table, file.path(tmp_root, "figure_summary.tsv"), sep = "\t")

caption <- paste(
  "Outcome-locked fibrosis-NAS molecular remodeling map. The primary score is",
  "the L1-weighted mean of within-cohort gene-standardized expression for each",
  "of 117 frozen single-cell Hotspot programs; no histology-based feature",
  "selection, program rediscovery, or score fitting was performed. Each cohort",
  "models sex-adjusted fibrosis and NAS effects jointly, followed by REML",
  "random-effects meta-analysis with Knapp-Hartung inference across four cohorts",
  "(n = 469 biological samples). Of 113 testable programs, 27 were BH-supported",
  "for fibrosis only, four for NAS only, none for both axes, and 82 were",
  "unsupported. Both effect vectors were positive in all four held-out cohorts",
  "(Holm P = 0.00020 for fibrosis and 0.0042 for NAS). Unsupported is not evidence",
  "of no association, and the assay-native axis magnitudes are not directly",
  "comparable. The locked primary paired-biopsy validation failed (median cosine",
  sprintf("%.3f, 95%% CI %.3f to %.3f; permutation P = %.3f), so the claim is limited",
          paired_primary$estimate, paired_primary$ci_lower,
          paired_primary$ci_upper, paired_primary$empirical_p),
  "to reproducible cross-sectional stage-associated remodeling, not longitudinal",
  "progression or prediction."
)
writeLines(caption, file.path(tmp_root, "figure_caption.txt"))

render_spec <- data.table(
  field = c(
    "key_message", "input", "registry_n", "testable_n", "untestable_n",
    "labelled_n", "map_canvas_inches", "method_canvas_inches", "font", "unsupported_color",
    "fibrosis_color", "nas_color", "technical_detail_location"
  ),
  value = c(
    "Outcome-locked programs show broader conditional fibrosis support than NAS support, with no testable program supported on both axes",
    input_path, nrow(effect), sum(effect$testable), sum(!effect$testable),
    sum(effect$robust_display), "3.6 x 3.0", "5.2 x 1.15", "Helvetica 6 pt plain",
    group_colors[["Unsupported"]], group_colors[["Fibrosis"]],
    group_colors[["NAS"]], "source table and supplement"
  )
)
fwrite(render_spec, file.path(tmp_root, "render_spec.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(tmp_root, "sessionInfo.txt"))

sha256 <- function(paths) {
  output <- system2("sha256sum", paths, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) fail("sha256sum failed: ", paste(output, collapse = "\n"))
  sub("  .*", "", output)
}
artifacts <- c(
  "fig3_program_effect_map_impact.pdf", "fig3_program_score_method.pdf",
  "source_table.tsv", "figure_summary.tsv", "figure_caption.txt",
  "render_spec.tsv", "sessionInfo.txt"
)
artifact_paths <- file.path(tmp_root, artifacts)
manifest <- data.table(
  artifact = artifacts,
  sha256 = sha256(artifact_paths),
  size_bytes = file.info(artifact_paths)$size
)
fwrite(manifest, file.path(tmp_root, "artifact_manifest.tsv"), sep = "\t")

if (!file.rename(tmp_root, output_root)) fail("Atomic publication failed: ", output_root)
published <- TRUE
message("Published compact program map: ", output_root)
