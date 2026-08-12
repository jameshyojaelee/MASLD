#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))
source(file.path(project_root, "scripts/figures/publication_theme.R"))
grDevices::pdf.options(useDingbats = FALSE)

if (!file.exists(file.path(discovery_root, "DISCOVERY_READY.json"))) fail("Discovery release missing")
if (!file.exists(file.path(holdout_root, "HOLDOUT_READY.json"))) fail("Holdout release missing")
if (dir.exists(figure_root)) fail("Refusing to overwrite figure release: ", figure_root)
tmp <- atomic_dir(figure_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

effect <- fread(file.path(discovery_root, "program_map.tsv"))
cohort <- fread(file.path(discovery_root, "cohort_effects.tsv"))
meta <- fread(file.path(discovery_root, "meta_effects.tsv"))
loco <- fread(file.path(discovery_root, "validation_summary.tsv"))
scores <- fread(file.path(discovery_root, "feature_scores.tsv.gz"))
paired <- fread(file.path(holdout_root, "paired_validation_summary.tsv"))
participants <- fread(file.path(holdout_root, "paired_participant_validation.tsv"))
diagnostics <- fread(file.path(holdout_root, "paired_diagnostics.tsv"))
transport <- fread(file.path(holdout_root, "transport_validation.tsv"))
paired_nulls <- readRDS(file.path(holdout_root, "paired_permutation_nulls.rds"))

lineage_colors <- c(
  hepatocytes = "#C9265E",
  macrophages = "#1565C0",
  fibroblasts = "#6A51A3",
  cholangiocytes = "#E69F00",
  `T cells` = "#00695C",
  t_cells = "#00695C",
  tcells = "#00695C"
)
effect[, lineage := as.character(cell_type)]
effect[!lineage %in% names(lineage_colors), lineage := "Other"]
lineage_colors <- c(lineage_colors, Other = "#757575")

map_data <- effect[testable == TRUE & is.finite(beta_meta_fibrosis) & is.finite(beta_meta_nas)]
untestable_text <- paste(
  effect[testable == FALSE, paste0(cell_type, ": ", module_name)],
  collapse = "; "
)
x_limits <- range(map_data$beta_meta_fibrosis, finite = TRUE)
y_limits <- range(map_data$beta_meta_nas, finite = TRUE)
x_limits <- x_limits + c(-1, 1) * diff(x_limits) * 0.08
y_limits <- y_limits + c(-1, 1) * diff(y_limits) * 0.08

p_map <- ggplot(map_data, aes(beta_meta_fibrosis, beta_meta_nas, color = lineage)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(aes(shape = evidence_state), size = 1.8, alpha = 0.88) +
  geom_text(
    data = map_data[robust_display == TRUE],
    aes(label = module_name), color = "black", size = 6 / ggplot2::.pt,
    hjust = -0.08, vjust = -0.45, check_overlap = TRUE, show.legend = FALSE
  ) +
  annotate(
    "text", x = x_limits[1], y = y_limits[2], hjust = 0, vjust = 1,
    label = sprintf("117 frozen programs; %d testable; %d untestable",
                    sum(effect$testable), sum(!effect$testable)),
    size = 6 / ggplot2::.pt, color = "black"
  ) +
  scale_color_manual(values = lineage_colors, name = "Lineage") +
  scale_shape_manual(values = c(
    both_axes_same_direction = 16,
    both_axes_opposing = 17,
    fibrosis_supported_nas_unsupported = 15,
    nas_supported_fibrosis_unsupported = 18,
    unsupported = 1
  ), name = "Evidence state") +
  coord_cartesian(xlim = x_limits, ylim = y_limits) +
  labs(
    x = "Conditional fibrosis association\n(program-score SD per fibrosis unit)",
    y = "Conditional NAS association\n(program-score SD per NAS unit)",
    caption = paste0("Untestable (<80% frozen weight): ", untestable_text)
  ) +
  theme_masld(base_size = 6) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 5))
ggsave(file.path(tmp, "fig3_program_effect_map.pdf"), p_map,
       width = 4.7, height = 3.6, device = cairo_pdf)

primary_pair <- paired[view == "hotspot_program" & score_definition == "weighted_mean_z"]
null <- paired_nulls[["hotspot_program::weighted_mean_z"]]
null_dt <- data.table(value = null)
p_pair <- ggplot(null_dt, aes(value)) +
  geom_density(fill = "grey85", color = "grey35", linewidth = 0.35) +
  geom_vline(xintercept = primary_pair$estimate, color = "#C9265E", linewidth = 0.7) +
  annotate(
    "text", x = primary_pair$estimate, y = Inf, vjust = 1.3, hjust = -0.05,
    label = sprintf("Observed median = %.2f\n95%% CI %.2f to %.2f\nPermutation P = %.4g",
                    primary_pair$estimate, primary_pair$ci_lower,
                    primary_pair$ci_upper, primary_pair$empirical_p),
    size = 6 / ggplot2::.pt, color = "black"
  ) +
  labs(x = "Median participant cosine under permuted histologic changes", y = "Density") +
  theme_masld(base_size = 6)
ggsave(file.path(tmp, "fig3_paired_holdout_validation.pdf"), p_pair,
       width = 3.5, height = 2.5, device = cairo_pdf)

support <- unique(scores[
  view == "hotspot_program" & score_definition == "weighted_mean_z" &
    dataset %in% discovery_cohorts & !is.na(fibrosis_stage) & !is.na(nas_score),
  .(sample_id, dataset, fibrosis_stage, nas_score)
])
grid <- support[, .N, by = .(dataset, fibrosis_stage, nas_score)]
p_grid <- ggplot(grid, aes(factor(fibrosis_stage), factor(nas_score), fill = N)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = N), size = 6 / ggplot2::.pt) +
  scale_fill_gradient(low = "#EAF2FB", high = "#1565C0", name = "Samples") +
  facet_wrap(~ dataset) +
  labs(x = "Fibrosis stage", y = "NAS") +
  theme_masld(base_size = 6) +
  theme(panel.grid = element_blank())
ggsave(file.path(tmp, "figS_program_map_support_grid.pdf"), p_grid,
       width = 6.2, height = 5.0, device = cairo_pdf)

folds <- loco[split == "leave_one_cohort_out"]
p_loco <- ggplot(folds, aes(cohort, estimate, color = axis)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.35), size = 2) +
  scale_color_manual(values = c(fibrosis = "#6A51A3", nas = "#E69F00"),
                     labels = c(fibrosis = "Fibrosis", nas = "NAS"), name = NULL) +
  labs(x = "Held-out cohort", y = "Training vs held-out effect-vector Spearman rho") +
  theme_masld(base_size = 6) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")
ggsave(file.path(tmp, "figS_program_map_loco.pdf"), p_loco,
       width = 4.1, height = 2.6, device = cairo_pdf)

comparator <- meta[
  model_kind == "primary" &
    ((view == "hallmark" & score_definition == "unweighted_mean_z") |
     (view == "published_signature" & score_definition == "directional_mean_z") |
     (view == "composition" & score_definition == "clr_1e-6")) &
    estimable == TRUE
]
comp_wide <- dcast(comparator, view + feature_id ~ axis, value.var = "beta_meta")
p_comp <- ggplot(comp_wide, aes(fibrosis, nas)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.25) +
  geom_point(size = 1.5, color = "#1565C0") +
  facet_wrap(~ view, scales = "free") +
  labs(x = "Conditional fibrosis association", y = "Conditional NAS association") +
  theme_masld(base_size = 6)
ggsave(file.path(tmp, "figS_fixed_comparator_maps.pdf"), p_comp,
       width = 4.6, height = 2.5, device = cairo_pdf)

robust_ids <- effect[robust_display == TRUE, feature_id]
forest <- cohort[
  view == "hotspot_program" & score_definition == "weighted_mean_z" &
    model_kind == "primary" & feature_id %in% robust_ids
]
forest <- merge(forest, effect[, .(feature_id, module_name)], by = "feature_id")
p_forest <- ggplot(forest, aes(beta, cohort, color = axis)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.12,
                 position = position_dodge(width = 0.35), linewidth = 0.35) +
  geom_point(position = position_dodge(width = 0.35), size = 1.5) +
  facet_wrap(~ module_name, scales = "free_x") +
  scale_color_manual(values = c(fibrosis = "#6A51A3", nas = "#E69F00"), name = NULL) +
  labs(x = "Cohort coefficient with HC3 95% CI", y = NULL) +
  theme_masld(base_size = 6) +
  theme(legend.position = "top")
ggsave(file.path(tmp, "figS_robust_display_program_forests.pdf"), p_forest,
       width = 5.4, height = 2.8, device = cairo_pdf)

primary_meta <- meta[
  view == "hotspot_program" & score_definition == "weighted_mean_z" &
    model_kind == "primary" & estimable == TRUE,
  .(feature_id, axis, beta_meta)
]
all_cohort <- merge(
  cohort[
    view == "hotspot_program" & score_definition == "weighted_mean_z" &
      model_kind == "primary" & estimable == TRUE,
    .(feature_id, cohort, axis, beta)
  ],
  primary_meta, by = c("feature_id", "axis")
)
p_all_cohort <- ggplot(all_cohort, aes(beta_meta, beta)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.25) +
  geom_point(size = 0.7, alpha = 0.55, color = "#1565C0") +
  facet_grid(axis ~ cohort, scales = "free") +
  labs(x = "Four-cohort meta coefficient", y = "Cohort coefficient") +
  theme_masld(base_size = 6)
ggsave(file.path(tmp, "figS_all_program_cohort_effects.pdf"), p_all_cohort,
       width = 6.8, height = 3.4, device = cairo_pdf)

score_sensitivity <- meta[
  view == "hotspot_program" & model_kind == "primary" & estimable == TRUE &
    score_definition %in% c("weighted_mean_z", "unweighted_mean_z",
                            "weighted_within_sample_rank"),
  .(feature_id, axis, score_definition, beta_meta)
]
score_sensitivity <- dcast(
  score_sensitivity, feature_id + axis ~ score_definition, value.var = "beta_meta"
)
score_sensitivity <- melt(
  score_sensitivity,
  id.vars = c("feature_id", "axis", "weighted_mean_z"),
  measure.vars = c("unweighted_mean_z", "weighted_within_sample_rank"),
  variable.name = "alternative", value.name = "alternative_beta"
)
p_sensitivity <- ggplot(score_sensitivity, aes(weighted_mean_z, alternative_beta)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55", linewidth = 0.3) +
  geom_point(size = 0.8, alpha = 0.6, color = "#6A51A3") +
  facet_grid(axis ~ alternative, scales = "free") +
  labs(x = "Primary weighted mean-z coefficient", y = "Alternative-score coefficient") +
  theme_masld(base_size = 6)
ggsave(file.path(tmp, "figS_program_score_sensitivities.pdf"), p_sensitivity,
       width = 5.5, height = 3.8, device = cairo_pdf)

interactions <- meta[
  view == "hotspot_program" & score_definition == "weighted_mean_z" &
    model_kind == "interaction" & estimable == TRUE
]
interactions[, feature_order := frank(beta_meta, ties.method = "first")]
p_interaction <- ggplot(interactions, aes(feature_order, beta_meta, color = q_value < 0.05)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_point(size = 1) +
  scale_color_manual(values = c(`TRUE` = "#C9265E", `FALSE` = "#757575"),
                     labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "Unsupported"),
                     name = NULL) +
  labs(x = "Programs ordered by interaction coefficient",
       y = "Fibrosis × NAS interaction coefficient") +
  theme_masld(base_size = 6) +
  theme(legend.position = "top")
ggsave(file.path(tmp, "figS_program_interactions.pdf"), p_interaction,
       width = 4.2, height = 2.5, device = cairo_pdf)

p_diagnostics <- ggplot(diagnostics, aes(split, estimate, color = metric)) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = NULL, y = "Estimate") +
  theme_masld(base_size = 6) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
ggsave(file.path(tmp, "figS_paired_diagnostics.pdf"), p_diagnostics,
       width = 4.8, height = 2.5, device = cairo_pdf)

p_transport <- ggplot(transport, aes(cohort, estimate, color = view)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_point(position = position_dodge(width = 0.45), size = 1.5) +
  facet_wrap(~ axis, scales = "free_x") +
  labs(x = "Transport cohort", y = "Discovery vs transport effect-vector Spearman rho") +
  theme_masld(base_size = 6) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")
ggsave(file.path(tmp, "figS_native_transport.pdf"), p_transport,
       width = 6.2, height = 3.2, device = cairo_pdf)

global <- loco[split == "leave_one_cohort_out_global"]
two_axis_loco <- nrow(global) == 2L && all(global$p_holm < 0.05) && all(global$positive_folds >= 3L)
paired_pass <- nrow(primary_pair) == 1L && primary_pair$empirical_p < 0.05 && primary_pair$ci_lower > 0
named <- effect[robust_display == TRUE]
named_claims_robust <- all(
  (!is.finite(named$q_value_fibrosis) | named$q_value_fibrosis >= 0.05 | named$robust_fibrosis) &
  (!is.finite(named$q_value_nas) | named$q_value_nas >= 0.05 | named$robust_nas)
)
decision <- data.table(
  criterion = c("both_loco_axes", "paired_holdout", "predeclared_named_programs", "two_axis_main_figure"),
  passed = c(two_axis_loco, paired_pass, named_claims_robust,
             two_axis_loco && paired_pass && named_claims_robust),
  consequence = c(
    "Two conditional axes replicate across cohorts",
    "Within-participant molecular change aligns with histologic change",
    "Named claims meet score and cohort direction rules",
    if (two_axis_loco && paired_pass && named_claims_robust)
      "Eligible for scientific adjudication; not automatically promoted"
    else "Keep failed component supplementary and narrow claims"
  )
)
write_tsv(decision, file.path(tmp, "promotion_decision.tsv"))

caption <- c(
  "Fibrosis-NAS program map: conditional random-effects coefficients from four independent cohorts.",
  "Programs were frozen before this analysis; all 117 remain represented in the source table.",
  "Fibrosis and NAS coefficients retain different clinical units and their magnitudes are not compared.",
  "Cross-sectional associations are not longitudinal progression.",
  "Paired validation uses 54 biological participants from GSE193066 and no discovery refit."
)
writeLines(caption, file.path(tmp, "figure_captions.txt"))

artifacts <- list.files(tmp, recursive = TRUE, full.names = TRUE)
artifacts <- artifacts[file.info(artifacts)$isdir == FALSE]
write_tsv(data.table(
  artifact = substring(artifacts, nchar(tmp) + 2L),
  sha256 = vapply(artifacts, sha256_file, character(1)),
  size_bytes = file.info(artifacts)$size
), file.path(tmp, "artifact_manifest.tsv"))

publish_dir(tmp, figure_root)
Sys.chmod(list.files(figure_root, recursive = TRUE, full.names = TRUE), mode = "0440")
message("Published candidate figures: ", figure_root)
