#!/usr/bin/env Rscript
# KEY MESSAGE: All five Resource cohorts show the continuum-associated molecular patterns, while formal inference remains restricted to the two competitor-independent cohorts.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

options(digits = 17, scipen = 999)
grDevices::pdf.options(useDingbats = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 64_render_five_cohort_context.R ANALYSIS_CANDIDATE OUTPUT_CANDIDATE",
       call. = FALSE)
}
analysis_candidate <- normalizePath(args[[1L]], mustWork = TRUE)
output_candidate <- args[[2L]]
stopifnot(!file.exists(output_candidate))

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lib_molecular_layers.R"))
source(file.path(script_dir, "10_bulk_lib.R"))

contract <- ml_read_contract()
source_root <- ml_source_root(contract)
all_cohorts <- c(contract$source_overlap_cohorts, contract$evaluation_cohorts)
evaluation_cohorts <- contract$evaluation_cohorts
stage_cohorts <- setdiff(all_cohorts, "GSE126848")
ml_assert(length(all_cohorts) == 5L, "Expected five Resource cohorts")

panel_dir <- file.path(output_candidate, "panels")
source_dir <- file.path(output_candidate, "source_tables")
provenance_dir <- file.path(output_candidate, "provenance")
for (path in c(panel_dir, source_dir, provenance_dir)) ml_ensure_dir(path)

font_family <- "Helvetica"
text_size <- 6
cohort_colors <- c(
  GSE126848 = "#BDBDBD", GSE130970 = "#7A7A7A", GSE135251 = "#4D4D4D",
  GSE162694 = "#C9265E", GSE213621 = "#1565C0"
)
cohort_roles <- c(
  GSE126848 = "Source-overlap, descriptive; fibrosis unavailable",
  GSE130970 = "Source-overlap, descriptive",
  GSE135251 = "Source-overlap, descriptive",
  GSE162694 = "Independent evaluation",
  GSE213621 = "Independent evaluation"
)
cohort_labels <- paste0(names(cohort_roles), " (", unname(cohort_roles), ")")
names(cohort_labels) <- names(cohort_roles)
cohort_display <- c(
  GSE126848 = "GSE126848†", GSE130970 = "GSE130970†",
  GSE135251 = "GSE135251†", GSE162694 = "GSE162694",
  GSE213621 = "GSE213621"
)
cohort_display_plain <- setNames(all_cohorts, all_cohorts)
cohort_linetypes <- c(
  GSE126848 = "22", GSE130970 = "22", GSE135251 = "22",
  GSE162694 = "solid", GSE213621 = "solid"
)
stage_colors <- c(
  `0` = "#9E9E9E", `1` = "#F4A674", `2` = "#ED8795",
  `3` = "#C9265E", `4` = "#741816"
)
program_colors <- c("Ductular injury" = "#1565C0", "Stromal ECM" = "#C9265E")

theme_context <- function() {
  theme_classic(base_size = text_size, base_family = font_family) %+replace%
    theme(
      text = element_text(size = text_size, face = "plain"),
      axis.text = element_text(size = text_size, color = "black", face = "plain"),
      axis.title = element_text(size = text_size, face = "plain"),
      legend.text = element_text(size = text_size, face = "plain"),
      legend.title = element_text(size = text_size, face = "plain"),
      strip.text = element_text(size = text_size, face = "plain"),
      strip.background = element_blank(), plot.title = element_blank(),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.key.size = grid::unit(2.7, "mm"),
      panel.spacing = grid::unit(1.1, "mm"),
      plot.caption = element_text(size = text_size, face = "plain", hjust = 0),
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

axes <- ml_load_axes(contract)[
  axis_id == "fixed_projection" & dataset %in% all_cohorts,
  .(sample_id, dataset, axis_percentile, axis_raw)
]
ml_assert(setequal(unique(axes$dataset), all_cohorts),
          "Fixed projection is not available in all five cohorts")

window_summary <- function(scores, feature_columns, score_column) {
  feature_columns <- as.character(feature_columns)
  required <- c("sample_id", "dataset", feature_columns, score_column)
  ml_assert(all(required %in% names(scores)), "Window-score schema is incomplete")
  value <- merge(
    scores[dataset %in% all_cohorts, ..required], axes,
    by = c("sample_id", "dataset"), all = FALSE
  )
  value[, score_z := ml_standardize(get(score_column)),
        by = c("dataset", feature_columns)]
  centers <- as.numeric(contract$window_centers)
  width <- as.numeric(contract$window_width)
  result <- rbindlist(lapply(seq_along(centers), function(index) {
    center <- centers[[index]]
    lower <- center - width / 2
    upper <- center + width / 2
    keep <- if (index == length(centers)) {
      value$axis_percentile >= lower & value$axis_percentile <= upper
    } else {
      value$axis_percentile >= lower & value$axis_percentile < upper
    }
    value[keep, .(
      n_participants = uniqueN(sample_id), mean_score = mean(score_z, na.rm = TRUE),
      se_score = sd(score_z, na.rm = TRUE) / sqrt(sum(is.finite(score_z)))
    ), by = c("dataset", feature_columns)][, `:=`(
      window_id = index, window_center = center,
      window_lower = lower, window_upper = upper,
      visualization_only = TRUE
    )]
  }), use.names = TRUE)
  expected <- uniqueN(value[, c("dataset", feature_columns), with = FALSE]) * length(centers)
  ml_assert(nrow(result) == expected, "A five-cohort fixed window is missing")
  result
}

stage_summary <- function(scores, feature_columns, score_column) {
  feature_columns <- as.character(feature_columns)
  required <- c("sample_id", "dataset", "fibrosis_stage", feature_columns, score_column)
  value <- scores[
    dataset %in% stage_cohorts & !is.na(fibrosis_stage), ..required
  ]
  value[, score_z := ml_standardize(get(score_column)),
        by = c("dataset", feature_columns)]
  value[, .(
    n_participants = uniqueN(sample_id), mean_score = mean(score_z, na.rm = TRUE),
    se_score = sd(score_z, na.rm = TRUE) / sqrt(sum(is.finite(score_z)))
  ), by = c("dataset", "fibrosis_stage", feature_columns)]
}

cohort_scales <- function(labels = cohort_display) {
  list(
    scale_color_manual(values = cohort_colors, breaks = all_cohorts,
                       labels = labels[all_cohorts], name = NULL),
    scale_linetype_manual(values = cohort_linetypes, breaks = all_cohorts,
                          labels = labels[all_cohorts], name = NULL)
  )
}

make_window_plot <- function(data, facet_column, facet_levels, y_label,
                             ncol = 2L, show_error = FALSE,
                             legend_labels = cohort_display,
                             show_participant_size = TRUE,
                             caption_text = paste0(
                               "† Source-overlap cohort, descriptive only (dashed gray). ",
                               "Solid color: independent evaluation; GSE126848 fibrosis unavailable."
                             ),
                             x_label = "Fixed-projection window center") {
  data[, facet_value := factor(get(facet_column), levels = facet_levels)]
  plot <- ggplot(
    data,
    aes(window_center, mean_score, color = dataset, linetype = dataset, group = dataset)
  ) +
    geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
    geom_line(linewidth = 0.52) +
    facet_wrap(~facet_value, ncol = ncol, scales = "free_y", drop = FALSE) +
    cohort_scales(legend_labels) +
    scale_x_continuous(
      breaks = c(0.1, 0.5, 0.9), labels = c("10%", "50%", "90%"),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = x_label, y = y_label, caption = caption_text
    ) +
    theme_context() +
    theme(legend.position = "bottom", legend.box = "vertical")
  if (show_participant_size) {
    plot <- plot +
      geom_point(aes(size = n_participants), stroke = 0, alpha = 0.85) +
      scale_size_continuous(range = c(0.5, 1.35), name = "Participants")
  } else {
    plot <- plot + geom_point(size = 0.85, stroke = 0, alpha = 0.85)
  }
  if (show_error) {
    plot <- plot + geom_errorbar(
      aes(ymin = mean_score - 1.96 * se_score, ymax = mean_score + 1.96 * se_score),
      width = 0, linewidth = 0.22, alpha = 0.65
    )
  }
  plot
}

make_stage_plot <- function(data, facet_column, facet_levels, y_label, ncol = 2L) {
  data[, facet_value := factor(get(facet_column), levels = facet_levels)]
  data[, stage_factor := factor(as.integer(fibrosis_stage), levels = 0:4)]
  ggplot(data, aes(stage_factor, mean_score, color = dataset, group = dataset,
                   linetype = dataset)) +
    geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
    geom_errorbar(aes(ymin = mean_score - 1.96 * se_score,
                      ymax = mean_score + 1.96 * se_score),
                  width = 0, linewidth = 0.22) +
    geom_line(linewidth = 0.5) + geom_point(aes(size = n_participants), stroke = 0) +
    facet_wrap(~facet_value, ncol = ncol, scales = "free_y", drop = FALSE) +
    cohort_scales() + scale_size_continuous(range = c(0.55, 1.4), name = "Participants") +
    labs(
      x = "Recorded fibrosis stage", y = y_label,
      caption = paste0(
        "† Source-overlap cohort, descriptive only (dashed gray). ",
        "GSE126848 omitted because fibrosis stage is unavailable."
      )
    ) +
    theme_context() + theme(legend.position = "bottom", legend.box = "vertical")
}

message("Focal Hotspot programs")
hotspot_root <- file.path(analysis_candidate, "programs", "hotspot")
hotspot <- fread(file.path(hotspot_root, "hotspot_participant_scores.tsv.gz"))
focal_lookup <- data.table(
  focal_program = names(contract$focal_programs),
  program_uid = unname(unlist(contract$focal_programs))
)
hotspot <- merge(hotspot, focal_lookup, by = "program_uid", all.x = TRUE)
focal_scores <- hotspot[!is.na(focal_program)]
focal_windows <- window_summary(focal_scores, c("program_uid", "focal_program"), "outcome_z")
focal_stage <- stage_summary(focal_scores, c("program_uid", "focal_program"), "outcome_z")
ml_write_tsv_once(focal_windows, file.path(source_dir, "focal_hotspot_five_cohort_windows.tsv"))
ml_write_tsv_once(focal_stage, file.path(source_dir, "focal_hotspot_four_cohort_stage.tsv"))

trajectory_main <- make_window_plot(
  focal_windows, "focal_program", names(contract$focal_programs),
  "Standardized program score", ncol = 2L, show_error = FALSE,
  legend_labels = cohort_display_plain, show_participant_size = FALSE,
  caption_text = NULL, x_label = "Continuum percentile"
) + theme(legend.position = "bottom", legend.box = "horizontal",
          legend.margin = margin(0, 0, 0, 0)) +
  guides(linetype = "none", color = guide_legend(nrow = 1, byrow = TRUE))
figure4f <- trajectory_main
fig4f_path <- save_panel(
  figure4f, "fig4f_continuum_program_trajectories.pdf", 5.5, 2.55
)
focal_window_path <- save_panel(
  make_window_plot(focal_windows, "focal_program", names(contract$focal_programs),
                   "Cohort-standardized program score", 2L, TRUE),
  "s4_focal_hotspot_by_continuum_windows.pdf", 5.5, 3.1
)
focal_stage_path <- save_panel(
  make_stage_plot(focal_stage, "focal_program", names(contract$focal_programs),
                  "Cohort-standardized program score", 2L),
  "s4_focal_hotspot_by_fibrosis_stage.pdf", 5.5, 3.0
)

message("Complete Hotspot family")
membership <- fread(file.path(hotspot_root, "hotspot_continuum_membership.tsv"))
registry <- fread(ml_resolve(contract$program_registry))[, .(
  program_uid, registry_order = .I, cell_type, module, module_name
)]
hotspot_all <- hotspot[, focal_program := NULL]
hotspot_windows <- window_summary(
  hotspot_all, c("program_uid", "cell_type", "module", "module_name"), "outcome_z"
)
programs <- merge(registry, membership[, .(program_uid, continuum_associated)],
                  by = "program_uid", all.x = TRUE, sort = FALSE)
cell_labels <- c(cholangiocytes = "Cholangiocytes", fibroblasts = "Fibroblasts",
                 hepatocytes = "Hepatocytes", macrophages = "Macrophages",
                 tcells = "T cells")
programs[, program_label := sprintf(
  "%s %02d · %s", cell_labels[cell_type], as.integer(module), module_name
)]
hotspot_windows <- merge(
  hotspot_windows,
  programs[, .(program_uid, registry_order, program_label, continuum_associated)],
  by = "program_uid", all.x = TRUE
)
ml_write_tsv_once(hotspot_windows,
                  file.path(source_dir, "all117_hotspot_five_cohort_windows.tsv.gz"))

render_hotspot_heatmap <- function(supported_only, filename, height) {
  d <- hotspot_windows[!supported_only | continuum_associated == TRUE]
  order_table <- unique(d[, .(program_uid, registry_order, program_label)])[order(registry_order)]
  d[, program_factor := factor(program_label, levels = rev(order_table$program_label))]
  d[, cohort_label := factor(cohort_display[dataset], levels = cohort_display[all_cohorts])]
  finite <- d[is.finite(mean_score), mean_score]
  clip <- as.numeric(quantile(abs(finite), 0.98, names = FALSE, type = 8))
  gg <- ggplot(d, aes(window_center, program_factor, fill = mean_score)) +
    geom_tile(color = "white", linewidth = 0.04) +
    facet_grid(. ~ cohort_label) +
    scale_fill_gradient2(
      low = "#1565C0", mid = "white", high = "#C9265E", midpoint = 0,
      limits = c(-clip, clip), oob = scales::squish,
      name = "Mean program score\n(cohort z)"
    ) +
    scale_x_continuous(breaks = c(0.1, 0.5, 0.9), labels = c("10%", "50%", "90%")) +
    labs(
      x = "Fixed-projection window center (overlapping 20-percentile windows)",
      y = "Frozen Hotspot program",
      caption = paste0(
        "† Source-overlap cohort, descriptive windows only. ",
        "Continuum gates and meta-analysis use GSE162694 and GSE213621 only."
      )
    ) +
    theme_context() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.spacing.x = grid::unit(1.2, "mm"), legend.position = "bottom")
  save_panel(gg, filename, 7.1, height)
}
hotspot_all_path <- render_hotspot_heatmap(
  FALSE, "s4_all117_hotspot_continuum_windows_labeled.pdf", 13.0
)
hotspot_supported_path <- render_hotspot_heatmap(
  TRUE, "s4_continuum_associated_hotspot_windows_labeled.pdf", 10.5
)

message("NMF axes")
nmf <- fread(file.path(analysis_candidate, "programs", "nmf", "nmf_participant_scores.tsv.gz"))[
  score_variant == "signature_excluded_fixed_w" & inference_eligible == TRUE
]
nmf[, axis_label := paste0("k", k, " ", program_code, ": ", program_label)]
nmf_order <- unique(nmf[order(k, program_code), axis_label])
nmf_windows <- window_summary(
  nmf, c("outcome_id", "k", "program_code", "program_label", "axis_label"), "outcome_z"
)
nmf_stage <- stage_summary(
  nmf, c("outcome_id", "k", "program_code", "program_label", "axis_label"), "outcome_z"
)
ml_write_tsv_once(nmf_windows, file.path(source_dir, "all10_nmf_five_cohort_windows.tsv"))
ml_write_tsv_once(nmf_stage, file.path(source_dir, "all10_nmf_four_cohort_stage.tsv"))
nmf_window_path <- save_panel(
  make_window_plot(nmf_windows, "axis_label", nmf_order, "NMF loading (cohort z)", 2L),
  "s4_all10_nmf_by_continuum_windows.pdf", 5.5, 8.8
)
nmf_stage_path <- save_panel(
  make_stage_plot(nmf_stage, "axis_label", nmf_order, "NMF loading (cohort z)", 2L),
  "s4_all10_nmf_by_fibrosis_stage.pdf", 5.5, 8.8
)

message("Hallmark pathways")
hallmark <- fread(file.path(
  analysis_candidate, "pathway_tf", "pathways", "hallmark", "donor_scores.tsv.gz"
))[set_id %in% contract$display_hallmarks]
hallmark_labels <- setNames(
  gsub("_", " ", sub("^HALLMARK_", "", contract$display_hallmarks)),
  contract$display_hallmarks
)
hallmark[, set_label := hallmark_labels[set_id]]
hallmark_windows <- window_summary(hallmark, c("set_id", "set_label"), "pathway_score")
manifest <- fread(ml_resolve(contract$sample_manifest), na.strings = c("", "NA"))
hallmark_stage_input <- merge(
  hallmark, manifest[, .(sample_id, dataset, fibrosis_stage)],
  by = c("sample_id", "dataset"), all.x = TRUE
)
hallmark_stage <- stage_summary(hallmark_stage_input, c("set_id", "set_label"), "pathway_score")
ml_write_tsv_once(hallmark_windows, file.path(source_dir, "six_hallmark_five_cohort_windows.tsv"))
ml_write_tsv_once(hallmark_stage, file.path(source_dir, "six_hallmark_four_cohort_stage.tsv"))
hallmark_window_path <- save_panel(
  make_window_plot(hallmark_windows, "set_label", unname(hallmark_labels),
                   "Pathway score (cohort z)", 3L),
  "s3_six_hallmarks_by_continuum_windows.pdf", 5.5, 5.8
)
hallmark_stage_path <- save_panel(
  make_stage_plot(hallmark_stage, "set_label", unname(hallmark_labels),
                  "Pathway score (cohort z)", 3L),
  "s3_six_hallmarks_by_fibrosis_stage.pdf", 5.5, 5.8
)

message("Regulon scores")
tf <- fread(file.path(
  analysis_candidate, "pathway_tf", "tf", "donor_regulon_scores.tsv.gz"
))[tf %in% contract$display_tfs]
tf_windows <- window_summary(tf, "tf", "regulon_score")
ml_write_tsv_once(tf_windows, file.path(source_dir, "four_tf_five_cohort_windows.tsv"))
tf_path <- save_panel(
  make_window_plot(tf_windows, "tf", contract$display_tfs,
                   "Signed regulon score (cohort z)", 4L),
  "s4_four_tf_regulon_continuum_windows.pdf", 5.5, 3.2
)

message("Continuum calibration against stage")
unsupervised <- fread(file.path(source_root, "unsupervised", "participant_scores.tsv"))
projection <- fread(file.path(source_root, "projection", "participant_scores.tsv"))
all_axis <- unique(rbindlist(list(unsupervised, projection), fill = TRUE),
                   by = c("sample_id", "dataset", "axis_id"))
all_axis <- merge(
  all_axis, manifest[, .(sample_id, dataset, fibrosis_stage)],
  by = c("sample_id", "dataset"), all.x = TRUE, suffixes = c("", ".manifest")
)
if ("fibrosis_stage.manifest" %in% names(all_axis)) {
  all_axis[is.na(fibrosis_stage), fibrosis_stage := fibrosis_stage.manifest]
}
axis_labels <- c(
  signature_pc1 = "145-gene cohort PC1", fixed_projection = "Fixed discovery projection",
  full_transcriptome_pc1 = "Full-transcriptome PC1", rf_fibrosis = "Discovery-trained RF",
  consensus_rank = "Consensus rank"
)
calibration <- all_axis[
  dataset %in% stage_cohorts & axis_id %in% names(axis_labels) &
    is.finite(axis_percentile) & !is.na(fibrosis_stage)
]
calibration[, `:=`(
  axis_label = factor(axis_labels[axis_id], levels = axis_labels),
  fibrosis_factor = factor(as.integer(fibrosis_stage), levels = 0:4),
  cohort_label = factor(cohort_display[dataset], levels = cohort_display[stage_cohorts])
)]
medians <- calibration[, .(median_percentile = median(axis_percentile), n_participants = .N),
                       by = .(dataset, cohort_label, axis_id, axis_label, fibrosis_factor)]
ml_write_tsv_once(calibration, file.path(source_dir, "all_scores_four_cohort_by_stage.tsv.gz"))
ml_write_tsv_once(medians, file.path(source_dir, "all_scores_four_cohort_by_stage_medians.tsv"))

calibration_primary <- calibration[axis_id %in% contract$co_primary_axes]
median_primary <- medians[axis_id %in% contract$co_primary_axes]
calibration_plot <- ggplot(
  calibration_primary,
  aes(fibrosis_factor, axis_percentile, color = fibrosis_factor)
) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = contract$seed),
             size = 0.45, alpha = 0.33, stroke = 0) +
  geom_point(data = median_primary, aes(y = median_percentile),
             shape = 95, size = 3.4, color = "black") +
  facet_grid(axis_label ~ cohort_label) +
  scale_color_manual(values = stage_colors, guide = "none", drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", "50", "100")) +
  labs(
    x = "Recorded fibrosis stage", y = "Within-cohort continuum percentile",
    caption = paste0(
      "† Source-overlap cohort, descriptive only. GSE126848 omitted because ",
      "fibrosis stage is unavailable."
    )
  ) +
  theme_context()
calibration_path <- save_panel(
  calibration_plot, "s3_continuum_calibration_by_stage.pdf", 7.1, 3.9
)

all_stage_plot <- ggplot(calibration, aes(fibrosis_factor, axis_percentile,
                                          color = fibrosis_factor)) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = contract$seed),
             size = 0.32, alpha = 0.25, stroke = 0) +
  geom_point(data = medians, aes(y = median_percentile),
             shape = 95, size = 2.9, color = "black") +
  facet_grid(axis_label ~ cohort_label) +
  scale_color_manual(values = stage_colors, guide = "none", drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", "50", "100")) +
  labs(
    x = "Recorded fibrosis stage", y = "Within-cohort score percentile",
    caption = paste0(
      "† Source-overlap cohort, descriptive only. GSE126848 omitted because ",
      "fibrosis stage is unavailable."
    )
  ) +
  theme_context()
all_stage_path <- save_panel(all_stage_plot, "s3_all_continuum_scores_by_stage.pdf", 7.1, 7.7)

message("Five fixed genes")
dge <- readRDS(ml_resolve(contract$dge_rds))
log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
rownames(log_cpm) <- ml_base_gene_id(rownames(log_cpm))
gene_annotation <- fread(ml_resolve(contract$gene_annotation))
gene_annotation[, gene_id_base := ml_base_gene_id(gene_id)]
gene_map <- unique(gene_annotation[gene_name %in% contract$display_genes,
                                   .(gene_name, gene_id_base)])
ml_assert(nrow(gene_map) == length(contract$display_genes) && !anyDuplicated(gene_map$gene_name),
          "Five-gene display roster does not map one-to-one")
projection_loadings <- fread(file.path(source_root, "projection", "fixed_projection_loadings.tsv"))
gene_rows <- rbindlist(lapply(all_cohorts, function(cohort) {
  normalized <- readRDS(file.path(
    source_root, "unsupervised", "normalized_expression_by_cohort", paste0(cohort, ".rds")
  ))
  full_axis <- axes[dataset == cohort, .(sample_id, axis_raw)]
  rbindlist(lapply(contract$display_genes, function(symbol) {
    gene_id <- gene_map[gene_name == symbol, gene_id_base]
    axis <- if (symbol %in% c("CYP2C19", "IGFBP7")) {
      bulk_loo_axis("fixed_projection", gene_id, normalized, data.table(), projection_loadings)
    } else {
      setNames(full_axis$axis_raw, full_axis$sample_id)
    }
    data.table(
      sample_id = intersect(colnames(log_cpm), names(axis)), dataset = cohort,
      gene_name = symbol, gene_id_base = gene_id,
      expression = as.numeric(log_cpm[gene_id, intersect(colnames(log_cpm), names(axis))]),
      loo_axis_percentile = ml_percentile(axis[intersect(colnames(log_cpm), names(axis))]),
      axis_role = if (symbol %in% c("CYP2C19", "IGFBP7"))
        "target_specific_leave_one_gene_out_fixed_projection" else "frozen_fixed_projection"
    )
  }))
}))
gene_rows[, expression_z := ml_standardize(expression), by = .(dataset, gene_name)]
gene_axes <- gene_rows[, .(sample_id, dataset, gene_name, axis_percentile = loo_axis_percentile)]
gene_windows <- rbindlist(lapply(contract$display_genes, function(symbol) {
  d <- gene_rows[gene_name == symbol]
  d <- merge(d[, .(sample_id, dataset, gene_name, expression_z)],
             gene_axes[gene_name == symbol, .(sample_id, dataset, axis_percentile)],
             by = c("sample_id", "dataset"))
  centers <- as.numeric(contract$window_centers)
  width <- as.numeric(contract$window_width)
  rbindlist(lapply(seq_along(centers), function(index) {
    lower <- centers[[index]] - width / 2
    upper <- centers[[index]] + width / 2
    keep <- if (index == length(centers)) d$axis_percentile >= lower & d$axis_percentile <= upper else
      d$axis_percentile >= lower & d$axis_percentile < upper
    d[keep, .(n_participants = uniqueN(sample_id), mean_score = mean(expression_z),
              se_score = sd(expression_z) / sqrt(.N)), by = .(dataset, gene_name)][,
      `:=`(window_id = index, window_center = centers[[index]],
           window_lower = lower, window_upper = upper, visualization_only = TRUE)]
  }))
}))
ml_write_tsv_once(gene_windows, file.path(source_dir, "five_gene_five_cohort_windows.tsv"))
gene_path <- save_panel(
  make_window_plot(gene_windows, "gene_name", contract$display_genes,
                   "Expression (cohort z)", 3L),
  "s3_fixed_gene_continuum_trajectories.pdf", 5.5, 4.5
)

figure_manifest <- data.table(
  callout = c("4F", "S3I", "S3K", "S3P", "S3T", "S3U", "S4F", "S4G",
              "S4H", "S4I", "S4M", "S4N", "S4O"),
  filename = c(
    basename(fig4f_path), basename(calibration_path), basename(all_stage_path),
    basename(gene_path), basename(hallmark_window_path), basename(hallmark_stage_path),
    basename(focal_window_path), basename(focal_stage_path), basename(hotspot_all_path),
    basename(hotspot_supported_path), basename(nmf_window_path), basename(nmf_stage_path),
    basename(tf_path)
  ),
  path = c(
    fig4f_path, calibration_path, all_stage_path, gene_path, hallmark_window_path,
    hallmark_stage_path, focal_window_path, focal_stage_path, hotspot_all_path,
    hotspot_supported_path, nmf_window_path, nmf_stage_path, tf_path
  ),
  descriptive_cohorts = c(5L, 4L, 4L, 5L, 5L, 4L, 5L, 4L, 5L, 5L, 5L, 4L, 5L),
  primary_inference_cohorts = c(2L, rep(NA_integer_, 12L)),
  source_overlap_descriptive_only = TRUE,
  embedded_explanatory_text = c(FALSE, rep(NA, 12L)),
  forest_panel_displayed = c(FALSE, rep(NA, 12L)),
  participant_size_legend_displayed = c(FALSE, rep(NA, 12L))
)
figure_manifest[, sha256 := vapply(path, ml_sha256, character(1))]
ml_write_tsv_once(figure_manifest, file.path(output_candidate, "figure_manifest.tsv"))

input_paths <- c(
  contract = file.path(script_dir, "00_contract.json"),
  analysis_candidate = analysis_candidate,
  axes = file.path(source_root, "projection", "participant_scores.tsv"),
  manifest = ml_resolve(contract$sample_manifest), dge = ml_resolve(contract$dge_rds)
)
input_manifest <- data.table(
  input = names(input_paths), path = input_paths,
  sha256 = vapply(input_paths, function(path) if (file.exists(path) && !dir.exists(path))
    ml_sha256(path) else NA_character_, character(1))
)
ml_write_tsv_once(input_manifest, file.path(provenance_dir, "input_manifest.tsv"))
ml_write_session_info(file.path(provenance_dir, "sessionInfo.txt"))
writeLines(c(
  "# Five-cohort continuum figure context candidate", "",
  "All descriptive continuum-window panels show the five Resource cohorts. The three cohorts overlapping continuum construction are dashed and descriptive only; GSE162694 and GSE213621 are solid and remain the only cohorts used for stage/sex-adjusted continuum inference and meta-analysis.", "",
  "Stage displays contain GSE130970, GSE135251, GSE162694, and GSE213621. GSE126848 is absent only because the frozen Resource metadata contain no recorded fibrosis stage.", "",
  "No inferential model, membership, NMF basis, Hotspot weight, pathway set, or continuum score was changed."
), file.path(output_candidate, "README.md"))

message("FIVE_COHORT_CONTEXT_COMPLETE: ", output_candidate)
