#!/usr/bin/env Rscript
# KEY MESSAGE: The externally anchored continuum resolves coordinated multicellular remodeling beyond recorded fibrosis stage.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "lib_molecular_layers.R"))
source(file.path(script_dir, "10_bulk_lib.R"))

contract <- ml_read_contract()
candidate <- ml_out_root()
source_root <- ml_source_root(contract)
figure_out <- file.path(candidate, "figures")
panel_out <- file.path(figure_out, "panels")
source_out <- file.path(figure_out, "source_tables")
review_out <- file.path(figure_out, "review_assemblies")
for (path in c(figure_out, panel_out, source_out, review_out)) {
  ml_assert(!dir.exists(path) || identical(path, figure_out),
            paste0("Refusing to overwrite figure namespace: ", path))
  ml_ensure_dir(path)
}

grDevices::pdf.options(useDingbats = FALSE)
font_family <- "Helvetica"
text_size <- 6
program_colors <- c("Ductular injury" = "#1565C0", "Stromal ECM" = "#C9265E")
cohort_linetypes <- c("GSE162694" = "solid", "GSE213621" = "22")
stage_colors <- c(F0 = "#9E9E9E", F1 = "#F4A674", F2 = "#ED8795",
                  F3 = "#C9265E", F4 = "#741816")

theme_continuum <- function() {
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
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.key.size = grid::unit(3, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

save_panel <- function(plot, name, width, height) {
  path <- file.path(panel_out, name)
  ml_assert(!file.exists(path), paste0("Refusing to overwrite figure: ", path))
  ggsave(path, plot, width = width, height = height, device = cairo_pdf,
         units = "in", bg = "white")
  path
}

focal_lookup <- data.table(
  focal_program = names(contract$focal_programs),
  program_uid = unname(unlist(contract$focal_programs))
)
program_result_root <- file.path(candidate, "programs", "hotspot")
program_scores <- fread(file.path(program_result_root, "hotspot_participant_scores.tsv.gz"))
setnames(program_scores, "score_raw", "program_score")
program_scores <- merge(program_scores, focal_lookup, by = "program_uid")
axes <- ml_load_axes(contract)[axis_id == "fixed_projection" &
                                dataset %in% contract$evaluation_cohorts]
plot_data <- merge(
  program_scores[dataset %in% contract$evaluation_cohorts], axes,
  by = c("sample_id", "dataset"), all = FALSE
)
plot_data <- plot_data[
  is.finite(program_score) & is.finite(axis_raw) & !is.na(fibrosis_stage) &
    !is.na(inferred_sex)
]
plot_data[, axis_z := ml_standardize(axis_raw), by = dataset]
plot_data[, stage_factor := factor(fibrosis_stage)]
plot_data[, sex_factor := factor(inferred_sex)]
plot_data[, fibrosis_label := factor(paste0("F", fibrosis_stage), levels = paste0("F", 0:4))]

marginal_predictions <- function(data) {
  # Match the rank guard the analysis scripts use (20_pathway_continuum.R:147,
  # 62_render_stage_continuum_comparison.R:182). An aliased stage x sex cell
  # leaves an NA in coef() while vcov() drops the column, and the quadratic form
  # below then fails with "non-conformable arguments" rather than a clear error.
  design_check <- model.matrix(~ axis_z + stage_factor + sex_factor, data = data)
  ml_assert(qr(design_check)$rank == ncol(design_check),
            "Rank-deficient marginal-prediction model (aliased stage x sex cell)")
  fit <- lm(program_score ~ axis_z + stage_factor + sex_factor, data = data)
  percentile <- seq(0, 1, length.out = 101L)
  axis_grid <- as.numeric(quantile(data$axis_z, probs = percentile, type = 8, na.rm = TRUE))
  terms_no_response <- delete.response(terms(fit))
  covariance <- vcov(fit)
  coefficient <- coef(fit)
  rbindlist(lapply(seq_along(percentile), function(index) {
    new_data <- copy(data)
    new_data[, axis_z := axis_grid[[index]]]
    design <- model.matrix(terms_no_response, new_data)
    average_design <- colMeans(design)[names(coefficient)]
    estimate <- sum(average_design * coefficient)
    se <- sqrt(as.numeric(t(average_design) %*% covariance %*% average_design))
    data.table(
      continuum_percentile = percentile[[index]], estimate = estimate,
      se = se, ci_low = estimate - qnorm(0.975) * se,
      ci_high = estimate + qnorm(0.975) * se,
      model = "linear_stage_sex_adjusted_marginal_prediction",
      n_participants = nrow(data)
    )
  }))
}

prediction <- plot_data[, marginal_predictions(.SD), by = .(dataset, focal_program)]
ml_write_tsv_once(prediction, file.path(source_out, "fig4f_adjusted_trajectories.tsv"))

rug_data <- unique(plot_data[, .(
  sample_id, dataset, focal_program, continuum_percentile = axis_percentile,
  fibrosis_label
)])
ml_write_tsv_once(
  rug_data, file.path(source_out, "fig4f_participant_density_stage.tsv")
)
trajectory_plot <- ggplot(
  prediction,
  aes(continuum_percentile, estimate, color = focal_program, linetype = dataset)
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = focal_program),
              alpha = 0.12, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.65) +
  geom_rug(
    data = rug_data,
    aes(x = continuum_percentile, color = fibrosis_label, linetype = NULL),
    sides = "b", alpha = 0.32, linewidth = 0.3, inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  facet_wrap(~focal_program, nrow = 1, scales = "free_y") +
  scale_color_manual(
    values = c(program_colors, stage_colors),
    breaks = c(names(program_colors), names(stage_colors)),
    name = "Program / fibrosis stage"
  ) +
  scale_fill_manual(values = program_colors, guide = "none") +
  scale_linetype_manual(values = cohort_linetypes) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Fixed-projection percentile", y = "Adjusted program score (z)",
       linetype = NULL) +
  theme_continuum() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.box = "vertical")

meta <- fread(file.path(program_result_root, "hotspot_meta_analysis.tsv"))
setnames(
  meta,
  c("fixed_beta", "fixed_ci_low", "fixed_ci_high", "fixed_q_value"),
  c("beta", "ci_low", "ci_high", "q_value")
)
meta <- merge(meta[
  program_uid %in% focal_lookup$program_uid & axis_id %in% contract$co_primary_axes
], focal_lookup, by = "program_uid")
meta[, score_label := factor(
  fifelse(axis_id == "signature_pc1", "Cohort PC1", "Fixed projection"),
  levels = c("Fixed projection", "Cohort PC1")
)]
meta[, focal_program := factor(focal_program, levels = names(contract$focal_programs))]
ml_write_tsv_once(meta, file.path(source_out, "fig4f_adjusted_meta_estimates.tsv"))

forest_plot <- ggplot(meta, aes(beta, score_label, color = focal_program)) +
  geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.45) +
  geom_point(size = 1.5) +
  facet_wrap(~focal_program, ncol = 1) +
  scale_color_manual(values = program_colors, guide = "none") +
  labs(x = "Adjusted meta-analysis beta", y = NULL) +
  theme_continuum() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())

figure4f <- trajectory_plot + forest_plot + plot_layout(widths = c(3.9, 1.6))
figure4f_path <- save_panel(figure4f, "fig4f_continuum_program_trajectories.pdf", 5.5, 2.55)
compact_path <- save_panel(
  trajectory_plot + theme(legend.position = "none"),
  "fig4f_continuum_program_trajectories_compact.pdf", 2.71, 2.42
)

message("Rendering fixed-window descriptive focal trajectories")
windows <- fread(file.path(program_result_root, "hotspot_fixed_windows.tsv.gz"))
windows <- merge(windows[
  program_uid %in% focal_lookup$program_uid & axis_id == "fixed_projection" &
    dataset %in% contract$evaluation_cohorts
], focal_lookup, by = "program_uid")
ml_write_tsv_once(windows, file.path(source_out, "s4_focal_fixed_windows.tsv"))
window_plot <- ggplot(
  windows,
  aes(window_center, mean_score, color = focal_program, linetype = dataset)
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_errorbar(aes(ymin = mean_score - 1.96 * se_score,
                    ymax = mean_score + 1.96 * se_score),
                width = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.55) +
  geom_point(aes(size = n_participants), shape = 21, fill = "white", stroke = 0.35) +
  facet_wrap(~focal_program, nrow = 1, scales = "free_y") +
  scale_color_manual(values = program_colors) +
  scale_linetype_manual(values = cohort_linetypes) +
  scale_size_continuous(range = c(0.7, 1.8)) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     breaks = c(0.1, 0.5, 0.9)) +
  labs(x = "Fixed-projection window center", y = "Raw program score (z)",
       color = NULL, linetype = NULL, size = "Participants") +
  theme_continuum() +
  theme(legend.position = "bottom")
raw_window_path <- save_panel(window_plot, "s4_focal_fixed_windows.pdf", 5.1, 2.25)

message("Rendering all-117 adjusted effect heatmap")
program_meta <- fread(file.path(program_result_root, "hotspot_meta_analysis.tsv"))
setnames(
  program_meta,
  c("fixed_beta", "fixed_ci_low", "fixed_ci_high", "fixed_q_value"),
  c("beta", "ci_low", "ci_high", "q_value")
)
program_meta <- program_meta[axis_id %in% contract$co_primary_axes]
program_meta[, display_label := paste0(module_name, " [", cell_type, "]")]
program_order <- program_meta[, .(mean_beta = mean(beta, na.rm = TRUE)), by = display_label][
  order(mean_beta), display_label
]
program_meta[, display_label := factor(display_label, levels = program_order)]
program_meta[, axis_label := factor(
  fifelse(axis_id == "signature_pc1", "Cohort PC1", "Fixed projection"),
  levels = c("Cohort PC1", "Fixed projection")
)]
program_meta[, supported := q_value < 0.05 & direction_concordant]
ml_write_tsv_once(program_meta, file.path(source_out, "s4_all117_adjusted_effects.tsv"))
all117_plot <- ggplot(program_meta, aes(axis_label, display_label, fill = beta)) +
  geom_tile(color = "white", linewidth = 0.05) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C9265E",
                       midpoint = 0, name = "Beta") +
  labs(x = NULL, y = NULL) +
  theme_continuum() +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), legend.position = "bottom")
all117_path <- save_panel(all117_plot, "s4_all117_adjusted_effects.pdf", 4.2, 8.1)

paired_path <- file.path(candidate, "paired", "paired_program_results.tsv")
paired_figure <- NA_character_
if (file.exists(paired_path)) {
  paired <- fread(paired_path)
  paired <- merge(paired[feature_id %in% focal_lookup$program_uid], focal_lookup,
                  by.x = "feature_id", by.y = "program_uid")
  paired[, axis_label := fifelse(axis_id == "signature_pc1", "Cohort PC1", "Fixed projection")]
  ml_write_tsv_once(paired, file.path(source_out, "s4_paired_focal_programs.tsv"))
  paired_plot <- ggplot(paired, aes(beta, axis_label, color = focal_program)) +
    geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.4) +
    geom_point(size = 1.5) +
    facet_wrap(~focal_program, nrow = 1) +
    scale_color_manual(values = program_colors, guide = "none") +
    labs(x = "Within-participant delta beta", y = NULL) +
    theme_continuum() +
    theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
  paired_figure <- save_panel(paired_plot, "s4_paired_focal_programs.pdf", 4.1, 1.8)
}

message("Rendering bulk transcript and held-out validation panels")
bulk_root <- file.path(candidate, "bulk")
bulk_meta <- fread(file.path(bulk_root, "continuum_gene_meta.tsv.gz"))
bulk_loo <- fread(file.path(bulk_root, "signature_loo_gene_meta.tsv"))
display_genes <- contract$display_genes
ordinary_gene_forest <- bulk_meta[
  gene_name %in% display_genes & !axis_component,
  .(gene_name, axis_id, beta, ci_low, ci_high, q_value = bh_q_value,
    evidence_role)
]
loo_gene_forest <- bulk_loo[
  gene_name %in% display_genes,
  .(gene_name, axis_id, beta, ci_low, ci_high, q_value = bh_q_value,
    evidence_role)
]
gene_forest_data <- rbindlist(list(ordinary_gene_forest, loo_gene_forest), fill = TRUE)
gene_forest_data <- unique(gene_forest_data, by = c("gene_name", "axis_id"))
gene_forest_data[, axis_label := fifelse(
  axis_id == "signature_pc1", "Cohort PC1", "Fixed projection"
)]
gene_forest_data[, gene_name := factor(gene_name, levels = rev(display_genes))]
ml_write_tsv_once(gene_forest_data, file.path(source_out, "s3_fixed_gene_forest.tsv"))
gene_forest_plot <- ggplot(gene_forest_data, aes(beta, gene_name, color = axis_label)) +
  geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0,
                 linewidth = 0.4, position = position_dodge(width = 0.35)) +
  geom_point(size = 1.4, position = position_dodge(width = 0.35)) +
  scale_color_manual(values = c("Cohort PC1" = "#C9265E", "Fixed projection" = "#1565C0")) +
  labs(x = "Adjusted continuum beta", y = NULL, color = NULL) +
  theme_continuum() +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")
gene_forest_path <- save_panel(gene_forest_plot, "s3_fixed_gene_adjusted_forest.pdf", 3.5, 2.25)

message("Rendering five fixed-gene adjusted trajectories")
dge <- readRDS(ml_resolve(contract$dge_rds))
manifest <- fread(ml_resolve(contract$sample_manifest), na.strings = c("", "NA"))
manifest <- manifest[dataset %in% contract$evaluation_cohorts]
gene_annotation <- fread(ml_resolve(contract$gene_annotation))
gene_annotation[, gene_id_base := ml_base_gene_id(gene_id)]
gene_annotation <- unique(gene_annotation[, .(gene_id_base, gene_name)])
display_mapping <- gene_annotation[gene_name %in% display_genes]
ml_assert(setequal(display_mapping$gene_name, display_genes) &&
            !anyDuplicated(display_mapping$gene_name),
          "The five-gene display roster does not map one-to-one under GENCODE v49")
gene_id_lookup <- setNames(display_mapping$gene_id_base, display_mapping$gene_name)
log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
rownames(log_cpm) <- ml_base_gene_id(rownames(log_cpm))
projection_loadings <- fread(file.path(
  source_root, "projection", "fixed_projection_loadings.tsv"
))
fixed_axes <- ml_load_axes(contract)[
  axis_id == "fixed_projection" & dataset %in% contract$evaluation_cohorts,
  .(sample_id, dataset, full_axis_raw = axis_raw)
]

gene_trajectory_data <- rbindlist(lapply(contract$evaluation_cohorts, function(cohort) {
  cohort_manifest <- manifest[dataset == cohort]
  normalized_expression <- readRDS(file.path(
    source_root, "unsupervised", "normalized_expression_by_cohort",
    paste0(cohort, ".rds")
  ))
  cohort_full_axis <- fixed_axes[dataset == cohort]
  rbindlist(lapply(display_genes, function(gene_symbol) {
    gene_id <- unname(gene_id_lookup[[gene_symbol]])
    axis_component <- gene_symbol %in% c("CYP2C19", "IGFBP7")
    if (axis_component) {
      axis <- bulk_loo_axis(
        "fixed_projection", gene_id, normalized_expression,
        data.table(), projection_loadings
      )
      axis_data <- data.table(sample_id = names(axis), axis_raw = as.numeric(axis))
      axis_role <- "target_specific_leave_one_gene_out_fixed_projection"
    } else {
      axis_data <- cohort_full_axis[, .(sample_id, axis_raw = full_axis_raw)]
      axis_role <- "frozen_fixed_projection"
    }
    value <- data.table(
      sample_id = colnames(log_cpm),
      expression_raw = as.numeric(log_cpm[gene_id, ])
    )
    value <- merge(value, axis_data, by = "sample_id")
    value <- merge(
      value,
      cohort_manifest[, .(sample_id, dataset, fibrosis_stage, inferred_sex)],
      by = "sample_id"
    )
    value[, `:=`(
      gene_name = gene_symbol,
      gene_id_base = gene_id,
      axis_role = axis_role,
      outcome_z = ml_standardize(expression_raw),
      axis_z = ml_standardize(axis_raw),
      stage_factor = factor(fibrosis_stage),
      sex_factor = factor(inferred_sex)
    )]
    value
  }))
}))

gene_prediction <- gene_trajectory_data[, {
  fit <- lm(outcome_z ~ axis_z + stage_factor + sex_factor, data = .SD)
  percentile <- seq(0, 1, length.out = 101L)
  axis_grid <- as.numeric(quantile(axis_z, probs = percentile, type = 8, na.rm = TRUE))
  terms_no_response <- delete.response(terms(fit))
  covariance <- vcov(fit)
  coefficient <- coef(fit)
  rbindlist(lapply(seq_along(percentile), function(index) {
    new_data <- copy(.SD)
    new_data[, axis_z := axis_grid[[index]]]
    design <- model.matrix(terms_no_response, new_data)
    average_design <- colMeans(design)[names(coefficient)]
    estimate <- sum(average_design * coefficient)
    se <- sqrt(as.numeric(t(average_design) %*% covariance %*% average_design))
    data.table(
      continuum_percentile = percentile[[index]], estimate = estimate,
      se = se, ci_low = estimate - qnorm(0.975) * se,
      ci_high = estimate + qnorm(0.975) * se,
      n_participants = .N,
      model = "logCPM_z ~ fixed_projection_z + factor(fibrosis_stage) + inferred_sex"
    )
  }))
}, by = .(dataset, gene_name, gene_id_base, axis_role)]
ml_write_tsv_once(
  gene_prediction, file.path(source_out, "s3_fixed_gene_adjusted_trajectories.tsv")
)
gene_prediction[, gene_name := factor(gene_name, levels = display_genes)]
gene_trajectory_plot <- ggplot(
  gene_prediction,
  aes(continuum_percentile, estimate, color = dataset, linetype = dataset)
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = dataset),
              alpha = 0.12, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.55) +
  facet_wrap(~gene_name, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("GSE162694" = "#C9265E", "GSE213621" = "#1565C0")) +
  scale_fill_manual(values = c("GSE162694" = "#C9265E", "GSE213621" = "#1565C0")) +
  scale_linetype_manual(values = cohort_linetypes) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     breaks = c(0, 0.5, 1)) +
  labs(x = "Fixed-projection percentile", y = "Adjusted expression (z)",
       color = NULL, linetype = NULL) +
  theme_continuum() +
  theme(strip.text = element_text(face = "italic"), legend.position = "bottom")
gene_trajectory_path <- save_panel(
  gene_trajectory_plot, "s3_fixed_gene_adjusted_trajectories.pdf", 5.2, 3.5
)

canonical <- fread(file.path(bulk_root, "canonical_deg_continuum.tsv.gz"))
canonical_plot <- canonical[inference_eligible == TRUE & is.finite(beta) &
                              is.finite(canonical_logFC)]
canonical_plot[, axis_label := fifelse(
  axis_id == "signature_pc1", "Cohort PC1", "Fixed projection"
)]
ml_write_tsv_once(canonical_plot, file.path(source_out, "s3_canonical_deg_concordance.tsv.gz"))
canonical_concordance_plot <- ggplot(
  canonical_plot, aes(canonical_logFC, beta)
) +
  stat_bin_2d(bins = 55) +
  geom_hline(yintercept = 0, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "white", linewidth = 0.2) +
  facet_wrap(~axis_label, nrow = 1) +
  scale_fill_gradient(low = "#E6E6E5", high = "#C9265E", name = "Genes") +
  labs(x = "Canonical disease-control log2 fold change",
       y = "Stage-adjusted continuum beta") +
  theme_continuum() +
  theme(legend.position = "bottom")
canonical_concordance_path <- save_panel(
  canonical_concordance_plot, "s3_all_gene_disease_continuum_concordance.pdf", 4.7, 2.15
)

heldout <- fread(file.path(bulk_root, "heldout_score_validation.tsv"))
heldout[, dataset_label := factor(dataset, levels = rev(unique(dataset)))]
ml_write_tsv_once(heldout, file.path(source_out, "s3_heldout_signature_validation.tsv"))
heldout_plot <- ggplot(heldout, aes(beta, dataset_label)) +
  geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.4) +
  geom_point(size = 1.5, color = "#00695C") +
  labs(x = "Held-out disease-signature effect", y = NULL) +
  theme_continuum() +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
heldout_path <- save_panel(heldout_plot, "s3_heldout_disease_signature_validation.pdf", 3.15, 1.75)

message("Rendering Hallmark and all-collection pathway panels")
pathway_root <- file.path(candidate, "pathway_tf", "pathways")
hallmark_meta <- fread(file.path(pathway_root, "hallmark", "meta_analysis.tsv"))
hallmark_meta[, axis_label := fifelse(
  axis_id == "signature_pc1", "Cohort PC1", "Fixed projection"
)]
hallmark_order <- hallmark_meta[, .(mean_beta = mean(fixed_beta, na.rm = TRUE)), by = set_id][
  order(mean_beta), set_id
]
hallmark_meta[, set_label := factor(sub("^HALLMARK_", "", set_id),
                                    levels = sub("^HALLMARK_", "", hallmark_order))]
ml_write_tsv_once(hallmark_meta, file.path(source_out, "s3_all50_hallmark_meta.tsv"))
hallmark_heatmap <- ggplot(hallmark_meta, aes(axis_label, set_label, fill = fixed_beta)) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C9265E",
                       midpoint = 0, name = "Beta") +
  labs(x = NULL, y = NULL) +
  theme_continuum() +
  theme(axis.text.y = element_text(size = text_size), axis.ticks = element_blank(),
        axis.line = element_blank(), legend.position = "bottom")
hallmark_heatmap_path <- save_panel(hallmark_heatmap, "s3_all50_hallmark_heatmap.pdf", 4.1, 6.8)

hallmark_windows <- fread(file.path(pathway_root, "hallmark", "fixed_display_windows.tsv"))
hallmark_display_labels <- gsub(
  "_", " ", sub("^HALLMARK_", "", contract$display_hallmarks)
)
hallmark_windows[, set_label := factor(
  gsub("_", " ", sub("^HALLMARK_", "", set_id)),
  levels = hallmark_display_labels
)]
ml_write_tsv_once(hallmark_windows, file.path(source_out, "s3_six_hallmark_windows.tsv"))
hallmark_testability <- fread(file.path(pathway_root, "hallmark", "testability.tsv"))[
  set_id %in% contract$display_hallmarks
]
ml_write_tsv_once(
  hallmark_testability, file.path(source_out, "s3_six_hallmark_testability.tsv")
)
missing_hallmark_labels <- setdiff(
  hallmark_display_labels, as.character(unique(hallmark_windows$set_label))
)
hallmark_placeholders <- if (length(missing_hallmark_labels)) {
  data.table(
    set_label = factor(missing_hallmark_labels, levels = hallmark_display_labels),
    center = rep(0.5, length(missing_hallmark_labels)),
    label = rep("Untestable after\nsignature exclusion", length(missing_hallmark_labels))
  )
} else {
  data.table(
    set_label = factor(character(), levels = hallmark_display_labels),
    center = numeric(), label = character()
  )
}
hallmark_trajectory_plot <- ggplot(
  hallmark_windows,
  aes(center, mean_score, color = dataset, linetype = dataset)
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_errorbar(aes(ymin = mean_score - 1.96 * se_score,
                    ymax = mean_score + 1.96 * se_score), width = 0, linewidth = 0.25) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.8) +
  geom_text(
    data = hallmark_placeholders,
    aes(center, 0, label = label), inherit.aes = FALSE,
    color = "#666666", size = 2
  ) +
  facet_wrap(~set_label, ncol = 3, scales = "free_y", drop = FALSE) +
  scale_color_manual(values = c("GSE162694" = "#C9265E", "GSE213621" = "#1565C0")) +
  scale_linetype_manual(values = cohort_linetypes) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     breaks = c(0.1, 0.5, 0.9)) +
  labs(x = "Fixed-projection window center", y = "Pathway score (z)",
       color = NULL, linetype = NULL) +
  theme_continuum() +
  theme(legend.position = "bottom")
hallmark_trajectory_path <- save_panel(
  hallmark_trajectory_plot, "s3_six_fixed_hallmark_trajectories.pdf", 5.4, 4.2
)

collection_meta <- rbindlist(lapply(names(contract$pathway_collections), function(collection) {
  value <- fread(file.path(pathway_root, collection, "meta_analysis.tsv"))
  value[, collection := collection]
  value
}), use.names = TRUE, fill = TRUE)
collection_meta[, supported := fixed_q_value < 0.05 & direction_concordant]
collection_meta[, plot_q := pmax(fixed_q_value, .Machine$double.xmin)]
ml_write_tsv_once(collection_meta, file.path(source_out, "s3_all_pathway_meta.tsv.gz"))
pathway_skyline <- ggplot(
  collection_meta,
  aes(fixed_beta, -log10(plot_q), color = supported)
) +
  geom_point(size = 0.35, alpha = 0.5) +
  geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.25) +
  geom_hline(yintercept = -log10(0.05), color = "#9E9E9E", linewidth = 0.25,
             linetype = "22") +
  facet_grid(collection ~ axis_id, scales = "free_y") +
  scale_color_manual(values = c(`TRUE` = "#C9265E", `FALSE` = "#9E9E9E"), guide = "none") +
  labs(x = "Adjusted meta-analysis beta", y = "-log10 complete-family q") +
  theme_continuum()
pathway_skyline_path <- save_panel(pathway_skyline, "s3_all_pathway_family_skylines.pdf", 5.2, 6.4)

message("Rendering NMF and TF supplementary panels")
nmf_root <- file.path(candidate, "programs", "nmf")
nmf_meta <- fread(file.path(nmf_root, "nmf_meta_analysis.tsv"))
preferred_variant <- if (any(nmf_meta$score_variant == "signature_excluded_fixed_w" &
                             nmf_meta$inference_eligible == TRUE)) {
  "signature_excluded_fixed_w"
} else {
  "native_released_loading"
}
nmf_display <- nmf_meta[score_variant == preferred_variant]
nmf_display[, axis_label := fifelse(
  axis_id == "signature_pc1", "Cohort PC1", "Fixed projection"
)]
nmf_display[, outcome_label := paste0("k", k, " ", program_code, ": ", program_label)]
ml_write_tsv_once(nmf_display, file.path(source_out, "s4_all10_nmf_meta.tsv"))
nmf_forest <- ggplot(nmf_display, aes(fixed_beta, reorder(outcome_label, fixed_beta),
                                      color = axis_label)) +
  geom_vline(xintercept = 0, color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = fixed_ci_low, xmax = fixed_ci_high), height = 0,
                 linewidth = 0.35, position = position_dodge(width = 0.35)) +
  geom_point(size = 1.2, position = position_dodge(width = 0.35)) +
  scale_color_manual(values = c("Cohort PC1" = "#C9265E", "Fixed projection" = "#1565C0")) +
  labs(x = "Adjusted meta-analysis beta", y = NULL, color = NULL) +
  theme_continuum() +
  theme(legend.position = "bottom")
nmf_forest_path <- save_panel(nmf_forest, "s4_all10_nmf_forest.pdf", 4.4, 3.2)

nmf_windows <- fread(file.path(nmf_root, "nmf_fixed_windows.tsv.gz"))
nmf_windows <- nmf_windows[score_variant == preferred_variant]
nmf_windows[, outcome_label := paste0("k", k, " ", program_code, ": ", program_label)]
ml_write_tsv_once(nmf_windows, file.path(source_out, "s4_all10_nmf_windows.tsv"))
nmf_trajectory <- ggplot(nmf_windows, aes(window_center, mean_score,
                                         color = dataset, linetype = dataset)) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_line(linewidth = 0.5) + geom_point(size = 0.7) +
  facet_wrap(~outcome_label, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("GSE162694" = "#C9265E", "GSE213621" = "#1565C0")) +
  scale_linetype_manual(values = cohort_linetypes) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     breaks = c(0.1, 0.5, 0.9)) +
  labs(x = "Fixed-projection window center", y = "NMF loading (z)",
       color = NULL, linetype = NULL) +
  theme_continuum() + theme(legend.position = "bottom")
nmf_trajectory_path <- save_panel(nmf_trajectory, "s4_all10_nmf_trajectories.pdf", 5.2, 7.2)

tf_root <- file.path(candidate, "pathway_tf", "tf")
tf_windows <- fread(file.path(tf_root, "fixed_display_regulon_windows.tsv"))
tf_windows[, tf := factor(tf, levels = contract$display_tfs)]
ml_write_tsv_once(tf_windows, file.path(source_out, "s4_four_tf_regulon_windows.tsv"))
tf_testability <- fread(file.path(tf_root, "regulon_testability.tsv"))[
  tf %in% contract$display_tfs & dataset %in% contract$evaluation_cohorts
]
ml_write_tsv_once(
  tf_testability, file.path(source_out, "s4_four_tf_regulon_testability.tsv")
)
missing_tf_labels <- setdiff(contract$display_tfs, as.character(unique(tf_windows$tf)))
tf_placeholders <- if (length(missing_tf_labels)) {
  data.table(
    tf = factor(missing_tf_labels, levels = contract$display_tfs),
    center = rep(0.5, length(missing_tf_labels)),
    label = rep("Untestable after\nsignature exclusion", length(missing_tf_labels))
  )
} else {
  data.table(
    tf = factor(character(), levels = contract$display_tfs),
    center = numeric(), label = character()
  )
}
tf_trajectory <- ggplot(tf_windows, aes(center, mean_score,
                                       color = dataset, linetype = dataset)) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_line(linewidth = 0.5) + geom_point(size = 0.8) +
  geom_text(
    data = tf_placeholders,
    aes(center, 0, label = label), inherit.aes = FALSE,
    color = "#666666", size = 2
  ) +
  facet_wrap(~tf, nrow = 1, scales = "free_y", labeller = label_value,
             drop = FALSE) +
  scale_color_manual(values = c("GSE162694" = "#C9265E", "GSE213621" = "#1565C0")) +
  scale_linetype_manual(values = cohort_linetypes) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     breaks = c(0.1, 0.5, 0.9)) +
  labs(x = "Fixed-projection window center", y = "Signed regulon score (z)",
       color = NULL, linetype = NULL) +
  theme_continuum() +
  theme(strip.text = element_text(face = "italic"), legend.position = "bottom")
tf_trajectory_path <- save_panel(tf_trajectory, "s4_four_tf_regulon_trajectories.pdf", 5.2, 2.1)

figure_manifest <- data.table(
  figure_id = c("figure4f_candidate", "figure4f_compact", "s4_focal_windows", "s4_all117",
                "s4_paired_focal"),
  path = c(figure4f_path, compact_path, raw_window_path, all117_path, paired_figure),
  role = c("main_candidate", "review_only_compact", "supplementary_descriptive",
           "supplementary_complete_family", "supplementary_within_person"),
  inference_source = c("stage_sex_adjusted_models", "stage_sex_adjusted_models",
                       "fixed_windows_visualization_only", "stage_sex_adjusted_models",
                       "participant_delta_models"),
  current_figure3_modified = FALSE
)
figure_manifest <- rbindlist(list(
  figure_manifest,
  data.table(
    figure_id = c(
      "s3_fixed_gene_forest", "s3_all_gene_concordance", "s3_heldout_signature",
      "s3_fixed_gene_trajectories",
      "s3_all50_hallmark", "s3_six_hallmarks", "s3_pathway_skylines",
      "s4_all10_nmf_forest", "s4_all10_nmf_trajectories", "s4_four_tf_trajectories"
    ),
    path = c(
      gene_forest_path, canonical_concordance_path, heldout_path,
      gene_trajectory_path,
      hallmark_heatmap_path, hallmark_trajectory_path, pathway_skyline_path,
      nmf_forest_path, nmf_trajectory_path, tf_trajectory_path
    ),
    role = c(
      "supplementary_fixed_roster", "supplementary_same_substrate",
      "supplementary_heldout", "supplementary_fixed_roster",
      "supplementary_complete_family",
      "supplementary_fixed_roster", "supplementary_complete_families",
      "supplementary_complete_family", "supplementary_complete_family",
      "supplementary_fixed_roster"
    ),
    inference_source = c(
      "stage_sex_adjusted_models", "same_substrate_concordance",
      "source_overlap_signature_heldout_evaluation", "stage_sex_adjusted_models",
      "stage_sex_adjusted_models",
      "fixed_windows_visualization_only", "stage_sex_adjusted_models",
      "stage_sex_adjusted_models", "fixed_windows_visualization_only",
      "fixed_windows_visualization_only"
    ),
    current_figure3_modified = FALSE
  )
), use.names = TRUE)
figure_manifest[, exists := !is.na(path) & file.exists(path)]

message("Rendering three review-only Figure 4E/4F assemblies")
current_fig4e <- ml_resolve(contract$current_fig4e)
render_pdf_bitmap <- function(path, dpi = 180) {
  prefix <- tempfile(pattern = "hac_pdf_page_")
  status <- system2(
    "pdftoppm",
    c("-f", "1", "-l", "1", "-singlefile", "-png", "-r", dpi, path, prefix),
    stdout = TRUE, stderr = TRUE
  )
  exit_status <- attr(status, "status")
  ml_assert(is.null(exit_status) || identical(exit_status, 0L),
            paste0("pdftoppm failed for ", path, ": ", paste(status, collapse = " ")))
  png_path <- paste0(prefix, ".png")
  ml_assert(file.exists(png_path), paste0("pdftoppm did not create ", png_path))
  image <- png::readPNG(png_path)
  unlink(png_path)
  image
}
draw_bitmap <- function(bitmap, x, y, width, height) {
  grid::grid.raster(bitmap, x = x, y = y, width = width, height = height,
                    interpolate = TRUE)
}
save_review <- function(path, width, height, draw) {
  ml_assert(!file.exists(path), paste0("Refusing to overwrite review proof: ", path))
  grDevices::cairo_pdf(path, width = width, height = height, family = font_family)
  grid::grid.newpage()
  draw()
  grDevices::dev.off()
  path
}
bitmap_4e <- render_pdf_bitmap(current_fig4e)
bitmap_4f <- render_pdf_bitmap(figure4f_path)
bitmap_4f_compact <- render_pdf_bitmap(compact_path)

additive_path <- save_review(
  file.path(review_out, "review_additive_fig4e_plus_fullwidth_fig4f.pdf"), 6.2, 5.4,
  function() {
    draw_bitmap(bitmap_4e, 0.25, 0.72, 0.42, 0.47)
    draw_bitmap(bitmap_4f, 0.50, 0.25, 0.92, 0.43)
    grid::grid.text("current 4E", x = 0.04, y = 0.97,
                    gp = grid::gpar(fontsize = text_size, fontfamily = font_family))
    grid::grid.text("additive candidate 4F", x = 0.08, y = 0.48,
                    gp = grid::gpar(fontsize = text_size, fontfamily = font_family))
  }
)
replacement_path <- save_review(
  file.path(review_out, "review_space_neutral_fig4f_replaces_fig4e.pdf"), 3.2, 3.0,
  function() {
    draw_bitmap(bitmap_4f_compact, 0.50, 0.52, 0.92, 0.82)
    grid::grid.text("space-neutral candidate 4F; current 4E routes to S4",
                    x = 0.50, y = 0.05,
                    gp = grid::gpar(fontsize = text_size, fontfamily = font_family))
  }
)
combined_path <- save_review(
  file.path(review_out, "review_combined_fig4e_fig4f_row.pdf"), 7.0, 3.2,
  function() {
    draw_bitmap(bitmap_4e, 0.20, 0.52, 0.36, 0.82)
    draw_bitmap(bitmap_4f, 0.70, 0.52, 0.58, 0.82)
    grid::grid.text("current 4E", x = 0.03, y = 0.96,
                    gp = grid::gpar(fontsize = text_size, fontfamily = font_family))
    grid::grid.text("candidate 4F", x = 0.43, y = 0.96,
                    gp = grid::gpar(fontsize = text_size, fontfamily = font_family))
  }
)
assembly_manifest <- data.table(
  assembly = c("additive", "space_neutral_replacement", "combined_row"),
  path = c(additive_path, replacement_path, combined_path),
  review_only = TRUE,
  canonical_panels_remain_individual_vector_pdfs = TRUE
)
ml_write_tsv_once(assembly_manifest, file.path(figure_out, "review_assembly_manifest.tsv"))
ml_write_tsv_once(figure_manifest, file.path(figure_out, "figure_manifest.tsv"))
ml_write_session_info(file.path(figure_out, "sessionInfo.txt"))
message("MOLECULAR_LAYER_FIGURES_COMPLETE: ", figure_out)
