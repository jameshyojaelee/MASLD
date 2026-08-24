#!/usr/bin/env Rscript
# KEY MESSAGE: Fibrosis stage and the molecular continuum recover aligned remodeling, while the continuum retains substantial within-stage molecular structure.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 62_render_stage_continuum_comparison.R <analysis_candidate> <figure_candidate>")
}

analysis_candidate <- normalizePath(args[[1]], mustWork = TRUE)
figure_candidate <- args[[2]]
if (dir.exists(figure_candidate)) {
  stop("Refusing to overwrite existing figure candidate: ", figure_candidate)
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lib_molecular_layers.R"))

contract <- ml_read_contract()
panel_dir <- file.path(figure_candidate, "panels")
source_dir <- file.path(figure_candidate, "source_tables")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
grDevices::pdf.options(useDingbats = FALSE)

cohort_colors <- c("GSE162694" = "#C9265E", "GSE213621" = "#1565C0")
support_colors <- c(
  "Neither" = "#9E9E9E",
  "Stage only" = "#F4A674",
  "Continuum only" = "#7B1FA2",
  "Both" = "#00695C"
)

theme_compare <- function() {
  theme_classic(base_size = 6, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(size = 6, face = "plain"),
      axis.text = element_text(size = 6, color = "black", face = "plain"),
      axis.title = element_text(size = 6, face = "plain"),
      legend.text = element_text(size = 6, face = "plain"),
      legend.title = element_text(size = 6, face = "plain"),
      strip.text = element_text(size = 6, face = "plain"),
      strip.background = element_blank(),
      plot.title = element_blank(),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.key.size = grid::unit(3, "mm"),
      legend.position = "bottom",
      plot.margin = margin(2, 2, 2, 2)
    )
}

save_panel <- function(plot, filename, width, height) {
  path <- file.path(panel_dir, filename)
  stopifnot(!file.exists(path))
  ggsave(path, plot, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white")
  path
}

fixed_meta_table <- function(cohort_table, feature_columns, family_size) {
  feature_columns <- as.character(feature_columns)
  result <- cohort_table[, {
    meta <- ml_fixed_meta(beta, se)
    direction_concordant <- sum(is.finite(beta)) == 2L &&
      length(unique(sign(beta[is.finite(beta)]))) == 1L
    cbind(meta, direction_concordant = direction_concordant)
  }, by = feature_columns]
  result[, q_value := ml_complete_bh(p_value, family_size)]
  result[, bh_family_size := as.integer(family_size)]
  result
}

fit_stage_groups <- function(data, feature_columns, outcome_column = "outcome_z") {
  feature_columns <- as.character(feature_columns)
  data[
    is.finite(get(outcome_column)) & !is.na(fibrosis_stage),
    .(
      n_participants = .N,
      mean_score = mean(get(outcome_column)),
      se_score = sd(get(outcome_column)) / sqrt(.N)
    ),
    by = c(feature_columns, "dataset", "fibrosis_stage")
  ]
}

fit_stage_lm_groups <- function(data, feature_columns, outcome_column = "outcome_z") {
  feature_columns <- as.character(feature_columns)
  data[
    is.finite(get(outcome_column)) & !is.na(fibrosis_stage) & !is.na(inferred_sex),
    {
      stage_z <- ml_standardize(fibrosis_stage)
      fit <- tryCatch(lm(get(outcome_column) ~ stage_z + factor(inferred_sex)),
                      error = function(e) NULL)
      if (is.null(fit) || !"stage_z" %in% rownames(coef(summary(fit)))) {
        .(n = .N, beta = NA_real_, se = NA_real_, p_value = NA_real_)
      } else {
        coefficient <- coef(summary(fit))["stage_z", ]
        .(n = .N, beta = unname(coefficient[["Estimate"]]),
          se = unname(coefficient[["Std. Error"]]),
          p_value = unname(coefficient[["Pr(>|t|)"]]))
      }
    },
    by = c(feature_columns, "dataset")
  ]
}

make_trajectory_pair <- function(stage_data, window_data, outcome_column,
                                 row_order, stage_file, continuum_file,
                                 height, window_center = "window_center") {
  stage_data[, outcome_label := factor(get(outcome_column), levels = row_order)]
  window_data[, outcome_label := factor(get(outcome_column), levels = row_order)]
  stage_plot <- ggplot(
    stage_data,
    aes(fibrosis_stage, mean_score, color = dataset, linetype = dataset, group = dataset)
  ) +
    geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
    geom_errorbar(aes(ymin = mean_score - 1.96 * se_score,
                      ymax = mean_score + 1.96 * se_score),
                  width = 0, linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1) +
    facet_wrap(~outcome_label, ncol = 1, scales = "free_y") +
    scale_color_manual(values = cohort_colors) +
    scale_linetype_manual(values = c("GSE162694" = "solid", "GSE213621" = "22")) +
    scale_x_continuous(breaks = 0:4, labels = paste0("F", 0:4)) +
    labs(x = "Fibrosis stage", y = "Mean score (z)", color = NULL, linetype = NULL) +
    theme_compare() +
    theme(legend.position = "none")

  continuum_plot <- ggplot(
    window_data,
    aes(.data[[window_center]], mean_score, color = dataset,
        linetype = dataset, group = dataset)
  ) +
    geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
    geom_errorbar(aes(ymin = mean_score - 1.96 * se_score,
                      ymax = mean_score + 1.96 * se_score),
                  width = 0, linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1) +
    facet_wrap(~outcome_label, ncol = 1, scales = "free_y") +
    scale_color_manual(values = cohort_colors) +
    scale_linetype_manual(values = c("GSE162694" = "solid", "GSE213621" = "22")) +
    scale_x_continuous(breaks = c(0.1, 0.5, 0.9),
                       labels = c("10%", "50%", "90%")) +
    labs(x = "Fixed-projection window", y = "Mean score (z)",
         color = NULL, linetype = NULL) +
    theme_compare()

  save_panel(stage_plot, stage_file, 3.2, height)
  save_panel(continuum_plot, continuum_file, 3.2, height)
}

message("Fitting same-evaluation-cohort fibrosis-stage transcript trends")
manifest <- fread(ml_resolve(contract$sample_manifest), na.strings = c("", "NA"))
dge <- readRDS(ml_resolve(contract$dge_rds))
annotation <- fread(ml_resolve(contract$gene_annotation), select = c("gene_id", "gene_name"))
annotation[, gene_id_base := ml_base_gene_id(gene_id)]
annotation <- unique(annotation[, .(gene_id_base, gene_name)], by = "gene_id_base")
gene_ids <- ml_base_gene_id(rownames(dge$counts))
stopifnot(length(gene_ids) == contract$gene_family_size, !anyDuplicated(gene_ids))

stage_gene_cohort <- rbindlist(lapply(contract$evaluation_cohorts, function(cohort) {
  d <- manifest[
    dataset == cohort & !is.na(fibrosis_stage) & !is.na(inferred_sex)
  ]
  d[, `:=`(
    stage_z = ml_standardize(fibrosis_stage),
    sex_factor = factor(inferred_sex)
  )]
  design <- model.matrix(~ stage_z + sex_factor, data = d)
  stopifnot(qr(design)$rank == ncol(design))
  subset <- dge[, d$sample_id, keep.lib.sizes = FALSE]
  subset <- calcNormFactors(subset)
  voom <- voomWithQualityWeights(subset, design = design, plot = FALSE)
  fit <- eBayes(lmFit(voom, design), robust = TRUE)
  coefficient <- match("stage_z", colnames(design))
  moderated_se <- fit$stdev.unscaled[, coefficient] * sqrt(fit$s2.post)
  data.table(
    gene_id_base = gene_ids,
    dataset = cohort,
    beta = as.numeric(fit$coefficients[, coefficient]),
    se = as.numeric(moderated_se),
    p_value = as.numeric(fit$p.value[, coefficient]),
    n_participants = nrow(d),
    model = "expression ~ fibrosis_stage_z + inferred_sex"
  )
}))
stage_gene_meta <- fixed_meta_table(
  stage_gene_cohort, "gene_id_base", contract$gene_family_size
)
stage_gene_meta[, stage_supported := estimable & direction_concordant & q_value < 0.05]

continuum_gene <- fread(file.path(analysis_candidate, "bulk", "continuum_gene_meta.tsv.gz"))[
  axis_id == "fixed_projection" & inference_eligible == TRUE,
  .(gene_id_base, continuum_beta = beta, continuum_q = bh_q_value)
]
membership <- fread(file.path(analysis_candidate, "bulk", "continuum_membership.tsv.gz"))[
  axis_component == FALSE,
  .(gene_id_base, gene_name, continuum_associated)
]
transcript_compare <- Reduce(
  function(x, y) merge(x, y, by = "gene_id_base", all = FALSE),
  list(stage_gene_meta, continuum_gene, membership)
)
transcript_compare[, support_class := fifelse(
  stage_supported & continuum_associated, "Both",
  fifelse(stage_supported, "Stage only",
          fifelse(continuum_associated, "Continuum only", "Neither"))
)]
transcript_compare[, support_class := factor(support_class, levels = names(support_colors))]
transcript_rho <- cor(
  transcript_compare$beta,
  transcript_compare$continuum_beta,
  method = "spearman",
  use = "complete.obs"
)
fwrite(stage_gene_cohort, file.path(source_dir, "transcript_stage_trend_cohort.tsv.gz"), sep = "\t")
fwrite(transcript_compare, file.path(source_dir, "transcript_stage_vs_continuum.tsv.gz"), sep = "\t")

transcript_plot <- ggplot(
  transcript_compare,
  aes(beta, continuum_beta, color = support_class)
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_point(size = 0.32, alpha = 0.38, stroke = 0) +
  scale_color_manual(values = support_colors, drop = FALSE, name = NULL) +
  labs(x = "Fibrosis-stage trend beta",
       y = "Stage-adjusted fixed-projection beta") +
  theme_compare() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(size = 1.4, alpha = 1)))
transcript_path <- save_panel(
  transcript_plot, "s3_transcript_stage_vs_continuum_effects.pdf", 3.7, 3.55
)

message("Rendering adjacent-stage transcript contrast comparison")
adjacent <- fread(ml_resolve(contract$stage_results))[
  axis == "fibrosis_adjacent",
  .(gene_id_base = ml_base_gene_id(gene_id_base), gene_name, contrast,
    stage_logFC = logFC, stage_fdr = FDR)
]
adjacent <- merge(adjacent, continuum_gene, by = "gene_id_base", all = FALSE)
adjacent <- merge(adjacent, membership, by = "gene_id_base", all = FALSE,
                  suffixes = c("_stage", ""))
adjacent[, stage_deg := is.finite(stage_fdr) & stage_fdr < 0.05 & abs(stage_logFC) > 0.5]
adjacent[, support_class := fifelse(
  stage_deg & continuum_associated, "Both",
  fifelse(stage_deg, "Stage only",
          fifelse(continuum_associated, "Continuum only", "Neither"))
)]
adjacent[, support_class := factor(support_class, levels = names(support_colors))]
adjacent[, contrast_label := factor(
  contrast,
  levels = c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"),
  labels = c("F0 to F1", "F1 to F2", "F2 to F3", "F3 to F4")
)]
fwrite(adjacent, file.path(source_dir, "adjacent_stage_deg_vs_continuum.tsv.gz"), sep = "\t")
adjacent_plot <- ggplot(adjacent, aes(stage_logFC, continuum_beta, color = support_class)) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_point(size = 0.24, alpha = 0.28, stroke = 0) +
  facet_wrap(~contrast_label, nrow = 1) +
  scale_color_manual(values = support_colors, drop = FALSE, name = NULL) +
  labs(x = "Adjacent-stage log2 fold change",
       y = "Stage-adjusted fixed-projection beta") +
  theme_compare() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(size = 1.4, alpha = 1)))
adjacent_path <- save_panel(
  adjacent_plot, "s3_adjacent_stage_degs_vs_continuum.pdf", 6.6, 2.35
)

message("Comparing all 117 frozen Hotspot programs")
hotspot_root <- file.path(analysis_candidate, "programs", "hotspot")
hotspot_scores <- fread(file.path(hotspot_root, "hotspot_participant_scores.tsv.gz"))[
  dataset %in% contract$evaluation_cohorts
]
hotspot_stage_cohort <- fit_stage_lm_groups(hotspot_scores, "program_uid")
hotspot_stage_meta <- fixed_meta_table(
  hotspot_stage_cohort, "program_uid", contract$program_family_size
)
hotspot_continuum <- fread(file.path(hotspot_root, "hotspot_meta_analysis.tsv"))[
  axis_id == "fixed_projection",
  .(program_uid, continuum_beta = fixed_beta, continuum_q = fixed_q_value,
    module_name, cell_type)
]
hotspot_membership <- fread(file.path(hotspot_root, "hotspot_continuum_membership.tsv"))[
  , .(program_uid, continuum_associated)
]
hotspot_compare <- Reduce(
  function(x, y) merge(x, y, by = "program_uid", all = TRUE),
  list(hotspot_stage_meta, hotspot_continuum, hotspot_membership)
)
hotspot_compare[, stage_supported := estimable & direction_concordant & q_value < 0.05]
hotspot_compare[, support_class := fifelse(
  stage_supported & continuum_associated, "Both",
  fifelse(stage_supported, "Stage only",
          fifelse(continuum_associated, "Continuum only", "Neither"))
)]
hotspot_compare[, support_class := factor(support_class, levels = names(support_colors))]
fwrite(hotspot_compare, file.path(source_dir, "hotspot_stage_vs_continuum.tsv"), sep = "\t")
hotspot_effect_plot <- ggplot(
  hotspot_compare[is.finite(beta) & is.finite(continuum_beta)],
  aes(beta, continuum_beta, color = support_class)
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_point(size = 1, alpha = 0.75) +
  scale_color_manual(values = support_colors, drop = FALSE, name = NULL) +
  labs(x = "Fibrosis-stage trend beta",
       y = "Stage-adjusted fixed-projection beta") +
  theme_compare()
hotspot_effect_path <- save_panel(
  hotspot_effect_plot, "s4_all117_hotspot_stage_vs_continuum_effects.pdf", 3.5, 3.35
)

focal_lookup <- data.table(
  program_uid = unname(unlist(contract$focal_programs)),
  focal_program = names(contract$focal_programs)
)
hotspot_stage <- merge(hotspot_scores, focal_lookup, by = "program_uid")
hotspot_stage <- fit_stage_groups(hotspot_stage, c("program_uid", "focal_program"))
hotspot_windows <- merge(
  fread(file.path(hotspot_root, "hotspot_fixed_windows.tsv.gz")),
  focal_lookup,
  by = "program_uid"
)[axis_id == "fixed_projection"]
fwrite(hotspot_stage, file.path(source_dir, "focal_hotspot_stage_groups.tsv"), sep = "\t")
fwrite(hotspot_windows, file.path(source_dir, "focal_hotspot_continuum_windows.tsv"), sep = "\t")
make_trajectory_pair(
  hotspot_stage, hotspot_windows, "focal_program",
  names(contract$focal_programs),
  "s4_focal_hotspot_by_fibrosis_stage.pdf",
  "s4_focal_hotspot_by_continuum_windows.pdf",
  2.55
)

message("Comparing six frozen pathway collections")
pathway_root <- file.path(analysis_candidate, "pathway_tf", "pathways")
integrated_registry <- fread(file.path(
  analysis_candidate, "decision", "continuum_membership_registry.tsv.gz"
))
pathway_compare_rows <- list()
pathway_stage_group_rows <- list()
for (collection in names(contract$pathway_collections)) {
  family_size <- as.integer(contract$pathway_collections[[collection]])
  scores <- fread(file.path(pathway_root, collection, "donor_scores.tsv.gz"))
  scores <- merge(
    scores,
    manifest[, .(sample_id, dataset, inferred_sex, fibrosis_stage)],
    by = c("sample_id", "dataset"), all.x = TRUE
  )[dataset %in% contract$evaluation_cohorts]
  scores[, outcome_z := ml_standardize(pathway_score), by = .(dataset, set_id)]
  stage_cohort <- fit_stage_lm_groups(scores, "set_id")
  stage_meta <- fixed_meta_table(stage_cohort, "set_id", family_size)
  stage_meta[, stage_supported := estimable & direction_concordant & q_value < 0.05]
  continuum <- fread(file.path(pathway_root, collection, "meta_analysis.tsv"))[
    axis_id == "fixed_projection",
    .(set_id, continuum_beta = fixed_beta, continuum_q = fixed_q_value)
  ]
  continuum_membership <- integrated_registry[
    molecular_layer == paste0("pathway_", collection),
    .(set_id = feature_id, continuum_supported = continuum_associated)
  ]
  comparison <- Reduce(
    function(x, y) merge(x, y, by = "set_id", all = TRUE),
    list(stage_meta, continuum, continuum_membership)
  )
  comparison[, support_class := fifelse(
    stage_supported & continuum_supported, "Both",
    fifelse(stage_supported, "Stage only",
            fifelse(continuum_supported, "Continuum only", "Neither"))
  )]
  comparison[, `:=`(
    collection = collection,
    family_size = family_size,
    support_class = factor(support_class, levels = names(support_colors))
  )]
  pathway_compare_rows[[collection]] <- comparison
  if (collection == "hallmark") {
    pathway_stage_group_rows[[collection]] <- fit_stage_groups(
      scores[set_id %in% contract$display_hallmarks], "set_id"
    )
  }
}
pathway_compare <- rbindlist(pathway_compare_rows, use.names = TRUE, fill = TRUE)
pathway_compare[, collection_label := factor(
  collection,
  levels = c("hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc"),
  labels = c("Hallmark", "KEGG", "Reactome", "GO BP", "GO MF", "GO CC")
)]
fwrite(pathway_compare, file.path(source_dir, "pathway_stage_vs_continuum.tsv.gz"), sep = "\t")
pathway_effect_plot <- ggplot(
  pathway_compare[is.finite(beta) & is.finite(continuum_beta)],
  aes(beta, continuum_beta, color = support_class)
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_point(size = 0.3, alpha = 0.32, stroke = 0) +
  facet_wrap(~collection_label, ncol = 3, scales = "free") +
  scale_color_manual(values = support_colors, drop = FALSE, name = NULL) +
  labs(x = "Fibrosis-stage trend beta",
       y = "Stage-adjusted fixed-projection beta") +
  theme_compare() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(size = 1.4, alpha = 1)))
pathway_effect_path <- save_panel(
  pathway_effect_plot, "s3_pathways_stage_vs_continuum_effects.pdf", 5.2, 4.35
)

hallmark_stage <- pathway_stage_group_rows$hallmark
hallmark_windows <- fread(file.path(pathway_root, "hallmark", "fixed_display_windows.tsv"))
fwrite(hallmark_stage, file.path(source_dir, "six_hallmark_stage_groups.tsv"), sep = "\t")
fwrite(hallmark_windows, file.path(source_dir, "six_hallmark_continuum_windows.tsv"), sep = "\t")
hallmark_labels <- setNames(
  gsub("_", " ", sub("^HALLMARK_", "", contract$display_hallmarks)),
  contract$display_hallmarks
)
hallmark_stage[, set_label := hallmark_labels[set_id]]
hallmark_windows[, set_label := hallmark_labels[set_id]]
make_trajectory_pair(
  hallmark_stage, hallmark_windows, "set_label",
  unname(hallmark_labels[contract$display_hallmarks]),
  "s3_six_hallmarks_by_fibrosis_stage.pdf",
  "s3_six_hallmarks_by_continuum_windows.pdf",
  5.65,
  window_center = "center"
)

message("Comparing all ten fixed NMF axes")
nmf_root <- file.path(analysis_candidate, "programs", "nmf")
nmf_scores <- fread(file.path(nmf_root, "nmf_participant_scores.tsv.gz"))[
  score_variant == "signature_excluded_fixed_w" & inference_eligible == TRUE &
    dataset %in% contract$evaluation_cohorts
]
nmf_stage_cohort <- fit_stage_lm_groups(
  nmf_scores, c("outcome_id", "k", "program_code", "program_label", "score_variant")
)
nmf_stage_meta <- fixed_meta_table(
  nmf_stage_cohort,
  c("outcome_id", "k", "program_code", "program_label", "score_variant"),
  contract$nmf_family_size
)
nmf_stage_meta[, stage_supported := estimable & direction_concordant & q_value < 0.05]
nmf_continuum <- fread(file.path(nmf_root, "nmf_meta_analysis.tsv"))[
  score_variant == "signature_excluded_fixed_w" & axis_id == "fixed_projection",
  .(outcome_id, continuum_beta = fixed_beta, continuum_q = fixed_q_value)
]
nmf_membership <- fread(file.path(nmf_root, "nmf_continuum_membership.tsv"))[
  , .(outcome_id, continuum_associated)
]
nmf_compare <- Reduce(
  function(x, y) merge(x, y, by = "outcome_id", all = TRUE),
  list(nmf_stage_meta, nmf_continuum, nmf_membership)
)
nmf_compare[, support_class := fifelse(
  stage_supported & continuum_associated, "Both",
  fifelse(stage_supported, "Stage only",
          fifelse(continuum_associated, "Continuum only", "Neither"))
)]
nmf_compare[, support_class := factor(support_class, levels = names(support_colors))]
nmf_compare[, axis_label := paste0("k", k, " ", program_code, ": ", program_label)]
fwrite(nmf_compare, file.path(source_dir, "nmf_stage_vs_continuum.tsv"), sep = "\t")
nmf_effect_plot <- ggplot(
  nmf_compare[is.finite(beta) & is.finite(continuum_beta)],
  aes(beta, continuum_beta, color = support_class)
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_point(size = 1.2) +
  geom_text(aes(label = paste0("k", k, " ", program_code)),
            size = 6 / ggplot2::.pt, nudge_y = 0.035,
            check_overlap = TRUE, show.legend = FALSE) +
  scale_color_manual(values = support_colors, drop = FALSE, name = NULL) +
  labs(x = "Fibrosis-stage trend beta",
       y = "Stage-adjusted fixed-projection beta") +
  theme_compare()
nmf_effect_path <- save_panel(
  nmf_effect_plot, "s4_all10_nmf_stage_vs_continuum_effects.pdf", 3.5, 3.35
)

nmf_stage <- fit_stage_groups(
  nmf_scores, c("outcome_id", "k", "program_code", "program_label", "score_variant")
)
nmf_stage[, axis_label := paste0("k", k, " ", program_code, ": ", program_label)]
nmf_windows <- fread(file.path(nmf_root, "nmf_fixed_windows.tsv.gz"))[
  score_variant == "signature_excluded_fixed_w" & axis_id == "fixed_projection"
]
nmf_windows[, axis_label := paste0("k", k, " ", program_code, ": ", program_label)]
nmf_order <- nmf_compare[order(k, program_code), axis_label]
fwrite(nmf_stage, file.path(source_dir, "all10_nmf_stage_groups.tsv"), sep = "\t")
fwrite(nmf_windows, file.path(source_dir, "all10_nmf_continuum_windows.tsv"), sep = "\t")
make_trajectory_pair(
  nmf_stage, nmf_windows, "axis_label", nmf_order,
  "s4_all10_nmf_by_fibrosis_stage.pdf",
  "s4_all10_nmf_by_continuum_windows.pdf",
  8.7
)

summary <- rbindlist(list(
  data.table(
    layer = "transcript",
    family_size = contract$gene_family_size,
    # na.rm matches every other layer below: stage_supported is NA whenever
    # estimable is TRUE but q_value is NA, which would turn these counts into NA.
    n_stage_supported = sum(transcript_compare$stage_supported, na.rm = TRUE),
    n_continuum_supported = sum(transcript_compare$continuum_associated, na.rm = TRUE),
    n_both = sum(transcript_compare$support_class == "Both", na.rm = TRUE),
    spearman = transcript_rho
  ),
  data.table(
    layer = "hotspot",
    family_size = contract$program_family_size,
    n_stage_supported = sum(hotspot_compare$stage_supported, na.rm = TRUE),
    n_continuum_supported = sum(hotspot_compare$continuum_associated, na.rm = TRUE),
    n_both = sum(hotspot_compare$support_class == "Both", na.rm = TRUE),
    spearman = cor(hotspot_compare$beta, hotspot_compare$continuum_beta,
                   method = "spearman", use = "complete.obs")
  ),
  pathway_compare[, .(
    layer = paste0("pathway_", collection[[1]]),
    family_size = family_size[[1]],
    n_stage_supported = sum(stage_supported, na.rm = TRUE),
    n_continuum_supported = sum(continuum_supported, na.rm = TRUE),
    n_both = sum(support_class == "Both", na.rm = TRUE),
    spearman = cor(beta, continuum_beta, method = "spearman", use = "complete.obs")
  ), by = collection][, collection := NULL],
  data.table(
    layer = "nmf",
    family_size = contract$nmf_family_size,
    n_stage_supported = sum(nmf_compare$stage_supported, na.rm = TRUE),
    n_continuum_supported = sum(nmf_compare$continuum_associated, na.rm = TRUE),
    n_both = sum(nmf_compare$support_class == "Both", na.rm = TRUE),
    spearman = cor(nmf_compare$beta, nmf_compare$continuum_beta,
                   method = "spearman", use = "complete.obs")
  )
), use.names = TRUE, fill = TRUE)
summary[, interpretation := "same_expression_substrate_method_comparison_not_independent_validation"]
fwrite(summary, file.path(figure_candidate, "comparison_summary.tsv"), sep = "\t")

files <- c(
  transcript_path, adjacent_path, hotspot_effect_path,
  pathway_effect_path, nmf_effect_path,
  list.files(panel_dir, pattern = "by_(fibrosis_stage|continuum_windows)\\.pdf$",
             full.names = TRUE)
)
stopifnot(length(files) == 11L, all(file.exists(files)))
writeLines(capture.output(sessionInfo()), file.path(figure_candidate, "sessionInfo.txt"))
manifest_out <- data.table(
  artifact = basename(files),
  path = normalizePath(files),
  sha256 = vapply(files, ml_sha256, character(1)),
  claim_boundary = "same-substrate comparison; continuum windows descriptive only"
)
fwrite(manifest_out, file.path(figure_candidate, "figure_manifest.tsv"), sep = "\t")
