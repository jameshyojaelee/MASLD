#!/usr/bin/env Rscript
# Figure 3D endpoint PCA and synchronized endpoint-DE sensitivity.
# The five-cohort, 844-participant disease-vs-control model remains canonical.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(matrixStats)
})
grDevices::pdf.options(useDingbats = FALSE)
set.seed(42)

base <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
candidate_root <- normalizePath(Sys.getenv("FIGURE_CANDIDATE_ROOT", ""), mustWork = TRUE)
release_root <- normalizePath(Sys.getenv("STAGE_RELEASE_ROOT", ""), mustWork = TRUE)
stage_extension_root <- normalizePath(Sys.getenv("STAGE_EXTENSION_ROOT", ""), mustWork = TRUE)
source(file.path(base, "scripts/figures/publication_theme.R"))

panel_dir <- file.path(candidate_root, "figure3", "panels")
supp_dir <- file.path(candidate_root, "supplementary", "figureS3", "panels")
analysis_dir <- file.path(candidate_root, "analysis", "endpoint_validation")
source_dir <- file.path(candidate_root, "source_tables")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(file.path(release_root, "model_input_manifest.tsv"))
stopifnot(nrow(manifest) == 1L, manifest$n_genes == 23370L,
          manifest$n_samples == 844L, manifest$n_cohorts == 5L,
          manifest$n_control == 157L, manifest$n_disease == 687L)
dge_all <- readRDS(manifest$dge_path)
meta_all <- as.data.table(readRDS(manifest$meta_path))
meta <- meta_all[match(colnames(dge_all), sample_id)]
stopifnot(nrow(meta) == 844L, !anyNA(meta$sample_id),
          identical(meta$sample_id, colnames(dge_all)), uniqueN(meta$dataset) == 5L)

# Prespecified endpoint definitions. Advanced histology takes precedence for
# the single source-label conflict. Cohorts must contribute >=5 to both arms.
meta[, advanced := (!is.na(nas_score) & nas_score >= 5) |
                   (!is.na(fibrosis_stage) & fibrosis_stage >= 3)]
meta[, strict_control := group_binary == "Control" &
                         (is.na(nas_score) | nas_score <= 1) &
                         (is.na(fibrosis_stage) | fibrosis_stage <= 0)]
meta[, endpoint_code := fifelse(advanced, "Advanced",
                         fifelse(strict_control, "StrictControl", NA_character_))]
eligible_counts <- meta[!is.na(endpoint_code), .N, by = .(dataset, endpoint_code)]
eligible_wide <- dcast(eligible_counts, dataset ~ endpoint_code, value.var = "N", fill = 0L)
eligible_cohorts <- eligible_wide[Advanced >= 5L & StrictControl >= 5L, dataset]
endpoint_meta <- copy(meta[dataset %in% eligible_cohorts & !is.na(endpoint_code)])
endpoint_meta[, endpoint := factor(endpoint_code,
  levels = c("StrictControl", "Advanced"), labels = c("Strict control", "Advanced disease"))]
setorder(endpoint_meta, dataset, sample_id)
stopifnot(nrow(endpoint_meta) == 418L, uniqueN(endpoint_meta$dataset) == 4L,
          sum(endpoint_meta$endpoint_code == "Advanced") == 303L,
          sum(endpoint_meta$endpoint_code == "StrictControl") == 115L,
          !anyDuplicated(endpoint_meta$sample_id))

endpoint_counts <- endpoint_meta[, .(
  n = .N, n_fibrosis = sum(!is.na(fibrosis_stage)), n_nas = sum(!is.na(nas_score)),
  source_label_conflicts = sum((endpoint_code == "Advanced" & group_binary == "Control") |
                               (endpoint_code == "StrictControl" & group_binary == "Disease"))
), by = .(dataset, endpoint)]
fwrite(endpoint_meta[, .(sample_id, dataset, inferred_sex, source_group = group_binary,
                         endpoint, endpoint_code, fibrosis_stage, nas_score,
                         diagnosis_harmonized, condition)],
       file.path(analysis_dir, "endpoint_roster.tsv"), sep = "\t")
fwrite(endpoint_counts, file.path(analysis_dir, "endpoint_cohort_census.tsv"), sep = "\t")

idx <- match(endpoint_meta$sample_id, colnames(dge_all))
dge <- calcNormFactors(dge_all[, idx])
info <- data.frame(
  dataset = droplevels(factor(endpoint_meta$dataset)),
  inferred_sex = droplevels(factor(endpoint_meta$inferred_sex)),
  endpoint = factor(endpoint_meta$endpoint_code, levels = c("StrictControl", "Advanced")),
  row.names = endpoint_meta$sample_id
)
design <- model.matrix(~ dataset + inferred_sex + endpoint, info)
stopifnot(qr(design)$rank == ncol(design), "endpointAdvanced" %in% colnames(design))

v <- voomWithQualityWeights(dge, design, plot = FALSE)
fit0 <- lmFit(v, design)
fit <- eBayes(fit0)
tt <- topTable(fit, coef = "endpointAdvanced", number = Inf, sort.by = "none")
genes <- rownames(dge)
endpoint_de <- data.table(
  gene = genes, logFC = tt[genes, "logFC"], SE = abs(tt[genes, "logFC"] / tt[genes, "t"]),
  t = tt[genes, "t"], P.Value = tt[genes, "P.Value"],
  padj = p.adjust(tt[genes, "P.Value"], method = "BH"), AveExpr = tt[genes, "AveExpr"]
)
treat_fit <- treat(fit0, lfc = 0.25)
treat_tt <- topTreat(treat_fit, coef = "endpointAdvanced", number = Inf, sort.by = "none")
endpoint_de[, `:=`(
  treat_lfc = 0.25,
  treat_p = treat_tt[genes, "P.Value"],
  treat_fdr = p.adjust(treat_tt[genes, "P.Value"], method = "BH")
)]
canonical <- fread(file.path(release_root, "deg_results.csv"))
stopifnot(nrow(canonical) == 23370L, setequal(canonical$gene, endpoint_de$gene))
endpoint_de[, symbol := canonical$symbol[match(gene, canonical$gene)]]
fwrite(endpoint_de, file.path(analysis_dir, "endpoint_de_results.tsv.gz"), sep = "\t")

comparison <- merge(
  canonical[, .(gene, symbol, canonical_logFC = logFC, canonical_padj = padj,
                canonical_treat_fdr = treat_fdr)],
  endpoint_de[, .(gene, endpoint_logFC = logFC, endpoint_padj = padj,
                  endpoint_treat_fdr = treat_fdr)], by = "gene", all = FALSE
)
comparison[, canonical_deg := canonical_padj < 0.05 & abs(canonical_logFC) > 0.5]
comparison[, endpoint_deg := endpoint_padj < 0.05 & abs(endpoint_logFC) > 0.5]
comparison[, same_direction := sign(canonical_logFC) == sign(endpoint_logFC)]
comparison[, status := fcase(
  canonical_deg & endpoint_deg & same_direction, "Recovered",
  canonical_deg & !endpoint_deg & same_direction, "Subthreshold",
  canonical_deg & !same_direction, "Discordant",
  !canonical_deg & endpoint_deg, "Endpoint only",
  default = "Neither"
)]
n_can <- sum(comparison$canonical_deg)
n_end <- sum(comparison$endpoint_deg)
n_overlap <- sum(comparison$canonical_deg & comparison$endpoint_deg)
metrics <- data.table(
  metric = c("genes_tested", "canonical_degs", "endpoint_degs", "overlap_degs",
             "endpoint_only_degs", "logFC_spearman_all_genes", "direction_agreement_all_genes",
             "canonical_recovery", "jaccard", "canonical_treat_degs", "endpoint_treat_degs",
             "canonical_treat_recovery"),
  value = c(nrow(comparison), n_can, n_end, n_overlap,
            sum(!comparison$canonical_deg & comparison$endpoint_deg),
            cor(comparison$canonical_logFC, comparison$endpoint_logFC, method = "spearman"),
            mean(comparison$same_direction), n_overlap / n_can,
            n_overlap / sum(comparison$canonical_deg | comparison$endpoint_deg),
            sum(comparison$canonical_treat_fdr < 0.05),
            sum(comparison$endpoint_treat_fdr < 0.05),
            sum(comparison$canonical_treat_fdr < 0.05 & comparison$endpoint_treat_fdr < 0.05) /
              sum(comparison$canonical_treat_fdr < 0.05))
)
fwrite(comparison, file.path(analysis_dir, "endpoint_vs_canonical.tsv.gz"), sep = "\t")
fwrite(metrics, file.path(analysis_dir, "endpoint_vs_canonical_metrics.tsv"), sep = "\t")

# Figure 3C: direct five-cohort direction and leave-one-cohort-out stability.
# Per-cohort effects come from the synchronized stage-extension candidate; LOO
# models are recomputed here from the exact 844-participant substrate.
cohort_de <- fread(file.path(stage_extension_root, "cohort_disease_all_gene_results.tsv"))
stopifnot(nrow(cohort_de) == 5L * 23370L,
          all(cohort_de[, uniqueN(gene_id_versioned), by = dataset]$V1 == 23370L))
canonical_effects <- canonical[, .(gene, full_logFC = logFC, full_padj = padj,
                                   canonical_deg = padj < 0.05 & abs(logFC) > 0.5)]
cohort_join <- merge(cohort_de, canonical_effects,
                     by.x = "gene_id_versioned", by.y = "gene", all.x = TRUE)
stopifnot(!anyNA(cohort_join$full_logFC))
cohort_direction <- cohort_join[canonical_deg == TRUE, .(
  canonical_genes = .N,
  same_direction_n = sum(sign(logFC) == sign(full_logFC)),
  same_direction_fraction = mean(sign(logFC) == sign(full_logFC)),
  cohort_significant_n = sum(FDR < 0.05 & abs(logFC) > 0.5),
  n_control = unique(n_control), n_disease = unique(n_disease)
), by = dataset]

fit_loo <- function(held_out) {
  keep <- meta$dataset != held_out
  dge_loo <- calcNormFactors(dge_all[, keep])
  loo_meta <- meta[keep]
  loo_info <- data.frame(
    dataset = droplevels(factor(loo_meta$dataset)),
    inferred_sex = droplevels(factor(loo_meta$inferred_sex)),
    group = factor(loo_meta$group_binary, levels = c("Control", "Disease")),
    row.names = loo_meta$sample_id
  )
  loo_design <- model.matrix(~ dataset + inferred_sex + group, loo_info)
  stopifnot(qr(loo_design)$rank == ncol(loo_design), "groupDisease" %in% colnames(loo_design))
  loo_v <- voomWithQualityWeights(dge_loo, loo_design, plot = FALSE)
  loo_fit <- eBayes(lmFit(loo_v, loo_design))
  loo_tt <- topTable(loo_fit, coef = "groupDisease", number = Inf, sort.by = "none")
  out <- data.table(
    held_out_cohort = held_out, gene = rownames(dge_loo),
    loo_logFC = loo_tt[rownames(dge_loo), "logFC"],
    loo_padj = p.adjust(loo_tt[rownames(dge_loo), "P.Value"], method = "BH")
  )
  out[, loo_deg := loo_padj < 0.05 & abs(loo_logFC) > 0.5]
  out
}
loo_results <- rbindlist(lapply(sort(unique(meta$dataset)), fit_loo))
loo_join <- merge(loo_results, canonical_effects, by = "gene", all.x = TRUE)
loo_summary <- loo_join[, .(
  n_training_participants = 844L - sum(meta$dataset == unique(held_out_cohort)),
  all_gene_spearman = cor(loo_logFC, full_logFC, method = "spearman"),
  canonical_same_direction_fraction = mean(sign(loo_logFC[canonical_deg]) ==
                                             sign(full_logFC[canonical_deg])),
  canonical_gate_retained_n = sum(loo_deg & canonical_deg),
  canonical_gate_retained_fraction = sum(loo_deg & canonical_deg) / sum(canonical_deg)
), by = held_out_cohort]
fwrite(cohort_join, file.path(analysis_dir, "fig3c_per_cohort_effects.tsv.gz"), sep = "\t")
fwrite(loo_results, file.path(analysis_dir, "fig3c_loo_all_gene_results.tsv.gz"), sep = "\t")
fwrite(cohort_direction, file.path(analysis_dir, "fig3c_cohort_direction_summary.tsv"), sep = "\t")
fwrite(loo_summary, file.path(analysis_dir, "fig3c_loo_summary.tsv"), sep = "\t")
fig3c_source <- merge(cohort_direction, loo_summary,
                      by.x = "dataset", by.y = "held_out_cohort", all = TRUE)
fwrite(fig3c_source, file.path(source_dir, "fig3c_cohort_robustness.tsv"), sep = "\t")

model_audit <- data.table(
  model = "Advanced disease versus strict control sensitivity",
  n_participants = nrow(endpoint_meta), n_cohorts = uniqueN(endpoint_meta$dataset),
  n_strict_control = sum(endpoint_meta$endpoint_code == "StrictControl"),
  n_advanced = sum(endpoint_meta$endpoint_code == "Advanced"),
  covariates = "fixed cohort + inferred sex", bh_family_size = nrow(endpoint_de),
  canonical_write = FALSE, seed = 42L
)
fwrite(model_audit, file.path(analysis_dir, "endpoint_model_audit.tsv"), sep = "\t")

# Group-blind PCA: cohort and inferred sex are removed without preserving the
# endpoint label. HVGs are selected after this adjustment.
sex_cov <- model.matrix(~ inferred_sex, info)[, -1L, drop = FALSE]
adjusted <- removeBatchEffect(v$E, batch = info$dataset, covariates = sex_cov)
hvg <- order(rowVars(adjusted), decreasing = TRUE)[seq_len(2000L)]
pca <- prcomp(t(adjusted[hvg, , drop = FALSE]), center = TRUE, scale. = FALSE)
if (mean(pca$x[info$endpoint == "Advanced", 1]) <
    mean(pca$x[info$endpoint == "StrictControl", 1])) pca$x[, 1] <- -pca$x[, 1]
pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
coords <- cbind(endpoint_meta, data.table(PC1 = pca$x[, 1], PC2 = pca$x[, 2]))
fwrite(coords[, .(sample_id, dataset, inferred_sex, source_group = group_binary,
                  endpoint, fibrosis_stage, nas_score, PC1, PC2)],
       file.path(source_dir, "fig3d_endpoint_pca_coordinates.tsv"), sep = "\t")
fwrite(endpoint_counts, file.path(source_dir, "fig3d_endpoint_census.tsv"), sep = "\t")
fwrite(metrics, file.path(source_dir, "figs3_endpoint_deg_concordance_metrics.tsv"), sep = "\t")

source(file.path(base, "scripts/figures/render_fig3d_endpoint_pca.R"))
render_endpoint_pca(coords, pve,
  file.path(panel_dir, "fig3d_pca_fibrosis_gradient.pdf"))

source(file.path(base, "scripts/figures/render_fig3c_robustness.R"))

plot_data <- comparison[canonical_deg | endpoint_deg]
plot_data[, status := factor(status,
  levels = c("Recovered", "Subthreshold", "Discordant", "Endpoint only"))]
status_colors <- c(Recovered = masld_colors$mash, Subthreshold = red_gradient[[2L]],
                   Discordant = "#4C78A8", `Endpoint only` = "#BDBDBD")
lim <- max(abs(c(plot_data$canonical_logFC, plot_data$endpoint_logFC)), na.rm = TRUE) * 1.03
concordance_plot <- ggplot(plot_data, aes(canonical_logFC, endpoint_logFC, color = status)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linewidth = 0.3, linetype = "dashed") +
  geom_point(size = 0.45, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = status_colors, name = NULL) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "Five-cohort disease vs control log2FC",
       y = "Advanced disease vs strict control log2FC") +
  theme_masld(base_size = 6) + theme_pub() +
  theme(legend.position = "bottom")
ggsave(file.path(supp_dir, "figs3_endpoint_deg_concordance.pdf"), concordance_plot,
       width = 3.25, height = 3.0, device = cairo_pdf)

writeLines(capture.output(sessionInfo()), file.path(candidate_root, "sessionInfo.txt"))
cat(sprintf("ENDPOINT_VALIDATION_COMPLETE n=%d controls=%d advanced=%d cohorts=%d\n",
            nrow(endpoint_meta), sum(endpoint_meta$endpoint_code == "StrictControl"),
            sum(endpoint_meta$endpoint_code == "Advanced"), uniqueN(endpoint_meta$dataset)))
print(metrics)
