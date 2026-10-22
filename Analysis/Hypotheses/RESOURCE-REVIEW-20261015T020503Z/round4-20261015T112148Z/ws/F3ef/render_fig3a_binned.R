# KEY MESSAGE: GSE213621 fibrosis is recorded only as source bins (F0-F1, F2, F3-F4) and does not enter the exact-stage contrasts of Fig 3E/3F.
# Candidate Fig 3A (F3ef). Adapted from scripts/figures/render_pi_bulk_opening_panels.R (3A block; copied as
# render_pi_bulk_opening_panels.orig.R). Change: a fourth tile state "Binned" for GSE213621 fibrosis.
# The installed panel shows the Sex/NAS/Fibrosis columns and a Units bar, so this candidate does too.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
set.seed(20260923)
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
source(file.path(base, "scripts/figures/publication_theme.R"))
grDevices::pdf.options(useDingbats = FALSE)
A <- file.path(base, "figures/candidates/fig3e-s3-774-20261011T165944Z/analysis/stage_extensions")
out <- file.path(base, "figures/candidates/review-r3-fixes-20261015T112148Z/F3ef")
samples <- fread(file.path(A, "five_cohort_sample_manifest.tsv"))
summary <- fread(file.path(A, "five_cohort_summary.tsv"))
stopifnot(nrow(samples) == 774)
cohorts <- summary[order(-n_participants), dataset]
binned <- "GSE213621"   # GEO fibrosis is Control / F0-F1 / F2 / F3-F4; stored as 0/1/2/3
availability <- rbindlist(lapply(cohorts, function(cohort) {
  d <- samples[dataset == cohort]
  data.table(dataset = cohort, field = c("Sex", "NAS", "Fibrosis"),
             n_observed = c(sum(!is.na(d$inferred_sex)), sum(!is.na(d$nas_score)), sum(!is.na(d$fibrosis_stage))),
             n_units = nrow(d))
}))
availability[, fraction := n_observed / n_units]
availability[, state := as.character(cut(fraction, breaks = c(-Inf, 0, 0.99, Inf), labels = c("Absent", "Partial", "Available")))]
availability[dataset == binned & field == "Fibrosis", state := "Binned"]
availability[, state := factor(state, levels = c("Absent", "Partial", "Binned", "Available"))]
availability[, dataset := factor(dataset, levels = rev(cohorts))]
availability[, field := factor(field, levels = c("Sex", "NAS", "Fibrosis"))]
p_meta <- ggplot(availability, aes(field, dataset, fill = state)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(values = c(Absent = "#ECECEC", Partial = cat_palette[[5]], Binned = cat_palette[[6]],
                               Available = cat_palette[[1]]), name = NULL, drop = FALSE) +
  labs(x = NULL, y = NULL) + theme_masld_compact() +
  theme(axis.text.x = element_text(size = 6, face = "plain", angle = 35, hjust = 1),
        axis.text.y = element_text(size = 6, face = "plain"), legend.position = "bottom")
composition <- melt(summary[, .(dataset, Control = n_control, Disease = n_disease)], id.vars = "dataset",
                    variable.name = "group", value.name = "n")
composition[, dataset := factor(dataset, levels = rev(cohorts))]
totals <- summary[, .(dataset = factor(dataset, levels = rev(cohorts)), n_participants)]
p_bar <- ggplot(composition, aes(dataset, n, fill = group)) +
  geom_col(width = 0.72) +
  geom_text(data = totals, aes(dataset, n_participants, label = n_participants), inherit.aes = FALSE,
            hjust = -0.2, size = 6 / .pt, color = "black") +
  coord_flip(clip = "off") + scale_y_continuous(limits = c(0, 400), breaks = c(0, 200)) +
  scale_fill_manual(values = c(Control = "#9E9E9E", Disease = masld_colors$mash), name = NULL) +
  labs(x = NULL, y = "Units") + theme_masld_compact() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
panel_a <- (p_meta | p_bar) + plot_layout(widths = c(1.35, 0.8), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))
ggsave(file.path(out, "fig3a_cohort_metadata_matrix_binned.pdf"), panel_a, width = 3.0, height = 2.25, device = cairo_pdf)
fwrite(availability[order(dataset, field)], file.path(out, "fig3a_metadata_availability_binned.tsv"), sep = "\t")
stage_tab <- samples[, .N, by = .(dataset, group_binary, fibrosis_stage)][order(dataset, group_binary, fibrosis_stage)]
fwrite(stage_tab, file.path(out, "fig3a_fibrosis_codes_by_cohort.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_fig3a.txt"))
