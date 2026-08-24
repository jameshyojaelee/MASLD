#!/usr/bin/env Rscript
# Render Figure 3C from synchronized per-cohort and leave-one-cohort-out summaries.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
grDevices::pdf.options(useDingbats = FALSE)

base <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
candidate_root <- normalizePath(Sys.getenv("FIGURE_CANDIDATE_ROOT", ""), mustWork = TRUE)
if (!exists("theme_masld") || !exists("cat_palette")) {
  source(file.path(base, "scripts/figures/publication_theme.R"))
}
source_path <- file.path(candidate_root, "source_tables", "figs3c_cohort_robustness.tsv")
out_path <- file.path(candidate_root, "supplementary", "figureS3", "panels", "figs3c_cohort_robustness.pdf")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
fig3c_source <- fread(source_path)
stopifnot(nrow(fig3c_source) == 5L, all(fig3c_source$canonical_genes == 1347L))

cohort_order <- fig3c_source[order(-(n_control + n_disease)), dataset]
fig3c_source[, dataset_factor := factor(dataset, levels = rev(cohort_order))]
p_direction <- ggplot(fig3c_source,
                      aes(100 * same_direction_fraction, dataset_factor, color = dataset)) +
  geom_segment(aes(x = 50, xend = 100 * same_direction_fraction,
                   yend = dataset_factor), color = "grey82", linewidth = 0.45) +
  geom_point(size = 1.4) +
  scale_x_continuous(limits = c(50, 100), breaks = c(50, 75, 100)) +
  scale_color_manual(values = setNames(cat_palette[seq_along(cohort_order)], cohort_order),
                     guide = "none") +
  labs(x = "Same direction (%)", y = NULL) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(axis.text = element_text(size = 6, face = "plain"), panel.grid.major.y = element_blank())
p_loo <- ggplot(fig3c_source,
                aes(all_gene_spearman, dataset_factor, color = dataset)) +
  geom_segment(aes(x = 0.8, xend = all_gene_spearman, yend = dataset_factor),
               color = "grey82", linewidth = 0.45) +
  geom_point(size = 1.4) +
  scale_x_continuous(limits = c(0.8, 1.0), breaks = c(0.8, 0.9, 1.0)) +
  scale_color_manual(values = setNames(cat_palette[seq_along(cohort_order)], cohort_order),
                     guide = "none") +
  labs(x = "LOO effect ρ", y = NULL) +
  theme_masld(base_size = 6) + theme_pub() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 6, face = "plain"), panel.grid.major.y = element_blank())
fig3c_panel <- p_direction | p_loo
ggsave(out_path, fig3c_panel, width = 2.65, height = 2.15, device = cairo_pdf)
cat("FIGS3C_RENDER_COMPLETE", out_path, "\n")
