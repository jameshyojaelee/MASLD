# KEY MESSAGE (panel 1): the F2-F4 versus F0 counts change little with the F0 reference; the F1 count falls from 1,338 to 783-814.
# KEY MESSAGE (panel 2): the F1-versus-F0 count depends on GSE130970 (14 of 1,338 genes without it).
# F3ef repair round. Adapted from ws/F3ef/render_s3_stage_sensitivity.R.
# Changes: panel 1 prints every count and gives n for all three F0 references;
# panel 2 reads the repair refit (job 7550) and drops the library-complexity bars,
# which belong to the H82 panel (figS3_h82_complexity_retention.pdf).
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
set.seed(20260923)
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
source(file.path(base, "scripts/figures/publication_theme.R"))
grDevices::pdf.options(useDingbats = FALSE)
out <- file.path(base, "figures/candidates/review-r3-fixes-20261015T112148Z/F3ef/repair")
refit <- file.path(base, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F3ef/refit-20261015T171823Z")

arms <- fread(file.path(base, "figures/candidates/f0-reference-arms-r8-20261008T193309Z/analysis/f0_arms/f0_arm_summary.tsv"))
arms <- arms[contrast %in% c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0")]
arms[, stage := factor(sub("_vs_F0", "", contrast), levels = c("F1", "F2", "F3", "F4"))]
arms[, arm_label := factor(arm, levels = c("A_all_F0", "B_disease_only_F0", "C_all_F0_plus_term"),
                           labels = c("All F0 (n=103)", "Disease-only F0 (n=89)", "All F0 + control-status term (n=103)"))]
stopifnot(arms[arm == "A_all_F0", n_BH_sig] == c(1338, 4266, 6190, 7043),
          arms[arm == "B_disease_only_F0", n_BH_sig] == c(783, 3946, 6459, 6929),
          arms[arm == "C_all_F0_plus_term", n_BH_sig] == c(814, 3679, 6166, 6640),
          arms[arm == "C_all_F0_plus_term", unique(n_reference)] == 103,
          arms[arm == "B_disease_only_F0", unique(n_reference)] == 89)
arm_cols <- c("All F0 (n=103)" = "#252525", "Disease-only F0 (n=89)" = "#C9265E",
              "All F0 + control-status term (n=103)" = "#4aa2c2")
# One printed table of counts under the points, so every legend number is readable on the panel.
lab <- dcast(arms, arm_label ~ stage, value.var = "n_BH_sig")
tab <- melt(lab, id.vars = "arm_label", variable.name = "stage", value.name = "n")
tab[, y := c(-900, -1500, -2100)[as.integer(arm_label)]]
p1 <- ggplot(arms, aes(stage, n_BH_sig, color = arm_label, shape = arm_label, group = arm_label)) +
  geom_line(linewidth = 0.3) + geom_point(size = 1.2) +
  geom_text(data = tab, aes(stage, y, label = scales::comma(n), color = NULL, shape = NULL, group = NULL),
            size = 6 / .pt, color = "black", inherit.aes = FALSE) +
  geom_point(data = tab, aes(x = 0.45, y = y, color = arm_label, shape = arm_label), size = 1.0, inherit.aes = FALSE) +
  scale_color_manual(values = arm_cols, name = NULL) + scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  scale_y_continuous(limits = c(-2400, 7500), breaks = seq(0, 7500, 2500), labels = scales::comma) +
  coord_cartesian(xlim = c(0.4, 4.4), clip = "off") +
  guides(y = guide_axis(cap = "both")) +
  labs(x = "Fibrosis stage versus F0", y = "Genes, BH q<0.05") + theme_masld_compact() +
  theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(out, "figs3_f0_reference_arms.pdf"), p1, width = 2.2, height = 2.7, device = cairo_pdf)
fwrite(arms[, .(arm, contrast, n_reference, n_comparison, n_ref_control, n_ref_disease, covariates, n_BH_sig, n_up, n_down)],
       file.path(out, "figs3_f0_reference_arms.tsv"), sep = "\t")

fc <- fread(file.path(refit, "fit_counts.tsv"))[contrast == "F0_to_F1"]
fc[, family := fifelse(fit == "base", "Released model",
               fifelse(grepl("^drop_", fit), "Leave one cohort out", "One cohort alone"))]
lab_map <- c(base = "Three cohorts (released)", drop_GSE130970 = "Without GSE130970",
             drop_GSE135251 = "Without GSE135251", drop_GSE162694 = "Without GSE162694",
             alone_GSE130970 = "GSE130970 alone", alone_GSE135251 = "GSE135251 alone", alone_GSE162694 = "GSE162694 alone")
fc[, label := sprintf("%s (n=%d)", lab_map[fit], n_units)]
fc[, family := factor(family, levels = c("Released model", "Leave one cohort out", "One cohort alone"))]
setorder(fc, family, -n_bh05); fc[, label := factor(label, levels = rev(label))]
stopifnot(fc[fit == "base", n_bh05] == 1338, nrow(fc) == 7L)
p2 <- ggplot(fc, aes(n_bh05, label)) +
  geom_col(width = 0.7, fill = "#C9265E") +
  geom_text(aes(label = scales::comma(n_bh05)), hjust = -0.15, size = 6 / .pt, color = "black") +
  facet_grid(family ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1600), breaks = c(0, 500, 1000, 1500), labels = scales::comma) +
  labs(x = "F1 versus F0 genes, BH q<0.05", y = NULL) + theme_masld_compact() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 6, face = "plain"),
        strip.background = element_blank())
ggsave(file.path(out, "figs3_f1_vs_f0_cohort_omission.pdf"), p2, width = 3.4, height = 1.6, device = cairo_pdf)
fwrite(fc, file.path(out, "figs3_f1_vs_f0_cohort_omission.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_s3_repair.txt"))
