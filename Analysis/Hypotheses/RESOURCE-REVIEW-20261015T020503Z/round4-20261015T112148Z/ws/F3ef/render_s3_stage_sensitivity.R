# KEY MESSAGE (panel 1): the F2-F4 versus F0 counts change little with the F0 reference; the F1 count falls from 1,338 to 783-814.
# KEY MESSAGE (panel 2): the F1-versus-F0 count depends on GSE130970 and on library-complexity covariates.
# Candidate Figure S3 panels (F3ef). Reads saved fits only; no model is refitted here.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
set.seed(20260923)
base <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
source(file.path(base, "scripts/figures/publication_theme.R"))
grDevices::pdf.options(useDingbats = FALSE)
out <- file.path(base, "figures/candidates/review-r3-fixes-20261015T112148Z/F3ef")
RUN <- file.path(base, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z")

arms <- fread(file.path(base, "figures/candidates/f0-reference-arms-r8-20261008T193309Z/analysis/f0_arms/f0_arm_summary.tsv"))
arms <- arms[contrast %in% c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0")]
arms[, stage := factor(sub("_vs_F0", "", contrast), levels = c("F1", "F2", "F3", "F4"))]
arms[, arm_label := factor(arm, levels = c("A_all_F0", "B_disease_only_F0", "C_all_F0_plus_term"),
                           labels = c("All F0 (n=103)", "Disease-only F0 (n=89)", "All F0 + control-status term"))]
stopifnot(arms[arm == "A_all_F0", n_BH_sig] == c(1338, 4266, 6190, 7043))
arm_cols <- c("All F0 (n=103)" = "#252525", "Disease-only F0 (n=89)" = "#C9265E", "All F0 + control-status term" = "#4aa2c2")
p1 <- ggplot(arms, aes(stage, n_BH_sig, color = arm_label, shape = arm_label, group = arm_label)) +
  geom_line(linewidth = 0.3) + geom_point(size = 1.2) +
  scale_color_manual(values = arm_cols, name = NULL) + scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  scale_y_continuous(limits = c(0, 7500), breaks = seq(0, 7500, 2500), labels = scales::comma) +
  labs(x = "Fibrosis stage versus F0", y = "Genes, BH q<0.05") + theme_masld_compact() +
  theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(out, "figs3_f0_reference_arms.pdf"), p1, width = 2.0, height = 2.3, device = cairo_pdf)
fwrite(arms[, .(arm, contrast, n_reference, n_comparison, n_ref_control, n_ref_disease, covariates, n_BH_sig, n_up, n_down)],
       file.path(out, "figs3_f0_reference_arms.tsv"), sep = "\t")

loco <- fread(file.path(RUN, "jobs/4191/f1_vs_f0_counts.tsv"))
cx <- fread(file.path(RUN, "jobs/4211/fit_summary.tsv"))[contrast == "F0_to_F1" & label %in% c("lib", "gini", "nf", "joint")]
f1 <- rbind(
  loco[, .(fit = label, family = fifelse(label == "pooled_all_three", "Released",
                                  fifelse(grepl("^loco", label), "Leave one cohort out", "One cohort alone")),
           n_units, n_bh05 = bh_q05)],
  cx[, .(fit = paste0("complexity_", label), family = "Complexity covariate", n_units, n_bh05)])
lab_map <- c(pooled_all_three = "Three cohorts (released)", loco_drop_GSE130970 = "Without GSE130970",
             loco_drop_GSE135251 = "Without GSE135251", loco_drop_GSE162694 = "Without GSE162694",
             alone_GSE130970 = "GSE130970 alone", alone_GSE135251 = "GSE135251 alone", alone_GSE162694 = "GSE162694 alone",
             complexity_lib = "+ log2 library size", complexity_gini = "+ Gini", complexity_nf = "+ TMM norm factor",
             complexity_joint = "+ all three")
f1[, label := sprintf("%s (n=%d)", lab_map[fit], n_units)]
f1[, family := factor(family, levels = c("Released", "Leave one cohort out", "One cohort alone", "Complexity covariate"))]
setorder(f1, family, -n_bh05); f1[, label := factor(label, levels = rev(label))]
stopifnot(f1[fit == "pooled_all_three", n_bh05] == 1338, f1[fit == "loco_drop_GSE130970", n_bh05] == 14)
p2 <- ggplot(f1, aes(n_bh05, label)) +
  geom_col(width = 0.7, fill = "#C9265E") +
  geom_text(aes(label = scales::comma(n_bh05)), hjust = -0.15, size = 6 / .pt, color = "black") +
  facet_grid(family ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1600), breaks = c(0, 500, 1000, 1500), labels = scales::comma) +
  labs(x = "F1 versus F0 genes, BH q<0.05", y = NULL) + theme_masld_compact() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 6, face = "plain"),
        strip.background = element_blank())
ggsave(file.path(out, "figs3_f1_vs_f0_cohort_complexity.pdf"), p2, width = 3.4, height = 2.1, device = cairo_pdf)
fwrite(f1, file.path(out, "figs3_f1_vs_f0_cohort_complexity.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_s3.txt"))
