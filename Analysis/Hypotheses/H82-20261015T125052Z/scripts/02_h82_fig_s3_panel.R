#!/usr/bin/env Rscript
# KEY MESSAGE: For each Figure 3E/S3 stage contrast, the share of released DEG calls kept under
# library-complexity adjustment, set against the share kept when the same covariates are permuted within cohort.
# H82 Supplementary Figure S3 panel, drawn only from the output tables of 01_h82_refits.R.
# Unit: bulk analysis unit (genetic-identity cluster); the retained fraction summarizes a unit-level refit over genes.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- normalizePath(args[[match("--out", args) + 1L]])
v <- jsonlite::fromJSON(file.path(OUT, "verdict.json"))
if (!v$status %in% c("completed", "untestable_S5")) {
  message("H82 panel not drawn: run status ", v$status); quit(save = "no", status = 0L)
}
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))
FS <- 6
d1 <- fread(file.path(OUT, "d1_arm_table.tsv"))[arm %in% c("lib", "gini", "nf", "joint")]
ns <- fread(file.path(OUT, "null_summary.tsv"))
lev <- c("F1_vs_F0", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0", "F1_to_F2", "F2_to_F3", "F3_to_F4",
         "NAS1_2_vs_NAS0", "NAS3_4_vs_NAS0", "NAS5_8_vs_NAS0", "NAS5_8_vs_NAS3_4", "NAS3_4_vs_NAS1_2")
lab <- c("F1 vs F0", "F2 vs F0", "F3 vs F0", "F4 vs F0", "F1 to F2", "F2 to F3", "F3 to F4",
         "NAS1-2 vs NAS0", "NAS3-4 vs NAS0", "NAS5-8 vs NAS0", "NAS5-8 vs NAS3-4", "NAS3-4 vs NAS1-2")
stopifnot(setequal(d1$contrast, lev), all(ns$contrast %in% lev))
arm_lev <- c("lib", "gini", "nf", "joint")
arm_lab <- c(lib = "log2 library size", gini = "Gini coefficient", nf = "stored norm factor", joint = "all three (joint)")
off <- c(lib = -0.24, gini = -0.08, nf = 0.08, joint = 0.24)
d1[, x := match(contrast, lev) + off[arm]]
d1[, arm := factor(arm, levels = arm_lev)]
ns[, x := match(contrast, lev)]
cols <- c(lib = "#40b499", gini = "#518dc9", nf = "#F4A674", joint = "#C9265E")
shp <- c(lib = 17, gini = 15, nf = 18, joint = 16)
p <- ggplot() +
  geom_vline(xintercept = 7.5, linewidth = 0.25, linetype = "dashed", colour = "#9E9E9E") +
  geom_rect(data = ns, aes(xmin = x - 0.36, xmax = x + 0.36, ymin = q025_R, ymax = q975_R, fill = "permuted-covariate null, 2.5-97.5%"),
            colour = NA, alpha = 0.55) +
  geom_segment(data = ns, aes(x = x - 0.36, xend = x + 0.36, y = median_R, yend = median_R), colour = "#9E9E9E", linewidth = 0.5) +
  geom_point(data = d1, aes(x = x, y = retained_frac, colour = arm, shape = arm), size = 1.3, stroke = 0.3) +
  scale_fill_manual(values = c("permuted-covariate null, 2.5-97.5%" = "#9E9E9E"), name = NULL) +
  scale_colour_manual(values = cols, labels = arm_lab, name = "covariates added") +
  scale_shape_manual(values = shp, labels = arm_lab, name = "covariates added") +
  scale_x_continuous(breaks = seq_along(lev), labels = lab, expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "retained fraction of released BH calls") +
  theme_classic(base_size = FS) +
  theme(text = element_text(size = FS, colour = "black", face = "plain"),
        axis.text = element_text(size = FS, colour = "black"), axis.title = element_text(size = FS, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = FS, colour = "black"), legend.title = element_text(size = FS, colour = "black"),
        legend.key.size = unit(0.25, "cm"), legend.position = "right",
        axis.line = element_line(linewidth = 0.3, colour = "black"), axis.ticks = element_line(linewidth = 0.3, colour = "black"))
FIG <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE)
f <- file.path(FIG, "figS3_h82_complexity_retention.pdf")
if (file.exists(f)) stop("Refusing overwrite: ", f)
grDevices::pdf(f, width = 125 / 25.4, height = 62 / 25.4, useDingbats = FALSE)
print(p); invisible(grDevices::dev.off())
pr <- fread(file.path(OUT, "primary_tests.tsv"))
f1 <- pr[contrast == "F1_vs_F0"]
cap <- sprintf(paste0(
  "Retained fraction of released stage DEG calls (BH < 0.05 over 24,115 genes) when within-cohort z-scores of log2 library size, ",
  "Gini coefficient and stored normalization factor are added singly or jointly to ~ dataset + inferred sex + group (limma voomWithQualityWeights). ",
  "Gray bar: 2.5-97.5 percentile of the joint retained fraction when the three covariates are permuted together within cohort (%d permutations); ",
  "gray line: null median. F1 vs F0 joint retains %d of %d released calls (%.3f; null median %.3f; one-sided permutation P = %.4f). ",
  "Units: bulk analysis units (fibrosis roster 387, NAS roster 395; three cohorts). NAS3-4 vs NAS1-2 has 3 released calls and no null. ",
  "Adjusted counts are a sensitivity, not a replacement, because Gini and norm factor can partly reflect disease-related RNA composition."),
  as.integer(fread(file.path(OUT, "null_summary.tsv"))$B[1]), f1$joint_retained, f1$released, f1$R, f1$null_median, f1$p)
message("CAPTION: ", cap)
cat(c("", "## Figure", "", paste0("Panel: ", f), "", paste0("Caption (draft): ", cap)), file = file.path(OUT, "RESULT_H82.md"), sep = "\n", append = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_fig.txt"))
