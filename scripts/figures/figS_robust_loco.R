#!/usr/bin/env Rscript
##############################################################################
# figS_robust_loco.R  (rewritten 2026-05-08; simplified per PI request)
# Pillar B condensed: leave-one-cohort-out (LOCO) cross-validation.
#
# KEY MESSAGE: When dream is trained on 4 cohorts and the resulting signature
#              is used to predict disease vs control on the held-out 5th
#              cohort, AUROC stays well above the random-label null and the
#              0.8 "good performance" benchmark for every cohort.
#
# Single panel: per-cohort AUROC bar (dot-product method; ssGSEA shown in
# companion CSV but omitted from the figure for clarity), with null markers
# and 0.5 / 0.8 reference lines.
##############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
AUDIT   <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
OUT_DIR <- file.path(FIG_SUPP, "figS_methods_validation/robustness")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(OUT_DIR, "figS_robust_loco.pdf")

# Theme matched to other figS_robust_* figures.
theme_robust <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    axis.line   = element_line(colour = "black", linewidth = 0.3),
    axis.ticks  = element_line(colour = "black", linewidth = 0.3),
    legend.position = "top",
    legend.title  = element_blank(),
    legend.margin = margin(b = -3),
    legend.text   = element_text(size = 9),
    plot.title    = element_text(face = "bold", size = 11,
                                 margin = margin(b = 2)),
    plot.subtitle = element_text(size = 9, colour = "grey30",
                                 margin = margin(b = 6)),
    plot.margin   = margin(8, 10, 8, 10)
  )

# Cohort first-author labels (match other supplementary figures)
COHORT_LABELS <- c(
  GSE126848 = "Suppli",  GSE130970 = "Hoang",
  GSE135251 = "Govaere", GSE162694 = "Bril",
  GSE213621 = "Chen"
)

dt <- fread(file.path(AUDIT, "pillar_B_loco_prediction.csv"))
# Use dot-product only (cleanest method, ssGSEA is similar)
dt <- dt[method == "dotprod"]
dt[, cohort_lbl := unname(COHORT_LABELS[held_out])]
setorder(dt, auroc)
dt[, cohort_lbl := factor(cohort_lbl, levels = cohort_lbl)]

p <- ggplot(dt, aes(y = cohort_lbl)) +
  # Dumbbell: connecting segment from null mean to real AUROC
  geom_segment(aes(x = null_mean, xend = auroc, yend = cohort_lbl),
               colour = "grey60", linewidth = 0.6) +
  # Null reference: grey diamond + horizontal SD whiskers
  geom_errorbarh(aes(xmin = null_mean - null_sd,
                     xmax = null_mean + null_sd,
                     y = cohort_lbl),
                 height = 0.18, colour = "grey45", linewidth = 0.4) +
  geom_point(aes(x = null_mean,
                 colour = "Null (shuffled labels)",
                 shape  = "Null (shuffled labels)"),
             size = 2.6, stroke = 0.6) +
  # Real AUROC: blue filled circle
  geom_point(aes(x = auroc,
                 colour = "Real DEG signature",
                 shape  = "Real DEG signature"),
             size = 3.4, stroke = 0) +
  geom_text(aes(x = auroc, label = sprintf("%.3f", auroc)),
            hjust = -0.30, size = 3, colour = "grey15") +
  scale_colour_manual(values = c("Null (shuffled labels)" = "grey45",
                                 "Real DEG signature"     = "#0072B2"),
                      breaks = c("Real DEG signature",
                                 "Null (shuffled labels)"),
                      name = NULL) +
  scale_shape_manual(values = c("Null (shuffled labels)" = 18,
                                "Real DEG signature"     = 16),
                     breaks = c("Real DEG signature",
                                "Null (shuffled labels)"),
                     name = NULL) +
  guides(
    colour = guide_legend(order = 1,
                          override.aes = list(shape = c(16, 18),
                                              size  = c(3.4, 2.6))),
    shape  = "none"
  ) +
  scale_x_continuous(breaks = seq(0.4, 1, 0.1)) +
  coord_cartesian(xlim = c(0.4, 1.05)) +
  labs(title = "Leave-one-cohort-out cross-validation",
       subtitle = paste0("Train dream on 4 cohorts, score held-out 5th. ",
                         "Whiskers = +/- SD over 100 random-label permutations."),
       x = "AUROC (held-out cohort)", y = NULL) +
  theme_robust

ggsave(OUT_PDF, p, width = 7.5, height = 4, device = cairo_pdf)
cat("Saved:", OUT_PDF, "\n")
