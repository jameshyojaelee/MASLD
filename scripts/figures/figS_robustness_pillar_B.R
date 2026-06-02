#!/usr/bin/env Rscript
# figS_robustness_pillar_B.R
# Pillar B — Cross-cohort biological transferability (LOCO disease prediction).
# Panels:
#   B1: AUROC per held-out cohort × method, with random-label null mean ± SD
#   B2: AUROC vs train n (transferability NOT just driven by training-set size)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

AUDIT <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
OUT_DIR <- file.path(BASE, "figures/supplementary/robustness")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dt <- fread(file.path(AUDIT, "pillar_B_loco_prediction.csv"))
dt[, method_lbl := factor(method, levels = c("dotprod", "signed_ssgsea"),
                          labels = c("Dot-product", "Signed ssGSEA"))]

# B1 — AUROC per cohort + null
p1 <- ggplot(dt, aes(x = held_out, y = auroc, fill = method_lbl)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  geom_errorbar(aes(ymin = null_mean - null_sd, ymax = null_mean + null_sd,
                    group = method_lbl),
                position = position_dodge(width = 0.85), width = 0.3,
                colour = "grey30") +
  geom_point(aes(y = null_mean, group = method_lbl),
             position = position_dodge(width = 0.85), shape = 4, size = 2,
             colour = "grey30") +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey50") +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey30") +
  scale_y_continuous(limits = c(0.4, 1.02), breaks = seq(0.4, 1, 0.1)) +
  labs(title = "B1 — Cross-cohort disease prediction",
       subtitle = "Train: dream LOO on 4 cohorts; Test: held-out 5th cohort. ✕ = random-label null mean (B = 100)",
       x = "Held-out cohort", y = "AUROC", fill = "") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

# B2 — AUROC vs test n (sanity: not just sample-size driven)
p2 <- ggplot(dt, aes(x = n_test, y = auroc, colour = method_lbl, label = held_out)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey30") +
  scale_x_continuous(labels = comma) +
  labs(title = "B2 — AUROC vs held-out cohort n",
       x = "Held-out cohort sample size", y = "AUROC", colour = "") +
  theme_pub()

combined <- p1 / p2 + plot_layout(heights = c(1.2, 1))
out_pdf <- file.path(OUT_DIR, "figS_robustness_pillar_B.pdf")
ggsave(out_pdf, combined, width = 9, height = 8, device = cairo_pdf)
cat("Saved:", out_pdf, "\n")
