#!/usr/bin/env Rscript
# figS_multimethod_power.R  — Panel B (simulation-based power) for the 3-method comparison.
#   panelB_power_curves.pdf  - power (TPR@FDR<0.05) vs n/group, facet tau2 x true_lfc, line per method
#   panelB2_fdr_control.pdf  - observed FDR per method across all 48 cells vs nominal 0.05
# NB simulation seeded from real data; K=5 cohorts; own-universe metrics.
# metafor uses the REML no-HKSJ primary (so it is powered, not the HKSJ ~0 floor).
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
P   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/multimethod_validation/power")
OUT <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

s <- fread(file.path(P, "power_summary.csv"))
mcol   <- c(dream = "#1b9e77", deseq2 = "#66a61e", metafor = "#1f78b4")
pretty <- c(dream = "dream", deseq2 = "DESeq2", metafor = "metafor")
s[, method := factor(method, levels = c("dream","deseq2","metafor"))]
s[, tau_lab := factor(paste0("tau2=", tau2), levels = paste0("tau2=", sort(unique(tau2))))]
s[, lfc_lab := factor(paste0("LFC=", true_lfc), levels = paste0("LFC=", sort(unique(true_lfc))))]

# ---- Panel B1: power curves ----
p1 <- ggplot(s, aes(n_per_group, power_own_mean, colour = method)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.1) +
  facet_grid(tau_lab ~ lfc_lab) +
  scale_colour_manual(values = mcol, labels = pretty, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = c(10, 25, 50, 100)) +
  labs(x = "samples per group (per cohort)", y = "power (TPR @ FDR<0.05)",
       title = "Simulation power: dream >= DESeq2 > metafor",
       subtitle = "K=5 cohorts; rows = between-cohort heterogeneity tau2, cols = true effect size; metafor uses REML no-HKSJ primary") +
  theme_masld(base_size = 7) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 5.8, colour = "grey40"),
        strip.text = element_text(size = 6), panel.spacing = unit(0.25, "lines"))
ggsave(file.path(OUT, "panelB_power_curves.pdf"), p1, width = 6.6, height = 4.6, useDingbats = FALSE)
cat("Wrote panelB_power_curves.pdf\n")

# ---- Panel B2: FDR control ----
nfail <- s[, .(n_gt = sum(fdr_own_mean > 0.05), n = .N,
               mean_fdr = round(mean(fdr_own_mean), 3)), by = method]
lab <- nfail[, sprintf("%s\n%d/%d cells >0.05", pretty[as.character(method)], n_gt, n)]
p2 <- ggplot(s, aes(method, fdr_own_mean, fill = method)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#D6604D", linewidth = 0.4) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.5, colour = "grey25") +
  scale_fill_manual(values = mcol, guide = "none") +
  scale_x_discrete(labels = setNames(lab, nfail$method)) +
  labs(x = NULL, y = "observed FDR (per grid cell)",
       title = "FDR control across 48 grid cells",
       subtitle = "nominal 0.05 (dashed); DESeq2 anti-conservative most often; metafor best controlled") +
  theme_masld(base_size = 7) +
  theme(plot.subtitle = element_text(size = 6, colour = "grey40"),
        axis.text.x = element_text(size = 6.5))
ggsave(file.path(OUT, "panelB2_fdr_control.pdf"), p2, width = 4.4, height = 3.4, useDingbats = FALSE)
cat("Wrote panelB2_fdr_control.pdf\n")
print(nfail)
cat("\nmean power/auroc by method:\n")
print(s[, .(power=round(mean(power_own_mean),3), auroc=round(mean(auroc_own_mean),3)), by=method])
