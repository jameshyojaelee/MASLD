#!/usr/bin/env Rscript
# figS_multimethod_power.R  — simulation-based power for the cohort-resample
# multimethod_validation harness. Plots WHATEVER DE engines are present in
# power_summary.csv (no hard-coded method filter); readable engine names +
# colours come from method_correction_labels.R. NB simulation seeded from real
# data; K=5 cohorts; own-universe metrics. metafor uses the REML no-HKSJ primary
# (so it is powered, not pinned at the HKSJ ~0 floor).
# Output (distinct from the degx-battery power_curves.pdf so both can coexist):
#   simulation_power_by_method.pdf  (was panelB_power_curves.pdf)
#   simulation_fdr_by_method.pdf    (was panelB2_fdr_control.pdf)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/method_correction_labels.R"))
P   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/multimethod_validation/power")
OUT <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

s <- fread(file.path(P, "power_summary.csv"))
# Plot whatever methods are present; order by engine family then readable name.
present <- unique(as.character(s$method))
present <- present[order(match(engine_family(present), ENGINE_FAMILY_LEVELS),
                         engine_label(present))]
mcol   <- setNames(ENGINE_FAMILY_PALETTE[engine_family(present)], present)
pretty <- setNames(engine_label(present), present)
s[, method := factor(method, levels = present)]
s[, tau_lab := factor(paste0("tau2=", tau2), levels = paste0("tau2=", sort(unique(tau2))))]
s[, lfc_lab := factor(paste0("LFC=", true_lfc), levels = paste0("LFC=", sort(unique(true_lfc))))]

# ---- power curves ----
p1 <- ggplot(s, aes(n_per_group, power_own_mean, colour = method)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.1) +
  facet_grid(tau_lab ~ lfc_lab) +
  scale_colour_manual(values = mcol, labels = pretty, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = c(10, 25, 50, 100)) +
  labs(x = "samples per group (per cohort)", y = "power (TPR @ FDR<0.05)",
       title = "Simulation power by DE method",
       subtitle = "K=5 cohorts; rows = between-cohort heterogeneity tau2, cols = true effect size; metafor uses REML no-HKSJ primary") +
  theme_masld(base_size = 7) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 5.8, colour = "grey40"),
        strip.text = element_text(size = 6), panel.spacing = unit(0.25, "lines"))
ggsave(file.path(OUT, "simulation_power_by_method.pdf"), p1, width = 6.6, height = 4.6, useDingbats = FALSE)
cat("Wrote simulation_power_by_method.pdf\n")

# ---- FDR control ----
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
       title = sprintf("FDR control across %d grid cells", nrow(s) / max(1L, length(present))),
       subtitle = "nominal 0.05 (dashed); boxes summarise observed FDR per simulation grid cell") +
  theme_masld(base_size = 7) +
  theme(plot.subtitle = element_text(size = 6, colour = "grey40"),
        axis.text.x = element_text(size = 6.5))
ggsave(file.path(OUT, "simulation_fdr_by_method.pdf"), p2, width = 4.4, height = 3.4, useDingbats = FALSE)
cat("Wrote simulation_fdr_by_method.pdf\n")
print(nfail)
cat("\nmean power/auroc by method:\n")
print(s[, .(power=round(mean(power_own_mean),3), auroc=round(mean(auroc_own_mean),3)), by=method])
