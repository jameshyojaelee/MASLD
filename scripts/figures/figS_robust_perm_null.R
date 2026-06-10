#!/usr/bin/env Rscript
##############################################################################
# figS_robust_perm_null.R  (new 2026-05-08)
# Pillar C condensed: within-cohort case/control label permutation null.
# KEY MESSAGE: Observed DEG count (1,885) is far above any permuted DEG count
#              (max = 2 across B = 100 permutations); empirical FDR < 0.001.
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
OUT_PDF <- file.path(OUT_DIR, "figS_robust_perm_null.pdf")

theme_robust <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
    axis.line   = element_line(colour = "black", linewidth = 0.3),
    axis.ticks  = element_line(colour = "black", linewidth = 0.3),
    plot.title    = element_text(face = "bold", size = 11,
                                 margin = margin(b = 2)),
    plot.subtitle = element_text(size = 9, colour = "grey30",
                                 margin = margin(b = 6)),
    plot.margin = margin(8, 10, 8, 10)
  )

count_dt <- fread(file.path(AUDIT, "pillar_C_perm_count_dist.csv"))
sum_dt   <- fread(file.path(AUDIT, "pillar_C_summary.csv"))
n_obs   <- sum_dt$n_deg_observed[1]
emp_fdr <- sum_dt$empirical_fdr[1]
n_perm  <- nrow(count_dt)
perm_max <- max(count_dt$n_deg_perm)
perm_q95 <- sum_dt$perm_count_q95[1]

# Plot the permutation histogram and put the observed value as an inset arrow
# on the right (since 1885 vs <= 2 would compress the histogram to nothing).
p <- ggplot(count_dt, aes(x = n_deg_perm)) +
  geom_histogram(fill = "grey80", colour = "grey45",
                 linewidth = 0.25, bins = 4) +
  geom_vline(xintercept = perm_q95, linetype = "dashed",
             colour = "grey25", linewidth = 0.4) +
  annotate("text", x = perm_q95, y = Inf, hjust = -0.15, vjust = 1.6,
           label = sprintf(" 95th percentile = %d", perm_q95),
           size = 3, colour = "grey25") +
  annotate("text", x = perm_max * 0.5, y = Inf, hjust = 0.5, vjust = 4.0,
           label = sprintf("Observed = %s\nEmpirical FDR < %.3f\n(observed >> perm null)",
                           comma(n_obs), max(emp_fdr, 0.001)),
           size = 3.6, colour = "#D55E00", fontface = "bold",
           lineheight = 1.05) +
  scale_x_continuous(breaks = pretty_breaks(n = 4),
                     limits = c(-0.5, perm_max + 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(title = "Within-cohort permutation null",
       subtitle = sprintf(paste0("B = %d label permutations preserving ",
                                  "per-cohort case/control counts"), n_perm),
       x = "Permuted DEG count (padj < 0.05, |LFC| > 0.3)",
       y = "Permutations") +
  theme_robust

ggsave(OUT_PDF, p, width = 7, height = 4, device = cairo_pdf)
cat("Saved:", OUT_PDF, "\n")
