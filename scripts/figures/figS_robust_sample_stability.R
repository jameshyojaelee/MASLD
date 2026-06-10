#!/usr/bin/env Rscript
##############################################################################
# figS_robust_sample_stability.R  (new 2026-05-08)
# Pillar A condensed: CPSS pi-hat (left) + bootstrap selection frequency (right).
# KEY MESSAGE: 98% of mega DEGs hit CPSS pi-hat >= 0.7; bootstrap frequency
#              for canonical DEGs is sharply right-shifted vs non-DEGs.
##############################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
AUDIT   <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
OUT_DIR <- file.path(FIG_SUPP, "figS_methods_validation/robustness")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(OUT_DIR, "figS_robust_sample_stability.pdf")

theme_robust <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
    axis.line     = element_line(colour = "black", linewidth = 0.3),
    axis.ticks    = element_line(colour = "black", linewidth = 0.3),
    legend.position = "top",
    legend.title  = element_blank(),
    legend.margin = margin(b = -3),
    plot.title    = element_text(face = "bold", size = 11,
                                 margin = margin(b = 2)),
    plot.subtitle = element_text(size = 9, colour = "grey30",
                                 margin = margin(b = 6)),
    plot.tag      = element_text(face = "bold", size = 12),
    plot.margin   = margin(8, 10, 8, 10)
  )

PAL <- c("Canonical DEG" = "#D55E00", "Non-DEG" = "#999999")

stab <- fread(file.path(AUDIT, "pillar_A_stability.csv"))
stab[, deg_class := fifelse(is_canonical_DEG, "Canonical DEG", "Non-DEG")]
stab[, deg_class := factor(deg_class, levels = c("Non-DEG", "Canonical DEG"))]

pct_pi <- round(100 * mean(stab[is_canonical_DEG == TRUE, cpss_pi_hat] >= 0.7,
                           na.rm = TRUE), 1)
pct_boot <- round(100 * mean(stab[is_canonical_DEG == TRUE, bootstrap_freq] >= 0.9,
                             na.rm = TRUE), 1)

# Panel A: CPSS pi-hat density
p_a <- ggplot(stab[!is.na(cpss_pi_hat)],
              aes(x = cpss_pi_hat, fill = deg_class, colour = deg_class)) +
  geom_density(alpha = 0.5, linewidth = 0.4) +
  geom_vline(xintercept = 0.7, linetype = "dashed",
             colour = "grey25", linewidth = 0.4) +
  annotate("text", x = 0.71, y = Inf, label = " stability\n cutoff = 0.7",
           hjust = 0, vjust = 1.4, size = 3, colour = "grey25") +
  annotate("text", x = 0.99, y = Inf, hjust = 1, vjust = 3,
           label = sprintf("%.1f%% of mega DEGs\npass cutoff", pct_pi),
           size = 3.2, colour = "#D55E00", fontface = "bold") +
  scale_fill_manual(values = PAL) +
  scale_colour_manual(values = PAL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(tag = "A",
       title = "CPSS selection probability",
       subtitle = "100 paired half-samples; pi-hat = max selection over pairs",
       x = expression(hat(pi)),
       y = "Density") +
  theme_robust

# Panel B: Bootstrap selection-frequency ECDF
p_b <- ggplot(stab[!is.na(bootstrap_freq)],
              aes(x = bootstrap_freq, colour = deg_class)) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  geom_vline(xintercept = 0.9, linetype = "dashed",
             colour = "grey25", linewidth = 0.4) +
  annotate("text", x = 0.91, y = 0.05, label = " freq >= 0.9",
           hjust = 0, vjust = 0, size = 3, colour = "grey25") +
  annotate("text", x = 0.05, y = 0.55, hjust = 0, vjust = 1,
           label = sprintf("%.1f%% of mega DEGs\nselected in >= 90%%\nof bootstraps",
                           pct_boot),
           size = 3.2, colour = "#D55E00", fontface = "bold") +
  scale_colour_manual(values = PAL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(tag = "B",
       title = "Bootstrap selection frequency",
       subtitle = "B = 1000 cohort-stratified subsamples (80% of cases + controls)",
       x = "P(selected as DEG)",
       y = "ECDF") +
  theme_robust

composite <- p_a + p_b + plot_layout(ncol = 2)
ggsave(OUT_PDF, composite, width = 9, height = 4, device = cairo_pdf)
cat("Saved:", OUT_PDF, "\n")
