#!/usr/bin/env Rscript
# figS_robustness_pillar_C.R
# Pillar C — Empirical null calibration via within-cohort permutation.
# Panels:
#   C1: Histogram of permuted DEG counts (B=1000) with observed count overlay
#   C2: Per-gene |perm_z| ECDF, canonical DEG vs non-DEG
#   C3: Volcano: observed dream t vs perm_emp_p (canonical DEGs flagged)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

AUDIT <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
OUT_DIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/robustness")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

count_dt <- fread(file.path(AUDIT, "pillar_C_perm_count_dist.csv"))
gene_dt  <- fread(file.path(AUDIT, "pillar_C_empirical_null.csv"))
sum_dt   <- fread(file.path(AUDIT, "pillar_C_summary.csv"))
atlas    <- fread(file.path(AUDIT, "dream_robustness_atlas.csv"))

n_obs <- sum_dt$n_deg_observed[1]
emp_fdr <- sum_dt$empirical_fdr[1]

# C1 — permuted DEG count histogram
p1 <- ggplot(count_dt, aes(x = n_deg_perm)) +
  geom_histogram(fill = "grey80", colour = "grey30", bins = 40) +
  geom_vline(xintercept = n_obs, colour = "red", linewidth = 1) +
  annotate("text", x = n_obs * 0.95, y = Inf, vjust = 1.5, hjust = 1,
           label = sprintf("Observed = %s", comma(n_obs)),
           colour = "red") +
  annotate("text", x = max(count_dt$n_deg_perm), y = Inf, vjust = 3, hjust = 1,
           label = sprintf("Empirical FDR = %.4f", emp_fdr)) +
  scale_x_continuous(labels = comma) +
  labs(title = "C1 — DEG count under within-cohort permutation",
       subtitle = sprintf("B = %d permutations; preserves per-cohort case/control counts",
                          nrow(count_dt)),
       x = "Permuted DEG count (padj < 0.05 & |LFC| > 0.3)", y = "Permutations") +
  theme_pub()

# C2 — |perm z| ECDF
gene_dt <- merge(gene_dt,
                 atlas[, .(gene, is_canonical_DEG)], by = "gene", all.x = TRUE)
gene_dt[, deg_class := fifelse(is_canonical_DEG, "Canonical DEG", "Non-DEG")]
p2 <- ggplot(gene_dt[!is.na(perm_z)], aes(x = abs(perm_z), colour = deg_class)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  geom_vline(xintercept = 3, linetype = "dashed", colour = "grey30") +
  scale_x_continuous(limits = c(0, quantile(abs(gene_dt$perm_z), 0.995, na.rm = TRUE)),
                     oob = scales::squish) +
  labs(title = "C2 — |permutation z| distribution",
       x = "|z| = (t_observed − μ_perm) / σ_perm", y = "ECDF",
       colour = "") + theme_pub()

# C3 — observed t vs empirical p (volcano-ish)
p3 <- ggplot(gene_dt[!is.na(perm_emp_p)], aes(x = t_obs, y = -log10(perm_emp_p + 1e-5),
                                              colour = deg_class)) +
  geom_point(alpha = 0.3, size = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey30") +
  labs(title = "C3 — observed t vs empirical permutation p",
       x = "dream t-statistic (observed)", y = "-log10(empirical p)",
       colour = "") + theme_pub()

combined <- (p1 / (p2 | p3)) + plot_layout(heights = c(1, 1))
out_pdf <- file.path(OUT_DIR, "figS_robustness_pillar_C.pdf")
ggsave(out_pdf, combined, width = 12, height = 9, device = cairo_pdf)
cat("Saved:", out_pdf, "\n")
