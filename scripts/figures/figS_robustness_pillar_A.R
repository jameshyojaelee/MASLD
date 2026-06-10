#!/usr/bin/env Rscript
# figS_robustness_pillar_A.R
# Pillar A — Sample-level stability (CPSS + bootstrap + k-fold).
# Panels:
#   A1: CPSS π̂ histogram, canonical DEG vs non-DEG overlay (target: bimodal)
#   A2: Bootstrap selection-frequency ECDF, stratified by canonical-DEG status
#   A3: K-fold recurrence stacked bar (canonical DEGs only) for K=10 and K=50
#   A4: Per-iter ρ (Spearman) box plot, by source

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

AUDIT  <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
OUT_DIR <- file.path(BASE, "figures/supplementary/figS_methods_validation/robustness")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stab <- fread(file.path(AUDIT, "pillar_A_stability.csv"))
iter <- fread(file.path(AUDIT, "pillar_A_iter_summary.csv"))

stab[, deg_class := fifelse(is_canonical_DEG, "Canonical DEG", "Non-DEG")]

# A1 — CPSS pi-hat histogram
p1 <- ggplot(stab[!is.na(cpss_pi_hat)],
             aes(x = cpss_pi_hat, fill = deg_class)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 40) +
  geom_vline(xintercept = 0.7, linetype = "dashed", colour = "grey30") +
  scale_y_continuous(labels = comma) +
  labs(title = "A1 — CPSS selection probability",
       x = expression(hat(pi)~" (max selection prob over 100 paired half-samples)"),
       y = "Genes", fill = "") +
  theme_pub()

# A2 — Bootstrap ECDF
p2 <- ggplot(stab[!is.na(bootstrap_freq)],
             aes(x = bootstrap_freq, colour = deg_class)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  geom_vline(xintercept = 0.9, linetype = "dashed", colour = "grey30") +
  labs(title = "A2 — Bootstrap selection frequency (B = 1000)",
       x = "P(selected) over 1000 cohort-stratified subsamples",
       y = "ECDF", colour = "") + theme_pub()

# A3 — K-fold recurrence (canonical DEGs only)
deg <- stab[is_canonical_DEG == TRUE]
kf_long <- rbind(
  deg[, .(K = "K = 10", recur = kfold10_recur)],
  deg[, .(K = "K = 50", recur = kfold50_recur)])
kf_long[, recur_bin := cut(recur, breaks = c(-1, 0, 5, 10, 20, 50),
                           labels = c("0", "1-5", "6-10", "11-20", "21-50"))]
p3 <- ggplot(kf_long, aes(x = recur_bin, fill = K)) +
  geom_bar(position = "dodge") +
  labs(title = "A3 — K-fold recurrence of canonical DEGs",
       x = "Folds where DEG was re-selected", y = "DEGs", fill = "") +
  theme_pub()

# A4 — per-iter Spearman rho across sources
p4 <- ggplot(iter[!is.na(rho)], aes(x = source, y = rho)) +
  geom_boxplot(outlier.size = 0.4, fill = "grey90") +
  geom_hline(yintercept = 0.85, linetype = "dashed", colour = "grey30") +
  coord_flip() +
  labs(title = "A4 — Per-iter logFC concordance",
       x = "", y = "Spearman ρ vs full dream") + theme_pub()

combined <- (p1 | p2) / (p3 | p4) + plot_layout(heights = c(1, 1))
out_pdf <- file.path(OUT_DIR, "figS_robustness_pillar_A.pdf")
ggsave(out_pdf, combined, width = 12, height = 9, device = cairo_pdf)
cat("Saved:", out_pdf, "\n")
