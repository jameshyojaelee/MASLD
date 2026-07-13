#!/usr/bin/env Rscript
##############################################################################
# figS_robustness.R  (new 2026-05-08)
# Consolidated robustness supplement -- replaces 5 separate per-pillar PDFs
# (figS_robustness_pillar_{A,B,C,D}.pdf + figS_robustness_atlas.pdf) with a
# single 6-panel composite that pulls the strongest panel from each pillar.
#
# KEY MESSAGE: The 1,885 canonical mega DEGs survive 4 orthogonal stress tests
# -- sample-level resampling (A), cross-cohort biological transfer (B),
# empirical null calibration (C), and batch-RE handling (D, E) -- and 92% pass
# all three per-gene pillars jointly (F).
#
# Panels:
#   A: CPSS pi-hat distribution (Pillar A1)        -- sample-level stability
#   B: LOCO disease-prediction AUROC (Pillar B1)   -- biological transfer
#   C: Within-cohort permutation null (Pillar C1)  -- DEG count vs null
#   D: Variance partition raw vs dream (Pillar D1) -- batch handled by RE
#   E: dream vs dream+SVA logFC (Pillar D2)        -- SVA adds ~nothing
#   F: Joint pillar-pass distribution (Atlas E1)   -- per-gene robustness summary
#
# Output:
#   figures/supplementary/figS_methods_validation/robustness/figS_robustness.pdf
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

AUDIT   <- file.path(BASE, "RNA-seq/results/audit_sensitivity")
OUT_DIR <- FIGS_ROBUST_DIR  # consolidated under figS_methods_validation/ (2026-06-04)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(OUT_DIR, "figS_robustness.pdf")

# Guideline-conformant publication theme. base_size 9pt = readable at the
# 9 x 6 in canvas below. Direct labels, no minor grid, thin axis lines.
theme_robust <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.3),
    axis.ticks       = element_line(colour = "black", linewidth = 0.3),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.position   = "top",
    legend.margin     = margin(b = -3),
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 8),
    plot.title        = element_text(face = "plain", size = 10,
                                     margin = margin(b = 2)),
    plot.subtitle     = element_text(size = 8, colour = "grey30",
                                     margin = margin(b = 4)),
    plot.tag          = element_text(face = "plain", size = 11),
    plot.margin       = margin(6, 8, 6, 8)
  )

# Colorblind-safe categorical palette (Okabe-Ito-ish via Set2)
PAL_DEG  <- c("Canonical DEG" = "#D55E00", "Non-DEG" = "#999999")
PAL_PART <- c("raw expression" = "#999999", "dream residuals" = "#0072B2")

# -----------------------------------------------------------------------------
# A: CPSS pi-hat distribution (Pillar A)
# -----------------------------------------------------------------------------
stab <- fread(file.path(AUDIT, "pillar_A_stability.csv"))
stab[, deg_class := fifelse(is_canonical_DEG, "Canonical DEG", "Non-DEG")]
stab[, deg_class := factor(deg_class, levels = c("Non-DEG", "Canonical DEG"))]

n_deg <- stab[is_canonical_DEG == TRUE, .N]
pct_pass_A <- round(100 * mean(stab[is_canonical_DEG == TRUE, cpss_pi_hat] >= 0.7,
                               na.rm = TRUE), 0)

p_a <- ggplot(stab[!is.na(cpss_pi_hat)],
              aes(x = cpss_pi_hat, fill = deg_class, colour = deg_class)) +
  geom_density(alpha = 0.5, linewidth = 0.4) +
  geom_vline(xintercept = 0.7, linetype = "dashed",
             colour = "grey20", linewidth = 0.4) +
  annotate("text", x = 0.7, y = Inf, label = " stability\n cutoff", hjust = 0,
           vjust = 1.4, size = 2.6, colour = "grey20") +
  scale_fill_manual(values = PAL_DEG, name = NULL) +
  scale_colour_manual(values = PAL_DEG, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(tag = "A",
       title = "Sample-level stability",
       subtitle = sprintf("%d%% of mega DEGs hit pi-hat >= 0.7 (CPSS, 100 paired half-samples)",
                          pct_pass_A),
       x = expression("CPSS selection probability " * hat(pi)),
       y = "Density") +
  theme_robust

# -----------------------------------------------------------------------------
# B: LOCO disease-prediction AUROC (Pillar B)
# -----------------------------------------------------------------------------
loco <- fread(file.path(AUDIT, "pillar_B_loco_prediction.csv"))
loco[, method_lbl := factor(method, levels = c("dotprod", "signed_ssgsea"),
                            labels = c("Dot-product", "Signed ssGSEA"))]
loco[, held_out := factor(held_out, levels = unique(held_out)[order(loco[
  method == "dotprod", auroc])])]

p_b <- ggplot(loco, aes(x = auroc, y = held_out)) +
  geom_segment(aes(x = null_mean, xend = auroc, yend = held_out,
                   colour = method_lbl),
               position = position_dodge(width = 0.55),
               linewidth = 0.4, alpha = 0.6, show.legend = FALSE) +
  geom_point(aes(x = null_mean, group = method_lbl),
             position = position_dodge(width = 0.55),
             shape = 4, size = 1.6, colour = "grey50", show.legend = FALSE) +
  geom_point(aes(colour = method_lbl, shape = method_lbl),
             position = position_dodge(width = 0.55), size = 2.2) +
  geom_vline(xintercept = 0.5, linetype = "dotted",
             colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0.8, linetype = "dashed",
             colour = "grey20", linewidth = 0.3) +
  scale_colour_manual(values = c("Dot-product" = "#0072B2",
                                 "Signed ssGSEA" = "#009E73"),
                      name = NULL) +
  scale_shape_manual(values = c("Dot-product" = 16,
                                "Signed ssGSEA" = 17),
                     name = NULL) +
  scale_x_continuous(limits = c(0.4, 1.02), breaks = seq(0.4, 1, 0.1)) +
  labs(tag = "B",
       title = "Cross-cohort transfer",
       subtitle = "Train on 4 cohorts, test on held-out 5th. x = random-label null.",
       x = "AUROC", y = NULL) +
  theme_robust

# -----------------------------------------------------------------------------
# C: Within-cohort permutation null (Pillar C)
# -----------------------------------------------------------------------------
count_dt <- fread(file.path(AUDIT, "pillar_C_perm_count_dist.csv"))
sum_dt   <- fread(file.path(AUDIT, "pillar_C_summary.csv"))
n_obs <- sum_dt$n_deg_observed[1]
emp_fdr <- sum_dt$empirical_fdr[1]
n_perm  <- nrow(count_dt)
perm_q95 <- sum_dt$perm_count_q95[1]

# Display: histogram of permutation counts with observed annotated separately
# because n_obs is far above the perm distribution and would compress the bins.
p_c <- ggplot(count_dt, aes(x = n_deg_perm)) +
  geom_histogram(fill = "grey85", colour = "grey40", linewidth = 0.25,
                 bins = 25) +
  geom_vline(xintercept = perm_q95, linetype = "dashed",
             colour = "grey20", linewidth = 0.3) +
  annotate("text", x = perm_q95, y = Inf, hjust = -0.1, vjust = 1.4,
           label = sprintf(" 95th = %d", perm_q95),
           size = 2.6, colour = "grey20") +
  annotate("text", x = max(count_dt$n_deg_perm) * 0.6, y = Inf,
           hjust = 0, vjust = 3.0,
           label = sprintf("Observed = %s\nEmpirical FDR = %.3f",
                           comma(n_obs), emp_fdr),
           size = 3.0, colour = "#D55E00", fontface = "plain") +
  scale_x_continuous(labels = comma,
                     limits = c(0, max(perm_q95 + 2, max(count_dt$n_deg_perm)))) +
  labs(tag = "C",
       title = "Empirical null DEG count",
       subtitle = sprintf("B = %d within-cohort case/control permutations",
                          n_perm),
       x = "Permuted DEG count (padj<0.05, |LFC|>0.3)", y = "Permutations") +
  theme_robust

# -----------------------------------------------------------------------------
# D: Variance partition raw vs dream residual (Pillar D)
# -----------------------------------------------------------------------------
INT  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
raw  <- fread(file.path(INT,   "variance_partition.csv"))
resd <- fread(file.path(AUDIT, "pillar_D_residual_varpart.csv"))

melt_long <- function(dt, lbl) {
  cols <- intersect(c("dataset", "condition", "group_binary",
                      "sex", "inferred_sex", "Residuals"), names(dt))
  m <- melt(dt[, c("gene", cols), with = FALSE],
            id.vars = "gene", variable.name = "component",
            value.name = "var_frac")
  m[, partition := lbl]
  m[]
}
raw_l  <- melt_long(raw,  "raw expression")
resd_l <- melt_long(resd, "dream residuals")
both <- rbind(raw_l, resd_l, fill = TRUE)
both[, component_lbl := fcase(
  component %in% c("dataset"), "Cohort",
  component %in% c("condition", "group_binary"), "Disease",
  component %in% c("sex", "inferred_sex"), "Sex",
  component == "Residuals", "Residual",
  default = as.character(component))]
agg <- both[, .(median_frac = median(var_frac, na.rm = TRUE)),
            by = .(partition, component_lbl)]
agg[, component_lbl := factor(component_lbl,
                              levels = c("Cohort", "Disease", "Sex", "Residual"))]
agg[, partition := factor(partition,
                          levels = c("raw expression", "dream residuals"))]

p_d <- ggplot(agg, aes(x = component_lbl, y = median_frac,
                       colour = partition, group = partition)) +
  geom_line(linewidth = 0.4, alpha = 0.6,
            position = position_dodge(width = 0.35)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = PAL_PART, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, NA)) +
  labs(tag = "D",
       title = "Variance partition",
       subtitle = "Cohort variance collapses after dream RE",
       x = NULL, y = "Median variance fraction") +
  theme_robust

# -----------------------------------------------------------------------------
# E: dream vs dream+SVA logFC (Pillar D)
# -----------------------------------------------------------------------------
sva <- fread(file.path(AUDIT, "pillar_D_sva_concordance.csv"))
sva_summary <- fread(file.path(AUDIT, "pillar_D_sva_summary.csv"))

# Restrict scatter to canonical DEGs and a sampled background for clarity.
sva[, panel_class := fifelse(sig_primary | sig_sva, "DEG", "Non-DEG")]
set.seed(42)
sva_plot <- rbind(
  sva[panel_class == "DEG"],
  sva[panel_class == "Non-DEG"][sample(.N, min(2000, .N))]
)
sva_plot[, panel_class := factor(panel_class, levels = c("Non-DEG", "DEG"))]

p_e <- ggplot(sva_plot, aes(x = logFC_primary, y = logFC_sva,
                            colour = panel_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.45, size = 0.4) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
           label = sprintf("rho = %.3f\nJaccard = %.3f\nn_SV = %d",
                           sva_summary$rho_logFC[1],
                           sva_summary$jaccard[1],
                           sva_summary$n_sv[1]),
           size = 2.7, colour = "grey20") +
  scale_colour_manual(values = c("Non-DEG" = "grey75",
                                 "DEG"     = "#0072B2"),
                      name = NULL) +
  labs(tag = "E",
       title = "dream vs dream+SVA",
       subtitle = "SVA adds essentially no new effect-size information",
       x = "log2FC (dream)", y = "log2FC (dream + SVA)") +
  theme_robust +
  theme(legend.position = "top") +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1)))

# -----------------------------------------------------------------------------
# F: Joint pillar-pass distribution (Atlas E1)
# -----------------------------------------------------------------------------
atlas <- fread(file.path(AUDIT, "dream_robustness_atlas.csv"))
deg <- atlas[is_canonical_DEG == TRUE]
e1_dt <- deg[, .N, by = n_pillars_passed][order(n_pillars_passed)]
e1_dt[, pct := 100 * N / sum(N)]
e1_dt[, n_pillars_passed := factor(n_pillars_passed,
                                   levels = sort(unique(n_pillars_passed)))]

p_f <- ggplot(e1_dt, aes(x = n_pillars_passed, y = N)) +
  geom_segment(aes(xend = n_pillars_passed, y = 0, yend = N),
               colour = "grey60", linewidth = 0.4) +
  geom_point(aes(colour = n_pillars_passed), size = 4) +
  geom_text(aes(label = sprintf("%s (%.0f%%)",
                                comma(N), pct)),
            vjust = -0.9, size = 2.6, colour = "grey15") +
  scale_colour_manual(values = c("0" = "#999999", "1" = "#F0E442",
                                 "2" = "#56B4E9", "3" = "#009E73"),
                      guide = "none") +
  scale_y_continuous(labels = comma,
                     expand = expansion(mult = c(0, 0.20))) +
  labs(tag = "F",
       title = "Joint per-gene robustness",
       subtitle = sprintf("Pillars A + C + D, mega DEGs (n = %s)",
                          comma(nrow(deg))),
       x = "# pillars passed (out of 3)", y = "Mega DEGs") +
  theme_robust

# -----------------------------------------------------------------------------
# Compose 2 x 3
# -----------------------------------------------------------------------------
composite <- (p_a | p_b | p_c) / (p_d | p_e | p_f)

ggsave(OUT_PDF, composite, width = 9.5, height = 6.0, device = cairo_pdf)
cat("Saved: ", OUT_PDF, "\n", sep = "")
cat("\nDone.\n")
