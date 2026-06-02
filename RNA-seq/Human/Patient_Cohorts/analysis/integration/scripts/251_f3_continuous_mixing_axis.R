# 251_f3_continuous_mixing_axis.R
# F3 sub-state (bulk) — Continuous F3a-F3b mixing-fraction axis test.
#
# Phase 1 result: pooled F3 distribution is unimodal on the F0→F4 PCA principal
# axis (PC2). That ruled out interpretation (1) "discrete bimodal sub-states in
# bulk." This script tests interpretation (2): continuous mixing-fraction
# variation across biopsies.
#
# Hypothesis: F3 biopsies vary continuously in their F3a-vs-F3b proportion.
# If so, the F3a-F3b ssGSEA score axis (independent of CRN) should:
#   (a) Show negative correlation between F3a and F3b scores (mixing fraction
#       trade-off), OR positive correlation (co-expression — sub-states co-occur)
#   (b) Span a wide range relative to within-class variance
#   (c) Optionally bimodal on its own axis (separate from the F0→F4 PC)
#
# Inputs:
#   results/granular_staging/f3_substate_signature_scores.csv (from 243)
#   results/integration/meta_matched.rds
#   results/granular_staging/continuous_fibrosis_index.csv (existing PC analysis)
#
# Outputs (results/granular_staging/):
#   f3_substate_bulk_scatter.csv      — per-sample F3a/F3b scores + mixing axis
#   f3_substate_bulk_axis_diag.csv — bimodality + correlation diagnostics
#   figures/supplementary/figS_granular_staging/figS_granular_f3a_vs_f3b_bulk.pdf  — diagnostic plots

suppressPackageStartupMessages({
  library(ggplot2)
  library(matrixStats)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
FIG_DIR <- file.path(PROJECT_ROOT, "figures/supplementary/figS_granular_staging")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== F3 sub-state (bulk) — F3a-F3b continuous mixing-fraction axis ==\n")

# ---------------------------------------------------------------------------
# Load existing F3 signature scores from Phase 1.3
# ---------------------------------------------------------------------------
ssgsea_path <- file.path(OUT_DIR, "f3_substate_signature_scores.csv")
ssgsea <- read.csv(ssgsea_path, row.names = 1, check.names = FALSE,
                    stringsAsFactors = FALSE)
cat(sprintf("ssGSEA scores: %d samples × %d signatures\n",
            nrow(ssgsea), ncol(ssgsea)))
print(colnames(ssgsea))

# Verify the Li et al. anchor signatures exist
stopifnot("LI_F3A_HYPOTHESIS" %in% colnames(ssgsea),
          "LI_F3B_HYPOTHESIS" %in% colnames(ssgsea))

# Load metadata to get cohort/stage info
meta <- readRDS(file.path(INT_DIR, "meta_matched.rds"))
meta_match <- meta[match(rownames(ssgsea), meta$sample_id), ]
stopifnot(all(!is.na(meta_match$dataset)))

cat(sprintf("F3 samples after metadata join: %d\n", nrow(meta_match)))

# Add additional anchors as comparators
F3a_score <- ssgsea[, "LI_F3A_HYPOTHESIS"]
F3b_score <- ssgsea[, "LI_F3B_HYPOTHESIS"]
ECM_score <- if ("REACTOME_ECM_PROTEOGLYCANS" %in% colnames(ssgsea))
  ssgsea[, "REACTOME_ECM_PROTEOGLYCANS"] else NA
COL_score <- if ("REACTOME_COLLAGEN_FORMATION" %in% colnames(ssgsea))
  ssgsea[, "REACTOME_COLLAGEN_FORMATION"] else NA
SEN_score <- if ("SENMAYO_SENESCENCE" %in% colnames(ssgsea))
  ssgsea[, "SENMAYO_SENESCENCE"] else NA
UPR_score <- if ("HALLMARK_UNFOLDED_PROTEIN_RESPONSE" %in% colnames(ssgsea))
  ssgsea[, "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"] else NA
SREBP_score <- if ("REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF" %in% colnames(ssgsea))
  ssgsea[, "REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF"] else NA

# Cohort-correct (regress out cohort effect on each score) for fair across-cohort comparison
adjust_for_cohort <- function(x) {
  fit <- lm(x ~ factor(meta_match$dataset))
  resid(fit) + mean(x, na.rm = TRUE)
}
F3a_adj <- adjust_for_cohort(F3a_score)
F3b_adj <- adjust_for_cohort(F3b_score)

# ---------------------------------------------------------------------------
# (a) Correlation between F3a and F3b scores
# ---------------------------------------------------------------------------
cat("\n--- (a) F3a vs F3b score correlation ---\n")
corr_raw_pearson <- cor(F3a_score, F3b_score, method = "pearson")
corr_raw_spearman <- cor(F3a_score, F3b_score, method = "spearman")
corr_adj_pearson <- cor(F3a_adj, F3b_adj, method = "pearson")
corr_adj_spearman <- cor(F3a_adj, F3b_adj, method = "spearman")

cat(sprintf("Raw      Pearson:  %+.3f  Spearman: %+.3f\n",
            corr_raw_pearson, corr_raw_spearman))
cat(sprintf("CohortAdj Pearson: %+.3f  Spearman: %+.3f\n",
            corr_adj_pearson, corr_adj_spearman))

interpret <- if (abs(corr_adj_pearson) < 0.2) {
  "WEAK correlation — sub-state activities are mostly independent in bulk; supports trade-off interpretation but at small effect size."
} else if (corr_adj_pearson > 0.2) {
  "POSITIVE correlation — F3a and F3b activate together (co-expression). NOT a trade-off; sub-states co-occur in F3 biopsies."
} else {
  "NEGATIVE correlation — F3a-rich biopsies are F3b-poor and vice versa. Supports continuous mixing-fraction interpretation."
}
cat("\nInterpretation:", interpret, "\n")

# ---------------------------------------------------------------------------
# (b) Mixing-axis (F3b - F3a, normalized) — test bimodality
# ---------------------------------------------------------------------------
cat("\n--- (b) Mixing-fraction axis bimodality ---\n")
mixing_axis <- F3b_adj - F3a_adj
mixing_axis_z <- scale(mixing_axis)[, 1]

# Same BIC + dip MC tests as 242
fit_gmm_2 <- function(x, max_iter = 200, tol = 1e-6) {
  n <- length(x)
  ord <- sort(x)
  mu1 <- mean(ord[1:floor(n/2)]); mu2 <- mean(ord[(floor(n/2)+1):n])
  sd1 <- max(sd(ord[1:floor(n/2)]), 1e-3); sd2 <- max(sd(ord[(floor(n/2)+1):n]), 1e-3)
  p1 <- 0.5; ll_old <- -Inf
  for (iter in 1:max_iter) {
    d1 <- p1 * dnorm(x, mu1, sd1); d2 <- (1 - p1) * dnorm(x, mu2, sd2)
    g1 <- d1 / (d1 + d2 + 1e-300)
    p1 <- mean(g1)
    mu1 <- sum(g1 * x) / sum(g1); mu2 <- sum((1 - g1) * x) / sum(1 - g1)
    sd1 <- max(sqrt(sum(g1 * (x - mu1)^2) / sum(g1)), 1e-3)
    sd2 <- max(sqrt(sum((1 - g1) * (x - mu2)^2) / sum(1 - g1)), 1e-3)
    ll <- sum(log(p1 * dnorm(x, mu1, sd1) + (1 - p1) * dnorm(x, mu2, sd2) + 1e-300))
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, p1 = p1, ll = ll)
}
n <- length(mixing_axis_z)
ll1 <- sum(dnorm(mixing_axis_z, 0, 1, log = TRUE))
bic1 <- -2 * ll1 + 2 * log(n)
fit2 <- fit_gmm_2(mixing_axis_z)
bic2 <- -2 * fit2$ll + 5 * log(n)
bic_diff <- bic1 - bic2
sep_sigma <- abs(fit2$mu1 - fit2$mu2) / sqrt((fit2$sd1^2 + fit2$sd2^2) / 2)

# Dip MC
hartigan_dip_pvalue <- function(x, n_mc = 999) {
  obs_dip <- {
    fx <- ecdf(x); grid <- seq(min(x), max(x), length.out = 200)
    e <- fx(grid); g <- pnorm(grid, mean(x), sd(x))
    max(abs(e - g))
  }
  null_dips <- replicate(n_mc, {
    y <- rnorm(length(x), mean(x), sd(x))
    fy <- ecdf(y); grid <- seq(min(y), max(y), length.out = 200)
    e <- fy(grid); g <- pnorm(grid, mean(y), sd(y))
    max(abs(e - g))
  })
  list(dip_obs = obs_dip,
       p_value = (sum(null_dips >= obs_dip) + 1) / (n_mc + 1))
}
dip_res <- hartigan_dip_pvalue(mixing_axis_z, n_mc = 999)

cat(sprintf("Mixing axis bimodality:\n"))
cat(sprintf("  ΔBIC (k=1 vs k=2): %.1f (positive = mixture preferred)\n", bic_diff))
cat(sprintf("  GMM separation:    %.2fσ\n", sep_sigma))
cat(sprintf("  GMM mixing prop:   p1 = %.2f\n", fit2$p1))
cat(sprintf("  Dip-MC p-value:    %.3f\n", dip_res$p_value))

# ---------------------------------------------------------------------------
# (c) Range / variance characterization
# ---------------------------------------------------------------------------
cat("\n--- (c) Score-axis dynamic range vs noise ---\n")
F3a_range_z <- (max(F3a_adj) - min(F3a_adj)) / sd(F3a_adj)
F3b_range_z <- (max(F3b_adj) - min(F3b_adj)) / sd(F3b_adj)
mix_range_z <- (max(mixing_axis_z) - min(mixing_axis_z))
cat(sprintf("  F3a score dynamic range / σ: %.2f\n", F3a_range_z))
cat(sprintf("  F3b score dynamic range / σ: %.2f\n", F3b_range_z))
cat(sprintf("  Mixing axis dynamic range:    %.2f σ\n", mix_range_z))

# ---------------------------------------------------------------------------
# (d) Per-cohort cohort-resilience
# ---------------------------------------------------------------------------
cat("\n--- (d) Per-cohort F3a-F3b correlation (cohort-resilience) ---\n")
per_cohort_corr <- list()
for (cohort in unique(meta_match$dataset)) {
  c_idx <- meta_match$dataset == cohort
  if (sum(c_idx) < 6) next
  c_corr_p <- cor(F3a_score[c_idx], F3b_score[c_idx], method = "pearson")
  c_corr_s <- cor(F3a_score[c_idx], F3b_score[c_idx], method = "spearman")
  per_cohort_corr[[cohort]] <- data.frame(
    cohort = cohort, n = sum(c_idx),
    pearson = c_corr_p, spearman = c_corr_s,
    stringsAsFactors = FALSE
  )
}
per_cohort_corr_df <- do.call(rbind, per_cohort_corr)
print(per_cohort_corr_df)

# ---------------------------------------------------------------------------
# Diagnostics CSV
# ---------------------------------------------------------------------------
diag_df <- data.frame(
  metric = c("F3a_F3b_pearson_raw", "F3a_F3b_spearman_raw",
             "F3a_F3b_pearson_cohortadj", "F3a_F3b_spearman_cohortadj",
             "mixing_axis_BIC_diff_2vs1", "mixing_axis_GMM_sep_sigma",
             "mixing_axis_GMM_p1", "mixing_axis_dipMC_p",
             "F3a_dynamic_range_sigma", "F3b_dynamic_range_sigma",
             "mixing_axis_dynamic_range_sigma"),
  value = c(corr_raw_pearson, corr_raw_spearman,
            corr_adj_pearson, corr_adj_spearman,
            bic_diff, sep_sigma, fit2$p1, dip_res$p_value,
            F3a_range_z, F3b_range_z, mix_range_z),
  stringsAsFactors = FALSE
)
write.csv(diag_df, file.path(OUT_DIR, "f3_substate_bulk_axis_diag.csv"),
          row.names = FALSE)

# Per-sample scatter table
scatter_df <- data.frame(
  sample_id = rownames(ssgsea),
  dataset = meta_match$dataset,
  fibrosis_stage = meta_match$fibrosis_stage,
  F3a_score = F3a_score,
  F3b_score = F3b_score,
  F3a_adj_score = F3a_adj,
  F3b_adj_score = F3b_adj,
  mixing_axis = mixing_axis,
  mixing_axis_z = mixing_axis_z,
  ECM_score = ECM_score,
  COL_score = COL_score,
  SEN_score = SEN_score,
  UPR_score = UPR_score,
  SREBP_score = SREBP_score,
  stringsAsFactors = FALSE
)
write.csv(scatter_df, file.path(OUT_DIR, "f3_substate_bulk_scatter.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Diagnostic plots — publication theme
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source(file.path(PROJECT_ROOT, "scripts/figures/publication_theme.R"))
source(file.path(PROJECT_ROOT, "scripts/figures/load_figure_data.R"))

scatter_plot_df <- data.frame(
  F3a_raw = F3a_score, F3b_raw = F3b_score,
  F3a_adj = F3a_adj,    F3b_adj = F3b_adj,
  mixing_z = mixing_axis_z,
  cohort = meta_match$dataset,
  stringsAsFactors = FALSE
)

# Panel A: F3a vs F3b raw (across cohorts)
panel_A <- ggplot(scatter_plot_df,
                  aes(x = F3a_raw, y = F3b_raw, color = cohort)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "gray25", fill = "gray80", linewidth = 0.4) +
  scale_color_brewer(palette = "Set2", guide = "none") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Pearson r = %+.2f\nSpearman ρ = %+.2f",
                           corr_raw_pearson, corr_raw_spearman),
           hjust = -0.1, vjust = 1.4, size = PUB_GEOM_TEXT, color = "gray25") +
  labs(title = "F3a vs F3b ssGSEA (raw)",
       subtitle = "Per F3 sample (cohort-coloured points)",
       x = "Li F3a (UPR / SREBP / lipid)",
       y = "Li F3b (ECM / senescence / Ig)") +
  theme_masld() + theme_pub()

# Panel B: cohort-adjusted scatter
panel_B <- ggplot(scatter_plot_df,
                  aes(x = F3a_adj, y = F3b_adj)) +
  geom_point(size = 0.7, alpha = 0.55, color = masld_colors$mash) +
  geom_smooth(method = "lm", se = TRUE, color = "gray25",
              fill = "gray80", linewidth = 0.4) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Pearson r = %+.2f (cohort-adj)\nessentially independent",
                           corr_adj_pearson),
           hjust = -0.1, vjust = 1.4, size = PUB_GEOM_TEXT, color = "gray25") +
  labs(title = "Cohort-adjusted scatter",
       subtitle = "F3a / F3b vary independently in bulk",
       x = "Li F3a (cohort-adj.)", y = "Li F3b (cohort-adj.)") +
  theme_masld() + theme_pub()

# Panel C: mixing axis density
mix_df <- data.frame(x = mixing_axis_z)
panel_C <- ggplot(mix_df, aes(x = x)) +
  geom_density(fill = masld_colors$fibrosis, color = masld_colors$fibrosis,
               alpha = 0.25, linewidth = 0.5) +
  geom_rug(alpha = 0.25, length = unit(2, "mm"), color = masld_colors$fibrosis) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("n = %d\nΔBIC = %.1f\nDip-MC p = %.3f\nsep = %.2fσ",
                           n, bic_diff, dip_res$p_value, sep_sigma),
           hjust = -0.05, vjust = 1.4, size = PUB_GEOM_TEXT, color = "gray25") +
  labs(title = "Mixing axis (F3b − F3a)",
       subtitle = "Unimodal at biopsy aggregation",
       x = "z-scored mixing axis", y = "Density") +
  theme_masld() + theme_pub()

# Panel D: per-cohort correlation forest
if (!is.null(per_cohort_corr_df)) {
  pc_df <- per_cohort_corr_df
  pc_df <- pc_df[order(pc_df$pearson), ]
  pc_df$cohort <- factor(pc_df$cohort, levels = pc_df$cohort)
  pc_df$direction <- ifelse(pc_df$pearson > 0, "positive", "negative")

  panel_D <- ggplot(pc_df, aes(y = cohort, x = pearson, color = direction)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "gray60") +
    geom_segment(aes(x = 0, xend = pearson, y = cohort, yend = cohort),
                 linewidth = 0.5) +
    geom_point(size = 2.0) +
    geom_text(aes(label = sprintf("n=%d", n)),
              hjust = -0.3, size = PUB_GEOM_TEXT, color = "gray35") +
    scale_color_manual(values = c(negative = masld_colors$mash,
                                   positive = masld_colors$control),
                       guide = "none") +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(title = "Per-cohort F3a-F3b correlation",
         subtitle = "Pearson r within each cohort",
         x = "Pearson r", y = NULL) +
    theme_masld() + theme_pub()
} else {
  panel_D <- ggplot() + theme_void()
}

combined <- (panel_A | panel_B) / (panel_C | panel_D) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 9)))

pdf_path <- file.path(FIGS_GRANULAR_DIR, "figS_granular_f3a_vs_f3b_bulk.pdf")
ggsave(pdf_path, combined, width = 7.5, height = 6.5, units = "in",
       device = cairo_pdf)
cat("Wrote", pdf_path, "\n")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n== Decision summary ==\n")
cat(sprintf("F3a vs F3b cohort-adj Pearson: %+.3f\n", corr_adj_pearson))
cat(sprintf("Mixing axis bimodality: ΔBIC=%.1f, dip-MC p=%.3f, sep=%.2fσ\n",
            bic_diff, dip_res$p_value, sep_sigma))
cat(sprintf("Per-cohort sign concordance: %d / %d cohorts have same correlation sign as pooled\n",
            sum(sign(per_cohort_corr_df$pearson) == sign(corr_adj_pearson)),
            nrow(per_cohort_corr_df)))

cat("\n== F3 sub-state (bulk) complete. ==\n")
