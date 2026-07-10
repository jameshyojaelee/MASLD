#!/usr/bin/env Rscript
# ===========================================================================
# figS_multimethod_diagnostics.R
# Cross-method DE statistical-diagnostic panels (light render; reads only the
# CSV intermediates from diag_panels_LMNO_compute.R).
#   panelL_variance_partition.pdf    — shared raw variance partition (846 set)
#   panelM_residual_pvca.pdf         — per-method residual PVCA (before/after)
#   panelN_genomic_inflation_qq.pdf  — lambda_GC + QQ, per method
#   panelO_pvalue_distribution.pdf   — p-value histograms + Storey pi1, per method
# Output: figures/.../multimethod_validation/panels/de_method_diagnostics/
# Conventions: theme_masld; Control/neutral #9E9E9E; cairo_pdf; individual PDFs.
# ===========================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
DIAG <- file.path(INT, "multimethod_validation/diagnostics")
OUT  <- file.path(BASE, "figures/supplementary/figS_methods_validation/multimethod_validation/panels/de_method_diagnostics")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

rd <- function(f) { p <- file.path(DIAG, f); if (file.exists(p)) fread(p) else { cat("SKIP (missing):", f, "\n"); NULL } }
sav <- function(p, file, w, h) { ggsave(file.path(OUT, file), p, width = w, height = h, device = cairo_pdf); cat("Wrote", file, "\n") }
wrap <- function(s, w = 100) paste(strwrap(s, width = w), collapse = "\n")   # avoid subtitle/caption clipping

method_levels <- c("limma-voom QW (canonical)", "dream", "DESeq2", "metafor (RE)", "edgeR-QLF")
method_pal <- c("limma-voom QW (canonical)" = "#C9265E", "dream" = "#40b499",
                "DESeq2" = "#4aa2c2", "metafor (RE)" = "#9b75d6", "edgeR-QLF" = "#E8A33D")

# ===========================================================================
# Panel L — shared raw variance partition (846 mega set)
# ===========================================================================
vp <- rd("variance_partition_mega846.csv")
if (!is.null(vp)) {
  comp_cols <- setdiff(names(vp), "gene")
  long <- melt(vp, id.vars = "gene", measure.vars = comp_cols,
               variable.name = "component", value.name = "frac")
  meds <- long[, .(med = median(frac, na.rm = TRUE)), by = component][order(-med)]
  long[, component := factor(component, levels = meds$component)]
  comp_pal <- c(dataset = "#4aa2c2", group_binary = "#C9265E",
                inferred_sex = "#9b75d6", Residuals = "#9E9E9E")
  lab <- c(dataset = "dataset\n(cohort/batch)", group_binary = "group_binary\n(disease)",
           inferred_sex = "inferred_sex", Residuals = "Residuals")
  meds[, ytxt := med + 0.04]
  pL <- ggplot(long, aes(component, frac, fill = component)) +
    geom_violin(scale = "width", linewidth = 0.2, colour = "grey30") +
    geom_boxplot(width = 0.12, outlier.size = 0.15, outlier.alpha = 0.2,
                 linewidth = 0.25, fill = "white") +
    geom_text(data = meds, aes(component, ytxt, label = sprintf("%.1f%%", 100 * med)),
              inherit.aes = FALSE, size = GEOM_TEXT_6PT, fontface = "plain", colour = "black") +
    scale_fill_manual(values = comp_pal, guide = "none") +
    scale_x_discrete(labels = lab) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = NULL, y = "Variance fraction per gene",
         caption = "Source: variance_partition_mega846.csv (846 mega subset; NOT the stale 1,260-sample variance_partition.csv).") +
    theme_masld(base_size = 6) +
    theme(plot.caption = element_text(size = 6, colour = "black", hjust = 0))
  message(sprintf("[caption] Variance partition of the pooled expression data (shared across methods): all five methods condition on dataset+sex, so this raw decomposition is identical for every method. n=%s genes x 846 samples; model ~ (1|group_binary)+(1|inferred_sex)+(1|dataset). Cohort dominates (median %.0f%%, agreeing with the cohort-driven PC1 in panelJ); disease is a minor per-gene axis (median %.1f%%) yet drives >1,400 Tier-1 DEGs - small variance fraction != weak signal at n=846. Per-method differences appear in panel M.",
                  format(nrow(vp), big.mark = ","),
                  100 * meds[component == "dataset", med], 100 * meds[component == "group_binary", med]))
  sav(pL, "variance_partition.pdf", 6.6, 3.8)
  fwrite(meds[, .(component, median_frac = med)], file.path(OUT, "variance_partition_data.csv"))
}

# ===========================================================================
# Panel M — per-method residual PVCA (before vs after batch correction)
# ===========================================================================
pv <- rd("pvca_before_after.csv")
if (!is.null(pv)) {
  cov_lab <- c(dataset = "dataset\n(cohort/batch)", group_binary = "group_binary\n(disease)",
               inferred_sex = "inferred_sex", Residuals = "Residuals")
  pv[, covariate := factor(covariate, levels = rev(c("dataset", "group_binary", "inferred_sex", "Residuals")))]
  corr_levels <- c("Raw", "Fixed (limma-QW/DESeq2/edgeR)", "Random (dream)")
  pv[, correction := factor(correction, levels = corr_levels)]
  corr_pal <- setNames(c("#9E9E9E", "#4aa2c2", "#40b499"), corr_levels)
  # pillar_D validation note (dream all-gene residual cohort variance ~ 0)
  pdsum <- file.path(BASE, "RNA-seq/results/audit_sensitivity/pillar_D_residual_varpart_summary.csv")
  pd_note <- if (file.exists(pdsum)) {
    s <- fread(pdsum); dv <- s[component == "dataset"]
    sprintf("Validation: dream all-gene mixed-residual cohort variance ~ 0%% (median %.0f%%, mean %.3f%%; pillar_D, %s genes).",
            100 * dv$median_var, 100 * dv$mean_var, format(dv$n_genes, big.mark = ","))
  } else ""
  pM <- ggplot(pv, aes(pct_variance, covariate, fill = correction)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72, linewidth = 0) +
    scale_fill_manual(values = corr_pal, name = "Batch correction") +
    scale_y_discrete(labels = cov_lab) +
    labs(x = "% variance (eigenvalue-weighted PVCA)", y = NULL) +
    theme_masld(base_size = 6) +
    theme(legend.position = "top")
  message(sprintf("[caption] Residual variance after each method's batch model: eigenvalue-weighted PVCA on top-2,000 within-cohort-HVG log-CPM. Fixed = limma-voom-QW / DESeq2 / edgeR-QLF (removeBatchEffect, ~dataset); Random = dream ((1|dataset) shrunken BLUPs). Both collapse cohort variance to ~0 (27.5%% -> 0.7%%) while disease is preserved (re-proportioned upward); fixed ~ random (r~0.9999 PCA equivalence). metafor (RE): n/a - never pools, so a joint residual is undefined. %s",
                  pd_note))
  sav(pM, "residual_pvca.pdf", 6.8, 4.0)
  fwrite(pv, file.path(OUT, "residual_pvca_data.csv"))
}

# ===========================================================================
# Panel N — genomic inflation lambda_GC + QQ, per method
# ===========================================================================
qq <- rd("qq_data.csv"); lam <- rd("lambda_table.csv")
if (!is.null(qq) && !is.null(lam)) {
  qq[, method := factor(method, levels = method_levels)]
  lam[, method := factor(method, levels = method_levels)]
  band <- unique(qq[method == method_levels[1], .(exp, lo, hi)])[order(exp)]
  lam_lab <- setNames(sprintf("%s  (λ=%.2f)", lam$method, lam$lambda_common), as.character(lam$method))
  n_common <- max(lam$n_common, na.rm = TRUE)

  qqp <- ggplot() +
    geom_ribbon(data = band, aes(exp, ymin = lo, ymax = hi), fill = "grey85", alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "grey45") +
    geom_line(data = qq, aes(exp, obs, colour = method), linewidth = 0.5) +
    scale_colour_manual(values = method_pal, labels = lam_lab, name = NULL) +
    labs(x = expression(Expected ~ -log[10](p)), y = expression(Observed ~ -log[10](p))) +
    theme_masld(base_size = 6) +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_rect(fill = scales::alpha("white", 0.6), colour = NA))
  message(sprintf("[caption] Genomic inflation (QQ), per DE method: common universe n=%s genes. With pi1~0.6-0.7 most genes are non-null, so the QQ departs the diagonal early and lambda_GC is mechanically large - this is real polygenic-scale signal, not miscalibration. The calibrated-null check is the flat p-value shelf in panel O.",
                  format(n_common, big.mark = ",")))

  lamp <- ggplot(lam, aes(lambda_common, method, colour = method)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey55") +
    geom_segment(aes(x = 1, xend = lambda_common, yend = method), linewidth = 0.5) +
    geom_point(size = 2) +
    geom_text(aes(label = sprintf("%.2f", lambda_common)), vjust = -0.9, size = GEOM_TEXT_6PT, colour = "black") +
    scale_colour_manual(values = method_pal, guide = "none") +
    scale_y_discrete(limits = rev(method_levels)) +
    labs(x = expression(lambda[GC] ~ "(common set)"), y = NULL) +
    theme_masld(base_size = 6)

  pN <- qqp + lamp + plot_layout(widths = c(2, 1)) +
    plot_annotation(
      caption = wrap("lambda_GC > 1 is EXPECTED, not miscalibration: an 846-sample contrast perturbs >1,400 Tier-1 genes (1,433 ashr / 1,853 raw). Because the median gene is non-null (pi1 > 0.5), the median-based lambda_GC is mechanically inflated by real biology. Read calibration from panel O (flat null shelf); type-I error is held at 0.05 under within-cohort permutation (degx panels G/G3). metafor (RE) is lowest, tracking its low pi1 (underpowered RE meta at high I-squared).", 150),
      theme = theme(plot.caption = element_text(size = 6, colour = "black", hjust = 0)))
  sav(pN, "genomic_inflation_qq.pdf", fig_full_width, 3.5)
  fwrite(lam, file.path(OUT, "genomic_inflation_qq_lambda_table.csv"))
}

# ===========================================================================
# Panel O — p-value distribution + Storey pi1, per method
# ===========================================================================
pl <- rd("pvalue_long.csv"); pi1 <- rd("pi1_table.csv")
if (!is.null(pl)) {
  pl[, method := factor(method, levels = method_levels)]
  NB <- 40
  expdf <- pl[, .(yint = .N / NB), by = method]                       # uniform expectation
  ann <- if (!is.null(pi1)) {
    pi1[, method := factor(method, levels = method_levels)]
    merge(pi1[, .(method, pi1_full)], pl[, .(n = .N), by = method], by = "method")
  } else NULL
  pO <- ggplot(pl, aes(P.Value, fill = method)) +
    geom_histogram(bins = NB, boundary = 0, colour = "white", linewidth = 0.1) +
    geom_hline(data = expdf, aes(yintercept = yint), linetype = "dashed",
               linewidth = 0.3, colour = "grey45") +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_manual(values = method_pal, guide = "none") +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = "raw p-value", y = "genes") +
    theme_masld(base_size = 6)
  message("[caption] P-value distribution per DE method: flat shelf toward p=1 (matching the dashed uniform reference) = calibrated null; near-zero spike = true-DEG mass. Storey pi1 = estimated non-null fraction (complement of pi0); pi1 and the panel-N lambda are two views of the same true signal. This shelf - not panel N's QQ - is the calibration read.")
  if (!is.null(ann))
    pO <- pO + geom_text(data = ann, aes(x = 0.97, y = Inf, label = sprintf("π₁=%.2f", pi1_full)),
                         inherit.aes = FALSE, hjust = 1, vjust = 1.6, size = GEOM_TEXT_6PT, colour = "black")
  sav(pO, "pvalue_distribution.pdf", fig_full_width, 2.2)
  if (!is.null(pi1)) fwrite(pi1, file.path(OUT, "pvalue_distribution_pi1_table.csv"))
}

cat("\n[figS_multimethod_diagnostics] done -> ", OUT, "\n")
