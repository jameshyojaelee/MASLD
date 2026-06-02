#!/usr/bin/env Rscript
# 341_ploidy_disease_comparison.R
# Phase B + A.5 cross-validation: donor-level scPloidy vs disease + signature correlation.
# Inputs:
#   Analysis/SingleCell/results_gpu_v2/ploidy/scploidy_donor_summary.csv     (Phase A output)
#   Analysis/SingleCell/results_gpu_v2/ploidy/signature_scores_donor_summary.csv (Phase A.5 output, optional)
#   Analysis/ATAC/Human_Multiome/results/sex_validation/donor_sex_inference.csv
# Outputs:
#   Analysis/SingleCell/results_gpu_v2/ploidy/ploidy_disease_stats.csv
#   Analysis/SingleCell/results_gpu_v2/ploidy/scploidy_vs_signature_correlation.csv (if Phase A.5 ran)
#   figures/supplementary/figS_ploidy_disease.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
setwd(ROOT)

OUT_DIR    <- "Analysis/SingleCell/results_gpu_v2/ploidy"
FIG_DIR    <- "figures/supplementary"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

SUMMARY    <- file.path(OUT_DIR, "scploidy_donor_summary.csv")
DONOR_META <- "Analysis/ATAC/Human_Multiome/results/sex_validation/donor_sex_inference.csv"
SIG_DONOR  <- file.path(OUT_DIR, "signature_scores_donor_summary.csv")

stopifnot(file.exists(SUMMARY), file.exists(DONOR_META))

ploidy <- fread(SUMMARY)
meta   <- fread(DONOR_META)[, .(donor_id, inferred_sex, condition_meta = condition)]
ploidy <- merge(ploidy, meta, by = "donor_id", all.x = TRUE)
# Prefer condition from donor_sex_inference (canonical) if present
if (!"condition" %in% names(ploidy) || any(is.na(ploidy$condition))) {
  ploidy[, condition := condition_meta]
}
ploidy[, condition_binary := ifelse(condition == "NORMAL", "NORMAL", "DISEASE")]
ploidy[, condition := factor(condition, levels = c("NORMAL", "MASL", "MASH"))]

# Keep only successfully fit donors
ok <- ploidy[status == "ok"]
cat("Donors successfully fit:", nrow(ok), "of", nrow(ploidy), "\n")
print(ok[, .(donor_id, condition, inferred_sex, n_cells_fit, frac_2N, frac_4N, frac_8N, polyploid_frac)])

# --- Phase B tests ----------------------------------------------------------
do_test <- function(dt, value, condition_var) {
  v <- dt[[value]]; g <- dt[[condition_var]]
  out <- list()
  # Wilcoxon NORMAL vs DISEASE
  if (length(unique(g)) >= 2 && all(c("NORMAL", "DISEASE") %in% as.character(g))) {
    w <- wilcox.test(v ~ g)
    out$wilcox_p <- w$p.value
    out$wilcox_W <- w$statistic
  } else {
    out$wilcox_p <- NA; out$wilcox_W <- NA
  }
  # Kruskal-Wallis across all levels (nonparam alternative to ANOVA)
  if (length(unique(g)) >= 3) {
    k <- kruskal.test(v ~ g)
    out$kw_p <- k$p.value
    out$kw_chisq <- k$statistic
  } else {
    out$kw_p <- NA; out$kw_chisq <- NA
  }
  out
}

tests <- list()
# Primary metric: continuous mean fractional ploidy (coverage-robust)
# Secondary: discrete calls + depth-filtered sensitivity
candidate_metrics <- c("mean_fractional_ploidy", "median_fractional_ploidy",
                       "polyploid_frac", "polyploid_frac_hi_depth",
                       "frac_8N", "frac_8N_hi_depth",
                       "frac_4N", "frac_2N")
metric_subset <- intersect(candidate_metrics, names(ok))
for (metric in metric_subset) {
  bin <- do_test(ok, metric, "condition_binary")
  tri <- do_test(ok, metric, "condition")
  # also report mean per group
  means_bin <- ok[, .(N = .N, mean = mean(.SD[[1]], na.rm = TRUE)),
                  by = condition_binary, .SDcols = metric]
  means_tri <- ok[, .(N = .N, mean = mean(.SD[[1]], na.rm = TRUE)),
                  by = condition,        .SDcols = metric]
  tests[[length(tests) + 1L]] <- data.table(
    metric = metric,
    n_NORMAL = means_bin[condition_binary == "NORMAL", N][1],
    n_DISEASE = means_bin[condition_binary == "DISEASE", N][1],
    mean_NORMAL = means_bin[condition_binary == "NORMAL", mean][1],
    mean_DISEASE = means_bin[condition_binary == "DISEASE", mean][1],
    diff_DISEASE_minus_NORMAL = means_bin[condition_binary == "DISEASE", mean][1] -
                                means_bin[condition_binary == "NORMAL", mean][1],
    wilcox_p_NORMAL_vs_DISEASE = bin$wilcox_p,
    kruskal_p_3lvl = tri$kw_p,
    mean_MASL = means_tri[condition == "MASL", mean][1],
    mean_MASH = means_tri[condition == "MASH", mean][1]
  )
}
disease_stats <- rbindlist(tests, fill = TRUE)
fwrite(disease_stats, file.path(OUT_DIR, "ploidy_disease_stats.csv"))
cat("\nDisease stats:\n"); print(disease_stats)

# Sex interaction model
if ("inferred_sex" %in% names(ok) && length(unique(ok$inferred_sex)) > 1) {
  sex_model <- lm(polyploid_frac ~ condition_binary * inferred_sex, data = ok)
  cat("\nSex-interaction model (polyploid_frac):\n")
  print(summary(sex_model)$coefficients)
  sink(file.path(OUT_DIR, "sex_interaction_model.txt"))
  print(summary(sex_model))
  sink()
}

# --- Phase A.5 cross-validation --------------------------------------------
DONOR_PAIRING <- "data/GSE244832/metadata/donor_pairing.csv"
if (file.exists(SIG_DONOR) && file.exists(DONOR_PAIRING)) {
  sig <- fread(SIG_DONOR)
  pairs <- fread(DONOR_PAIRING)
  # Expand pairs to one row per rna_srr (semicolon-separated in rna_srrs)
  pairs_long <- pairs[, .(rna_srr = unlist(strsplit(rna_srrs, ";"))),
                      by = .(donor_id, condition, atac_srr, confidence)]
  # Merge sig with pairs by sample(=SRR) -> add donor_id
  sig <- merge(sig, pairs_long[, .(sample = rna_srr, donor_id)], by = "sample")
  # Aggregate per donor (mean across the donor's SRRs)
  sig_cols <- grep("^sig_", names(sig), value = TRUE)
  agg_cols <- c(sig_cols, "S_score", "G2M_score", "log1p_total_counts", "frac_G1")
  agg_cols <- intersect(agg_cols, names(sig))
  sig_donor <- sig[, lapply(.SD, mean, na.rm = TRUE), .SDcols = agg_cols, by = donor_id]
  cat("\nPer-donor signature scores (n donors after merge):", nrow(sig_donor), "\n")
  fwrite(sig_donor, file.path(OUT_DIR, "signature_scores_perdonor_merged.csv"))
  m <- merge(ok, sig_donor, by = "donor_id", all.x = TRUE, suffixes = c("", ".sig"))
  corr_rows <- list()
  # Use mean_fractional_ploidy as the primary scPloidy metric (coverage-robust)
  primary_scploidy <- if ("mean_fractional_ploidy" %in% names(m)) "mean_fractional_ploidy" else "polyploid_frac"
  for (sc in sig_cols) {
    if (sum(!is.na(m[[sc]])) < 3) next
    cor_obj <- suppressWarnings(cor.test(m[[primary_scploidy]], m[[sc]], method = "spearman"))
    corr_rows[[length(corr_rows) + 1L]] <- data.table(
      signature = sc,
      scploidy_metric = primary_scploidy,
      spearman_rho = unname(cor_obj$estimate),
      spearman_p = cor_obj$p.value,
      n_paired = sum(!is.na(m[[sc]]))
    )
  }
  if (length(corr_rows)) {
    corr_dt <- rbindlist(corr_rows)
    fwrite(corr_dt, file.path(OUT_DIR, "scploidy_vs_signature_correlation.csv"))
    cat("\nscPloidy vs signature correlations (donor level):\n")
    print(corr_dt)
  } else {
    cat("\nNo signature columns matched between donor summaries.\n")
  }
} else {
  cat("\n[Phase A.5 not run yet — skipping signature cross-validation]\n")
}

# --- Figures ----------------------------------------------------------------
theme_set(theme_classic(base_size = 11))
pal3 <- c("NORMAL" = "#4d9221", "MASL" = "#fdb863", "MASH" = "#b2182b")

# Panel A: stacked proportions per donor
long <- melt(ok[, .(donor_id, condition, frac_2N, frac_4N, frac_8N)],
             id.vars = c("donor_id", "condition"),
             variable.name = "ploidy", value.name = "frac")
long[, ploidy := factor(ploidy, levels = c("frac_2N", "frac_4N", "frac_8N"),
                                labels = c("2N", "4N", "8N"))]
long[, donor_id := factor(donor_id, levels = ok[order(condition, -polyploid_frac)]$donor_id)]
pA <- ggplot(long, aes(donor_id, frac, fill = ploidy)) +
  geom_col() +
  scale_fill_manual(values = c("2N" = "#4d9221", "4N" = "#e7298a", "8N" = "#7570b3")) +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Ploidy fraction", title = "scPloidy hepatocyte ploidy by donor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Panel B: polyploid fraction boxplot
pB <- ggplot(ok, aes(condition, polyploid_frac, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.18, height = 0, size = 1.6, alpha = 0.9) +
  scale_fill_manual(values = pal3) +
  labs(x = NULL, y = "Polyploid fraction (4N+8N)", title = "Polyploid fraction by stage") +
  theme(legend.position = "none")

# Panel C: 8N fraction
pC <- ggplot(ok, aes(condition, frac_8N, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.18, height = 0, size = 1.6, alpha = 0.9) +
  scale_fill_manual(values = pal3) +
  labs(x = NULL, y = "8N fraction", title = "8N hepatocyte fraction") +
  theme(legend.position = "none")

# Panel D: sex stratification (if present)
if ("inferred_sex" %in% names(ok) && length(unique(ok$inferred_sex)) > 1) {
  pD <- ggplot(ok, aes(condition, polyploid_frac, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.18, height = 0, size = 1.6, alpha = 0.9) +
    scale_fill_manual(values = pal3) +
    facet_wrap(~ inferred_sex) +
    labs(x = NULL, y = "Polyploid fraction", title = "Polyploid fraction by sex × stage") +
    theme(legend.position = "none")
} else pD <- patchwork::plot_spacer()

# Panel E: signature correlation scatter (if Phase A.5 ran)
pE <- patchwork::plot_spacer()
if (exists("m") && !is.null(m) && exists("corr_dt")) {
  # plot the strongest correlation
  top_sig <- corr_dt[which.max(abs(spearman_rho)), signature]
  if (!is.na(top_sig) && top_sig %in% names(m)) {
    pE <- ggplot(m, aes(polyploid_frac, get(top_sig), color = condition)) +
      geom_point(size = 2.4) +
      ggrepel::geom_text_repel(aes(label = donor_id), size = 3, show.legend = FALSE) +
      scale_color_manual(values = pal3) +
      labs(title = paste0("scPloidy vs ", top_sig,
                          " (ρ=", round(corr_dt[signature == top_sig, spearman_rho], 2),
                          ", p=", signif(corr_dt[signature == top_sig, spearman_p], 2), ")"),
           x = "Polyploid fraction (scPloidy)",
           y = top_sig)
  }
}

# Compose
fig <- (pA / (pB | pC) / (pD | pE)) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG_DIR, "figS_ploidy_disease.pdf"), fig,
       width = 10, height = 11, device = cairo_pdf)
cat("\nWrote", file.path(FIG_DIR, "figS_ploidy_disease.pdf"), "\n")

# Decision gate summary - use the primary continuous metric (coverage-robust)
cat("\n========== DECISION GATE ==========\n")
primary <- if ("mean_fractional_ploidy" %in% disease_stats$metric) "mean_fractional_ploidy" else "polyploid_frac"
poly_test <- disease_stats[metric == primary]
cat("Primary metric: ", primary, "\n", sep = "")
cat(primary, ":\n", sep = "")
cat("  NORMAL mean =", round(poly_test$mean_NORMAL, 3),
    "(n=", poly_test$n_NORMAL, ")\n")
cat("  DISEASE mean =", round(poly_test$mean_DISEASE, 3),
    "(n=", poly_test$n_DISEASE, ")\n")
cat("  diff =", round(poly_test$diff_DISEASE_minus_NORMAL, 3), "\n")
cat("  Wilcoxon NORMAL vs DISEASE p =", signif(poly_test$wilcox_p_NORMAL_vs_DISEASE, 3), "\n")
cat("  Kruskal-Wallis 3-level p =", signif(poly_test$kruskal_p_3lvl, 3), "\n")
cat("  MASL mean =", round(poly_test$mean_MASL, 3), "\n")
cat("  MASH mean =", round(poly_test$mean_MASH, 3), "\n")
cat("\nIMPORTANT CAVEAT: scPloidy was validated on rat scATAC at ~28k fragments/cell.\n")
cat("Our 10x multiome ATAC has ~1.6k fragments/cell (17x sparser). At this depth,\n")
cat("the discrete 2N/4N/8N calls are biased toward higher ploidy. The continuous\n")
cat("`mean_fractional_ploidy` is more robust to coverage; relative comparisons\n")
cat("across donors are still meaningful if coverage is balanced.\n")
if (!is.na(poly_test$wilcox_p_NORMAL_vs_DISEASE) &&
    poly_test$wilcox_p_NORMAL_vs_DISEASE < 0.05 &&
    poly_test$diff_DISEASE_minus_NORMAL > 0) {
  cat(">>> STRONG SIGNAL: disease-up polyploid_frac, p<0.05. Proceed to Phase C.\n")
} else if (!is.na(poly_test$wilcox_p_NORMAL_vs_DISEASE) &&
           poly_test$wilcox_p_NORMAL_vs_DISEASE < 0.1) {
  cat(">>> TREND ONLY: nominal trend; report with caveats, consider replication before Phase C.\n")
} else {
  cat(">>> NULL: no evidence of disease-associated polyploidy shift in scATAC. Document and stop.\n")
}
