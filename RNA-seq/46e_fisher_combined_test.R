#!/usr/bin/env Rscript
##############################################################################
# 46e_fisher_combined_test.R
#
# Fisher's combined probability test across truly INDEPENDENT evidence sources:
#   1. Human bulk RNA-seq (bulk_padj)
#   2. Mouse bulk RNA-seq (mouse_meta_padj)
#   3. GWAS genetics (INTACT posterior -> p-value, or coloc_susie_best_pp4)
#   4. Spatial (spatial_morans_i via atlas spatial_sig or spatial_is_svg)
#
# Fisher's combined chi-sq = -2 * sum(ln(p_i)), df = 2k
# where k = number of non-NA sources per gene.
#
# Also computes:
#   - Spearman correlation between Fisher rank and 46d heuristic rank
#   - COLOC enzyme weighting sensitivity (Task 3): runs 46d scoring with
#     and without the 0.5x liver-enzyme penalty, reports rank correlation
#
# Output:
#   RNA-seq/results/multi_evidence/convergence_evidence_fisher.csv
#   RNA-seq/results/multi_evidence/coloc_enzyme_weighting_sensitivity.csv
#
# Reviewer response: R1 C9, R2 #6, R3 P0-3
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)
cat("=== 46e Fisher Combined Probability Test ===\n")

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ME   <- file.path(BASE, "RNA-seq/results/multi_evidence")

# ============================================================================
# 1. LOAD DATA
# ============================================================================
cat("\n--- 1. Load atlas + INTACT scores ---\n")

atlas <- fread(file.path(ME, "multi_evidence_atlas.csv"))
cat(sprintf("  Atlas: %d rows x %d cols\n", nrow(atlas), ncol(atlas)))
stopifnot(all(c("bulk_padj","bulk_logFC") %in% names(atlas)))

# Load 46d heuristic output for rank comparison
conv_path <- file.path(ME, "convergence_evidence.csv")
if (!file.exists(conv_path)) {
  stop("convergence_evidence.csv not found. Run 46d first.")
}
conv <- fread(conv_path)
cat(sprintf("  46d convergence evidence: %d rows\n", nrow(conv)))

# Genetic-causal score = COLOC PP.H4 (INTACT dropped 2026-06-19; mirrors 46d).
coloc_pp4 <- rep(NA_real_, nrow(atlas))
if ("coloc_susie_best_pp4" %in% names(atlas)) {
  coloc_pp4 <- atlas$coloc_susie_best_pp4
  cat(sprintf("  COLOC PP.H4 (genetic-causal): %d genes\n", sum(!is.na(coloc_pp4))))
}

# ============================================================================
# 2. ASSEMBLE INDEPENDENT P-VALUES
# ============================================================================
# ----------------------------------------------------------------------------
# P1 fix 2026-05-28: CALIBRATION CAVEAT (documentation only; computation unchanged)
# ----------------------------------------------------------------------------
# Two of the four input "p-values" below are NOT null-calibrated Uniform(0,1)
# p-values, they are posterior-predictive / heuristic approximations:
#   * Genetics (S3): p = 1 - INTACT_score (or 1 - coloc PP4 fallback). INTACT is
#     a posterior probability of causality, so 1 - posterior is NOT uniform
#     under the null of no causal effect.
#   * Spatial (S4): the SVG fallback uses p ~ exp(-5 * Moran's I) and a
#     hard-coded 0.01 for SVGs lacking Moran's I — these are bounded heuristics,
#     not a calibrated test statistic's tail probability.
# Because Fisher's method assumes each input p ~ Uniform(0,1) under the null,
# the combined chi-square p produced here is APPROXIMATE. It should be read as a
# convergence/ranking heuristic, NOT a calibrated frequentist combined test.
# Do NOT re-interpret fisher_padj as a strict family-wise/FDR-controlled p.
# ----------------------------------------------------------------------------
cat("\n--- 2. Assemble independent p-values ---\n")
message("CAVEAT: the genetics p (1-INTACT) and spatial heuristic p are ",
        "posterior-predictive approximations, not null-calibrated Uniform(0,1) ",
        "p-values; Fisher combined p is therefore approximate and should be read ",
        "as a ranking heuristic, not a calibrated frequentist test.")

n <- nrow(atlas)

# Source 1: Human bulk RNA-seq (bulk_padj)
p_human <- atlas$bulk_padj
p_human[p_human <= 0] <- NA_real_
p_human[p_human > 1]  <- 1
cat(sprintf("  S1 Human bulk: %d non-NA p-values\n", sum(!is.na(p_human))))

# Source 2: Mouse bulk RNA-seq (mouse_meta_padj)
p_mouse <- atlas$mouse_meta_padj
p_mouse[p_mouse <= 0] <- NA_real_
p_mouse[p_mouse > 1]  <- 1
cat(sprintf("  S2 Mouse bulk: %d non-NA p-values\n", sum(!is.na(p_mouse))))

# Source 3: GWAS genetics (COLOC PP.H4; INTACT dropped 2026-06-19, mirrors 46d).
# Strategy: convert the COLOC posterior PP.H4 to a one-sided p-value via the
# posterior-predictive approximation p = 1 - PP.H4. Primary = SuSiE-COLOC
# (coloc_susie_best_pp4); fall back to ABF-COLOC where SuSiE is absent.
# If both are NA, the gene gets NA for this source.
p_genetic <- rep(NA_real_, n)

# Primary: SuSiE-COLOC PP.H4
has_coloc <- !is.na(coloc_pp4)
p_genetic[has_coloc] <- 1 - coloc_pp4[has_coloc]

# Fallback: coloc_susie_best_pp4
if ("coloc_susie_best_pp4" %in% names(atlas)) {
  pp4 <- atlas$coloc_susie_best_pp4
  fallback_idx <- is.na(p_genetic) & !is.na(pp4) & pp4 > 0
  p_genetic[fallback_idx] <- 1 - pp4[fallback_idx]
}
# Fallback 2: coloc_abf_best_pp4
if ("coloc_abf_best_pp4" %in% names(atlas)) {
  pp4_abf <- atlas$coloc_abf_best_pp4
  fallback_idx2 <- is.na(p_genetic) & !is.na(pp4_abf) & pp4_abf > 0
  p_genetic[fallback_idx2] <- 1 - pp4_abf[fallback_idx2]
}

# Floor at machine epsilon to avoid log(0)
p_genetic[!is.na(p_genetic) & p_genetic <= 0] <- .Machine$double.eps
p_genetic[!is.na(p_genetic) & p_genetic > 1] <- 1
cat(sprintf("  S3 Genetics: %d non-NA p-values (COLOC-SuSiE=%d, COLOC-ABF fallback=%d)\n",
            sum(!is.na(p_genetic)),
            sum(has_coloc),
            sum(!is.na(p_genetic) & !has_coloc)))

# Source 4: Spatial
# Use spatial_is_svg (boolean) + spatial_morans_i to construct a p-value.
# Moran's I p-value is not directly in atlas. We use a conservative proxy:
# for SVG genes, compute p from Moran's I using the normal approximation
# (Z = I / SE(I), p = 2*pnorm(-|Z|)). Since we lack SE(I), we use the
# atlas boolean + spatial_sig column if present.
p_spatial <- rep(NA_real_, n)
if ("spatial_sig" %in% names(atlas)) {
  # spatial_sig is a p-value or boolean. Check type.
  sp_sig <- atlas$spatial_sig
  if (is.numeric(sp_sig)) {
    p_spatial <- sp_sig
    p_spatial[p_spatial <= 0] <- .Machine$double.eps
    p_spatial[p_spatial > 1] <- 1
  }
}
# If spatial_sig is not usable, construct from spatial_hep_wilcoxon_padj_bh
if (sum(!is.na(p_spatial)) < 100 && "spatial_hep_wilcoxon_padj_bh" %in% names(atlas)) {
  sp_wilcox <- atlas$spatial_hep_wilcoxon_padj_bh
  sp_ok <- !is.na(sp_wilcox) & sp_wilcox > 0 & sp_wilcox <= 1
  p_spatial[sp_ok] <- sp_wilcox[sp_ok]
}
# If still sparse, use the SVG boolean as a crude 0.05 / 1.0 proxy
if (sum(!is.na(p_spatial)) < 100 && "spatial_is_svg" %in% names(atlas)) {
  svg_flag <- atlas$spatial_is_svg %in% TRUE
  morans_i <- atlas$spatial_morans_i
  # For SVG genes with Moran's I, approximate p via exponential decay
  # p ~ exp(-5 * max(morans_i, 0)) [heuristic, bounded]
  has_svg <- svg_flag & !is.na(morans_i) & morans_i > 0
  p_spatial[has_svg] <- pmax(exp(-5 * morans_i[has_svg]), .Machine$double.eps)
  # SVG genes without Moran's I: set p = 0.01 (conservative proxy)
  svg_no_i <- svg_flag & (is.na(morans_i) | morans_i <= 0)
  p_spatial[svg_no_i] <- 0.01
  # Non-SVG genes with Moran's I present: set to 1 (non-significant)
  non_svg_with_i <- !svg_flag & !is.na(morans_i)
  p_spatial[non_svg_with_i] <- 1.0
}
cat(sprintf("  S4 Spatial: %d non-NA p-values\n", sum(!is.na(p_spatial))))

# ----------------------------------------------------------------------------
# NEW orthogonal axes (2026-06-01 modality expansion). Each is genuinely
# orthogonal to the existing 4 sources and disease-anchored. df = 2*n_sources
# auto-scales, so an uninformative axis (NA -> contributes 0 to chisq, +0 df)
# cannot inflate the FDR; a real small p boosts chisq but also raises df by 2.
# ----------------------------------------------------------------------------
# S5 metabolite/lipid-QTL COLOC: posterior -> p = 1 - PP4 (same posterior-
# predictive convention as the genetics axis above; orthogonal molecular trait).
p_mqtl <- rep(NA_real_, n)
if ("mqtl_best_pp4" %in% names(atlas)) {
  mq <- atlas$mqtl_best_pp4
  ok <- !is.na(mq) & mq > 0
  p_mqtl[ok] <- pmax(1 - mq[ok], .Machine$double.eps)
}
cat(sprintf("  S5 mQTL: %d non-NA p-values\n", sum(!is.na(p_mqtl))))

# S6 rare-variant burden (Verma 2025 MASLD/PDFF gene-based test): a TRUE
# null-calibrated p-value -- the ideal Fisher input. Orthogonal (rare variants).
p_burden <- rep(NA_real_, n)
if ("burden_masld_pval" %in% names(atlas)) {
  bp <- atlas$burden_masld_pval
  ok <- !is.na(bp) & bp > 0 & bp <= 1
  p_burden[ok] <- bp[ok]
}
cat(sprintf("  S6 rare-variant burden: %d non-NA p-values\n", sum(!is.na(p_burden))))

# S7 liver sQTL credible-set overlap: the raw GTEx nominal p is a molecular-QTL
# p (near-deterministic), NOT a disease-association p, so feeding it directly
# would conflate QTL strength with disease evidence. Use a CONSERVATIVE fixed
# proxy (mirrors the spatial-axis proxy convention in this script): 0.01 for a
# MASLD-disease-GWAS-driven sQTL hit, 0.05 for an enzyme-proxy-only hit.
p_sqtl <- rep(NA_real_, n)
if ("sqtl_credset_hit" %in% names(atlas)) {
  hit    <- atlas$sqtl_credset_hit %in% TRUE
  driven <- if ("sqtl_masld_gwas_driven" %in% names(atlas)) atlas$sqtl_masld_gwas_driven %in% TRUE else rep(FALSE, n)
  p_sqtl[hit]            <- 0.05
  p_sqtl[hit & driven]  <- 0.01
}
cat(sprintf("  S7 sQTL (conservative proxy): %d non-NA p-values\n", sum(!is.na(p_sqtl))))

# ============================================================================
# 3. FISHER'S COMBINED PROBABILITY TEST
# ============================================================================
cat("\n--- 3. Fisher's combined test ---\n")

# Build p-value matrix (n x 7: 4 original + 3 new orthogonal axes)
p_mat <- cbind(
  human_bulk  = p_human,
  mouse_bulk  = p_mouse,
  genetics    = p_genetic,
  spatial     = p_spatial,
  mqtl        = p_mqtl,
  burden      = p_burden,
  sqtl        = p_sqtl
)

# Floor all p-values to avoid -Inf in log
p_mat[!is.na(p_mat) & p_mat < .Machine$double.eps] <- .Machine$double.eps

# Count non-NA sources per gene
n_sources <- rowSums(!is.na(p_mat))

# Source mask (binary string indicating which sources are present)
source_mask <- apply(!is.na(p_mat), 1, function(x) {
  paste(as.integer(x), collapse = "")
})

# Fisher's chi-squared: -2 * sum(ln(p_i))
# Only for genes with at least 1 source
fisher_chisq <- rep(NA_real_, n)
fisher_pval  <- rep(NA_real_, n)

has_sources <- n_sources >= 1

log_p_mat <- log(p_mat)
log_p_mat[is.na(log_p_mat)] <- 0  # NA sources contribute 0

fisher_chisq[has_sources] <- -2 * rowSums(log_p_mat[has_sources, , drop = FALSE])

# Degrees of freedom = 2 * k (number of non-NA sources)
df <- 2 * n_sources

# P-value from chi-squared distribution
fisher_pval[has_sources] <- pchisq(fisher_chisq[has_sources],
                                    df = df[has_sources],
                                    lower.tail = FALSE)

# BH correction (only among genes with at least 1 source)
fisher_padj <- rep(NA_real_, n)
fisher_padj[has_sources] <- p.adjust(fisher_pval[has_sources], method = "BH")

cat(sprintf("  Genes with >=1 source: %d\n", sum(n_sources >= 1)))
cat(sprintf("  Genes with >=2 sources: %d\n", sum(n_sources >= 2)))
cat(sprintf("  Genes with >=3 sources: %d\n", sum(n_sources >= 3)))
cat(sprintf("  Genes with >=4 sources: %d (max observed: %d of %d axes)\n",
            sum(n_sources >= 4), max(n_sources), ncol(p_mat)))
cat(sprintf("  Fisher padj < 0.05: %d\n", sum(fisher_padj < 0.05, na.rm = TRUE)))
cat(sprintf("  Fisher padj < 0.01: %d\n", sum(fisher_padj < 0.01, na.rm = TRUE)))
cat(sprintf("  Fisher padj < 0.001: %d\n", sum(fisher_padj < 0.001, na.rm = TRUE)))

# ============================================================================
# 4. RANK COMPARISON WITH 46d HEURISTIC
# ============================================================================
cat("\n--- 4. Rank comparison with 46d heuristic ---\n")

# Fisher rank (lower p = better rank = rank 1)
fisher_rank <- rep(NA_integer_, n)
ok <- !is.na(fisher_pval)
fisher_rank[ok] <- rank(fisher_pval[ok], ties.method = "first")

# Match 46d ranks (by gene symbol)
heuristic_rank <- rep(NA_integer_, n)
m46d <- match(atlas$human_symbol, conv$human_symbol)
heuristic_rank <- conv$convergence_rank[m46d]
heuristic_score <- conv$convergence_score[m46d]

# Spearman correlation (both ranked, overlapping set)
both_ranked <- !is.na(fisher_rank) & !is.na(heuristic_rank)
if (sum(both_ranked) > 100) {
  rho <- cor(fisher_rank[both_ranked], heuristic_rank[both_ranked],
             method = "spearman", use = "complete.obs")
  rho_test <- cor.test(fisher_rank[both_ranked], heuristic_rank[both_ranked],
                        method = "spearman")
  cat(sprintf("  Spearman(Fisher rank, 46d heuristic rank): rho = %.4f (p = %.2e, n = %d)\n",
              rho, rho_test$p.value, sum(both_ranked)))
} else {
  rho <- NA_real_
  cat("  Too few overlapping ranked genes for correlation\n")
}

# Top-100 overlap
fisher_top100 <- atlas$human_symbol[which(fisher_rank <= 100)]
heur_top100   <- conv$human_symbol[which(conv$convergence_rank <= 100)]
overlap_100   <- length(intersect(fisher_top100, heur_top100))
cat(sprintf("  Top-100 overlap: %d / 100 (Jaccard = %.3f)\n",
            overlap_100, overlap_100 / length(union(fisher_top100, heur_top100))))

# Top-500 overlap
fisher_top500 <- atlas$human_symbol[which(fisher_rank <= 500)]
heur_top500   <- conv$human_symbol[which(conv$convergence_rank <= 500)]
overlap_500   <- length(intersect(fisher_top500, heur_top500))
cat(sprintf("  Top-500 overlap: %d / 500 (Jaccard = %.3f)\n",
            overlap_500, overlap_500 / length(union(fisher_top500, heur_top500))))

# ============================================================================
# 5. OUTPUT
# ============================================================================
cat("\n--- 5. Write output ---\n")

out <- data.table(
  gene           = atlas$human_symbol,
  ensembl_id     = atlas$ensembl_id,
  n_sources      = n_sources,
  fisher_chisq   = round(fisher_chisq, 4),
  fisher_pval    = fisher_pval,
  fisher_padj    = fisher_padj,
  fisher_rank    = fisher_rank,
  source_mask    = source_mask,
  p_human_bulk   = p_human,
  p_mouse_bulk   = p_mouse,
  p_genetics     = p_genetic,
  p_spatial      = p_spatial,
  heuristic_rank = heuristic_rank,
  heuristic_score = heuristic_score,
  spearman_rho_fisher_vs_heuristic = rho
)

out_path <- file.path(ME, "convergence_evidence_fisher.csv")
fwrite(out, out_path)
cat(sprintf("  Wrote %s (%d rows, %d cols)\n", out_path, nrow(out), ncol(out)))

# Top-20 by Fisher rank
cat("\nTop 20 by Fisher rank:\n")
top20 <- out[!is.na(fisher_rank)][order(fisher_rank)][1:20,
  .(fisher_rank, gene, n_sources, fisher_chisq, fisher_padj,
    source_mask, heuristic_rank)]
print(top20)

# Fisher-specific discoveries (top-100 Fisher but NOT top-100 heuristic)
fisher_only <- setdiff(fisher_top100, heur_top100)
cat(sprintf("\nFisher top-100 unique (not in heuristic top-100): %d genes\n",
            length(fisher_only)))
if (length(fisher_only) > 0) {
  fisher_only_dt <- out[gene %in% fisher_only][order(fisher_rank)][1:min(20, length(fisher_only)),
    .(fisher_rank, gene, n_sources, fisher_padj, heuristic_rank)]
  print(fisher_only_dt)
}

# ============================================================================
# 6. COLOC ENZYME WEIGHTING SENSITIVITY (Task 3)
# ============================================================================
cat("\n--- 6. COLOC enzyme weighting sensitivity ---\n")
# Compare convergence scores: (a) with 0.5x enzyme weighting vs (b) without.
# Since 46d was edited to remove the weighting (task 1), we simulate both
# configurations here by re-computing the COLOC channel with/without penalty.

# Identify COLOC columns
coloc_cols <- grep("_coloc_pp4$", names(atlas), value = TRUE)
extra_cols <- intersect(c("coloc_susie_best_pp4","coloc_abf_best_pp4",
                           "broadaway_coloc_pp4","sceqtl_coloc_best_pp4"),
                        names(atlas))
coloc_cols <- union(coloc_cols, extra_cols)

pp4_mat_raw <- as.matrix(atlas[, ..coloc_cols])

# (a) WITH 0.5x enzyme down-weighting (old behavior)
enzyme_pattern <- "(?i)(alt|ast|ggt|pdff)_coloc_pp4$"
col_weights_old <- ifelse(grepl(enzyme_pattern, coloc_cols, perl = TRUE), 0.5, 1.0)
pp4_weighted <- sweep(pp4_mat_raw, 2, col_weights_old, "*")
pp4_max_weighted <- apply(pp4_weighted, 1, function(x)
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))

# (b) WITHOUT enzyme down-weighting (new behavior = equal weights)
pp4_max_equal <- apply(pp4_mat_raw, 1, function(x)
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))

# COLOC BF conversion (same as 46d: coloc_posterior_upgrade)
COLOC_P12 <- 1e-5
COLOC_P2  <- 1e-4
COLOC_PRIOR_ODDS <- COLOC_P12 / COLOC_P2
LOG_BF_FLOOR <- log(0.1)
LOG_BF_CEIL  <- log(1e6)

coloc_bf <- function(pp4) {
  out <- rep(0, length(pp4))
  ok <- !is.na(pp4) & pp4 > 0 & pp4 < 1
  posterior_odds <- pp4[ok] / (1 - pp4[ok])
  out[ok] <- log(posterior_odds / COLOC_PRIOR_ODDS)
  out[!is.na(pp4) & pp4 >= 1] <- LOG_BF_CEIL
  pmin(pmax(out, LOG_BF_FLOOR), LOG_BF_CEIL)
}

bf_weighted <- coloc_bf(pp4_max_weighted)
bf_equal    <- coloc_bf(pp4_max_equal)

# How many genes are affected?
affected <- !is.na(pp4_max_weighted) & !is.na(pp4_max_equal) &
            abs(pp4_max_weighted - pp4_max_equal) > 1e-10
cat(sprintf("  Genes affected by enzyme weighting: %d / %d tested\n",
            sum(affected), sum(!is.na(pp4_max_equal))))
cat(sprintf("  Enzyme-weighted COLOC columns: %d / %d\n",
            sum(col_weights_old == 0.5), length(col_weights_old)))

# Rank correlation of the BF channel
both_bf <- is.finite(bf_weighted) & is.finite(bf_equal) &
           (bf_weighted != 0 | bf_equal != 0)
if (sum(both_bf) > 100) {
  rho_bf <- cor(bf_weighted[both_bf], bf_equal[both_bf],
                method = "spearman", use = "complete.obs")
  cat(sprintf("  Spearman(BF_weighted, BF_equal): rho = %.4f (n = %d)\n",
              rho_bf, sum(both_bf)))
} else {
  rho_bf <- NA_real_
}

# PP4 max correlation
both_pp4 <- !is.na(pp4_max_weighted) & !is.na(pp4_max_equal)
if (sum(both_pp4) > 100) {
  rho_pp4 <- cor(pp4_max_weighted[both_pp4], pp4_max_equal[both_pp4],
                 method = "spearman", use = "complete.obs")
  cat(sprintf("  Spearman(PP4_max_weighted, PP4_max_equal): rho = %.4f (n = %d)\n",
              rho_pp4, sum(both_pp4)))
} else {
  rho_pp4 <- NA_real_
}

# Per-gene comparison
sens_out <- data.table(
  gene              = atlas$human_symbol,
  ensembl_id        = atlas$ensembl_id,
  pp4_max_equal     = round(pp4_max_equal, 6),
  pp4_max_weighted  = round(pp4_max_weighted, 6),
  pp4_delta         = round(pp4_max_equal - pp4_max_weighted, 6),
  bf_equal          = round(bf_equal, 4),
  bf_weighted       = round(bf_weighted, 4),
  bf_delta          = round(bf_equal - bf_weighted, 4),
  affected          = affected
)

# Rank both configurations and compare
rank_eq <- rep(NA_integer_, n)
rank_wt <- rep(NA_integer_, n)
ok_eq <- bf_equal != 0 | !is.na(pp4_max_equal)
ok_wt <- bf_weighted != 0 | !is.na(pp4_max_weighted)
rank_eq[ok_eq] <- rank(-bf_equal[ok_eq], ties.method = "first")
rank_wt[ok_wt] <- rank(-bf_weighted[ok_wt], ties.method = "first")
sens_out[, rank_equal := rank_eq]
sens_out[, rank_weighted := rank_wt]
sens_out[, rank_delta := rank_wt - rank_eq]  # positive = promoted by removing weighting

# Summary statistics
summary_row <- data.table(
  metric = c("n_genes_affected", "n_coloc_cols_enzyme", "n_coloc_cols_total",
             "spearman_bf", "spearman_pp4_max",
             "mean_pp4_delta_affected", "max_pp4_delta",
             "n_promoted_top100", "n_demoted_top100"),
  value  = c(sum(affected),
             sum(col_weights_old == 0.5),
             length(col_weights_old),
             round(rho_bf, 4),
             round(rho_pp4, 4),
             round(mean(sens_out$pp4_delta[affected], na.rm = TRUE), 6),
             round(max(sens_out$pp4_delta, na.rm = TRUE), 6),
             sum(sens_out$rank_delta > 0 & rank_eq <= 100, na.rm = TRUE),
             sum(sens_out$rank_delta < 0 & rank_eq <= 100, na.rm = TRUE))
)

sens_path <- file.path(ME, "coloc_enzyme_weighting_sensitivity.csv")
fwrite(sens_out[order(-abs(pp4_delta))], sens_path)
cat(sprintf("  Wrote %s (%d rows)\n", sens_path, nrow(sens_out)))

# Also write summary
summ_path <- file.path(ME, "coloc_enzyme_weighting_sensitivity_summary.csv")
fwrite(summary_row, summ_path)
cat(sprintf("  Wrote %s\n", summ_path))

cat("\nSensitivity summary:\n")
print(summary_row)

# Top affected genes (largest PP4 change)
cat("\nTop 20 most affected genes by enzyme weighting removal:\n")
top_affected <- sens_out[affected == TRUE][order(-pp4_delta)][1:min(20, sum(affected)),
  .(gene, pp4_max_equal, pp4_max_weighted, pp4_delta, rank_equal, rank_weighted, rank_delta)]
print(top_affected)

cat("\nDone.\n")
