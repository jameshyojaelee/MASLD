#!/usr/bin/env Rscript
# sex_v3/04_dream_M2_extended.R
# ---------------------------------------------------------------------------
# Module 04 — limma-voom-quality-weighted (LVQW) interaction fit with extended
# covariate set, repeated over 5 mice age imputations, pooled via Rubin's rules.
#
# LVQW re-engineering 2026-06-08: dataset random->fixed.
#   ENGINE swap: voomWithDreamWeights + dream (variancePartition mixed model,
#   random slope) -> voomWithQualityWeights(dge, design) + eBayes(lmFit(v, design)).
#   The dataset random-effect structure `(1 + group_binary | dataset)` is
#   replaced by `dataset` as a FIXED effect. Everything else (interaction,
#   composition, MI age, SVA surrogate variables) is preserved EXACTLY.
#
# Per-imputation formula (design matrix):
#   ~ dataset                                                     (FIXED cohort effect)
#     + group_binary * inferred_sex                              (main effects + interaction)
#     + Hepatocytes + Macrophages + Endothelial + Cholangiocytes (composition)
#     + age_imputed                                              (MI age)
#     + SV1 + ... + SVk                                          (SVA hidden factors, k from Module 02)
#
# RAW vs MODERATED variance (match-the-original): the dream path consumed RAW
# (unmoderated) Wald SEs for Rubin pooling (moderation belongs downstream — mash
# provides its own shrinkage). We preserve that: per-gene Var(beta) and
# Cov(beta_F, beta_int) are built from the UNMODERATED sigma^2 (fit$sigma^2),
# NOT eBayes's moderated s2.post.
#
# Coefficient parameterization (level alphabetics — sex_F is reference):
#   coef = "group_binaryDisease"                       => β_F (disease-in-F)
#   coef = "group_binaryDisease:inferred_sexM"         => β_int (M minus F effect)
#   β_M  = β_F + β_int   (delta sum)
#
# Per-gene outputs (pooled across 5 MI fits with Rubin's rules):
#   β_F, se_F, β_M, se_M, β_int, se_int, t_int, padj_BH
#
# Reads:  intermediates/sex_v3_input.rds
#         intermediates/sva_factors.rds
#         intermediates/age_mi.rds
# Writes: dream_M2_v3.rds         (full pooled fit + vcov per gene; filename kept
#                                  for downstream-path compatibility with 05/06)
#         dream_topTable_v3.csv   (flat per-gene table)
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

# LVQW re-engineering 2026-06-08: dataset random->fixed.
# Engine is now limma voomWithQualityWeights + lmFit + eBayes (fixed-effects),
# so the lme4/reformulas/variancePartition stack and its namespace injection are
# no longer needed. limma carries voomWithQualityWeights + lmFit + eBayes.
suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(edgeR)
})

# Record limma version (LVQW engine) for downstream provenance.
limma_ver <- packageVersion("limma")
cat("limma version:", as.character(limma_ver), "\n")
if (limma_ver < "3.50") {
  stop("limma >= 3.50 required for voomWithQualityWeights + cov.coefficients; have ",
       as.character(limma_ver))
}

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RDIR   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                    "results/integration")
SEXV3  <- file.path(RDIR, "sex_v3")
IDIR   <- file.path(SEXV3, "intermediates")
GENE_META <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")

# Shared utilities (atomic writes + sessionInfo dump) -- R5 Issue 2 fix
source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

cat("============================================================\n")
cat("Module 04 — Dream M2 extended (random-slope, 5 MI imputations)\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")

# ---------------------------------------------------------------------------
# Parallel
# ---------------------------------------------------------------------------
# LVQW re-engineering 2026-06-08: dataset random->fixed.
# limma's voomWithQualityWeights + lmFit + eBayes are single-threaded matrix
# algebra (no BiocParallel dispatch like dream). We retain the cpu count for
# log parity but no longer build a BiocParallel param object.
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
cat("Using", ncpus, "CPU core(s) reported by SLURM (limma LVQW is single-threaded)\n")

# ---------------------------------------------------------------------------
# Load intermediates
# ---------------------------------------------------------------------------
inp     <- readRDS(file.path(IDIR, "sex_v3_input.rds"))
sva_o   <- readRDS(file.path(IDIR, "sva_factors.rds"))
age_mi  <- readRDS(file.path(IDIR, "age_mi.rds"))

dge   <- inp$dge
info0 <- inp$meta
sv_mat <- sva_o$sv
n_sv   <- ncol(sv_mat)
m      <- length(age_mi$age_list)
stopifnot(m == 5L)
stopifnot(nrow(sv_mat) == nrow(info0))
stopifnot(identical(rownames(sv_mat), rownames(info0)))
cat("Samples:", ncol(dge), "  Genes:", nrow(dge), "\n")
cat("SVs:", n_sv, "  MI imputations:", m, "\n")

# ---------------------------------------------------------------------------
# Confirm coefficient parameterization (R1 Issue 6 / Major fix)
# ---------------------------------------------------------------------------
# Enforce strict factor levels so β_F / β_M column assignment in Module 05
# (which hard-codes "F"/"M") is never silently swapped by upstream label drift.
if (!"inferred_sex" %in% names(info0)) {
  stop("info0$inferred_sex missing — Module 01 should produce this column.")
}
stopifnot(is.factor(info0$inferred_sex))
sex_levels <- levels(info0$inferred_sex)
cat("inferred_sex factor levels (ref =", sex_levels[1], "):",
    paste(sex_levels, collapse = ", "), "\n")
stopifnot(identical(sex_levels, c("F", "M")))
stopifnot(all(as.character(info0$inferred_sex) %in% c("F", "M")))
message("Confirmed: inferred_sex levels are c('F','M'); F is reference for β_F derivation")

# ---------------------------------------------------------------------------
# Build per-imputation info data frames
# ---------------------------------------------------------------------------
# Combine info + SV + age per imputation. SVs are constant across imputations;
# only age_imputed differs.
sv_df <- as.data.frame(sv_mat)
colnames(sv_df) <- paste0("SV", seq_len(n_sv))

build_info_k <- function(k) {
  out <- info0
  out$age_imputed <- age_mi$age_list[[k]]
  out <- cbind(out, sv_df)
  rownames(out) <- rownames(info0)
  out
}

# LVQW re-engineering 2026-06-08: dataset random->fixed.
# The dream-era random slope `(1 + group_binary | dataset)` (R1 Issue 4: over-
# parameterized for k=5 cohorts -> singular fits, hence the slope->intercept
# fallback ladder) is replaced by `dataset` as a leading FIXED effect. A fixed
# design is always full-rank here (5 cohorts -> 4 dataset dummy columns) so there
# is NO fallback ladder: a single design is built per imputation.
sv_terms <- paste(colnames(sv_df), collapse = " + ")
# diagnosis_harmonized REMOVED 2026-05-14: was collinear with group_binary
# (318/645 Disease samples have diagnosis_harmonized=="" which perfectly
# predicts group_binary=="Disease"; collinearity flipped β_F sign for any
# gene whose disease effect was sensitive to the coding). The within-Disease
# severity question (NAFL/Borderline/NASH) is a separate analysis that needs
# a properly coded variable and a non-collinear design.
form_str_fixed <- paste0(
  "~ dataset",
  " + group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed",
  " + ", sv_terms
)
form_fixed <- as.formula(form_str_fixed)
cat("\nFormula (LVQW fixed-effects design, dataset fixed):\n  ", form_str_fixed, "\n")

# Track which variant was used per imputation (for downstream provenance).
# LVQW has a single variant; kept as a vector for rds-schema parity.
variant_used <- rep("lvqw_fixed", m)

# ---------------------------------------------------------------------------
# Helper: fit LVQW for one imputation (single fixed-effects design, no fallback)
# LVQW re-engineering 2026-06-08: dataset random->fixed.
# ---------------------------------------------------------------------------
fit_one_imputation <- function(k) {
  cat("\n========== Imputation", k, "of", m, "==========\n")
  info_k <- build_info_k(k)

  # Build the model.matrix from the fixed-effects formula. dataset is a leading
  # fixed effect (4 dummy columns for 5 cohorts). All other terms are unchanged.
  design <- model.matrix(form_fixed, data = info_k)
  stopifnot(nrow(design) == ncol(dge))
  qrr <- qr(design)
  if (qrr$rank < ncol(design)) {
    stop("LVQW design is rank-deficient (rank ", qrr$rank, " < ncol ",
         ncol(design), "); inspect collinearity among dataset/SV/composition terms.")
  }

  cat("  voomWithQualityWeights (fixed-effects design)...\n")
  t0 <- Sys.time()
  v <- suppressWarnings(voomWithQualityWeights(dge, design))
  cat("    voom elapsed:", format(Sys.time() - t0), "\n")

  cat("  lmFit + eBayes...\n")
  t0 <- Sys.time()
  fit <- eBayes(lmFit(v, design))
  cat("    lmFit+eBayes elapsed:", format(Sys.time() - t0), "\n")

  variant <- "lvqw_fixed"
  variant_used[k] <<- variant
  cat("  variant used for imp", k, ":", variant, "\n")

  # Locate coefficient names. With factor(c("F","M")), expect:
  #   "group_binaryDisease"   -- main effect (= β_F)
  #   "group_binaryDisease:inferred_sexM" or "inferred_sexM:group_binaryDisease"
  all_coefs <- colnames(fit$coefficients)
  group_coef <- grep("^group_binary", all_coefs, value = TRUE)
  group_coef <- setdiff(group_coef, grep(":", group_coef, value = TRUE))
  if (length(group_coef) != 1) {
    stop("Expected exactly one main effect coef matching ^group_binary; got: ",
         paste(group_coef, collapse = ", "))
  }
  int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary",
                   all_coefs, value = TRUE)
  if (length(int_coef) == 0) {
    stop("No interaction coefficient found among: ", paste(all_coefs, collapse = ", "))
  }
  int_coef <- int_coef[1]
  cat("  group coef (= β_F):", group_coef, "\n")
  cat("  interaction coef    :", int_coef, "\n")

  # Extract β_F, β_int per gene
  beta_F   <- fit$coefficients[, group_coef]
  beta_int <- fit$coefficients[, int_coef]

  # ---- RAW (unmoderated) Wald SEs — match the dream-era choice (R1 Issue 2).
  # LVQW re-engineering 2026-06-08: dataset random->fixed.
  # lmFit's per-gene Var(beta) = (stdev.unscaled[gene, coef])^2 * sigma[gene]^2
  # where sigma is the UNMODERATED residual SD (we deliberately do NOT use the
  # eBayes-moderated s2.post, mirroring the dream path which fed RAW Wald SEs to
  # Rubin pooling — moderation belongs downstream, mash provides its own shrink).
  stopifnot("stdev.unscaled" %in% names(fit), "sigma" %in% names(fit))
  raw_se_mat <- fit$stdev.unscaled * fit$sigma   # genes x K, recycles sigma per row
  se_F   <- raw_se_mat[, group_coef]
  se_int <- raw_se_mat[, int_coef]
  names(se_F)   <- rownames(fit$coefficients)
  names(se_int) <- rownames(fit$coefficients)
  se_F   <- se_F[names(beta_F)]
  se_int <- se_int[names(beta_int)]

  # ---- Per-gene Cov(β_F, β_int) — R1 Issue 1 (omitting 2*Cov biases Var(β_M)).
  # LVQW re-engineering 2026-06-08: dataset random->fixed.
  # dream stored a FULL per-gene vcov in cov.coefficients.list; lmFit instead
  # returns a SINGLE shared UNSCALED coefficient covariance (fit$cov.coefficients,
  # a K x K named matrix = (X'WX)^-1) plus a per-gene fit$sigma. The per-gene
  # raw covariance is therefore:
  #   Cov_g(β_F, β_int) = cov.coefficients[group_coef, int_coef] * sigma[g]^2
  # This is exactly consistent with the diagonal SEs above, since
  #   stdev.unscaled[g, coef]^2 == cov.coefficients[coef, coef]   (shared across g)
  # so Var_g(β) = cov.coefficients[coef,coef] * sigma[g]^2. Under treatment-coded
  # interaction the unscaled off-diagonal is strongly NEGATIVE (see the toy QR
  # example: Cov(gD, gD:sM) = -1 of the unscaled matrix), so Cov_g is negative and
  # correctly shrinks Var(β_M) = Var(β_F) + Var(β_int) + 2*Cov.
  stopifnot("cov.coefficients" %in% names(fit))
  ccm <- fit$cov.coefficients
  if (is.null(ccm) || !is.matrix(ccm) ||
      !(group_coef %in% rownames(ccm)) || !(int_coef %in% colnames(ccm))) {
    stop("fit$cov.coefficients missing or coef names not found: ",
         paste(rownames(ccm), collapse = ", "))
  }
  cov_unscaled_F_int <- as.numeric(ccm[group_coef, int_coef])
  sigma2_g <- (fit$sigma)^2
  names(sigma2_g) <- rownames(fit$coefficients)
  sigma2_g <- sigma2_g[names(beta_F)]
  cov_F_int <- cov_unscaled_F_int * sigma2_g
  names(cov_F_int) <- names(beta_F)
  if (any(is.na(cov_F_int))) {
    stop(sum(is.na(cov_F_int)),
         " genes have missing Cov(β_F, β_int) — check coef name lookup.")
  }
  # Diagnostic on Cov-derived correlation (should be strongly negative for this design)
  rho_F_int <- cov_F_int / pmax(se_F * se_int, .Machine$double.eps)
  rho_F_int[!is.finite(rho_F_int)] <- NA_real_
  cat(sprintf("  Cov(β_F, β_int) summary: median rho = %.3f (q25=%.3f, q75=%.3f)\n",
              median(rho_F_int, na.rm = TRUE),
              quantile(rho_F_int, 0.25, na.rm = TRUE),
              quantile(rho_F_int, 0.75, na.rm = TRUE)))
  # Per-gene t_int + p_int from raw Wald (NOT moderated) so it is consistent
  # with the raw SEs we pass to Rubin's pooling.
  t_int_raw <- beta_int / pmax(se_int, .Machine$double.eps)
  p_int_raw <- 2 * pnorm(-abs(t_int_raw))

  list(
    beta_F     = beta_F,
    beta_int   = beta_int,
    se_F       = se_F,
    se_int     = se_int,
    cov_F_int  = cov_F_int,
    t_int      = t_int_raw,
    p_int      = p_int_raw,
    gene       = names(beta_F),
    group_coef = group_coef,
    int_coef   = int_coef,
    variant    = variant
  )
}

# ---------------------------------------------------------------------------
# Run 5 fits sequentially (each uses ncpus internally)
# ---------------------------------------------------------------------------
fits <- vector("list", m)
for (k in seq_len(m)) {
  fits[[k]] <- fit_one_imputation(k)
}

# Validate row alignment across imputations
gene_names <- fits[[1]]$gene
for (k in seq_len(m)) {
  stopifnot(identical(fits[[k]]$gene, gene_names))
}
n_genes <- length(gene_names)
cat("\nAligned across MI fits:", n_genes, "genes\n")

# ---------------------------------------------------------------------------
# Rubin's rules pooling
#   pooled β  = mean(within-imputation β)
#   within-var W = mean(within-imp var)
#   between-var B = var(within-imp β, ddof = m-1)
#   pooled var T = W + (1 + 1/m) * B
#
# CRITICAL invariant (R1 Issue 2 fix):
#   - Within-imp SE is the RAW Wald SE = sigma * stdev.unscaled (NOT topTable's
#     moderated logFC/t quotient). Rubin's pooling MUST consume frequentist
#     sampling variances; moderation belongs downstream (mash provides its own
#     shrinkage).
#
# For β_M = β_F + β_int, the per-imp variance is the linear-combination rule:
#   Var(β_F + β_int) = Var(β_F) + Var(β_int) + 2*Cov(β_F, β_int)
# where Cov(β_F, β_int) is extracted from fit$cov.coefficients.list[[g]]
# (R1 Issue 1 fix). The previous code assumed Cov = 0 which over-estimated
# Var(β_M) by ~30-60% under the treatment-coded interaction.
# ---------------------------------------------------------------------------
pool_rubin <- function(beta_mat, se_mat, m) {
  # beta_mat: n_genes x m; se_mat: same shape
  pooled_beta <- rowMeans(beta_mat)
  W <- rowMeans(se_mat^2)               # within-imp variance
  B <- apply(beta_mat, 1, var)          # between-imp variance (n-1 / ddof = m-1)
  Tvar <- W + (1 + 1 / m) * B
  list(beta = pooled_beta, var = Tvar, W = W, B = B)
}

beta_F_mat    <- sapply(fits, function(x) x$beta_F)
beta_int_mat  <- sapply(fits, function(x) x$beta_int)
se_F_mat      <- sapply(fits, function(x) x$se_F)
se_int_mat    <- sapply(fits, function(x) x$se_int)
cov_Fint_mat  <- sapply(fits, function(x) x$cov_F_int)   # n_genes x m

stopifnot(nrow(beta_F_mat)   == n_genes, ncol(beta_F_mat)   == m)
stopifnot(nrow(beta_int_mat) == n_genes, ncol(beta_int_mat) == m)
stopifnot(nrow(cov_Fint_mat) == n_genes, ncol(cov_Fint_mat) == m)

pool_F   <- pool_rubin(beta_F_mat,   se_F_mat,   m)
pool_int <- pool_rubin(beta_int_mat, se_int_mat, m)

# β_M = β_F + β_int (per imputation), pool the linear combination with
# proper within-imp variance Var(β_F + β_int) = SE_F^2 + SE_int^2 + 2 Cov.
# (R1 Issue 1 fix: previously this set Cov(β_F, β_int) = 0.)
beta_M_mat <- beta_F_mat + beta_int_mat
se_M_mat   <- sqrt(pmax(se_F_mat^2 + se_int_mat^2 + 2 * cov_Fint_mat, 0))
# Floor at machine eps to avoid sqrt(<0) numerical fallout when cov is very
# negative and dominates the variance (rare in practice, but guard regardless).
pool_M     <- pool_rubin(beta_M_mat, se_M_mat, m)

# Pooled Cov(β_F, β_int) across imputations. Each imputation contributes one
# per-gene covariance value; Rubin pooling on a scalar covariance reduces to
# the mean of the per-imp covariances (the between-imp variance of a single
# scalar across imputations is captured implicitly in pool_F$B / pool_int$B).
pool_cov_F_int  <- rowMeans(cov_Fint_mat)

# Test statistic for the interaction
t_int <- pool_int$beta / sqrt(pool_int$var)
# Approximate Barnard-Rubin DF; fall back to a large df for simple Wald-like
# p-value. (Use normal approximation — conservative since for large n the
# Barnard-Rubin DF is very large; for small effective DF this slightly
# inflates significance. Use simple normal Wald p-value here for
# downstream BH; calibration module re-checks with simulation.)
p_int <- 2 * pnorm(-abs(t_int))
padj_int <- p.adjust(p_int, method = "BH")

cat("\nInteraction p-value distribution:\n")
cat("  median p:", round(median(p_int), 4),
    "  min p:", signif(min(p_int), 3), "\n")
cat("  n padj < 0.05:", sum(padj_int < 0.05),
    "  n padj < 0.1:",  sum(padj_int < 0.1), "\n")

# ---------------------------------------------------------------------------
# Build per-gene 3x3 vcov (β_F, β_M, β_int) — R1 Issue 1 fix.
# Use the empirical Cov(β_F, β_int) pooled from fit$cov.coefficients.list
# (NOT zero as the old build_vcov assumed). Linear-combination rules:
#   Var(β_F)        = pool_F$var
#   Var(β_int)      = pool_int$var
#   Cov(β_F, β_int) = pool_cov_F_int  (typically strongly negative)
#   β_M = β_F + β_int  →
#     Var(β_M)        = Var(β_F) + Var(β_int) + 2*Cov(β_F, β_int)   = pool_M$var
#     Cov(β_F, β_M)   = Var(β_F) + Cov(β_F, β_int)
#     Cov(β_int, β_M) = Cov(β_F, β_int) + Var(β_int)
# Stored as per-gene 3x3 array indexed by (F, M, int).
# ---------------------------------------------------------------------------
build_vcov <- function(vF, vM, vInt, cFi) {
  c_F_M   <- vF   + cFi
  c_int_M <- cFi  + vInt
  array(c(
    vF,    c_F_M,   cFi,     # row 1: F
    c_F_M, vM,      c_int_M, # row 2: M
    cFi,   c_int_M, vInt     # row 3: int
  ), dim = c(3, 3))
}
vcov_arr <- array(NA_real_, dim = c(n_genes, 3, 3),
                  dimnames = list(gene_names,
                                  c("F","M","int"),
                                  c("F","M","int")))
for (i in seq_len(n_genes)) {
  vcov_arr[i, , ] <- build_vcov(pool_F$var[i],
                                pool_M$var[i],
                                pool_int$var[i],
                                pool_cov_F_int[i])
}

# ---------------------------------------------------------------------------
# Gene symbols (gencode_v49)
# ---------------------------------------------------------------------------
cat("\nLoading gene_id -> gene_name mapping...\n")
gene_meta <- fread(GENE_META)
sym_lookup <- setNames(gene_meta$gene_name, gene_meta$gene_id)
gene_symbol <- sym_lookup[gene_names]
if (any(is.na(gene_symbol))) {
  cat("  Warning:", sum(is.na(gene_symbol)),
      "genes without symbol — using gene_id as symbol\n")
  gene_symbol[is.na(gene_symbol)] <- gene_names[is.na(gene_symbol)]
}

# ---------------------------------------------------------------------------
# Flat output table
# ---------------------------------------------------------------------------
out_dt <- data.table(
  gene        = gene_names,
  gene_symbol = unname(gene_symbol),
  beta_F      = pool_F$beta,
  se_F        = sqrt(pool_F$var),
  beta_M      = pool_M$beta,
  se_M        = sqrt(pool_M$var),
  beta_int    = pool_int$beta,
  se_int      = sqrt(pool_int$var),
  cov_F_int   = pool_cov_F_int,
  t_int       = t_int,
  p_int       = p_int,
  padj_BH     = padj_int,
  W_F         = pool_F$W,
  B_F         = pool_F$B,
  W_int       = pool_int$W,
  B_int       = pool_int$B
)

out_csv <- file.path(SEXV3, "dream_topTable_v3.csv")
write_atomic_csv(out_dt, out_csv)
cat("\nSaved (atomic) flat table:", out_csv, " (", nrow(out_dt), "rows )\n")

# ---------------------------------------------------------------------------
# R3 Issue 2 fix: emit `sex_interaction_dream_v3.csv` — a v2-schema-compatible
# interaction stats file for 27a (line 981 sex_int_layer), 218c pathway
# gsea, and figS_sex_dimorphism panel E. The v2 file shipped columns
# {logFC, AveExpr, t, P.Value, padj, z.std, gene}; we mirror the
# downstream-consumed columns (gene, logFC, t, P.Value, padj) and also
# expose adj.P.Val (limma-style) so newer v3-aware consumers have both
# naming conventions available.
# ---------------------------------------------------------------------------
sex_int_v3 <- data.table(
  logFC     = pool_int$beta,
  AveExpr   = NA_real_,                     # v2 schema slot (not used downstream)
  t         = t_int,
  P.Value   = p_int,
  padj      = padj_int,                     # v2 schema name (read by 27a:988)
  adj.P.Val = padj_int,                     # limma alias for v3-aware consumers
  z.std     = NA_real_,                     # v2 schema slot
  gene      = gene_names,
  se_int    = sqrt(pool_int$var)            # bonus column; downstream is free to ignore
)
sex_int_v3_csv <- file.path(SEXV3, "sex_interaction_dream_v3.csv")
write_atomic_csv(sex_int_v3, sex_int_v3_csv)
cat("Saved (atomic) sex_interaction_dream_v3:", sex_int_v3_csv,
    "  rows:", nrow(sex_int_v3), "\n")
# Also stash a copy under intermediates/ so downstream agents that scan
# intermediates can pick it up without having to know the public SEXV3 path.
sex_int_v3_csv_idir <- file.path(IDIR, "sex_interaction_dream_v3.csv")
write_atomic_csv(sex_int_v3, sex_int_v3_csv_idir)
cat("Saved (atomic) sex_interaction_dream_v3 (intermediates):",
    sex_int_v3_csv_idir, "\n")

# ---------------------------------------------------------------------------
# Full pooled fit RDS
# ---------------------------------------------------------------------------
out_rds <- file.path(SEXV3, "dream_M2_v3.rds")
write_atomic_rds(list(
  pooled = list(
    gene_names = gene_names,
    gene_symbol = unname(gene_symbol),
    beta_F   = pool_F$beta,
    beta_M   = pool_M$beta,
    beta_int = pool_int$beta,
    var_F    = pool_F$var,
    var_M    = pool_M$var,
    var_int  = pool_int$var,
    W_F      = pool_F$W,
    B_F      = pool_F$B,
    W_int    = pool_int$W,
    B_int    = pool_int$B,
    t_int    = t_int,
    p_int    = p_int,
    padj_int = padj_int,
    cov_F_int = pool_cov_F_int,
    vcov     = vcov_arr
  ),
  per_imputation = list(
    beta_F   = beta_F_mat,
    beta_int = beta_int_mat,
    beta_M   = beta_M_mat,
    se_F     = se_F_mat,
    se_int   = se_int_mat,
    se_M     = se_M_mat,
    cov_F_int = cov_Fint_mat
  ),
  # LVQW re-engineering 2026-06-08: dataset random->fixed.
  # formula_slope / formula_intercept are retained as KEYS (downstream rds
  # readers that look them up by name keep working) but now both carry the
  # single LVQW fixed-effects design string; engine = "lvqw_fixed".
  formula_slope     = form_str_fixed,
  formula_intercept = form_str_fixed,
  formula_fixed     = form_str_fixed,
  engine            = "lvqw_fixed",
  variant_used      = variant_used,
  group_coef    = fits[[1]]$group_coef,
  int_coef      = fits[[1]]$int_coef,
  sex_reference = sex_levels[1],
  sex_levels    = sex_levels,
  n_sv          = n_sv,
  m_imputations = m,
  vp_version    = as.character(limma_ver)
), out_rds)
cat("Saved (atomic) pooled fit:", out_rds, "\n")

# ---------------------------------------------------------------------------
# Console headline
# ---------------------------------------------------------------------------
cat("\n===== SUMMARY =====\n")
cat("  n_genes              :", n_genes, "\n")
cat("  n_samples            :", ncol(dge), "\n")
cat("  n_SVs                :", n_sv, "\n")
cat("  MI imputations       :", m, "\n")
cat("  median |beta_F|      :", round(median(abs(pool_F$beta), na.rm = TRUE), 4), "\n")
cat("  median |beta_M|      :", round(median(abs(pool_M$beta), na.rm = TRUE), 4), "\n")
cat("  median |beta_int|    :", round(median(abs(pool_int$beta), na.rm = TRUE), 4), "\n")
cat("  beta_int padj < 0.05 :", sum(padj_int < 0.05), "\n")
cat("  beta_int padj < 0.1  :", sum(padj_int < 0.1), "\n")

# R5 Issue 2 fix: sessionInfo dump for per-module reproducibility audit
dump_session_info(IDIR, "04")

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
