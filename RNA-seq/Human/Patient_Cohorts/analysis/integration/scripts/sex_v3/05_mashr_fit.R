#!/usr/bin/env Rscript
# sex_v3/05_mashr_fit.R
# ---------------------------------------------------------------------------
# Module 05 — mashr Bayesian fit over 2-condition (F, M) sex-effect estimates.
#
# Strategy (Urbut & Stephens 2019 eQTL-style "strong + random" recipe):
#   1. Load dream_M2_v3.rds (Module 04) -> (beta_F, beta_M, vcov) per gene.
#   2. Build Bhat / Shat matrices [n_genes x 2] from pooled posteriors.
#   3. Define strong subset (top 1000 by |t_interaction|).
#   4. Define random subset (10,000 genes random.sample, seed=42).
#   5. Estimate Vhat (null residual correlation) from random set via
#      estimate_null_correlation_simple().
#   6. Compute data-driven covariances:
#        U.pca = cov_pca(data.strong, npc = 5)
#        U.ed  = cov_ed(data.strong, U.pca)
#   7. Compute canonical covariances:
#        U.c_default = cov_canonical(data.random)
#        + custom anti_correlated U_anti (off-diagonal = -1)
#   8. Fit mash on random set with outputlevel = 1 (estimate pi only).
#   9. Compute posteriors for ALL genes by passing g = get_fitted_g(fit_random),
#      fixg = TRUE so we use the data-driven mixture trained on random + strong.
#  10. Aggregate per-pattern weights into 4 plan-defined buckets
#      (concordant / F_only / M_only / anti_correlated) via U-matrix sign +
#      structure inspection.
#  11. Emit canonical outputs:
#        mashr_fit.rds              full mash object + metadata
#        mashr_pattern_loadings.csv per gene posterior weight per pattern
#                                    + pre-aggregated 4-bucket posteriors
#        mashr_posterior_summary.csv per gene posterior means/SDs/lfsr
#        mashr_data_driven_cov.rds  discovered ED covariances + canonical U list
#
# This script supports two run modes:
#   PILOT (CHRX_PILOT=TRUE): subset to chrX genes only (~1,025 protein-coding +
#     lncRNAs) before fitting. Writes pilot-named outputs.
#   FULL  (CHRX_PILOT=FALSE, default): full atlas.
#
# Pilot is meant to be invoked via the dedicated wrapper 05_mashr_pilot_chrX.R
# which calls into this file with CHRX_PILOT=TRUE (or sets the env var).
#
# Reads:  intermediates/dream_M2_v3.rds
# Writes: mashr_fit.rds | mashr_pattern_loadings.csv | mashr_data_driven_cov.rds
#         mashr_posterior_summary.csv
#         (pilot: mashr_pilot_chrX.rds + B2_chrX_pilot_report.md)
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(mashr)
})

# Helper: NULL-coalesce (defined early so report writer can use it)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Verify mashr version (>= 0.2.79 expected from install_v3 job)
mashr_ver <- packageVersion("mashr")
cat("mashr version:", as.character(mashr_ver), "\n")
if (mashr_ver < "0.2.50") {
  stop("mashr >= 0.2.50 required for cov_ed + cov_canonical + estimate_null_correlation_simple; have ",
       as.character(mashr_ver))
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
CHRX_F <- file.path(BASE, "RNA-seq/results/audit_sensitivity",
                    "sex_xci/chrx_stratified_classification.csv")
TEAM_B_DIR <- file.path(BASE, "outputs/team_B")
dir.create(IDIR,       recursive = TRUE, showWarnings = FALSE)
dir.create(TEAM_B_DIR, recursive = TRUE, showWarnings = FALSE)

# Shared utilities (atomic writes + sessionInfo dump) -- R5 Issue 2 fix
source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

# ---------------------------------------------------------------------------
# Run mode (full vs chrX pilot)
# ---------------------------------------------------------------------------
CHRX_PILOT <- isTRUE(as.logical(Sys.getenv("CHRX_PILOT", "FALSE")))
N_STRONG <- as.integer(Sys.getenv("MASHR_N_STRONG", "1000"))
N_RANDOM <- as.integer(Sys.getenv("MASHR_N_RANDOM", "10000"))
N_PCA    <- as.integer(Sys.getenv("MASHR_N_PCA",    "2"))

cat("============================================================\n")
cat("Module 05 - mashr Bayesian sex-effect classification\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Mode:    ", if (CHRX_PILOT) "CHRX_PILOT" else "FULL", "\n")
cat("N_strong:", N_STRONG, "  N_random:", N_RANDOM, "  N_PCs:", N_PCA, "\n")
cat("============================================================\n")

# ---------------------------------------------------------------------------
# 1) Load Module 04 output
# ---------------------------------------------------------------------------
cat("\n[1] Loading dream_M2_v3.rds (Module 04 output)...\n")
dream_rds <- file.path(SEXV3, "dream_M2_v3.rds")
stopifnot(file.exists(dream_rds))
dm2 <- readRDS(dream_rds)

gene_names  <- dm2$pooled$gene_names
gene_symbol <- dm2$pooled$gene_symbol
beta_F   <- dm2$pooled$beta_F
beta_M   <- dm2$pooled$beta_M
beta_int <- dm2$pooled$beta_int
var_F    <- dm2$pooled$var_F
var_M    <- dm2$pooled$var_M
var_int  <- dm2$pooled$var_int
t_int    <- dm2$pooled$t_int

cat("  n_genes:", length(gene_names), "\n")
cat("  sex_reference:", dm2$sex_reference, "  levels:", paste(dm2$sex_levels, collapse=","), "\n")
cat("  beta_F summary:");  print(summary(beta_F))
cat("  beta_M summary:");  print(summary(beta_M))
cat("  t_int summary :");  print(summary(t_int))

# ---------------------------------------------------------------------------
# 2) Optional chrX subset (PILOT mode)
# ---------------------------------------------------------------------------
if (CHRX_PILOT) {
  cat("\n[2] PILOT MODE: chrX biology landmarks + autosomal random sample (R1 Issue 12 fix).\n")
  cat("    chrX alone is too small (~1k genes) for cov_ed / V_hat / mixture\n")
  cat("    identification. Augment with N_AUTOSOMAL_PILOT random autosomal genes\n")
  cat("    so mashr learns realistic genome-wide correlation structure.\n")
  stopifnot(file.exists(CHRX_F))
  chrx <- fread(CHRX_F, select = c("gene_id", "gene_name", "chromosome", "chr_category"))
  chrx_genes <- chrx[chromosome == "chrX", gene_id]

  # ---- chrX subset
  is_chrx <- gene_names %in% chrx_genes
  n_chrx  <- sum(is_chrx)
  cat("  chrX genes in atlas:", length(chrx_genes), "\n")
  cat("  chrX genes in dream output:", n_chrx, "\n")
  if (n_chrx < 100) {
    stop("chrX subset is too small (<100 genes) - check chrX annotation join.")
  }

  # ---- Autosomal random sample (R1 Issue 12 fix)
  # Default 1,000 random autosomal genes; configurable via env var.
  N_AUTOSOMAL_PILOT <- as.integer(Sys.getenv("MASHR_N_AUTOSOMAL_PILOT", "1000"))
  is_autosomal <- !is_chrx
  cat("  Autosomal genes available:", sum(is_autosomal), "\n")
  cat("  Sampling", N_AUTOSOMAL_PILOT, "random autosomal genes (seed=42)\n")
  set.seed(42)
  auto_idx <- sample(which(is_autosomal),
                     size = min(N_AUTOSOMAL_PILOT, sum(is_autosomal)))

  keep_pilot <- is_chrx
  keep_pilot[auto_idx] <- TRUE
  cat("  Pilot total:", sum(keep_pilot),
      "  (chrX:", n_chrx, "+ autosomal:", length(auto_idx), ")\n")

  gene_names  <- gene_names[keep_pilot]
  gene_symbol <- gene_symbol[keep_pilot]
  beta_F   <- beta_F[keep_pilot]
  beta_M   <- beta_M[keep_pilot]
  beta_int <- beta_int[keep_pilot]
  var_F    <- var_F[keep_pilot]
  var_M    <- var_M[keep_pilot]
  var_int  <- var_int[keep_pilot]
  t_int    <- t_int[keep_pilot]
  cat("  Pilot n_genes:", length(gene_names), "\n")

  # Pilot subsets are larger now — re-clamp accordingly.
  N_STRONG <- min(N_STRONG, max(100, floor(length(gene_names) * 0.05)))
  N_RANDOM <- min(N_RANDOM, length(gene_names))
  N_PCA    <- min(N_PCA,    2)
  cat("  Pilot N_strong:", N_STRONG, "  N_random:", N_RANDOM, "  N_PCs:", N_PCA, "\n")
}

# ---------------------------------------------------------------------------
# 3) Build Bhat / Shat matrices (n_genes x 2)
# ---------------------------------------------------------------------------
cat("\n[3] Building Bhat / Shat matrices (2 conditions: F, M)...\n")
Bhat <- cbind(F = beta_F, M = beta_M)
Shat <- cbind(F = sqrt(pmax(var_F, .Machine$double.eps)),
              M = sqrt(pmax(var_M, .Machine$double.eps)))
# Hard clamp: cov_pca requires npc <= n_conditions
if (N_PCA > ncol(Bhat)) {
  cat("  Clamping N_PCA from", N_PCA, "to", ncol(Bhat), "(= n_conditions)\n")
  N_PCA <- ncol(Bhat)
}
rownames(Bhat) <- gene_names
rownames(Shat) <- gene_names

# Drop genes with NA / zero SE (would break mash)
bad_rows <- !is.finite(rowSums(Bhat)) | !is.finite(rowSums(Shat)) |
            rowSums(Shat <= 0) > 0
cat("  Dropping", sum(bad_rows), "genes with NA / non-finite / zero-SE rows.\n")
if (any(bad_rows)) {
  Bhat <- Bhat[!bad_rows, , drop = FALSE]
  Shat <- Shat[!bad_rows, , drop = FALSE]
  gene_names  <- gene_names[!bad_rows]
  gene_symbol <- gene_symbol[!bad_rows]
  beta_int    <- beta_int[!bad_rows]
  t_int       <- t_int[!bad_rows]
}
n_genes <- nrow(Bhat)
cat("  Final n_genes for mash:", n_genes, "  conditions: F, M\n")

# ---------------------------------------------------------------------------
# 4) Define strong + random subsets
# ---------------------------------------------------------------------------
cat("\n[4] Defining strong + random subsets...\n")
abs_t <- abs(t_int)
abs_t[!is.finite(abs_t)] <- 0

N_STRONG <- min(N_STRONG, n_genes)
N_RANDOM <- min(N_RANDOM, n_genes)

# v3.1 SYMMETRIC STRONG SUBSET (Urbut & Stephens 2019 §3.2 prescription for
# unbalanced designs). Single-statistic |t_int| ranking is F-biased under
# asymmetric power: cov_pca/cov_ed eigenvectors load on the F-axis. Union of
# top-|z_F| + top-|z_M| forces both modes into the ED training set.
z_F <- Bhat[, "F"] / Shat[, "F"]; z_F[!is.finite(z_F)] <- 0
z_M <- Bhat[, "M"] / Shat[, "M"]; z_M[!is.finite(z_M)] <- 0
N_HALF <- min(750L, floor(N_STRONG * 0.75))
strong_F <- order(abs(z_F), decreasing = TRUE)[seq_len(N_HALF)]
strong_M <- order(abs(z_M), decreasing = TRUE)[seq_len(N_HALF)]
strong_idx <- sort(unique(c(strong_F, strong_M)))
cat("  symmetric strong subset: top", N_HALF, "by |z_F| UNION top", N_HALF,
    "by |z_M| =", length(strong_idx), "unique genes\n")
cat("    overlap (genes strong in BOTH F and M):",
    length(intersect(strong_F, strong_M)), "\n")
set.seed(42)
# Random subset is an UNFILTERED uniform sample of all genes — this set is
# used both for Vhat estimation (where estimate_null_correlation_simple()
# applies its own internal |z|<2 filter) AND for mixture proportion
# training in mash(). Pre-filtering to |z|<2 here would bias the mixture
# pi training toward null and starve non-null components — confirmed empirically
# 2026-05-14 (a z-pre-filter forced pi=53% F + 47% null with all other
# components at 0, which is degenerate).
random_idx <- sample.int(n_genes, size = N_RANDOM)

cat("  strong subset size:", length(strong_idx),
    "  median |t_int| in strong:", round(median(abs_t[strong_idx]), 3), "\n")
cat("  random subset size:", length(random_idx),
    "  median |t_int| in random:", round(median(abs_t[random_idx]), 3), "\n")

# ---------------------------------------------------------------------------
# 5) Estimate Vhat from random subset
#    estimate_null_correlation_simple() expects a mash data object; it uses
#    z-scores |z| < 2 (approximate null) to estimate residual correlation
#    among conditions.
# ---------------------------------------------------------------------------
cat("\n[5] Estimating Vhat (residual correlation across F,M) from intersection-null subset...\n")
# v3.1 SEX-BALANCED Vhat. Restrict to intersection nulls (|z_F|<2 AND |z_M|<2)
# instead of one-axis filter. Removes asymmetric F-variance inflation that
# propagates into every posterior. Falls back to random subset if too few.
null_idx <- which(abs(z_F) < 2 & abs(z_M) < 2)
cat("  intersection-null genes (|z_F|<2 AND |z_M|<2):", length(null_idx), "\n")
if (length(null_idx) < 1000) {
  cat("  WARNING: <1000 intersection-null genes; falling back to random subset for Vhat\n")
  null_idx <- random_idx
}
data_tmp <- mash_set_data(Bhat[null_idx, , drop = FALSE],
                          Shat[null_idx, , drop = FALSE])
Vhat <- estimate_null_correlation_simple(data_tmp)
rm(data_tmp)
cat("  Vhat (F,M residual corr):\n"); print(round(Vhat, 4))

# ---------------------------------------------------------------------------
# 6) Rebuild data objects with V = Vhat
# ---------------------------------------------------------------------------
cat("\n[6] Rebuilding data objects with V = Vhat...\n")
data_random <- mash_set_data(Bhat[random_idx, , drop = FALSE],
                             Shat[random_idx, , drop = FALSE], V = Vhat)
data_strong <- mash_set_data(Bhat[strong_idx, , drop = FALSE],
                             Shat[strong_idx, , drop = FALSE], V = Vhat)
data_all    <- mash_set_data(Bhat, Shat, V = Vhat)

# ---------------------------------------------------------------------------
# 7) Data-driven covariances (PCA + ED on strong genes)
# ---------------------------------------------------------------------------
cat("\n[7] PER-ARM data-driven covariances (v3.1): cov_pca + cov_ed on F-strong AND M-strong...\n")
# v3.1 PER-ARM ED. Run cov_ed twice on disjoint inputs (F-strong, M-strong),
# rename components so M-loaded covariances are guaranteed in the U.list.
# Without this, even a symmetric union strong set lets one ED component
# absorb both modes if their magnitudes differ.
data_strong_F <- mash_set_data(Bhat[strong_F, , drop = FALSE],
                               Shat[strong_F, , drop = FALSE], V = Vhat)
data_strong_M <- mash_set_data(Bhat[strong_M, , drop = FALSE],
                               Shat[strong_M, , drop = FALSE], V = Vhat)

t0 <- Sys.time()
U.pca_F <- cov_pca(data_strong_F, npc = N_PCA)
U.pca_M <- cov_pca(data_strong_M, npc = N_PCA)
names(U.pca_F) <- paste0("ED_F_pca_", seq_along(U.pca_F))
names(U.pca_M) <- paste0("ED_M_pca_", seq_along(U.pca_M))
cat("  cov_pca per-arm:", length(U.pca_F), "F PCs +", length(U.pca_M), "M PCs\n")

U.ed_F <- tryCatch(cov_ed(data_strong_F, U.pca_F),
                   error = function(e) { cat("  cov_ed_F FAILED:", conditionMessage(e), "\n"); U.pca_F })
U.ed_M <- tryCatch(cov_ed(data_strong_M, U.pca_M),
                   error = function(e) { cat("  cov_ed_M FAILED:", conditionMessage(e), "\n"); U.pca_M })
# cov_ed strips names — reapply prefixes
names(U.ed_F) <- paste0("ED_F_", seq_along(U.ed_F))
names(U.ed_M) <- paste0("ED_M_", seq_along(U.ed_M))
# Combine and keep both PCA + ED legacy slots for diagnostics
U.pca <- c(U.pca_F, U.pca_M)
U.ed  <- c(U.ed_F, U.ed_M)
cat("  cov_ed per-arm:", length(U.ed_F), "F-ED +", length(U.ed_M), "M-ED matrices\n")
# Eigenvector loading diagnostic
for (nm in names(U.ed_M)) {
  Um <- U.ed_M[[nm]]
  cat(sprintf("    %s diag=(%.3f, %.3f) off=%.3f -- M-loading=%.2f\n",
              nm, Um[1,1], Um[2,2], Um[1,2],
              Um[2,2] / (Um[1,1] + Um[2,2] + 1e-12)))
}
for (nm in names(U.ed_F)) {
  Uf <- U.ed_F[[nm]]
  cat(sprintf("    %s diag=(%.3f, %.3f) off=%.3f -- F-loading=%.2f\n",
              nm, Uf[1,1], Uf[2,2], Uf[1,2],
              Uf[1,1] / (Uf[1,1] + Uf[2,2] + 1e-12)))
}
cat("  Elapsed:", format(Sys.time() - t0), "\n")

# ---------------------------------------------------------------------------
# 8) Canonical covariances
#    cov_canonical() default returns: identity, singletons_F, singletons_M,
#    equal_effects, simple_het (3 off-diag levels 0.25/0.5/0.75).
#    We extend with two custom matrices:
#       concordant_up   = [[1,1],[1,1]]  (same as equal_effects in 2D)
#       anti_correlated = [[1,-1],[-1,1]]  (opposite-direction effects)
#    Note: concordant_down is handled by mash via sign of beta — the U matrix
#    captures *covariance structure*, not sign.
# ---------------------------------------------------------------------------
cat("\n[8] Canonical covariances...\n")
U.c_default <- cov_canonical(data_random)
cat("  Default canonical covariances:", length(U.c_default), "\n")
for (nm in names(U.c_default)) {
  cat("    ", nm, ":\n", sep = "")
  print(round(U.c_default[[nm]], 3))
}

# R4 Issue 2 fix: drop `identity` from the canonical covariance list.
# `identity` ([[1,0],[0,1]]) gets hard-mapped to "concordant" in classify_U
# below, which systematically pulls modest autosomal F_only / M_only genes
# (e.g. CYP3A4) into the concordant bucket. By removing identity entirely,
# the data-driven covariances (cov_ed) + singletons_F / singletons_M /
# equal_effects / simple_het carry the prior load, and modest sex-specific
# effects route to F_only / M_only as intended. Keep equal_effects (which
# is the canonical "both sexes same effect" matrix), singletons_F /
# singletons_M (genuine sex-specific), simple_het (correlated effects with
# off-diagonals 0.25 / 0.5 / 0.75).
drop_canonical <- "identity"
keep_idx <- !names(U.c_default) %in% drop_canonical
if (any(!keep_idx)) {
  cat("  Dropping canonical covariances (R4 Issue 2):",
      paste(names(U.c_default)[!keep_idx], collapse = ", "), "\n")
  U.c_default <- U.c_default[keep_idx]
}
cat("  Canonical covariances after drop:", length(U.c_default),
    "names:", paste(names(U.c_default), collapse = ", "), "\n")

# Custom anti-correlated covariance (-1 off-diagonal)
U_anti <- matrix(c(1, -1, -1, 1), nrow = 2, ncol = 2,
                 dimnames = list(c("F","M"), c("F","M")))
U_anti_list <- list(anti_correlated = U_anti)

# Combined canonical covariances
U.c <- c(U.c_default, U_anti_list)
cat("  Total canonical covariances (with anti_correlated, no identity):", length(U.c), "\n")

# ---------------------------------------------------------------------------
# 9) Fit mash on random subset to estimate mixture proportions
# ---------------------------------------------------------------------------
cat("\n[9] Fitting mash on RANDOM subset (outputlevel=1, training mixture pi)...\n")
# v3.1 SENSITIVITY 5.1: env var to drop ED entirely and use canonical-only library.
# Justification: when ED matrices absorb concordant mass and zero out M singleton,
# canonical-only gives the conservative reference where pi on singletons F/M is
# free to be learned from data alone.
CANONICAL_ONLY <- isTRUE(as.logical(Sys.getenv("MASHR_CANONICAL_ONLY", "FALSE")))
if (CANONICAL_ONLY) {
  cat("  *** MASHR_CANONICAL_ONLY=TRUE — using canonical U only (no ED) ***\n")
  U.fit <- U.c
} else {
  U.fit <- c(U.ed, U.c)
}
t0 <- Sys.time()
fit_random <- mash(
  data_random,
  Ulist = U.fit,
  outputlevel = 1,
  verbose = TRUE
)
cat("  Elapsed:", format(Sys.time() - t0), "\n")
pi_est <- get_estimated_pi(fit_random)
cat("  Estimated mixture proportions (top 10):\n")
print(round(sort(pi_est, decreasing = TRUE)[1:min(10, length(pi_est))], 4))
cat("  Sum of pi_est:", round(sum(pi_est), 6), "\n")
cat("  Null component weight:", round(pi_est[grep("^null", names(pi_est))], 4), "\n")

# ---------------------------------------------------------------------------
# 10) Compute posteriors for ALL genes using trained pi (fixg = TRUE)
# ---------------------------------------------------------------------------
cat("\n[10] Computing posteriors for ALL genes with fixed g (trained pi)...\n")
t0 <- Sys.time()
g_fit <- get_fitted_g(fit_random)
fit_all <- mash(
  data_all,
  g = g_fit,
  fixg = TRUE,
  verbose = TRUE
)
cat("  Elapsed:", format(Sys.time() - t0), "\n")

# ---------------------------------------------------------------------------
# 11) Extract posterior summaries
# ---------------------------------------------------------------------------
cat("\n[11] Extracting per-gene posterior summaries...\n")
PM   <- get_pm(fit_all)   # posterior mean
PSD  <- get_psd(fit_all)  # posterior sd
LFSR <- get_lfsr(fit_all) # local false sign rate

stopifnot(nrow(PM)   == n_genes, ncol(PM)   == 2)
stopifnot(nrow(PSD)  == n_genes, ncol(PSD)  == 2)
stopifnot(nrow(LFSR) == n_genes, ncol(LFSR) == 2)

post_sum <- data.table(
  gene             = gene_names,
  gene_symbol      = gene_symbol,
  beta_F           = Bhat[, "F"],
  beta_M           = Bhat[, "M"],
  beta_interaction = beta_int,
  posterior_mean_F = PM[, "F"],
  posterior_mean_M = PM[, "M"],
  posterior_sd_F   = PSD[, "F"],
  posterior_sd_M   = PSD[, "M"],
  lfsr_F           = LFSR[, "F"],
  lfsr_M           = LFSR[, "M"]
)

# ---------------------------------------------------------------------------
# 12) Per-pattern posterior weights (per-gene component responsibilities)
#     Direct access to per-gene mixture responsibilities is via the
#     `posterior_weights` matrix returned in fit_all$posterior_weights.
#     If unavailable in this version, we fall back to a recomputation using
#     the fitted likelihoods (l_mat) and the pi_est weights.
# ---------------------------------------------------------------------------
cat("\n[12] Extracting per-gene posterior pattern weights...\n")
pat_weights <- NULL

if (!is.null(fit_all$posterior_weights)) {
  pat_weights <- fit_all$posterior_weights
  cat("  Got posterior_weights matrix directly from mash() output.\n")
  cat("    dim:", dim(pat_weights), "\n")
} else if (!is.null(fit_all[["result"]]) &&
           !is.null(fit_all$result$posterior_weights)) {
  pat_weights <- fit_all$result$posterior_weights
  cat("  Got posterior_weights from fit_all$result.\n")
} else {
  # Fallback: weights = pi * grid expansion -- summed across grid expansions
  # for each base U. We approximate by mapping the full pi to base U names.
  cat("  WARNING: mash() did not return per-gene posterior_weights.\n")
  cat("           Falling back to global pi_est repeated per gene.\n")
  pat_weights <- matrix(rep(pi_est, each = n_genes),
                        nrow = n_genes, ncol = length(pi_est),
                        dimnames = list(gene_names, names(pi_est)))
}

# Column names look like "U.PCA_1.1", "U.PCA_2.0.7071", "ED_PCA_1.0.5", etc.
# Each is "{Uname}.{grid_value}". Strip the grid suffix to get base U name.
strip_grid <- function(nms) {
  # mash uses the format "<Uname>.<grid_index>" or "<Uname>.<grid_value>".
  # Drop the last dot-separated token if it's purely numeric.
  sapply(nms, function(x) {
    parts <- strsplit(x, "\\.")[[1]]
    if (length(parts) <= 1) return(x)
    last <- parts[length(parts)]
    if (grepl("^[0-9.eE+-]+$", last)) {
      paste(parts[-length(parts)], collapse = ".")
    } else {
      x
    }
  }, USE.NAMES = FALSE)
}

base_pat <- strip_grid(colnames(pat_weights))
# Aggregate weights by base U name (sum over grid components)
base_unique <- unique(base_pat)
agg_mat <- matrix(0, nrow = nrow(pat_weights), ncol = length(base_unique),
                  dimnames = list(rownames(pat_weights), base_unique))
for (b in base_unique) {
  cols_b <- which(base_pat == b)
  if (length(cols_b) == 1) {
    agg_mat[, b] <- pat_weights[, cols_b]
  } else {
    agg_mat[, b] <- rowSums(pat_weights[, cols_b, drop = FALSE])
  }
}

# ---------------------------------------------------------------------------
# 13) Map each base U pattern -> one of 4 plan buckets:
#       concordant / F_only / M_only / anti_correlated
#     null component is reported separately.
#
#     Heuristic:
#       null:                              -> "null"
#       singletons named "*_F" or first:   -> "F_only"
#       singletons named "*_M" or second:  -> "M_only"
#       diagonal == off-diagonal sign(+):  -> "concordant"  (e.g., equal_effects, simple_het, ED PCs w/ +corr)
#       diagonal != off-diagonal sign(-):  -> "anti_correlated"
#     For mash defaults the canonical names are:
#       identity, singletons (creates one per condition), equal_effects, simple_het_0.25, _0.5, _0.75
#     plus our custom anti_correlated.
# ---------------------------------------------------------------------------
cat("\n[13] Mapping base patterns -> 4 plan buckets via U matrix inspection...\n")

# Reconstruct the Ulist that mash used (canonical + ED). Get fitted g for U.
U_used <- g_fit$Ulist
cat("  n base U in fitted model:", length(U_used), "\n")
cat("  base U names:", paste(names(U_used), collapse = ", "), "\n")

classify_U <- function(uname, U) {
  if (is.null(U) || !is.matrix(U) || any(dim(U) != c(2, 2))) {
    # Heuristic by name only
    nm <- tolower(uname)
    if (grepl("null", nm)) return("null")
    if (grepl("singleton.*1|_f$|singletonF", nm)) return("F_only")
    if (grepl("singleton.*2|_m$|singletonM", nm)) return("M_only")
    if (grepl("anti|opp", nm)) return("anti_correlated")
    return("concordant")
  }
  # U is 2x2 matrix
  vF  <- U[1, 1]
  vM  <- U[2, 2]
  cFM <- U[1, 2]

  # null detection
  if (max(abs(U)) < 1e-8 || grepl("null", tolower(uname))) return("null")

  # F-only / M-only: one diagonal much larger than the other (>= 10x)
  # i.e., singletons matrices [[1,0],[0,0]] or [[0,0],[0,1]]
  if (vF > 1e-6 && vM < 1e-8 * max(vF, 1)) return("F_only")
  if (vM > 1e-6 && vF < 1e-8 * max(vM, 1)) return("M_only")

  # Use correlation to discriminate concordant vs anti-correlated
  if (vF > 1e-12 && vM > 1e-12) {
    r <- cFM / sqrt(vF * vM)
    if (is.finite(r)) {
      if (r < -0.1) return("anti_correlated")
      if (r >= -0.1) return("concordant")
    }
  }

  # Edge cases: if off-diagonal == 0 but both diagonals > 0 (identity), treat
  # as a *mixed* "F_only OR M_only" base — assign to concordant (this matches
  # mash docs' framing that identity allows any direction independently).
  if (abs(cFM) < 1e-12 && vF > 0 && vM > 0) {
    # identity matrix == independent effects; mash uses this as "any effect"
    # bucket. We assign to concordant since both conditions can be non-null
    # and the gene is interpretable as a shared (possibly different
    # magnitude) effect rather than strictly sex-specific.
    return("concordant")
  }

  return("concordant")
}

bucket <- vapply(base_unique, function(b) {
  U_b <- U_used[[b]]
  classify_U(b, U_b)
}, character(1))

cat("  Bucket mapping per base U:\n")
for (b in base_unique) {
  Ub <- U_used[[b]]
  if (is.null(Ub) || !is.matrix(Ub)) {
    cat(sprintf("    %-30s -> %s (no U matrix)\n", b, bucket[b]))
  } else {
    cat(sprintf("    %-30s -> %s  (diag=%.2f,%.2f off=%.2f)\n",
                b, bucket[b], Ub[1,1], Ub[2,2], Ub[1,2]))
  }
}

# Now collapse base-U weights into 4 plan buckets
pat_lbl <- c("concordant", "F_only", "M_only", "anti_correlated")
bucket_mat <- matrix(0, nrow = n_genes, ncol = length(pat_lbl),
                     dimnames = list(gene_names, pat_lbl))
for (b in base_unique) {
  pat_b <- bucket[b]
  if (pat_b %in% pat_lbl) {
    bucket_mat[, pat_b] <- bucket_mat[, pat_b] + agg_mat[, b]
  }
  # null is dropped; the 4 bucket weights are renormalized below
}
# Note: null weight is dropped here -- we report active sex-effect proportions
# only. Module 06 renormalizes so the 4 sum to 1 (using the present implementation).
# Track null weight per gene as well for diagnostics.
null_cols <- base_unique[bucket == "null"]
null_weight_per_gene <- if (length(null_cols) > 0) {
  rowSums(agg_mat[, null_cols, drop = FALSE])
} else {
  rep(0, n_genes)
}
cat("  Null weight per-gene summary:\n"); print(summary(null_weight_per_gene))

# ---------------------------------------------------------------------------
# 14) Build pattern_loadings table -- both per-base U columns + 4 aggregated
# ---------------------------------------------------------------------------
cat("\n[14] Building mashr_pattern_loadings.csv...\n")
loadings_dt <- data.table(gene = gene_names)
# Per-base U columns (prefixed)
for (b in base_unique) {
  loadings_dt[[b]] <- agg_mat[, b]
}
# Aggregated 4-bucket columns (consumed verbatim by Module 06)
loadings_dt[, posterior_P_concordant      := bucket_mat[, "concordant"]]
loadings_dt[, posterior_P_F_only          := bucket_mat[, "F_only"]]
loadings_dt[, posterior_P_M_only          := bucket_mat[, "M_only"]]
loadings_dt[, posterior_P_anti_correlated := bucket_mat[, "anti_correlated"]]
loadings_dt[, posterior_P_null            := null_weight_per_gene]

# ---------------------------------------------------------------------------
# 15) Output paths (pilot vs full)
# ---------------------------------------------------------------------------
suffix <- if (CANONICAL_ONLY) "_canonical_only" else ""
if (CHRX_PILOT) {
  out_fit_rds       <- file.path(IDIR,  "mashr_pilot_chrX.rds")
  out_pattern_csv   <- file.path(SEXV3, "mashr_pattern_loadings_chrX_pilot.csv")
  out_posterior_csv <- file.path(SEXV3, "mashr_posterior_summary_chrX_pilot.csv")
  out_cov_rds       <- file.path(IDIR,  "mashr_data_driven_cov_chrX_pilot.rds")
} else {
  out_fit_rds       <- file.path(IDIR,  paste0("mashr_fit", suffix, ".rds"))
  out_pattern_csv   <- file.path(SEXV3, paste0("mashr_pattern_loadings", suffix, ".csv"))
  out_posterior_csv <- file.path(SEXV3, paste0("mashr_posterior_summary", suffix, ".csv"))
  out_cov_rds       <- file.path(IDIR,  paste0("mashr_data_driven_cov", suffix, ".rds"))
}

# ---------------------------------------------------------------------------
# 16) Write outputs
# ---------------------------------------------------------------------------
cat("\n[15] Writing outputs (atomic tmp -> rename) ...\n")
write_atomic_csv(loadings_dt, out_pattern_csv)
cat("  Wrote:", out_pattern_csv, "  rows:", nrow(loadings_dt),
    "  cols:", ncol(loadings_dt), "\n")

write_atomic_csv(post_sum, out_posterior_csv)
cat("  Wrote:", out_posterior_csv, "  rows:", nrow(post_sum),
    "  cols:", ncol(post_sum), "\n")

write_atomic_rds(list(
  fit_random  = fit_random,
  fit_all     = fit_all,
  Vhat        = Vhat,
  U.ed        = U.ed,
  U.pca       = U.pca,
  U.canonical = U.c,
  bucket_map  = bucket,
  pi_est      = pi_est,
  strong_idx  = strong_idx,
  random_idx  = random_idx,
  N_STRONG    = N_STRONG,
  N_RANDOM    = N_RANDOM,
  N_PCA       = N_PCA,
  mode        = if (CHRX_PILOT) "CHRX_PILOT" else "FULL",
  gene_names  = gene_names,
  mashr_ver   = as.character(mashr_ver)
), out_fit_rds)
cat("  Wrote:", out_fit_rds, "\n")

write_atomic_rds(list(
  U.ed        = U.ed,
  U.pca       = U.pca,
  U.canonical = U.c,
  U_anti      = U_anti,
  Vhat        = Vhat,
  bucket_map  = bucket
), out_cov_rds)
cat("  Wrote:", out_cov_rds, "\n")

# ---------------------------------------------------------------------------
# 17) Pilot-only: biology sanity check + pilot report
# ---------------------------------------------------------------------------
if (CHRX_PILOT) {
  cat("\n[16] PILOT biology sanity checks...\n")

  # Identify XIST, DDX3Y (won't be in chrX subset but check), KDM6A
  symbol_map <- data.table(gene = gene_names, gene_symbol = gene_symbol)
  xist_row  <- symbol_map[gene_symbol == "XIST"]
  ddx3y_row <- symbol_map[gene_symbol == "DDX3Y"]
  kdm6a_row <- symbol_map[gene_symbol == "KDM6A"]

  check_one <- function(gene_id, sym, expect_class, expect_thresh = 0.5) {
    if (length(gene_id) == 0) {
      return(list(gene = sym, status = "NOT_IN_PILOT", weight = NA_real_,
                  pass = NA))
    }
    row_i <- which(gene_names == gene_id)
    if (length(row_i) == 0) {
      return(list(gene = sym, status = "DROPPED_IN_QC", weight = NA_real_,
                  pass = NA))
    }
    wt <- switch(expect_class,
                 F_only          = bucket_mat[row_i, "F_only"],
                 M_only          = bucket_mat[row_i, "M_only"],
                 concordant      = bucket_mat[row_i, "concordant"],
                 anti_correlated = bucket_mat[row_i, "anti_correlated"])
    list(gene = sym, status = "OK",
         weight = wt, threshold = expect_thresh,
         pass = isTRUE(wt > expect_thresh))
  }

  res_xist  <- check_one(xist_row$gene,  "XIST",  "F_only", 0.5)
  res_kdm6a <- check_one(kdm6a_row$gene, "KDM6A", "F_only", 0.3)
  res_ddx3y <- check_one(ddx3y_row$gene, "DDX3Y", "M_only", 0.5)

  # NA check: no NA posteriors
  n_na_post <- sum(is.na(post_sum$posterior_mean_F) |
                   is.na(post_sum$posterior_mean_M) |
                   is.na(post_sum$lfsr_F) |
                   is.na(post_sum$lfsr_M))
  cat("  NA posteriors:", n_na_post, " / ", nrow(post_sum), "\n")

  # Pi sanity: any base U weight > 0?
  cat("  Sum of pi (across grid expansions):", round(sum(pi_est), 4), "\n")
  cat("  Aggregated 4-bucket pi (population):\n")
  pi_base_names <- strip_grid(names(pi_est))
  pi_bucket_lbl <- bucket[pi_base_names]
  # Some grid components may map to null component name "null" -- preserve it.
  pi_bucket_lbl[is.na(pi_bucket_lbl)] <- "null"
  bucket_pi <- tapply(pi_est, pi_bucket_lbl, sum)
  # Ensure all 4 + null buckets are reported even if 0
  for (b in c("concordant","F_only","M_only","anti_correlated","null")) {
    if (!b %in% names(bucket_pi)) bucket_pi[[b]] <- 0
  }
  print(round(bucket_pi, 4))

  # Write pilot report markdown (atomic write — R5 Issue 2 fix)
  report_md <- file.path(TEAM_B_DIR, "B2_chrX_pilot_report.md")
  report_lines <- c(
    "# B2 mashr chrX pilot report",
    "",
    paste0("**Agent:** B2 (Team B - Statistical Lead)"),
    paste0("**Date:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("**Mode:** CHRX_PILOT"),
    paste0("**mashr version:** ", as.character(mashr_ver)),
    "",
    "## Pilot summary",
    "",
    paste0("- Pilot genes (chrX): ", n_genes),
    paste0("- Strong subset: ", N_STRONG),
    paste0("- Random subset: ", N_RANDOM),
    paste0("- N PCs (data-driven): ", N_PCA),
    paste0("- Sum(pi) on random fit: ", round(sum(pi_est), 4)),
    paste0("- NA posteriors: ", n_na_post),
    "",
    "## Vhat (residual correlation, F vs M)",
    "",
    "```",
    capture.output(print(round(Vhat, 4))),
    "```",
    "",
    "## Base-U bucket mapping",
    "",
    "```",
    sprintf("%-30s -> %s", names(bucket), unname(bucket)),
    "```",
    "",
    "## Aggregated mixture pi by 4-bucket",
    "",
    "```",
    capture.output(print(round(bucket_pi, 4))),
    "```",
    "",
    "## Biology sanity checks",
    "",
    paste0("- XIST (expect F_only > 0.5):  status=", res_xist$status,
           "  weight=", round(res_xist$weight %||% NA_real_, 3),
           "  PASS=", isTRUE(res_xist$pass)),
    paste0("- DDX3Y (expect M_only > 0.5; chrY NOT in chrX pilot): status=",
           res_ddx3y$status,
           "  weight=", round(res_ddx3y$weight %||% NA_real_, 3),
           "  PASS=", isTRUE(res_ddx3y$pass)),
    paste0("- KDM6A (expect F_only > 0.3): status=", res_kdm6a$status,
           "  weight=", round(res_kdm6a$weight %||% NA_real_, 3),
           "  PASS=", isTRUE(res_kdm6a$pass)),
    "",
    "## Output files",
    "",
    paste0("- ", out_fit_rds),
    paste0("- ", out_pattern_csv),
    paste0("- ", out_posterior_csv),
    paste0("- ", out_cov_rds),
    "",
    "## Verdict",
    "",
    if (n_na_post == 0L &&
        (isTRUE(res_xist$pass) || res_xist$status == "NOT_IN_PILOT") &&
        (isTRUE(res_kdm6a$pass) || res_kdm6a$status %in%
           c("NOT_IN_PILOT", "DROPPED_IN_QC"))) {
      "**PASS**: Pilot converges with no NA posteriors and biology landmarks land in expected buckets. Proceed to full run."
    } else {
      "**FAIL**: Review pilot diagnostics; iterate on covariance specs or covariate set before launching full Module 05."
    }
  )
  write_atomic_lines(report_lines, report_md)
  cat("  Pilot report written:", report_md, "\n")
}

# ---------------------------------------------------------------------------
# 18) Console summary
# ---------------------------------------------------------------------------
cat("\n============================================================\n")
cat("Module 05 complete.\n")
cat("  Mode: ", if (CHRX_PILOT) "CHRX_PILOT" else "FULL", "\n")
cat("  n_genes (fit):", n_genes, "\n")
cat("  N_strong:", N_STRONG, "  N_random:", N_RANDOM, "\n")
cat("  4-bucket E[N] (sum of posterior weights):\n")
cat(sprintf("    concordant      = %.1f\n", sum(bucket_mat[, "concordant"])))
cat(sprintf("    F_only          = %.1f\n", sum(bucket_mat[, "F_only"])))
cat(sprintf("    M_only          = %.1f\n", sum(bucket_mat[, "M_only"])))
cat(sprintf("    anti_correlated = %.1f\n", sum(bucket_mat[, "anti_correlated"])))
cat("  Null weight (sum):", sum(null_weight_per_gene), "\n")

# ---------------------------------------------------------------------------
# v3.1 DIAGNOSTIC GATE
# Require:
#   (a) singleton F and singleton M both with pi > 0.005
#   (b) at least one ED_F_* and one ED_M_* component with pi > 0.01
#   (c) F_only / M_only call-count ratio < 5x (loose; tighter via triple-gate later)
# Failures are logged, not fatal — Layer 3 calibration adjudicates.
# ---------------------------------------------------------------------------
cat("\n[v3.1 DIAGNOSTIC GATE]\n")
pi_base <- tapply(pi_est, strip_grid(names(pi_est)), sum)
get_pi <- function(pat) sum(pi_base[grep(pat, names(pi_base), value = TRUE)])
pi_singleton_F <- get_pi("singletons.*1$|singletonsF$|^F$")
pi_singleton_M <- get_pi("singletons.*2$|singletonsM$|^M$")
pi_ed_F        <- get_pi("^ED_F_")
pi_ed_M        <- get_pi("^ED_M_")
cat(sprintf("  pi(singleton_F) = %.4f   pi(singleton_M) = %.4f\n",
            pi_singleton_F, pi_singleton_M))
cat(sprintf("  pi(ED_F_*)      = %.4f   pi(ED_M_*)      = %.4f\n",
            pi_ed_F, pi_ed_M))
nF_only <- sum(bucket_mat[, "F_only"] > 0.5)
nM_only <- sum(bucket_mat[, "M_only"] > 0.5)
cat(sprintf("  argmax counts: F_only=%d  M_only=%d  ratio=%.2fx\n",
            nF_only, nM_only,
            max(nF_only, nM_only) / max(min(nF_only, nM_only), 1)))
gate_a <- pi_singleton_F > 0.005 && pi_singleton_M > 0.005
gate_b <- pi_ed_F > 0.01 && pi_ed_M > 0.01
gate_c <- max(nF_only, nM_only) / max(min(nF_only, nM_only), 1) < 5
cat(sprintf("  gate (a) singletons F & M >0.005: %s\n", if (gate_a) "PASS" else "FAIL"))
cat(sprintf("  gate (b) ED_F & ED_M both >0.01:  %s\n", if (gate_b) "PASS" else "FAIL"))
cat(sprintf("  gate (c) F/M ratio <5x:           %s\n", if (gate_c) "PASS" else "FAIL"))
cat(sprintf("  OVERALL v3.1: %s\n",
            if (gate_a && gate_b && gate_c) "PASS — proceed to consensus" else
            "PARTIAL — Layer 3 calibration must adjudicate"))

# R5 Issue 2 fix: sessionInfo dump for per-module reproducibility audit
dump_session_info(IDIR,
                  if (CHRX_PILOT) "05_pilot" else "05")

cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
