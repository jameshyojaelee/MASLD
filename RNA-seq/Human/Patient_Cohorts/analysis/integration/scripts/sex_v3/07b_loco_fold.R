#!/usr/bin/env Rscript
# sex_v3/07b_loco_fold.R
# ---------------------------------------------------------------------------
# Module 07b -- ONE leave-one-cohort-out posterior-predictive fold.
# Controlled by env var LOCO_FOLD=1..5; the corresponding cohort is held
# out (alphabetic order: Bril, Chen, Govaere, Hoang, Suppli), dream M2 is
# fit on the remaining 4, mash with trained g predicts (beta_F, beta_M) on
# all genes, and the fold predictions vs the full-data posterior are
# stashed to `intermediates/loco_fold_${LOCO_FOLD}_${COHORT}.rds`. Module
# 07c aggregates all 5 folds + Test 1 into the final figure.
#
# Reads:
#   intermediates/sex_v3_input.rds       (Module 01)
#   intermediates/sva_factors.rds        (Module 02)
#   intermediates/age_mi.rds             (Module 03)
#   dream_M2_v3.rds                      (Module 04 / 04b)
#   intermediates/mashr_fit.rds          (Module 05)
#   mashr_posterior_summary.csv          (Module 05)
#
# Writes: intermediates/loco_fold_${LOCO_FOLD}_${COHORT}.rds
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
})

# lme4 namespace injection BEFORE variancePartition load
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({
      unlockBinding(fn, ns_lme4)
      assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
      lockBinding(fn, ns_lme4)
    }, silent = TRUE)
  }
}

suppressPackageStartupMessages({
  library(variancePartition)
  library(BiocParallel)
  library(mashr)
})

# ---------------------------------------------------------------------------
# Resolve LOCO fold + cohort
# ---------------------------------------------------------------------------
loco_fold <- as.integer(Sys.getenv("LOCO_FOLD", "0"))
if (is.na(loco_fold) || loco_fold < 1L || loco_fold > 5L) {
  stop("LOCO_FOLD env var must be an integer in [1,5]; got '",
       Sys.getenv("LOCO_FOLD"), "'")
}

# Deterministic cohort list -- alphabetic order over factor levels.
# 5-cohort canonical mega: GSE126848(Suppli) / GSE130970(Hoang) /
# GSE135251(Govaere) / GSE162694(Bril) / GSE213621(Chen).
COHORTS_SORTED <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
target_cohort <- COHORTS_SORTED[loco_fold]
cat("LOCO_FOLD:", loco_fold, "  target cohort (target HOLDOUT):", target_cohort, "\n")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")
SEXV3  <- file.path(RDIR, "sex_v3")
IDIR   <- file.path(SEXV3, "intermediates")

# Contrast routing (sex_v3_utils.R::contrast_paths) — overrides SEXV3 / IDIR
# when CONTRAST_NAME != "disease_vs_ctrl"; default preserves legacy layout.
if (!exists("contrast_paths")) source(file.path(
  Sys.getenv("MASLD_PROJECT_ROOT",
             "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3/sex_v3_utils.R"))
.cpaths <- contrast_paths()
SEXV3 <- .cpaths$sexv3
IDIR  <- .cpaths$idir
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

cat("============================================================\n")
cat("Module 07b -- LOCO fold ", loco_fold, " (", target_cohort, ")\n", sep = "")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")

# ---------------------------------------------------------------------------
# Parallel
# ---------------------------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(workers = ncpus, RNGseed = 42L) else SerialParam()

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
cat("\n[1] Loading inputs...\n")
inp     <- readRDS(file.path(IDIR, "sex_v3_input.rds"))
sva_o   <- readRDS(file.path(IDIR, "sva_factors.rds"))
age_mi  <- readRDS(file.path(IDIR, "age_mi.rds"))
mashr_o <- readRDS(file.path(IDIR, "mashr_fit.rds"))

dge   <- inp$dge
info0 <- inp$meta
sv_mat <- sva_o$sv
n_sv   <- ncol(sv_mat)

sv_df <- as.data.frame(sv_mat)
colnames(sv_df) <- paste0("SV", seq_len(n_sv))
sv_terms <- paste(colnames(sv_df), collapse = " + ")
form_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed",
  " + ", sv_terms,
  " + (1 + group_binary | dataset)"
)
form_full <- as.formula(form_str)

info_template <- info0
info_template$age_imputed <- age_mi$age_list[[1]]
info_template <- cbind(info_template, sv_df)
rownames(info_template) <- rownames(info0)

# Sanity-check cohort labels
all_cohorts <- levels(info_template$dataset)
cat("  All cohorts in info$dataset:", paste(all_cohorts, collapse = ", "), "\n")
if (!target_cohort %in% all_cohorts) {
  stop("Target cohort '", target_cohort,
       "' is not in info_template$dataset levels: ",
       paste(all_cohorts, collapse = ", "))
}

# ---------------------------------------------------------------------------
# Helpers (identical to legacy 07_calibration.R)
# ---------------------------------------------------------------------------
fit_dream_panel <- function(dge_panel, info_k, label = "fit") {
  cat("    [", label, "] voomWithDreamWeights...\n", sep = "")
  t0 <- Sys.time()
  v <- suppressWarnings(
    voomWithDreamWeights(dge_panel, form_full, info_k, BPPARAM = param,
                         useWeights = TRUE)
  )
  cat("      elapsed:", format(Sys.time() - t0), "\n")

  cat("    [", label, "] dream (random-slope)...\n", sep = "")
  t0 <- Sys.time()
  fit <- suppressWarnings(
    dream(v, form_full, info_k, BPPARAM = param, useWeights = TRUE)
  )
  cat("      elapsed:", format(Sys.time() - t0), "\n")

  all_coefs <- colnames(fit$coefficients)
  group_coef <- grep("^group_binary", all_coefs, value = TRUE)
  group_coef <- setdiff(group_coef, grep(":", group_coef, value = TRUE))
  int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary",
                   all_coefs, value = TRUE)[1]
  stopifnot(length(group_coef) == 1, length(int_coef) >= 1)

  tt_F   <- topTable(fit, coef = group_coef, number = Inf, sort.by = "none")
  tt_int <- topTable(fit, coef = int_coef,   number = Inf, sort.by = "none")
  stopifnot(identical(rownames(tt_F), rownames(tt_int)))

  beta_F  <- tt_F$logFC
  beta_int <- tt_int$logFC
  beta_M  <- beta_F + beta_int
  # R1-Issue-1/2 fix (re-applied per meta-review B0-4): raw Wald SE (not
  # moderated logFC/t which is sigma-shrunk) + per-gene Cov(β_F, β_int) from
  # fit$cov.coefficients.list. Omitting 2*Cov biases Var(β_M).
  raw_se_mat <- fit$stdev.unscaled * fit$sigma
  se_F   <- raw_se_mat[, group_coef]
  se_int <- raw_se_mat[, int_coef]
  if (!"cov.coefficients.list" %in% names(fit) ||
      length(fit$cov.coefficients.list) != length(beta_F)) {
    stop("fit$cov.coefficients.list missing or wrong length: ",
         length(fit$cov.coefficients.list), " vs ", length(beta_F))
  }
  ccl <- fit$cov.coefficients.list
  cov_F_int <- vapply(seq_along(ccl), function(g) {
    M <- ccl[[g]]
    if (group_coef %in% rownames(M) && int_coef %in% colnames(M))
      as.numeric(M[group_coef, int_coef])
    else NA_real_
  }, numeric(1))
  if (any(is.na(cov_F_int))) stop(sum(is.na(cov_F_int)),
                                  " genes have missing Cov(β_F, β_int)")
  se_M   <- sqrt(pmax(se_F^2 + se_int^2 + 2 * cov_F_int, .Machine$double.eps))

  list(
    gene = rownames(tt_F),
    beta_F = beta_F, beta_M = beta_M,
    se_F = se_F, se_M = se_M,
    beta_int = beta_int, se_int = se_int,
    cov_F_int = cov_F_int,
    t_int = beta_int / pmax(se_int, .Machine$double.eps)
  )
}

g_trained <- get_fitted_g(mashr_o$fit_random)
Vhat_trained <- mashr_o$Vhat

run_mash_with_trained_g <- function(beta_F, beta_M, se_F, se_M, gene_names) {
  Bhat <- cbind(F = beta_F, M = beta_M)
  Shat <- cbind(F = sqrt(pmax(se_F^2, .Machine$double.eps)),
                M = sqrt(pmax(se_M^2, .Machine$double.eps)))
  rownames(Bhat) <- gene_names
  rownames(Shat) <- gene_names

  bad <- !is.finite(rowSums(Bhat)) | !is.finite(rowSums(Shat)) |
         rowSums(Shat <= 0) > 0
  if (any(bad)) {
    Bhat <- Bhat[!bad, , drop = FALSE]
    Shat <- Shat[!bad, , drop = FALSE]
    gene_names <- gene_names[!bad]
  }
  if (nrow(Bhat) < 5) {
    stop("Too few genes (", nrow(Bhat), ") after QC drop for mash to run.")
  }
  d <- mash_set_data(Bhat, Shat, V = Vhat_trained)
  mash(d, g = g_trained, fixg = TRUE, verbose = FALSE)
}

# ---------------------------------------------------------------------------
# Hold out target cohort, fit dream + mash, compare to full-data posterior
# ---------------------------------------------------------------------------
post_sum_full <- fread(file.path(SEXV3, "mashr_posterior_summary.csv"))
setkey(post_sum_full, gene)

cat("\n========== HELD OUT:", target_cohort, "==========\n")
t0 <- Sys.time()
keep_samp <- info_template$dataset != target_cohort
n_hold <- sum(!keep_samp); n_keep <- sum(keep_samp)
cat("  samples held out:", n_hold, "  samples retained:", n_keep, "\n")

dge_k <- dge[, keep_samp]
info_k <- info_template[keep_samp, , drop = FALSE]
info_k$dataset <- droplevels(info_k$dataset)

fit_k <- fit_dream_panel(dge_k, info_k,
                         label = paste0("LOCO_", target_cohort))
mash_k <- run_mash_with_trained_g(fit_k$beta_F, fit_k$beta_M,
                                  fit_k$se_F,   fit_k$se_M,
                                  fit_k$gene)

PM_k <- get_pm(mash_k)
gene_k <- rownames(PM_k)
m_idx <- match(gene_k, post_sum_full$gene)
ok <- !is.na(m_idx)

loco_dt <- data.table(
  cohort_held_out      = target_cohort,
  gene                 = gene_k[ok],
  predicted_beta_F     = PM_k[ok, "F"],
  predicted_beta_M     = PM_k[ok, "M"],
  actual_beta_F        = post_sum_full$posterior_mean_F[m_idx[ok]],
  actual_beta_M        = post_sum_full$posterior_mean_M[m_idx[ok]]
)

cor_F <- suppressWarnings(cor(loco_dt$predicted_beta_F, loco_dt$actual_beta_F,
                              use = "complete.obs"))
cor_M <- suppressWarnings(cor(loco_dt$predicted_beta_M, loco_dt$actual_beta_M,
                              use = "complete.obs"))
r2_F <- if (is.finite(cor_F)) cor_F^2 else NA_real_
r2_M <- if (is.finite(cor_M)) cor_M^2 else NA_real_
r2_avg <- mean(c(r2_F, r2_M), na.rm = TRUE)
loco_dt[, predictive_r2_F := r2_F]
loco_dt[, predictive_r2_M := r2_M]
loco_dt[, predictive_r2   := r2_avg]
cat("  predictive r^2: F=", round(r2_F, 4), "  M=", round(r2_M, 4),
    "  avg=", round(r2_avg, 4), "\n")
cat("  LOCO elapsed for", target_cohort, ":", format(Sys.time() - t0), "\n")

# Persist per-fold output
out_rds <- file.path(IDIR, sprintf("loco_fold_%d_%s.rds",
                                   loco_fold, target_cohort))
write_atomic_rds(list(
  loco_fold      = loco_fold,
  cohort_held    = target_cohort,
  loco_dt        = loco_dt,
  predictive_r2_F = r2_F,
  predictive_r2_M = r2_M,
  predictive_r2   = r2_avg,
  n_samples_held = n_hold,
  n_samples_keep = n_keep,
  n_genes        = nrow(loco_dt)
), out_rds)
cat("Saved (atomic):", out_rds, "\n")

# R5 Issue 2 fix: sessionInfo dump (one per fold, namespaced)
dump_session_info(IDIR, sprintf("07b_fold%d", loco_fold))

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
