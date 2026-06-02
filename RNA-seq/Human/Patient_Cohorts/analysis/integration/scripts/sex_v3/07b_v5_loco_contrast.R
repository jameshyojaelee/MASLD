#!/usr/bin/env Rscript
# sex_v3/07b_v5_loco_contrast.R
# ---------------------------------------------------------------------------
# Contrast-aware port of `archive/pre_v6_2026-05-17/07b_v5_loco.R`.
# Performs ONE leave-one-cohort-out posterior-predictive fold for
# whichever contrast `CONTRAST_NAME` selects (default disease_vs_ctrl).
#
# Differences from the archived script:
#   - sources sex_v3_utils.R for contrast_paths() routing
#   - discovers cohort list dynamically from the contrast-specific input
#     (5 cohorts for disease_vs_ctrl, 7 for mash_vs_masl, etc.)
#   - writes to {SEXV3}/intermediates/loco_v5/fold_{LOCO_FOLD}.rds where
#     SEXV3 is contrast-aware
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR)
})
ns_lme4 <- asNamespace("lme4")
for (fn in c("findbars", "nobars", "subbars", "rebuildFormula")) {
  if (exists(fn, envir = ns_lme4)) {
    try({ unlockBinding(fn, ns_lme4)
          assign(fn, get(fn, asNamespace("reformulas")), envir = ns_lme4)
          lockBinding(fn, ns_lme4) }, silent = TRUE)
  }
}
suppressPackageStartupMessages({
  library(variancePartition); library(BiocParallel); library(limma)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

CPATHS <- contrast_paths()
SEXV3  <- CPATHS$sexv3
IDIR   <- CPATHS$idir
OUTDIR <- file.path(IDIR, "loco_v5")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Discover cohorts from the contrast-specific input (alphabetic order)
inp     <- readRDS(file.path(IDIR, "sex_v3_input.rds"))
sva_o   <- readRDS(file.path(IDIR, "sva_factors.rds"))
age_mi  <- readRDS(file.path(IDIR, "age_mi.rds"))

dge   <- inp$dge
info0 <- inp$meta
COHORTS_SORTED <- sort(levels(droplevels(as.factor(info0$dataset))))
n_folds <- length(COHORTS_SORTED)

loco_fold <- as.integer(Sys.getenv("LOCO_FOLD", "0"))
if (is.na(loco_fold) || loco_fold < 1L || loco_fold > n_folds) {
  stop(sprintf("LOCO_FOLD must be integer in [1,%d]; got '%s'",
               n_folds, Sys.getenv("LOCO_FOLD")))
}
target_cohort <- COHORTS_SORTED[loco_fold]
out_rds <- file.path(OUTDIR, sprintf("fold_%d.rds", loco_fold))

cat("============================================================\n")
cat("Module 07b v5 (contrast) -- fold ", loco_fold, " of ", n_folds,
    " (HOLDOUT=", target_cohort, ")\n", sep = "")
cat("Contrast: ", CPATHS$contrast, "\n", sep = "")
cat("SEXV3:    ", SEXV3, "\n", sep = "")
cat("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(workers = ncpus, RNGseed = 42L) else SerialParam()

sv_mat <- sva_o$sv
n_sv   <- ncol(sv_mat)
sv_df  <- as.data.frame(sv_mat)
colnames(sv_df) <- paste0("SV", seq_len(n_sv))
sv_terms <- paste(colnames(sv_df), collapse = " + ")

info_template <- info0
info_template$age_imputed <- age_mi$age_list[[1]]
info_template <- cbind(info_template, sv_df)
rownames(info_template) <- rownames(info0)

form_train_slope_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 + group_binary | dataset)"
)
form_train_intercept_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 | dataset)"
)
form_holdout_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms
)

all_cohorts <- levels(info_template$dataset)
cat("  cohorts:", paste(all_cohorts, collapse = ", "), "\n")
if (!target_cohort %in% all_cohorts) {
  stop("Target cohort '", target_cohort, "' not in dataset levels")
}

extract_int <- function(fit) {
  all_coefs <- colnames(fit$coefficients)
  int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary",
                   all_coefs, value = TRUE)[1]
  if (is.na(int_coef)) stop("Interaction coef not found in: ",
                            paste(all_coefs, collapse = ", "))
  tt_int <- topTable(fit, coef = int_coef, number = Inf, sort.by = "none")
  list(gene = rownames(tt_int),
       beta_int = tt_int$logFC,
       se_int   = abs(tt_int$logFC / tt_int$t),
       t_int    = tt_int$t,
       p_int    = tt_int$P.Value)
}

cat("\n[2] HELD OUT:", target_cohort, "\n")
keep_samp <- info_template$dataset != target_cohort
n_hold <- sum(!keep_samp); n_keep <- sum(keep_samp)
cat("  samples held out:", n_hold, "  samples retained:", n_keep, "\n")

dge_train <- dge[, keep_samp]
info_train <- info_template[keep_samp, , drop = FALSE]
info_train$dataset <- droplevels(info_train$dataset)
dge_train <- calcNormFactors(dge_train)

rank_def_fallback <- FALSE
dream_method <- "slope"
fit_train <- NULL

cat("  [training] voom + dream (random slope)...\n")
t0 <- Sys.time()
fit_train <- tryCatch({
  v <- suppressWarnings(voomWithDreamWeights(dge_train, as.formula(form_train_slope_str),
                                             info_train, BPPARAM = param,
                                             useWeights = TRUE))
  suppressWarnings(dream(v, as.formula(form_train_slope_str), info_train,
                         BPPARAM = param, useWeights = TRUE))
}, error = function(e) {
  cat("    slope dream FAILED:", conditionMessage(e), "\n")
  NULL
})
if (is.null(fit_train)) {
  rank_def_fallback <- TRUE
  dream_method <- "intercept"
  cat("  [training] dream intercept fallback...\n")
  fit_train <- tryCatch({
    v <- suppressWarnings(voomWithDreamWeights(dge_train,
                                               as.formula(form_train_intercept_str),
                                               info_train, BPPARAM = param,
                                               useWeights = TRUE))
    suppressWarnings(dream(v, as.formula(form_train_intercept_str), info_train,
                           BPPARAM = param, useWeights = TRUE))
  }, error = function(e) {
    cat("    intercept dream FAILED:", conditionMessage(e), "\n")
    NULL
  })
}
if (is.null(fit_train)) {
  dream_method <- "limma_fixed"
  cat("  [training] limma fixed-effects fallback...\n")
  form_lm_fb <- paste0(
    "~ group_binary * inferred_sex",
    " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
    " + age_imputed + ", sv_terms, " + dataset")
  v <- voom(dge_train, model.matrix(as.formula(form_lm_fb), info_train))
  fit_train <- eBayes(lmFit(v, model.matrix(as.formula(form_lm_fb), info_train)))
}
cat("  training fit done in", format(Sys.time() - t0), "  method:", dream_method, "\n")

train_extr <- extract_int(fit_train)
cat("  training b_int extracted for", length(train_extr$gene), "genes\n")

cat("\n[3] Fitting held-out cohort (", target_cohort, ") with limma...\n", sep = "")
dge_hold <- dge[, !keep_samp]
info_hold <- info_template[!keep_samp, , drop = FALSE]
info_hold$dataset <- droplevels(info_hold$dataset)
dge_hold <- calcNormFactors(dge_hold)

sv_cols <- colnames(sv_df)
sv_keep <- sv_cols[vapply(sv_cols, function(s)
                          var(info_hold[[s]], na.rm = TRUE) > 1e-10,
                          logical(1))]
if (length(sv_keep) < length(sv_cols)) {
  cat("  Dropping", length(sv_cols) - length(sv_keep),
      "constant SVs in held-out cohort\n")
}
sv_terms_hold <- if (length(sv_keep) > 0)
                   paste("+", paste(sv_keep, collapse = " + ")) else ""
form_hold_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed ", sv_terms_hold
)
form_hold_min_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes"
)
fit_hold <- NULL
holdout_method <- "limma_full"
fit_hold <- tryCatch({
  mm <- model.matrix(as.formula(form_hold_str), info_hold)
  v <- voom(dge_hold, mm)
  eBayes(lmFit(v, mm))
}, error = function(e) {
  cat("    limma FULL FAILED:", conditionMessage(e), "\n")
  NULL
})
if (is.null(fit_hold)) {
  holdout_method <- "limma_minimal"
  cat("  [holdout] minimal-form fallback...\n")
  mm <- model.matrix(as.formula(form_hold_min_str), info_hold)
  v <- voom(dge_hold, mm)
  fit_hold <- eBayes(lmFit(v, mm))
}
hold_extr <- extract_int(fit_hold)
cat("  holdout b_int extracted for", length(hold_extr$gene),
    "genes  method:", holdout_method, "\n")

cat("\n[4] Computing sign concordance per gene...\n")
common <- intersect(train_extr$gene, hold_extr$gene)
cat("  common genes:", length(common), "\n")
mt <- match(common, train_extr$gene); mh <- match(common, hold_extr$gene)

loco_dt <- data.table(
  gene             = common,
  beta_int_pred    = train_extr$beta_int[mt],
  se_int_pred      = train_extr$se_int[mt],
  p_int_pred       = train_extr$p_int[mt],
  beta_int_obs     = hold_extr$beta_int[mh],
  se_int_obs       = hold_extr$se_int[mh],
  p_int_obs        = hold_extr$p_int[mh]
)
loco_dt[, sign_match := sign(beta_int_pred) == sign(beta_int_obs) &
                       beta_int_pred != 0 & beta_int_obs != 0]
overall_sign_conc <- mean(loco_dt$sign_match, na.rm = TRUE)
cat("  overall sign concordance:", round(overall_sign_conc, 4), "\n")

top_idx <- order(loco_dt$p_int_pred)[seq_len(min(500, nrow(loco_dt)))]
top_conc <- mean(loco_dt$sign_match[top_idx], na.rm = TRUE)
cat("  top-500 by p_int_pred sign concordance:", round(top_conc, 4), "\n")

write_atomic_rds(list(
  fold                   = loco_fold,
  cohort_held            = target_cohort,
  dream_method           = dream_method,
  holdout_method         = holdout_method,
  rank_deficient_fallback = rank_def_fallback,
  n_held                 = n_hold,
  n_keep                 = n_keep,
  n_genes                = nrow(loco_dt),
  loco_dt                = loco_dt,
  overall_sign_conc      = overall_sign_conc,
  top500_sign_conc       = top_conc
), out_rds)
cat("Saved (atomic):", out_rds, "\n")

dump_session_info(IDIR, sprintf("07b_v5_fold%d_%s", loco_fold, CPATHS$contrast))

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
