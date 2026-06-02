#!/usr/bin/env Rscript
# sex_v3/13_bootstrap_dream.R
# ---------------------------------------------------------------------------
# Pillar 4 — Dream-M2 bootstrap (replaces v5 limma+lmFit bootstrap).
# Per rep: cohort × sex × group_binary stratified equal-N subsample (without
# replacement, Meinshausen-Bühlmann), fresh voom + dream M2 with the F2
# random-slope formula from P1. Per-gene classification via v5 decision tree.
#
# Forked from 04g_bootstrap_interaction.R; replaces limma::voom + lmFit + ashr
# with voomWithDreamWeights + dream + ashr per arm.
#
# Output: intermediates/boot_v6/rep_{REP}.csv with
#   gene, beta_F, beta_M, beta_int, p_int, q_int, lfsr_F, lfsr_M, class_v6_rep
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(reformulas); library(lme4); library(data.table); library(edgeR); library(ashr)
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
  library(variancePartition); library(BiocParallel)
})

REP <- as.integer(Sys.getenv("REP",
                  Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
stopifnot(is.integer(REP), REP >= 1)
set.seed(20260515L + REP * 7L)
cat("=== sex_v6 P4 bootstrap dream rep", REP, "===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SEXV3 <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                   "results/integration/sex_v3")
IDIR  <- file.path(SEXV3, "intermediates")

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
BOOT_DIR <- file.path(IDIR, "boot_v6")
dir.create(BOOT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PATH <- sprintf("%s/rep_%03d.csv", BOOT_DIR, REP)

emit_na_and_exit <- function(reason) {
  cat("FAILURE reason:", reason, "— writing NA placeholder and exiting 0\n")
  na_df <- data.table(gene = NA_character_, beta_F = NA_real_, beta_M = NA_real_,
                      beta_int = NA_real_, p_int = NA_real_, q_int = NA_real_,
                      lfsr_F = NA_real_, lfsr_M = NA_real_,
                      class_v6_rep = NA_character_, failure_reason = reason)
  fwrite(na_df, OUT_PATH); quit(save = "no", status = 0)
}

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
param <- if (ncpus > 1) MulticoreParam(workers = ncpus, RNGseed = 42L + REP) else SerialParam()
cat("CPU cores:", ncpus, "\n")

ok <- tryCatch({
  inp    <- readRDS(file.path(IDIR, "sex_v3_input.rds"))
  sva    <- readRDS(file.path(IDIR, "sva_factors.rds"))
  age_mi <- readRDS(file.path(IDIR, "age_mi.rds"))
  TRUE
}, error = function(e) { cat("Input load failed:", conditionMessage(e), "\n"); FALSE })
if (!ok) emit_na_and_exit("input_load_failed")

info <- as.data.table(inp$meta, keep.rownames = "sample_id")
sv_mat <- sva$sv
sv_df  <- as.data.frame(sv_mat[match(info$sample_id, rownames(sv_mat)), ,
                               drop = FALSE])
colnames(sv_df) <- paste0("SV", seq_len(ncol(sv_df)))
info[, age_imputed := age_mi$age_list[[1]][match(sample_id, age_mi$sample_ids)]]
info <- cbind(info, sv_df)
cat("SV count:", ncol(sv_df), "\n")

# Stratified subsample
keep_ids <- character()
for (ds in levels(info$dataset)) {
  for (gr in c("Control", "Disease")) {
    ids_M <- info[dataset == ds & inferred_sex == "M" & group_binary == gr, sample_id]
    ids_F <- info[dataset == ds & inferred_sex == "F" & group_binary == gr, sample_id]
    n_keep <- min(length(ids_M), length(ids_F))
    if (n_keep == 0) next
    if (length(ids_M) > n_keep) ids_M <- sample(ids_M, n_keep)
    if (length(ids_F) > n_keep) ids_F <- sample(ids_F, n_keep)
    keep_ids <- c(keep_ids, ids_M, ids_F)
  }
}
info_bs <- info[sample_id %in% keep_ids]
info_bs <- info_bs[match(intersect(colnames(inp$dge), info_bs$sample_id),
                         sample_id)]
dge_bs <- inp$dge[, info_bs$sample_id]
dge_bs <- calcNormFactors(dge_bs)
df_bs <- as.data.frame(info_bs)
rownames(df_bs) <- info_bs$sample_id
stopifnot(identical(rownames(df_bs), colnames(dge_bs)))

cat("Subsample sizes:\n"); print(table(info_bs$inferred_sex, info_bs$group_binary))
cat("Per-cohort sizes:\n"); print(table(info_bs$dataset, info_bs$inferred_sex,
                                          info_bs$group_binary))

sv_cols <- grep("^SV", colnames(info_bs), value = TRUE)
sv_terms <- paste(sv_cols, collapse = " + ")
form_F2_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 + group_binary + inferred_sex + group_binary:inferred_sex || dataset)"
)
form_F5_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 + group_binary | dataset)"
)
form_F2 <- as.formula(form_F2_str)
form_F5 <- as.formula(form_F5_str)

cat("[1] voomWithDreamWeights (F2)...\n")
t0 <- Sys.time()
v <- tryCatch(
  suppressWarnings(voomWithDreamWeights(dge_bs, form_F2, df_bs,
                                        BPPARAM = param, useWeights = TRUE)),
  error = function(e) { cat("  voom F2 failed:", conditionMessage(e), "\n"); NULL })
cat("  voom elapsed:", format(Sys.time() - t0), "\n")

variant_used <- "F2"
fit_bs <- NULL
if (!is.null(v)) {
  cat("[2] dream (F2)...\n")
  t0 <- Sys.time()
  fit_bs <- tryCatch(
    suppressWarnings(dream(v, form_F2, df_bs,
                           BPPARAM = param, useWeights = TRUE)),
    error = function(e) { cat("  dream F2 failed:", conditionMessage(e), "\n"); NULL })
  cat("  dream elapsed:", format(Sys.time() - t0), "\n")
}
if (is.null(fit_bs)) {
  cat("[2b] Fallback to F5...\n")
  variant_used <- "F5"
  v <- tryCatch(suppressWarnings(voomWithDreamWeights(dge_bs, form_F5, df_bs,
                                                       BPPARAM = param, useWeights = TRUE)),
                error = function(e) NULL)
  if (is.null(v)) emit_na_and_exit("voom_F5_failed")
  fit_bs <- tryCatch(suppressWarnings(dream(v, form_F5, df_bs,
                                             BPPARAM = param, useWeights = TRUE)),
                     error = function(e) NULL)
}
if (is.null(fit_bs)) emit_na_and_exit("dream_failed")
cat("Variant used:", variant_used, "\n")

# Extract coefficients
all_coefs <- colnames(fit_bs$coefficients)
group_coef <- setdiff(grep("^group_binary", all_coefs, value = TRUE),
                     grep(":", all_coefs, value = TRUE))[1]
int_coef <- grep(":", grep("^group_binary", all_coefs, value = TRUE), value = TRUE)[1]
if (is.na(group_coef) || is.na(int_coef)) emit_na_and_exit("coefs_not_found")
cat("group_coef:", group_coef, "  int_coef:", int_coef, "\n")

beta_F   <- fit_bs$coefficients[, group_coef]
beta_int <- fit_bs$coefficients[, int_coef]
beta_M   <- beta_F + beta_int
raw_se   <- fit_bs$sigma * fit_bs$stdev.unscaled
se_F     <- raw_se[, group_coef]
se_int   <- raw_se[, int_coef]
ccl <- fit_bs$cov.coefficients.list
cov_F_int <- vapply(seq_along(ccl), function(g) {
  M <- ccl[[g]]
  if (!is.null(M) && group_coef %in% rownames(M) && int_coef %in% colnames(M))
    as.numeric(M[group_coef, int_coef]) else NA_real_
}, numeric(1))
var_M <- se_int^2 + se_F^2 + 2 * cov_F_int
se_M  <- sqrt(pmax(var_M, .Machine$double.eps))

t_int <- beta_int / se_int
p_int <- 2 * pnorm(-abs(t_int))
q_int <- p.adjust(p_int, method = "BH")

ash_F <- tryCatch(ash(betahat = beta_F, sebetahat = se_F, mixcompdist = "normal"),
                  error = function(e) NULL)
ash_M <- tryCatch(ash(betahat = beta_M, sebetahat = se_M, mixcompdist = "normal"),
                  error = function(e) NULL)
if (is.null(ash_F) || is.null(ash_M)) emit_na_and_exit("ash_fit_failed")
lfsr_F <- ash_F$result$lfsr
lfsr_M <- ash_M$result$lfsr

# v5 decision tree, q_int<0.20 stand-in for padj_int_5k<0.20 (per-rep no Tier-1 univ.)
classify_v5_rep <- function(bF, bM, bI, qI, lfF, lfM) {
  n <- length(bF)
  cls <- rep(NA_character_, n)
  sig_int <- !is.na(qI) & qI < 0.20
  same_sign <- !is.na(bF) & !is.na(bM) & bF != 0 & bM != 0 & sign(bF) == sign(bM)
  brA <- !is.na(qI) & qI >= 0.20
  conc_both    <- brA & lfF < 0.05 & lfM < 0.05 & same_sign
  conc_singleA <- brA & ((lfF < 0.05) | (lfM < 0.05)) & !conc_both
  not_deg      <- brA & !conc_both & !conc_singleA
  cls[conc_both] <- "Concordant"
  cls[conc_singleA] <- "Concordant"
  cls[not_deg] <- "Not_DEG"
  brB <- sig_int
  is_div <- brB & lfF < 0.05 & lfM < 0.05 & !same_sign
  cls[is_div] <- "Divergent"
  is_div_one <- brB & !is_div & !same_sign &
                ((lfF < 0.05 & lfM >= 0.05) | (lfM < 0.05 & lfF >= 0.05))
  cls[is_div_one] <- "Divergent_one_sided"
  brB3 <- brB & same_sign & is.na(cls)
  ratio_F <- abs(bF) / pmax(abs(bM), .Machine$double.eps)
  ratio_M <- abs(bM) / pmax(abs(bF), .Machine$double.eps)
  is_F     <- brB3 & ratio_F >= 2 & lfF < 0.05 & lfM > 0.20
  is_M     <- brB3 & ratio_M >= 2 & lfM < 0.05 & lfF > 0.20
  is_F_pu  <- brB3 & ratio_F >= 2 & lfF < 0.05 & lfM >= 0.05 & lfM <= 0.20
  is_M_pu  <- brB3 & ratio_M >= 2 & lfM < 0.05 & lfF >= 0.05 & lfF <= 0.20
  cls[is_F]    <- "Female_biased"
  cls[is_M]    <- "Male_biased"
  cls[is_F_pu] <- "Female_biased_M_underpowered"
  cls[is_M_pu] <- "Male_biased_F_underpowered"
  is_modif <- brB & is.na(cls)
  cls[is_modif] <- "Sex_modifier"
  cls[is.na(cls)] <- "Not_DEG"
  cls
}
class_v6 <- classify_v5_rep(beta_F, beta_M, beta_int, q_int, lfsr_F, lfsr_M)

out <- data.table(gene = names(beta_F),
                  beta_F = beta_F, beta_M = beta_M, beta_int = beta_int,
                  p_int = p_int, q_int = q_int,
                  lfsr_F = lfsr_F, lfsr_M = lfsr_M,
                  class_v6_rep = class_v6, variant_used = variant_used)
.tmp_path <- paste0(OUT_PATH, ".tmp.", Sys.getpid())
fwrite(out, .tmp_path); file.rename(.tmp_path, OUT_PATH)
cat("Wrote:", OUT_PATH, "  rows:", nrow(out), "\n")
cat("Class distribution:\n"); print(table(out$class_v6_rep, useNA = "ifany"))
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
