#!/usr/bin/env Rscript
# sex_v3/11_perm_fdr.R
# ---------------------------------------------------------------------------
# Pillar 2 — Permutation FDR (per-rep worker).
# Strategy: cached-voom permutation. Reads P1's dream_M2_v3_random.rds for the
# cached voom output ($v_voom) and refits dream LMM (F2 random-slope formula)
# with `inferred_sex` shuffled within each (cohort × group_binary) stratum.
# Voom is invariant to sex permutation under this stratification (marginal
# counts preserved within strata) so caching voom is safe.
#
# Fallback: if v_voom is NULL (P1 used F5/F6) or load fails, recompute voom
# fresh per rep using the F2 formula on real data, then refit with permuted
# labels. Slower but always correct.
#
# Inputs (cached from P1):
#   dream_M2_v3_random.rds  — list with v_voom (TestResults / EList), info_template,
#                              re_variant_used
#
# Output (per rep):
#   intermediates/perm_v6/rep_{REP}.csv with columns:
#     gene, t_int_perm, p_int_perm
# ---------------------------------------------------------------------------
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
  library(variancePartition); library(BiocParallel)
})

REP <- as.integer(Sys.getenv("REP",
                  Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
stopifnot(is.integer(REP), REP >= 1)
set.seed(20260515L + REP)
cat("=== sex_v6 P2 permutation rep", REP, "===\n")
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
PERM_DIR <- file.path(IDIR, "perm_v6")
dir.create(PERM_DIR, showWarnings = FALSE, recursive = TRUE)
P1_RDS <- file.path(SEXV3, "dream_M2_v3_random.rds")
OUT_PATH <- sprintf("%s/rep_%03d.csv", PERM_DIR, REP)

emit_na_and_exit <- function(reason) {
  cat("FAILURE reason:", reason, "— writing NA placeholder and exiting 0\n")
  na_df <- data.table(gene = NA_character_, t_int_perm = NA_real_,
                      p_int_perm = NA_real_, failure_reason = reason)
  fwrite(na_df, OUT_PATH); quit(save = "no", status = 0)
}

ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
param <- if (ncpus > 1) MulticoreParam(workers = ncpus, RNGseed = 42L + REP) else SerialParam()
cat("CPU cores:", ncpus, "\n")

# ---- Load P1 cache ----
if (!file.exists(P1_RDS)) emit_na_and_exit("p1_rds_missing")
cat("[1] Loading P1 cache (may be large)...\n")
t0 <- Sys.time()
p1 <- tryCatch(readRDS(P1_RDS), error = function(e) NULL)
cat("  loaded in", format(Sys.time() - t0), "\n")
if (is.null(p1)) emit_na_and_exit("p1_load_failed")
info <- p1$info_template
re_variant <- p1$re_variant_used
cat("P1 re_variant_used:", re_variant, "\n")
stopifnot(!is.null(info))

# ---- Permute inferred_sex within (dataset × group_binary) ----
info_perm <- as.data.table(info, keep.rownames = "sample_id")
info_perm[, inferred_sex_perm := inferred_sex]
strata_keys <- unique(info_perm[, .(dataset, group_binary)])
for (i in seq_len(nrow(strata_keys))) {
  ds <- strata_keys$dataset[i]; gr <- strata_keys$group_binary[i]
  idx <- which(info_perm$dataset == ds & info_perm$group_binary == gr)
  if (length(idx) > 1) {
    info_perm$inferred_sex_perm[idx] <- sample(info_perm$inferred_sex[idx])
  }
}
# Sanity: marginal counts preserved within strata
chk <- info_perm[, .(orig_F = sum(inferred_sex == "F"),
                     perm_F = sum(inferred_sex_perm == "F")),
                 by = .(dataset, group_binary)]
if (!all(chk$orig_F == chk$perm_F)) emit_na_and_exit("stratum_count_mismatch")
info_perm[, inferred_sex := factor(inferred_sex_perm,
                                   levels = levels(info$inferred_sex))]
info_perm[, inferred_sex_perm := NULL]
info_perm_df <- as.data.frame(info_perm)
rownames(info_perm_df) <- info_perm$sample_id
info_perm_df$sample_id <- NULL

# ---- Build formulas — use whichever variant P1 actually fit, so the
# permutation null mirrors the parametric null. Earlier this script hardcoded
# F2 regardless of P1's fallback ladder; if P1 fell back to F5/F6, the
# perm-FDR fits would either silently fail or use a different model than the
# observed run.
sv_cols <- grep("^SV", colnames(info_perm_df), value = TRUE)
sv_terms <- paste(sv_cols, collapse = " + ")
form_F2_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 + group_binary + inferred_sex + group_binary:inferred_sex || dataset)"
)
form_F2 <- as.formula(form_F2_str)
form_F5_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms,
  " + (1 + group_binary | dataset)"
)
form_F5 <- as.formula(form_F5_str)

re_variant_p1 <- if (!is.null(p1$re_variant_used)) p1$re_variant_used else "F2"
form_primary <- if (!is.null(p1$form_used)) p1$form_used else
                switch(re_variant_p1, F2 = form_F2, F5 = form_F5, form_F2)
cat("[2] Primary formula = P1 variant '", re_variant_p1, "'\n", sep = "")

# ---- Resolve voom: cached vs per-rep fresh ----
v_voom <- p1$v_voom
if (is.null(v_voom)) {
  cat("[2] v_voom is NULL in P1 cache; recomputing voom on real data with primary formula...\n")
  # Need DGEList — load from canonical input
  inp <- tryCatch(readRDS(file.path(IDIR, "sex_v3_input.rds")),
                  error = function(e) NULL)
  if (is.null(inp)) emit_na_and_exit("dge_load_failed_no_voom_cache")
  info_real <- as.data.frame(info)
  rownames(info_real) <- rownames(info)
  t1 <- Sys.time()
  v_voom <- tryCatch(
    suppressWarnings(voomWithDreamWeights(inp$dge, form_primary, info_real,
                                          BPPARAM = param, useWeights = TRUE)),
    error = function(e) NULL)
  cat("  voom recompute elapsed:", format(Sys.time() - t1), "\n")
  if (is.null(v_voom)) emit_na_and_exit("voom_recompute_failed")
}
rm(p1); invisible(gc())

# ---- Refit dream with permuted labels (primary formula from P1; fallback F5) ----
cat("[3] dream with permuted labels (", re_variant_p1, ")...\n", sep = "")
t0 <- Sys.time()
fit_perm <- tryCatch(
  suppressWarnings(dream(v_voom, form_primary, info_perm_df,
                         BPPARAM = param, useWeights = TRUE)),
  error = function(e) { cat("  dream primary failed:", conditionMessage(e), "\n"); NULL }
)
cat("  dream primary elapsed:", format(Sys.time() - t0), "\n")
if (is.null(fit_perm) && re_variant_p1 != "F5") {
  cat("[3b] Falling back to F5 with permuted labels...\n")
  t0 <- Sys.time()
  fit_perm <- tryCatch(
    suppressWarnings(dream(v_voom, form_F5, info_perm_df,
                           BPPARAM = param, useWeights = TRUE)),
    error = function(e) NULL)
  cat("  dream F5 elapsed:", format(Sys.time() - t0), "\n")
}
if (is.null(fit_perm)) emit_na_and_exit("dream_permuted_failed")

# ---- Extract permuted t_int / p_int ----
all_coefs <- colnames(fit_perm$coefficients)
int_coef <- grep(":", grep("^group_binary", all_coefs, value = TRUE), value = TRUE)[1]
if (is.na(int_coef)) emit_na_and_exit("int_coef_not_found")
beta_int <- fit_perm$coefficients[, int_coef]
raw_se   <- fit_perm$sigma * fit_perm$stdev.unscaled
se_int   <- raw_se[, int_coef]
t_int    <- beta_int / se_int
p_int    <- 2 * pnorm(-abs(t_int))

out <- data.table(gene = names(beta_int),
                  t_int_perm = t_int,
                  p_int_perm = p_int)
# Atomic write so a partial per-rep CSV (truncated by node kill or oom-kill)
# cannot be misread by the aggregator as a "complete" rep with NA tails.
.tmp_path <- paste0(OUT_PATH, ".tmp.", Sys.getpid())
fwrite(out, .tmp_path); file.rename(.tmp_path, OUT_PATH)
cat("Wrote:", OUT_PATH, "  rows:", nrow(out), "\n")
cat("  median |t_perm|:", round(median(abs(t_int), na.rm = TRUE), 3), "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
