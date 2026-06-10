#!/usr/bin/env Rscript
# sex_v3/11_perm_fdr.R
# ---------------------------------------------------------------------------
# Pillar 2 — Permutation FDR (per-rep worker).
# LVQW re-engineering 2026-06-08: dataset random->fixed.
#   Strategy: cached-voom permutation. Reads P1's dream_M2_v3_random.rds for the
#   cached LVQW voom output ($v_voom from voomWithQualityWeights) and refits the
#   LVQW lmFit + eBayes with `inferred_sex` shuffled within each
#   (cohort × group_binary) stratum. Voom is approximately invariant to sex
#   permutation under this stratification (marginal counts preserved within
#   strata) so caching voom is safe; only the design matrix is rebuilt per rep.
#
# Fallback: if v_voom is NULL or load fails, recompute voom fresh per rep using
# the LVQW fixed-effects formula on real data, then refit with permuted labels.
#
# Inputs (cached from P1):
#   dream_M2_v3_random.rds  — list with v_voom (EList), info_template,
#                              re_variant_used, form_used, design
#
# Output (per rep):
#   intermediates/perm_v6/rep_{REP}.csv with columns:
#     gene, t_int_perm, p_int_perm
# ---------------------------------------------------------------------------
# LVQW re-engineering 2026-06-08: dataset random->fixed.
suppressPackageStartupMessages({
  library(limma); library(data.table); library(edgeR)
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
# LVQW re-engineering 2026-06-08: dataset random->fixed (limma is single-threaded).
cat("CPU cores:", ncpus, "(limma LVQW is single-threaded)\n")

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

# ---- Build the LVQW fixed-effects formula — mirror whatever P1 fit so the
# permutation null mirrors the parametric null.
# LVQW re-engineering 2026-06-08: dataset random->fixed.
sv_cols <- grep("^SV", colnames(info_perm_df), value = TRUE)
sv_terms <- paste(sv_cols, collapse = " + ")
form_fixed_str <- paste0(
  "~ dataset",
  " + group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed + ", sv_terms
)
form_fixed <- as.formula(form_fixed_str)

re_variant_p1 <- if (!is.null(p1$re_variant_used)) p1$re_variant_used else "lvqw_fixed"
form_primary <- if (!is.null(p1$form_used)) p1$form_used else form_fixed
cat("[2] Primary formula = P1 variant '", re_variant_p1, "' (LVQW fixed-effects)\n", sep = "")

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
  # LVQW re-engineering 2026-06-08: dataset random->fixed.
  # voomWithQualityWeights takes a model.matrix (built from the REAL labels) —
  # the voom weights are estimated under the real design and reused for the
  # permuted lmFit (weights are approx. invariant to within-stratum sex shuffle).
  design_real <- model.matrix(form_primary, data = info_real)
  t1 <- Sys.time()
  v_voom <- tryCatch(
    suppressWarnings(voomWithQualityWeights(inp$dge, design_real)),
    error = function(e) NULL)
  cat("  voom recompute elapsed:", format(Sys.time() - t1), "\n")
  if (is.null(v_voom)) emit_na_and_exit("voom_recompute_failed")
}
rm(p1); invisible(gc())

# ---- Refit LVQW (lmFit + eBayes) with permuted labels ----
# LVQW re-engineering 2026-06-08: dataset random->fixed.
# Build the permuted design (inferred_sex shuffled within strata) and refit on
# the cached voom EList. The permuted design's interaction column carries the
# null sex-by-disease signal.
cat("[3] lmFit + eBayes with permuted labels (", re_variant_p1, ")...\n", sep = "")
design_perm <- tryCatch(model.matrix(form_primary, data = info_perm_df),
                        error = function(e) NULL)
if (is.null(design_perm)) emit_na_and_exit("perm_design_build_failed")
if (qr(design_perm)$rank < ncol(design_perm))
  emit_na_and_exit("perm_design_rank_deficient")
# Align voom EList sample columns to the permuted design rows (defensive).
if (!is.null(rownames(design_perm)) && !is.null(colnames(v_voom)) &&
    !identical(rownames(design_perm), colnames(v_voom))) {
  common <- intersect(colnames(v_voom), rownames(design_perm))
  if (length(common) < 5) emit_na_and_exit("voom_design_sample_mismatch")
  v_voom <- v_voom[, common]
  design_perm <- design_perm[common, , drop = FALSE]
}
t0 <- Sys.time()
fit_perm <- tryCatch(
  suppressWarnings(eBayes(lmFit(v_voom, design_perm))),
  error = function(e) { cat("  lmFit permuted failed:", conditionMessage(e), "\n"); NULL }
)
cat("  lmFit+eBayes permuted elapsed:", format(Sys.time() - t0), "\n")
if (is.null(fit_perm)) emit_na_and_exit("lmfit_permuted_failed")

# ---- Extract permuted t_int / p_int (RAW Wald, matching the observed run) ----
all_coefs <- colnames(fit_perm$coefficients)
int_coef <- grep(":", grep("^group_binary", all_coefs, value = TRUE), value = TRUE)[1]
if (is.na(int_coef)) emit_na_and_exit("int_coef_not_found")
beta_int <- fit_perm$coefficients[, int_coef]
raw_se   <- fit_perm$stdev.unscaled * fit_perm$sigma
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
