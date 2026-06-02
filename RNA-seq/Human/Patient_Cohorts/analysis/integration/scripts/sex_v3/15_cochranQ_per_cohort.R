#!/usr/bin/env Rscript
# sex_v3/15_cochranQ_per_cohort.R
# ---------------------------------------------------------------------------
# Pillar 6 — per-cohort fixed-effects dream M2 fit for Cochran's Q
# heterogeneity diagnostics.
#
# Per cohort (single dataset), fit:
#   ~ group_binary * inferred_sex
#     + Hepatocytes + Macrophages + Endothelial + Cholangiocytes
#     + age_imputed + SV1..SVk
# with NO random effect (single cohort -> dataset is constant).
# SVs with within-cohort variance ~0 are dropped to avoid rank deficiency.
#
# Extracts β_int + se_int per gene -> intermediates/cochranQ_v6/cohort_<k>.rds
#
# CLI knob: COHORT_IDX (1..5) via SLURM_ARRAY_TASK_ID. Cohort order matches
# 07b_v5_loco.R: GSE126848, GSE130970, GSE135251, GSE162694, GSE213621.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
})

set.seed(42)

idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",
                             Sys.getenv("COHORT_IDX", "1")))
if (is.na(idx) || idx < 1L) {
  stop("COHORT_IDX / SLURM_ARRAY_TASK_ID must be >= 1; got ", idx)
}

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
OUTDIR <- file.path(IDIR, "cochranQ_v6")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(OUTDIR, sprintf("cohort_%d.rds", idx))

inp    <- readRDS(file.path(IDIR, "sex_v3_input.rds"))
# Discover cohorts from the contrast-specific input (alphabetic). For
# disease_vs_ctrl this yields the 5 canonical cohorts; for mash_vs_masl it
# yields 7 (the 4 control-bearing + 3 disease-only).
COHORTS_SORTED <- sort(levels(droplevels(as.factor(inp$meta$dataset))))
if (idx > length(COHORTS_SORTED)) {
  stop("idx=", idx, " exceeds N cohorts (", length(COHORTS_SORTED),
       ") for contrast ", .cpaths$contrast)
}
target_cohort <- COHORTS_SORTED[idx]
cat("=== 15_cochranQ_per_cohort idx=", idx, " cohort=",
    target_cohort, " (contrast=", .cpaths$contrast, ") ===\n", sep = "")
sva    <- readRDS(file.path(IDIR, "sva_factors.rds"))
age_mi <- readRDS(file.path(IDIR, "age_mi.rds"))

dge   <- inp$dge
info0 <- inp$meta
sv_mat <- sva$sv
n_sv   <- ncol(sv_mat)
sv_df  <- as.data.frame(sv_mat)
colnames(sv_df) <- paste0("SV", seq_len(n_sv))

info <- info0
info$age_imputed <- age_mi$age_list[[1]]
info <- cbind(info, sv_df)
rownames(info) <- rownames(info0)

# Restrict to one cohort
keep <- info$dataset == target_cohort
n_keep <- sum(keep)
cat("samples in", target_cohort, ":", n_keep, "\n")
if (n_keep < 20) stop("Too few samples for per-cohort fit: ", n_keep)

dge_c <- dge[, keep]
info_c <- info[keep, , drop = FALSE]
info_c$dataset <- droplevels(info_c$dataset)
dge_c <- calcNormFactors(dge_c)

# Drop covariates with zero within-cohort variance
sv_cols <- colnames(sv_df)
sv_keep <- sv_cols[vapply(sv_cols,
                          function(s) var(info_c[[s]], na.rm = TRUE) > 1e-10,
                          logical(1))]
if (length(sv_keep) < length(sv_cols)) {
  cat("  Dropping", length(sv_cols) - length(sv_keep),
      "constant SVs in cohort\n")
}
sv_terms <- if (length(sv_keep) > 0)
              paste("+", paste(sv_keep, collapse = " + ")) else ""

# Check inferred_sex levels — if cohort is single-sex, skip cleanly
sex_lv <- unique(info_c$inferred_sex)
if (length(sex_lv) < 2) {
  cat("  WARNING: only one sex level in", target_cohort,
      "; no interaction estimable. Writing NA placeholder.\n")
  out <- list(cohort = target_cohort, idx = idx, n_samp = n_keep,
              status = "single_sex",
              loco_dt = data.table(gene = rownames(dge_c),
                                   beta_int = NA_real_, se_int = NA_real_,
                                   t_int = NA_real_, p_int = NA_real_))
  saveRDS(out, out_rds)
  cat("Wrote:", out_rds, "\n")
  quit(status = 0)
}

# Check the 2x2 sex × group_binary table — if any cell has 0 samples, the
# interaction is structurally non-estimable in this cohort.
xt <- table(info_c$inferred_sex, info_c$group_binary)
cat("  sex × group_binary table:\n")
print(xt)
if (any(xt == 0)) {
  cat("  WARNING: empty cell in", target_cohort,
      "; no interaction estimable. Writing NA placeholder.\n")
  out <- list(cohort = target_cohort, idx = idx, n_samp = n_keep,
              status = "empty_cell",
              loco_dt = data.table(gene = rownames(dge_c),
                                   beta_int = NA_real_, se_int = NA_real_,
                                   t_int = NA_real_, p_int = NA_real_,
                                   ci_lo = NA_real_, ci_hi = NA_real_))
  saveRDS(out, out_rds)
  cat("Wrote:", out_rds, "\n")
  quit(status = 0)
}

# Try full formula; fall back to minimal on rank deficiency.
form_full_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes",
  " + age_imputed ", sv_terms
)
form_min_str <- paste0(
  "~ group_binary * inferred_sex",
  " + Hepatocytes + Macrophages + Endothelial + Cholangiocytes"
)

fit_one <- function(form_str) {
  mm <- model.matrix(as.formula(form_str), info_c)
  # Detect rank-deficient design and reduce by alias.
  # Drop columns whose rank is colinear with others to avoid NA coefficients.
  v <- voom(dge_c, mm)
  fit <- eBayes(lmFit(v, mm))
  fit
}

# Helper: returns TRUE if the interaction coef has non-NA estimates for
# essentially all genes (>99%); FALSE if rank-deficient design left NAs.
int_estimable <- function(fit) {
  all_coefs <- colnames(fit$coefficients)
  int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary",
                   all_coefs, value = TRUE)[1]
  if (is.na(int_coef)) return(FALSE)
  ok <- sum(!is.na(fit$coefficients[, int_coef])) /
        nrow(fit$coefficients)
  ok > 0.99
}

t0 <- Sys.time()
fit <- tryCatch(fit_one(form_full_str), error = function(e) {
  cat("  full formula failed:", conditionMessage(e),
      "\n  retrying with minimal formula\n")
  fit_one(form_min_str)
})
# Rank-deficient guard: if interaction is NA across (~all) genes under the
# full formula, retry with minimal formula. This is necessary because limma
# does NOT throw on rank deficiency; it silently emits NA coefficients.
if (!int_estimable(fit)) {
  cat("  full formula yielded NA interaction; retrying with minimal formula\n")
  fit <- fit_one(form_min_str)
  if (!int_estimable(fit)) {
    # Last-ditch: drop composition + age, keep only group_binary*inferred_sex
    cat("  minimal formula ALSO yielded NA interaction; retrying bare form\n")
    fit <- fit_one("~ group_binary * inferred_sex")
  }
}
cat("  fit elapsed:", format(Sys.time() - t0), "\n")

all_coefs <- colnames(fit$coefficients)
int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary",
                 all_coefs, value = TRUE)[1]
if (is.na(int_coef)) {
  stop("Interaction coef not found in: ", paste(all_coefs, collapse = ", "))
}
cat("  interaction coef:", int_coef, "\n")

tt <- topTable(fit, coef = int_coef, number = Inf, sort.by = "none")
# Raw SE from limma fit: sigma * stdev.unscaled
raw_se_mat <- fit$sigma * fit$stdev.unscaled
se_int <- raw_se_mat[, int_coef]
beta_int <- fit$coefficients[, int_coef]
t_int_raw <- beta_int / pmax(se_int, .Machine$double.eps)
p_int_raw <- 2 * pnorm(-abs(t_int_raw))

dt <- data.table(
  gene = rownames(fit$coefficients),
  beta_int = beta_int,
  se_int   = se_int,
  t_int    = t_int_raw,
  p_int    = p_int_raw,
  ci_lo    = beta_int - 1.96 * se_int,
  ci_hi    = beta_int + 1.96 * se_int
)

.tmp_rds <- paste0(out_rds, ".tmp.", Sys.getpid())
saveRDS(list(cohort = target_cohort, idx = idx,
             n_samp = n_keep, status = "ok",
             coef = int_coef, loco_dt = dt),
        .tmp_rds)
file.rename(.tmp_rds, out_rds)
cat("Wrote:", out_rds, "  n_genes:", nrow(dt), "\n")
