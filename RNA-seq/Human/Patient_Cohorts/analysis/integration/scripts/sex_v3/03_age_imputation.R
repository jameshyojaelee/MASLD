#!/usr/bin/env Rscript
# sex_v3/03_age_imputation.R
# ---------------------------------------------------------------------------
# Module 03 — Multiple imputation of missing age via mice::mice (m = 5,
# method = "pmm"), conditional on: dataset, group_binary, inferred_sex,
# diagnosis_harmonized, and the 4 composition fractions.
#
# Reads:  intermediates/sex_v3_input.rds
# Writes: intermediates/age_mi.rds
#         intermediates/age_mi_qc.csv  (per-imputation mean ± sd vs observed)
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(mice)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RDIR   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                    "results/integration")
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

# Shared utilities (atomic writes + sessionInfo dump) -- R5 Issue 2 fix
source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

cat("============================================================\n")
cat("Module 03 — Multiple imputation of age (mice, m=5, pmm)\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")

# ---------------------------------------------------------------------------
# Load input bundle
# ---------------------------------------------------------------------------
in_path <- file.path(IDIR, "sex_v3_input.rds")
inp <- readRDS(in_path)
info <- inp$meta
n <- nrow(info)
cat("Samples:", n, "\n")
cat("Observed age non-NA:", sum(!is.na(info$age)), "/", n,
    "  missingness:", sprintf("%.1f%%", 100 * mean(is.na(info$age))), "\n")
cat("Observed age summary:\n"); print(summary(info$age))
cat("Observed age mean by dataset:\n")
print(round(aggregate(age ~ dataset, data = info, FUN = function(x)
  c(mean = mean(x, na.rm = TRUE), n = sum(!is.na(x))))$age, 2))

# ---------------------------------------------------------------------------
# Build mice input data frame
# ---------------------------------------------------------------------------
mi_df <- data.frame(
  age            = info$age,
  dataset        = info$dataset,
  group_binary   = info$group_binary,
  inferred_sex   = info$inferred_sex,
  # 2026-05-14: diagnosis_harmonized dropped — was miscoded with "" for
  # 318 Disease samples and collinear with group_binary. Age MI now
  # conditions on dataset + group_binary + inferred_sex + composition only.
  Hepatocytes    = info$Hepatocytes,
  Macrophages    = info$Macrophages,
  Endothelial    = info$Endothelial,
  Cholangiocytes = info$Cholangiocytes
)
cat("\nmice input frame:", nrow(mi_df), "rows x", ncol(mi_df), "cols\n")

# Confirm: only `age` should have missingness (other covariates already
# filtered in Module 01). If any covariate is also NA, mice will impute it
# too — we want to control that. Set imputation method to pmm for `age` only,
# "" (skip) for everything else.
methods_vec <- mice::make.method(mi_df)
methods_vec[]            <- ""        # skip all
methods_vec["age"]       <- "pmm"
cat("\nmice methods vector:\n"); print(methods_vec)

# Predictor matrix: use all other columns to predict age; ensure age isn't
# used to predict itself.
pred_mat <- mice::make.predictorMatrix(mi_df)
pred_mat[]      <- 0
pred_mat["age", ] <- 1
pred_mat["age", "age"] <- 0
cat("\nPredictor matrix (row=target, col=predictor):\n")
print(pred_mat)

# ---------------------------------------------------------------------------
# Run mice (m = 5, seed = 42)
# ---------------------------------------------------------------------------
cat("\nRunning mice(m = 5, method = 'pmm', seed = 42)...\n")
imp <- mice(mi_df,
            m              = 5,
            method         = methods_vec,
            predictorMatrix = pred_mat,
            printFlag      = TRUE,
            seed           = 42)

# ---------------------------------------------------------------------------
# Extract 5 complete vectors
# ---------------------------------------------------------------------------
age_list <- lapply(seq_len(5), function(k) {
  v <- mice::complete(imp, action = k)$age
  stopifnot(length(v) == n)
  stopifnot(!any(is.na(v)))
  v
})
names(age_list) <- paste0("age_imputed_", seq_len(5))

# ---------------------------------------------------------------------------
# Sanity QC: per-imputation mean within ±5 years of observed mean per cohort
# ---------------------------------------------------------------------------
cat("\nSanity QC — imputed vs observed mean per cohort:\n")
qc <- data.table()
for (ds in levels(info$dataset)) {
  obs <- info$age[info$dataset == ds]
  obs_mean <- mean(obs, na.rm = TRUE)
  for (k in seq_len(5)) {
    imp_vec  <- age_list[[k]]
    sub_imp  <- imp_vec[info$dataset == ds]
    diff_m   <- mean(sub_imp) - obs_mean
    qc <- rbind(qc, data.table(
      dataset     = ds,
      imputation  = k,
      n           = length(obs),
      n_missing   = sum(is.na(obs)),
      obs_mean    = round(obs_mean, 2),
      imputed_mean= round(mean(sub_imp), 2),
      delta       = round(diff_m, 2),
      within_5yrs = abs(diff_m) <= 5
    ))
  }
}
print(qc)

write_atomic_csv(qc, file.path(IDIR, "age_mi_qc.csv"))
cat("\nSaved QC (atomic):", file.path(IDIR, "age_mi_qc.csv"), "\n")

n_fail <- sum(!qc$within_5yrs, na.rm = TRUE)
n_na <- sum(is.na(qc$within_5yrs))
if (n_na > 0) {
  cat("\n! ", n_na, "imputation-cohort cells had no observed age to compare (entire cohort missing age) — skipped from QC bound check.\n")
}
if (isTRUE(n_fail > 0)) {
  cat("\n! ", n_fail, "imputation-cohort cells exceeded the ±5yrs sanity bound.\n")
  cat("   See qc table; review before downstream modules consume.\n")
} else {
  cat("\nOK — all comparable imputations within ±5 yrs of observed cohort means.\n")
}

# ---------------------------------------------------------------------------
# Save imputed age list
# ---------------------------------------------------------------------------
out_path <- file.path(IDIR, "age_mi.rds")
write_atomic_rds(list(
  age_list = age_list,
  m        = 5,
  method   = "pmm",
  seed     = 42,
  sample_ids = rownames(info),
  qc       = qc,
  mice_obj_summary = list(method = imp$method,
                          predictorMatrix = imp$predictorMatrix,
                          nmis  = imp$nmis,
                          loggedEvents = imp$loggedEvents)
), out_path)
cat("Saved (atomic):", out_path, "\n")

# R5 Issue 2 fix: sessionInfo dump for per-module reproducibility audit
dump_session_info(IDIR, "03")

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
