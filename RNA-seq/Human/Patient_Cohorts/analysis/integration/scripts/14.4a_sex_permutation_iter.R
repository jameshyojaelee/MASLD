#!/usr/bin/env Rscript
# 14.4a_sex_permutation_iter.R
# ---------------------------------------------------------------------------
# ONE permutation iteration for sex-label shuffling null distribution.
#
# Shuffles inferred_sex labels within each dataset (preserving per-dataset
# sex ratio) then re-runs male-only and female-only dream analyses.
# Parameterized by env var PERM_ITER (int 1-50).
#
# Output (in OUT_DIR):
#   perm_class_i{PERM}.csv   — all genes: gene, sex_class, logFC_M/F, padj_M/F
#   perm_summary_i{PERM}.csv — 1 row: counts per sex_class
# ---------------------------------------------------------------------------

t0 <- proc.time()

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(yaml)
  library(edgeR)
})

# Force injection into lme4 namespace BEFORE loading variancePartition
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
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
RDIR <- file.path(INT, "results/integration")
OUT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_subsampling/permutations"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parallel setup
# ---------------------------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# ---------------------------------------------------------------------------
# Read permutation index
# ---------------------------------------------------------------------------
PERM <- as.integer(Sys.getenv("PERM_ITER", "1"))
PERM_LABEL <- sprintf("%02d", PERM)
cat("Permutation iteration:", PERM, "(label:", PERM_LABEL, ")\n")

# ============================================================================
# Load data (same as script 26)
# ============================================================================
cat("\n===== Loading data =====\n")

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Exclude datasets with irrecoverable confounds (same as script 05/26)
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Mega cohorts (yaml):", paste(mega_cohorts, collapse = ", "), "
")
cat("Samples for mega-analysis:", ncol(dge_mega), "\n")

# Filter to samples that passed the sex check (pass_sex == TRUE)
qc_report <- fread(file.path(INT, "qc/sample_qc_report.csv"))
sex_pass_ids <- qc_report[pass_sex == TRUE, sample_id]
sex_fail <- !colnames(dge_mega) %in% sex_pass_ids
if (any(sex_fail)) {
  cat("Removing", sum(sex_fail), "samples that failed sex check (pass_sex=FALSE)\n")
  dge_mega <- dge_mega[, !sex_fail]
}
cat("Samples after sex check filter:", ncol(dge_mega), "\n")

# Load metadata with inferred sex
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

# Build info data.frame
info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_mega)

# Drop samples with NA sex
na_sex <- is.na(info$inferred_sex)
cat("Samples with NA inferred_sex:", sum(na_sex), "\n")
if (any(na_sex)) {
  dge_mega <- dge_mega[, !na_sex]
  info     <- info[!na_sex, , drop = FALSE]
  cat("Samples remaining after dropping NA sex:", ncol(dge_mega), "\n")
}

cat("Original sex distribution:\n")
print(table(info$inferred_sex))
cat("Group x sex distribution (before shuffling):\n")
print(table(info$group_binary, info$inferred_sex))

# ============================================================================
# Shuffle sex labels within each dataset
# ============================================================================
cat("\n===== Shuffling sex labels (PERM =", PERM, ") =====\n")

set.seed(42 + PERM * 7919)
for (ds in unique(info$dataset)) {
  idx <- which(info$dataset == ds)
  info$inferred_sex[idx] <- sample(info$inferred_sex[idx])
}

cat("Shuffled sex distribution:\n")
print(table(info$inferred_sex))
cat("Group x shuffled sex distribution:\n")
print(table(info$group_binary, info$inferred_sex))

# ============================================================================
# Run dream for each sex stratum (same logic as script 26 Part 2)
# ============================================================================
cat("\n===== Running sex-stratified dream (shuffled labels) =====\n")

run_stratum_dream <- function(sex_label, dge_full, info_full, param) {
  cat("\n--- Stratum:", sex_label, "---\n")

  info_s <- info_full[info_full$inferred_sex == sex_label, , drop = FALSE]
  sample_ids_s <- rownames(info_s)
  dge_s <- dge_full[, colnames(dge_full) %in% sample_ids_s]

  cat("  Samples:", ncol(dge_s), "\n")
  cat("  Group distribution:\n")
  print(table(info_s$group_binary, info_s$dataset))

  # Require both groups present in at least one dataset
  n_ctrl <- sum(info_s$group_binary == "Control")
  n_dis  <- sum(info_s$group_binary == "Disease")
  if (n_ctrl == 0 || n_dis == 0) {
    warning("  Stratum ", sex_label, " missing controls or disease samples — skipping.")
    return(NULL)
  }

  # Drop datasets with only one group (would cause singular fit)
  ds_counts <- table(info_s$dataset, info_s$group_binary)
  keep_ds <- rownames(ds_counts)[apply(ds_counts, 1, function(x) all(x > 0))]
  if (length(keep_ds) < nrow(ds_counts)) {
    dropped <- setdiff(rownames(ds_counts), keep_ds)
    cat("  Dropping single-group datasets:", paste(dropped, collapse = ", "), "\n")
    info_s <- info_s[info_s$dataset %in% keep_ds, , drop = FALSE]
    info_s$dataset <- droplevels(info_s$dataset)
    dge_s  <- dge_s[, colnames(dge_s) %in% rownames(info_s)]
  }

  # Require at least 2 datasets for random effect
  n_ds <- nlevels(info_s$dataset)
  if (n_ds < 2) {
    warning("  Only ", n_ds, " dataset(s) in stratum ", sex_label,
            " — random effect not estimable; skipping.")
    return(NULL)
  }

  form_s <- ~ group_binary + (1|dataset)
  cat("  Formula:", deparse(form_s), "\n")

  # CRITICAL: Do NOT re-run filterByExpr() or calcNormFactors()
  v_s <- suppressWarnings(
    voomWithDreamWeights(dge_s, form_s, info_s, BPPARAM = param)
  )
  fit_s <- tryCatch(
    suppressWarnings(
      dream(v_s, form_s, info_s, BPPARAM = param)
    ),
    error = function(e) {
      cat("  ERROR in dream for stratum", sex_label, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(fit_s)) return(NULL)

  res_s <- topTable(fit_s, coef = "group_binaryDisease", number = Inf, sort.by = "none")
  res_s$gene <- rownames(res_s)

  res_s_dt <- as.data.table(res_s)
  setnames(res_s_dt, "adj.P.Val", "padj")

  n_sig <- nrow(res_s_dt[padj < 0.05 & abs(logFC) > 0.25])
  cat("  DEGs padj < 0.05 & |LFC| > 0.25:", n_sig, "\n")

  return(res_s_dt)
}

# Detect sex labels from data
sex_levels <- levels(info$inferred_sex)
cat("Sex factor levels:", paste(sex_levels, collapse = ", "), "\n")
male_label   <- sex_levels[grepl("^M", sex_levels)][1]
female_label <- sex_levels[grepl("^F", sex_levels)][1]
cat("Using male label:", male_label, ", female label:", female_label, "\n")

res_male   <- run_stratum_dream(male_label,   dge_mega, info, param)
res_female <- run_stratum_dream(female_label, dge_mega, info, param)

# ============================================================================
# Sex-DEG classification (same logic as script 26 lines 280-291)
# ============================================================================
cat("\n===== Sex-DEG classification =====\n")

PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.25

if (!is.null(res_male) && !is.null(res_female)) {
  cls_dt <- merge(
    res_male[,   .(gene, logFC_M = logFC, padj_M = padj)],
    res_female[, .(gene, logFC_F = logFC, padj_F = padj)],
    by = "gene", all = TRUE
  )

  cls_dt[, sig_M := (!is.na(padj_M) & padj_M < PADJ_THRESH &
                     !is.na(logFC_M) & abs(logFC_M) > LFC_THRESH)]
  cls_dt[, sig_F := (!is.na(padj_F) & padj_F < PADJ_THRESH &
                     !is.na(logFC_F) & abs(logFC_F) > LFC_THRESH)]
  cls_dt[, same_sign := !is.na(logFC_M) & !is.na(logFC_F) & sign(logFC_M) == sign(logFC_F)]

  cls_dt[, sex_class := fcase(
    sig_M & !sig_F,                          "Male_specific",
    !sig_M & sig_F,                          "Female_specific",
    sig_M & sig_F & same_sign,               "Shared",
    sig_M & sig_F & !same_sign,              "Divergent",
    default = "Not_significant"
  )]

  cat("Sex-class distribution (shuffled):\n")
  print(cls_dt[, .N, by = sex_class][order(-N)])

  # --- Save perm_class ---
  out_class <- cls_dt[, .(gene, sex_class, logFC_M, padj_M, logFC_F, padj_F)]
  fwrite(out_class, file.path(OUT_DIR, paste0("perm_class_i", PERM_LABEL, ".csv")))
  cat("Saved:", paste0("perm_class_i", PERM_LABEL, ".csv"), "\n")

  # --- Save perm_summary ---
  summary_dt <- cls_dt[, .N, by = sex_class]
  summary_wide <- dcast(summary_dt, . ~ sex_class, value.var = "N", fill = 0)
  summary_wide[, `.` := NULL]
  # Ensure all columns present
  for (col in c("Female_specific", "Male_specific", "Shared", "Divergent", "Not_significant")) {
    if (!col %in% names(summary_wide)) summary_wide[, (col) := 0L]
  }
  setnames(summary_wide,
           c("Female_specific", "Male_specific", "Shared", "Divergent", "Not_significant"),
           c("n_Female_specific", "n_Male_specific", "n_Shared", "n_Divergent", "n_Not_significant"))
  fwrite(summary_wide, file.path(OUT_DIR, paste0("perm_summary_i", PERM_LABEL, ".csv")))
  cat("Saved:", paste0("perm_summary_i", PERM_LABEL, ".csv"), "\n")

} else {
  cat("WARNING: one or both strata failed — classification skipped.\n")
  # Write empty summary so downstream aggregation knows this iteration failed
  summary_wide <- data.table(
    n_Female_specific = NA_integer_,
    n_Male_specific   = NA_integer_,
    n_Shared          = NA_integer_,
    n_Divergent       = NA_integer_,
    n_Not_significant = NA_integer_
  )
  fwrite(summary_wide, file.path(OUT_DIR, paste0("perm_summary_i", PERM_LABEL, ".csv")))
  cat("Saved empty summary for failed iteration\n")
}

# ============================================================================
# Summary
# ============================================================================
elapsed <- (proc.time() - t0)["elapsed"]
cat("\n===== PERMUTATION", PERM, "COMPLETE =====\n")
cat("Samples per shuffled sex:\n")
print(table(info$inferred_sex))
if (!is.null(res_male))   cat("Male DEGs (padj<0.05 & |LFC|>0.25):",   nrow(res_male[padj < 0.05 & abs(logFC) > 0.25]),   "\n")
if (!is.null(res_female)) cat("Female DEGs (padj<0.05 & |LFC|>0.25):", nrow(res_female[padj < 0.05 & abs(logFC) > 0.25]), "\n")
cat("Elapsed time:", round(elapsed / 60, 1), "minutes\n")
cat("Done:", as.character(Sys.time()), "\n")
