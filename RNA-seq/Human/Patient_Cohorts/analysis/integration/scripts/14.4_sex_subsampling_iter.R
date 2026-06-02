#!/usr/bin/env Rscript
# 14.4_sex_subsampling_iter.R
# ---------------------------------------------------------------------------
# ONE iteration of sex-stratified dream subsampling sensitivity analysis.
#
# Env vars:
#   SUBSAMPLE_FRAC  — float in [0.30, 1.0)  (default 0.70)
#   SUBSAMPLE_ITER  — integer >= 1           (default 1)
#
# For each iteration:
#   1. Stratified subsampling (dataset x group x sex) at SUBSAMPLE_FRAC
#   2. Validate drawn samples (>=2 datasets with both groups per sex)
#   3. Run dream for male and female strata
#   4. Classify sex_class
#   5. Save 3 CSVs per iteration
#
# Output: .../audit_sensitivity/sex_subsampling/iterations/
#   sub_male_f{FF}_i{II}.csv    — gene, logFC, padj, t
#   sub_female_f{FF}_i{II}.csv  — gene, logFC, padj, t
#   sub_class_f{FF}_i{II}.csv   — gene, sex_class, logFC_M, padj_M, logFC_F, padj_F
# ---------------------------------------------------------------------------

t0 <- proc.time()

# ---------------------------------------------------------------------------
# Library loading + reformulas injection (MUST match script 26 lines 25-50)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(edgeR)
  library(yaml)
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
BASE    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT     <- file.path(BASE, "analysis/integration")
RDIR    <- file.path(INT, "results/integration")
OUT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_subsampling/iterations"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parallel setup (from script 26)
# ---------------------------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# ---------------------------------------------------------------------------
# Read env vars
# ---------------------------------------------------------------------------
FRAC <- as.numeric(Sys.getenv("SUBSAMPLE_FRAC", "0.70"))
ITER <- as.integer(Sys.getenv("SUBSAMPLE_ITER", "1"))

cat("SUBSAMPLE_FRAC:", FRAC, "\n")
cat("SUBSAMPLE_ITER:", ITER, "\n")

# Validate
if (is.na(FRAC) || FRAC < 0.3 || FRAC >= 1.0) {
  stop("SUBSAMPLE_FRAC must be in [0.3, 1.0). Got: ", FRAC)
}
if (is.na(ITER) || ITER < 1) {
  stop("SUBSAMPLE_ITER must be >= 1. Got: ", ITER)
}

# File tag strings
frac_tag <- sprintf("%02d", as.integer(FRAC * 100))  # e.g. 0.70 -> "70"
iter_tag <- sprintf("%02d", ITER)                     # e.g. 1   -> "01"

# ============================================================================
# Load data (same as script 26)
# ============================================================================
cat("\n===== Loading data =====\n")

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

# Mega-analysis cohort selection: enforce config/human_datasets.yaml include_in_mega (same as 05/26).
ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Mega cohorts (yaml):", paste(mega_cohorts, collapse = ", "), "\n")
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

cat("Group x sex distribution:\n")
print(table(info$group_binary, info$inferred_sex))

# Detect sex labels
sex_levels <- levels(info$inferred_sex)
male_label   <- sex_levels[grepl("^M", sex_levels)][1]
female_label <- sex_levels[grepl("^F", sex_levels)][1]
cat("Using male label:", male_label, ", female label:", female_label, "\n")

# ============================================================================
# Stratified subsampling
# ============================================================================
cat("\n===== Stratified subsampling =====\n")

perform_subsampling <- function(info_df, frac, seed) {
  set.seed(seed)
  # Get unique strata: dataset x group_binary x inferred_sex
  info_df$row_idx <- seq_len(nrow(info_df))
  strata <- interaction(info_df$dataset, info_df$group_binary, info_df$inferred_sex, drop = TRUE)
  drawn_indices <- c()
  for (s in levels(strata)) {
    cell_idx <- info_df$row_idx[strata == s]
    n_cell <- length(cell_idx)
    n_draw <- max(floor(n_cell * frac), min(n_cell, 2))
    drawn_indices <- c(drawn_indices, sample(cell_idx, n_draw, replace = FALSE))
  }
  drawn_indices
}

validate_draw <- function(info_df, drawn_idx, male_lab, female_lab) {
  # For each sex stratum: check >= 2 datasets have both Control and Disease
  for (sex_lab in c(male_lab, female_lab)) {
    info_sex <- info_df[drawn_idx, , drop = FALSE]
    info_sex <- info_sex[info_sex$inferred_sex == sex_lab, , drop = FALSE]
    ds_counts <- table(info_sex$dataset, info_sex$group_binary)
    # Datasets with both groups
    both_groups <- rownames(ds_counts)[apply(ds_counts, 1, function(x) all(x > 0))]
    if (length(both_groups) < 2) {
      return(paste0("Sex stratum '", sex_lab, "' has only ",
                    length(both_groups), " dataset(s) with both groups"))
    }
  }
  return(NULL)  # NULL = pass
}

base_seed <- 42 + ITER * 1000 + as.integer(FRAC * 100)
drawn_idx <- NULL
max_attempts <- 10

for (attempt in seq_len(max_attempts)) {
  seed <- base_seed + (attempt - 1)
  candidate <- perform_subsampling(info, FRAC, seed)
  err <- validate_draw(info, candidate, male_label, female_label)
  if (is.null(err)) {
    drawn_idx <- candidate
    cat("Subsampling succeeded on attempt", attempt, "(seed =", seed, ")\n")
    break
  } else {
    cat("Attempt", attempt, "failed:", err, "\n")
  }
}

if (is.null(drawn_idx)) {
  err_file <- file.path(OUT_DIR, sprintf("ERROR_f%s_i%s.txt", frac_tag, iter_tag))
  writeLines(paste("All", max_attempts, "redraw attempts failed for FRAC =", FRAC,
                   "ITER =", ITER), err_file)
  stop("All ", max_attempts, " subsampling attempts failed. See ", err_file)
}

n_drawn <- length(drawn_idx)
cat("Samples drawn:", n_drawn, "of", nrow(info), "\n")

# Subset DGE and info to drawn samples
# CRITICAL: Do NOT re-run filterByExpr() or calcNormFactors() — preserves gene universe + normalization
info_sub <- info[drawn_idx, , drop = FALSE]
dge_sub  <- dge_mega[, drawn_idx]

cat("Subsampled group x sex distribution:\n")
print(table(info_sub$group_binary, info_sub$inferred_sex))

# ============================================================================
# Run dream for each sex stratum
# ============================================================================
cat("\n===== Running sex-stratified dream =====\n")

run_stratum_dream <- function(sex_label, dge_full, info_full, param) {
  cat("\n--- Stratum:", sex_label, "---\n")

  info_s <- info_full[info_full$inferred_sex == sex_label, , drop = FALSE]
  sample_ids_s <- rownames(info_s)
  dge_s <- dge_full[, colnames(dge_full) %in% sample_ids_s]

  cat("  Samples:", ncol(dge_s), "\n")
  cat("  Group distribution:\n")
  print(table(info_s$group_binary, info_s$dataset))

  # Require both groups present
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

  # Wrap in tryCatch — on error, log warning and return NULL

  res_dt <- tryCatch({
    v_s <- suppressWarnings(
      voomWithDreamWeights(dge_s, form_s, info_s, BPPARAM = param)
    )
    fit_s <- suppressWarnings(
      dream(v_s, form_s, info_s, BPPARAM = param)
    )

    res_s <- topTable(fit_s, coef = "group_binaryDisease", number = Inf, sort.by = "none")
    res_s$gene <- rownames(res_s)

    res_s_dt <- as.data.table(res_s)
    setnames(res_s_dt, "adj.P.Val", "padj")

    n_sig <- nrow(res_s_dt[padj < 0.05 & abs(logFC) > 0.25])
    cat("  DEGs padj < 0.05 & |LFC| > 0.25:", n_sig, "\n")

    res_s_dt
  }, error = function(e) {
    warning("  Dream failed for stratum ", sex_label, ": ", conditionMessage(e))
    return(NULL)
  })

  return(res_dt)
}

res_male   <- run_stratum_dream(male_label,   dge_sub, info_sub, param)
res_female <- run_stratum_dream(female_label, dge_sub, info_sub, param)

# ============================================================================
# Classify sex_class (same logic as script 26 lines 280-291)
# ============================================================================
cat("\n===== Sex-DEG classification =====\n")

PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.25

if (!is.null(res_male) && !is.null(res_female)) {
  cls_dt <- merge(
    res_male[,   .(gene, logFC_M = logFC, padj_M = padj, t_M = t)],
    res_female[, .(gene, logFC_F = logFC, padj_F = padj, t_F = t)],
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

  cat("Sex-class distribution:\n")
  print(cls_dt[, .N, by = sex_class][order(-N)])
} else {
  cat("WARNING: one or both strata failed — classification skipped.\n")
  cls_dt <- NULL
}

# ============================================================================
# Save outputs
# ============================================================================
cat("\n===== Saving outputs =====\n")

n_degs_male   <- 0
n_degs_female <- 0

if (!is.null(res_male)) {
  out_male <- res_male[, .(gene, logFC, padj, t)]
  fwrite(out_male, file.path(OUT_DIR, sprintf("sub_male_f%s_i%s.csv", frac_tag, iter_tag)))
  n_degs_male <- nrow(res_male[padj < PADJ_THRESH & abs(logFC) > LFC_THRESH])
  cat("Saved: sub_male_f", frac_tag, "_i", iter_tag, ".csv\n", sep = "")
} else {
  cat("WARNING: Male stratum failed — no male output.\n")
}

if (!is.null(res_female)) {
  out_female <- res_female[, .(gene, logFC, padj, t)]
  fwrite(out_female, file.path(OUT_DIR, sprintf("sub_female_f%s_i%s.csv", frac_tag, iter_tag)))
  n_degs_female <- nrow(res_female[padj < PADJ_THRESH & abs(logFC) > LFC_THRESH])
  cat("Saved: sub_female_f", frac_tag, "_i", iter_tag, ".csv\n", sep = "")
} else {
  cat("WARNING: Female stratum failed — no female output.\n")
}

if (!is.null(cls_dt)) {
  out_class <- cls_dt[, .(gene, sex_class, logFC_M, padj_M, logFC_F, padj_F)]
  fwrite(out_class, file.path(OUT_DIR, sprintf("sub_class_f%s_i%s.csv", frac_tag, iter_tag)))
  cat("Saved: sub_class_f", frac_tag, "_i", iter_tag, ".csv\n", sep = "")
} else {
  cat("WARNING: Classification failed — no class output.\n")
}

# ============================================================================
# Summary
# ============================================================================
elapsed <- (proc.time() - t0)["elapsed"]
cat("\n===== SUBSAMPLING ITERATION COMPLETE =====\n")
cat("  SUBSAMPLE_FRAC:", FRAC, "\n")
cat("  SUBSAMPLE_ITER:", ITER, "\n")
cat("  n_samples_drawn:", n_drawn, "\n")
cat("  n_degs_male (padj < 0.05 & |LFC| > 0.25):", n_degs_male, "\n")
cat("  n_degs_female (padj < 0.05 & |LFC| > 0.25):", n_degs_female, "\n")
cat("  Elapsed time:", round(elapsed / 60, 1), "minutes\n")
cat("Done:", as.character(Sys.time()), "\n")
