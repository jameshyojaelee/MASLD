#!/usr/bin/env Rscript
# 14.4c_sex_interaction_subsampling_iter.R
# ---------------------------------------------------------------------------
# ONE iteration of interaction-based sex subsampling sensitivity analysis.
#
# Tests stability of the interaction model (group_binary * inferred_sex)
# under stratified subsampling. This validates the v2 interaction-based
# sex classification from Script 26.
#
# Env vars:
#   SUBSAMPLE_FRAC  — float in [0.30, 1.0)  (default 0.70)
#   SUBSAMPLE_ITER  — integer >= 1           (default 1)
#
# For each iteration:
#   1. Stratified subsampling (dataset x group x sex) at SUBSAMPLE_FRAC
#   2. Validate drawn samples (>=2 datasets with both groups per sex)
#   3. Run interaction dream model (group_binary * inferred_sex + (1|dataset))
#   4. Run per-sex dream models for logFC characterization
#   5. Classify using interaction padj + per-sex logFCs
#   6. Save outputs
#
# Output: .../audit_sensitivity/sex_interaction_subsampling/iterations/
#   sub_interaction_f{FF}_i{II}.csv  — gene, interaction_logFC, interaction_padj
#   sub_class_v2_f{FF}_i{II}.csv    — gene, sex_class, logFC_M, logFC_F, etc.
# ---------------------------------------------------------------------------

t0 <- proc.time()

# ---------------------------------------------------------------------------
# Library loading + reformulas injection (MUST match script 26)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
  library(yaml)
  library(edgeR)
})

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
OUT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_interaction_subsampling/iterations"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parallel setup
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

if (is.na(FRAC) || FRAC < 0.3 || FRAC >= 1.0) {
  stop("SUBSAMPLE_FRAC must be in [0.3, 1.0). Got: ", FRAC)
}
if (is.na(ITER) || ITER < 1) {
  stop("SUBSAMPLE_ITER must be >= 1. Got: ", ITER)
}

frac_tag <- sprintf("%02d", as.integer(FRAC * 100))
iter_tag <- sprintf("%02d", ITER)

# ============================================================================
# Load data (same as script 26)
# ============================================================================
cat("\n===== Loading data =====\n")

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Mega cohorts (yaml):", paste(mega_cohorts, collapse = ", "), "
")
cat("Samples for mega-analysis:", ncol(dge_mega), "\n")

qc_report <- fread(file.path(INT, "qc/sample_qc_report.csv"))
sex_pass_ids <- qc_report[pass_sex == TRUE, sample_id]
sex_fail <- !colnames(dge_mega) %in% sex_pass_ids
if (any(sex_fail)) {
  cat("Removing", sum(sex_fail), "samples that failed sex check\n")
  dge_mega <- dge_mega[, !sex_fail]
}
cat("Samples after sex check filter:", ncol(dge_mega), "\n")

meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]

info <- data.frame(
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control", "Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(matched_sex),
  stringsAsFactors = FALSE
)
rownames(info) <- colnames(dge_mega)

na_sex <- is.na(info$inferred_sex)
if (any(na_sex)) {
  dge_mega <- dge_mega[, !na_sex]
  info     <- info[!na_sex, , drop = FALSE]
}

cat("Group x sex distribution:\n")
print(table(info$group_binary, info$inferred_sex))

sex_levels <- levels(info$inferred_sex)
male_label   <- sex_levels[grepl("^M", sex_levels)][1]
female_label <- sex_levels[grepl("^F", sex_levels)][1]

# ============================================================================
# Stratified subsampling (same logic as 14.4)
# ============================================================================
cat("\n===== Stratified subsampling =====\n")

perform_subsampling <- function(info_df, frac, seed) {
  set.seed(seed)
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
  for (sex_lab in c(male_lab, female_lab)) {
    info_sex <- info_df[drawn_idx, , drop = FALSE]
    info_sex <- info_sex[info_sex$inferred_sex == sex_lab, , drop = FALSE]
    ds_counts <- table(info_sex$dataset, info_sex$group_binary)
    both_groups <- rownames(ds_counts)[apply(ds_counts, 1, function(x) all(x > 0))]
    if (length(both_groups) < 2) {
      return(paste0("Sex stratum '", sex_lab, "' has only ",
                    length(both_groups), " dataset(s) with both groups"))
    }
  }
  return(NULL)
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
  writeLines(paste("All", max_attempts, "redraw attempts failed"), err_file)
  stop("All ", max_attempts, " subsampling attempts failed.")
}

n_drawn <- length(drawn_idx)
cat("Samples drawn:", n_drawn, "of", nrow(info), "\n")

info_sub <- info[drawn_idx, , drop = FALSE]
dge_sub  <- dge_mega[, drawn_idx]

cat("Subsampled group x sex distribution:\n")
print(table(info_sub$group_binary, info_sub$inferred_sex))

# ============================================================================
# Run interaction dream model on subsample
# ============================================================================
cat("\n===== Running interaction dream model =====\n")

form_int <- ~ group_binary * inferred_sex + (1|dataset)
cat("Formula:", deparse(form_int), "\n")

res_int_dt <- tryCatch({
  v_int <- suppressWarnings(
    voomWithDreamWeights(dge_sub, form_int, info_sub, BPPARAM = param)
  )
  fit_int <- suppressWarnings(
    dream(v_int, form_int, info_sub, BPPARAM = param)
  )

  all_coefs <- colnames(fit_int$coefficients)
  int_coef <- grep("group_binary.*inferred_sex|inferred_sex.*group_binary",
                   all_coefs, value = TRUE)
  if (length(int_coef) == 0) stop("No interaction coefficient found")
  if (length(int_coef) > 1) int_coef <- int_coef[1]
  cat("Using interaction coefficient:", int_coef, "\n")

  res_int <- topTable(fit_int, coef = int_coef, number = Inf, sort.by = "none")
  res_int$gene <- rownames(res_int)
  res_dt <- as.data.table(res_int)
  setnames(res_dt, "adj.P.Val", "padj")

  cat("Sex-differential (interaction padj < 0.05):", nrow(res_dt[padj < 0.05]), "\n")
  cat("Sex-differential (interaction padj < 0.1):", nrow(res_dt[padj < 0.1]), "\n")

  res_dt
}, error = function(e) {
  warning("Interaction dream failed: ", conditionMessage(e))
  return(NULL)
})

# ============================================================================
# Run per-sex dream models for logFC characterization
# ============================================================================
cat("\n===== Running per-sex dream models =====\n")

run_stratum_dream <- function(sex_label, dge_full, info_full, param) {
  cat("\n--- Stratum:", sex_label, "---\n")

  info_s <- info_full[info_full$inferred_sex == sex_label, , drop = FALSE]
  dge_s <- dge_full[, colnames(dge_full) %in% rownames(info_s)]

  cat("  Samples:", ncol(dge_s), "\n")

  n_ctrl <- sum(info_s$group_binary == "Control")
  n_dis  <- sum(info_s$group_binary == "Disease")
  if (n_ctrl == 0 || n_dis == 0) {
    warning("  Stratum ", sex_label, " missing controls or disease — skipping.")
    return(NULL)
  }

  ds_counts <- table(info_s$dataset, info_s$group_binary)
  keep_ds <- rownames(ds_counts)[apply(ds_counts, 1, function(x) all(x > 0))]
  if (length(keep_ds) < nrow(ds_counts)) {
    dropped <- setdiff(rownames(ds_counts), keep_ds)
    cat("  Dropping single-group datasets:", paste(dropped, collapse = ", "), "\n")
    info_s <- info_s[info_s$dataset %in% keep_ds, , drop = FALSE]
    info_s$dataset <- droplevels(info_s$dataset)
    dge_s  <- dge_s[, colnames(dge_s) %in% rownames(info_s)]
  }

  n_ds <- nlevels(info_s$dataset)
  if (n_ds < 2) {
    warning("  Only ", n_ds, " dataset(s) — random effect not estimable; skipping.")
    return(NULL)
  }

  form_s <- ~ group_binary + (1|dataset)

  tryCatch({
    v_s <- suppressWarnings(voomWithDreamWeights(dge_s, form_s, info_s, BPPARAM = param))
    fit_s <- suppressWarnings(dream(v_s, form_s, info_s, BPPARAM = param))
    res_s <- topTable(fit_s, coef = "group_binaryDisease", number = Inf, sort.by = "none")
    res_s$gene <- rownames(res_s)
    res_s_dt <- as.data.table(res_s)
    setnames(res_s_dt, "adj.P.Val", "padj")
    cat("  DEGs padj < 0.05 & |LFC| > 0.25:", nrow(res_s_dt[padj < 0.05 & abs(logFC) > 0.25]), "\n")
    res_s_dt
  }, error = function(e) {
    warning("  Dream failed for stratum ", sex_label, ": ", conditionMessage(e))
    return(NULL)
  })
}

res_male   <- run_stratum_dream(male_label,   dge_sub, info_sub, param)
res_female <- run_stratum_dream(female_label, dge_sub, info_sub, param)

# ============================================================================
# Interaction-based classification (same as script 26 v2 PART 3)
# ============================================================================
cat("\n===== Interaction-based classification =====\n")

INT_PADJ_THRESH   <- 0.05
DIVERGENT_MIN_LFC <- 0.1

if (!is.null(res_int_dt) && !is.null(res_male) && !is.null(res_female)) {
  cls_dt <- merge(
    res_male[,   .(gene, logFC_M = logFC, padj_M = padj)],
    res_female[, .(gene, logFC_F = logFC, padj_F = padj)],
    by = "gene", all = TRUE
  )

  cls_dt <- merge(cls_dt,
    res_int_dt[, .(gene, interaction_logFC = logFC, interaction_padj = padj)],
    by = "gene", all.x = TRUE
  )

  cls_dt[, sex_dimorphic := !is.na(interaction_padj) & interaction_padj < INT_PADJ_THRESH]

  cls_dt[, sex_class := fcase(
    !sex_dimorphic, "Concordant",
    sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
      sign(logFC_F) != sign(logFC_M) &
      abs(logFC_F) > DIVERGENT_MIN_LFC & abs(logFC_M) > DIVERGENT_MIN_LFC,
      "Divergent",
    sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
      abs(logFC_F) >= abs(logFC_M), "Female_biased",
    sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
      abs(logFC_M) > abs(logFC_F), "Male_biased",
    default = "Unclassified"
  )]

  cls_dt[, lfc_diff := abs(logFC_M - logFC_F)]

  cat("Sex-class distribution:\n")
  print(cls_dt[, .N, by = sex_class][order(-N)])
} else {
  cat("WARNING: one or more models failed — classification skipped.\n")
  cls_dt <- NULL
}

# ============================================================================
# Save outputs
# ============================================================================
cat("\n===== Saving outputs =====\n")

if (!is.null(res_int_dt)) {
  out_int <- res_int_dt[, .(gene, logFC, padj, t)]
  fwrite(out_int, file.path(OUT_DIR, sprintf("sub_interaction_f%s_i%s.csv", frac_tag, iter_tag)))
  cat("Saved: sub_interaction\n")
}

if (!is.null(cls_dt)) {
  out_cls <- cls_dt[, .(gene, sex_class, sex_dimorphic, interaction_logFC, interaction_padj,
                         logFC_M, padj_M, logFC_F, padj_F, lfc_diff)]
  fwrite(out_cls, file.path(OUT_DIR, sprintf("sub_class_v2_f%s_i%s.csv", frac_tag, iter_tag)))
  cat("Saved: sub_class_v2\n")
}

elapsed <- (proc.time() - t0)[3]
cat(sprintf("\nDone. Elapsed: %.1f min (frac=%s, iter=%s)\n",
            elapsed / 60, frac_tag, iter_tag))
