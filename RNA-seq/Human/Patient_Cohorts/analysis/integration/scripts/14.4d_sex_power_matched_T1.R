#!/usr/bin/env Rscript
# 14.4d_sex_power_matched_T1.R
# ---------------------------------------------------------------------------
# T1 (Team A audit): Power-matched F-control subsampling.
#
# Issue addressed (A1 audit Issue 6 + reviewer concern): the canonical Script 26
# fits ~ group_binary * inferred_sex on F controls n=101 vs M controls n=52.
# The 2:1 F:M asymmetry inflates female stratum power and drives the 2.9:1 F:M
# bias in sex_class assignment. T1 subsamples F controls to n=52 (match M count)
# WITHIN EACH COHORT (stratified by dataset × group_binary × inferred_sex,
# downsampling only the female-control cell), refits dream interaction +
# per-sex models, and reclassifies. Repeated B=50 times.
#
# Per-gene class-switching frequency reveals which sex_class calls are stable
# under power equalization (likely real biology) vs which collapse to
# Concordant/Male_biased under matched n (likely power artifacts).
#
# Env vars:
#   SUBSAMPLE_ITER  — integer >= 1 (default 1)  array task id
#   N_PER_CELL_M    — integer (default 52)      M-control target n (= observed M_ctrl)
#                                               also caps F-control draws
#
# Outputs (one per iteration):
#   .../audit_sensitivity/sex_power_matched/iterations/
#     pm_class_i{II}.csv     — gene, sex_class, logFC_M, logFC_F, interaction_padj
#     pm_metrics_i{II}.csv   — iteration metadata + n_dimorphic + n_per_class
# ---------------------------------------------------------------------------

t0 <- proc.time()

# Library loading + reformulas injection (MUST match script 26)
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
OUT_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_power_matched/iterations"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parallel + seed setup
# ---------------------------------------------------------------------------
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
seed_R_par <- as.integer(Sys.getenv("R_PARALLEL_SEED", "42"))
param <- if (ncpus > 1) MulticoreParam(workers = ncpus, RNGseed = seed_R_par) else SerialParam()

# ---------------------------------------------------------------------------
# Read env vars
# ---------------------------------------------------------------------------
ITER         <- as.integer(Sys.getenv("SUBSAMPLE_ITER", "1"))
N_PER_CELL_M <- as.integer(Sys.getenv("N_PER_CELL_M", "52"))   # target M control n

cat("SUBSAMPLE_ITER:", ITER, "\n")
cat("N_PER_CELL_M (target M control count):", N_PER_CELL_M, "\n")

if (is.na(ITER) || ITER < 1) stop("SUBSAMPLE_ITER must be >= 1. Got: ", ITER)

iter_tag <- sprintf("%02d", ITER)

# ============================================================================
# Load data (same as script 26 / 14.4c)
# ============================================================================
cat("\n===== Loading data =====\n")

dge <- readRDS(file.path(RDIR, "merged_dge.rds"))

ycfg <- yaml::read_yaml(file.path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design", "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
keep_samples <- dge$samples$dataset %in% mega_cohorts
dge_mega <- dge[, keep_samples]
cat("Mega cohorts (yaml):", paste(mega_cohorts, collapse = ", "), "\n")
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

cat("Full group x sex distribution:\n")
print(table(info$group_binary, info$inferred_sex))

sex_levels <- levels(info$inferred_sex)
male_label   <- sex_levels[grepl("^M", sex_levels)][1]
female_label <- sex_levels[grepl("^F", sex_levels)][1]

# ============================================================================
# Power-matched draw
#
# Strategy: keep ALL M samples, ALL F-disease samples; downsample F-control to
# match TOTAL M-control count (N_PER_CELL_M). The F-control downsample is
# stratified across cohorts proportionally to each cohort's F-control share, so
# we preserve the cohort composition of the F-control set rather than dropping
# whole cohorts.
# ============================================================================
draw_power_matched <- function(info_df, n_target_F_ctrl, seed) {
  set.seed(seed)
  fctrl_idx <- which(info_df$inferred_sex == female_label & info_df$group_binary == "Control")
  n_fctrl   <- length(fctrl_idx)
  if (n_target_F_ctrl >= n_fctrl) {
    return(seq_len(nrow(info_df)))  # nothing to drop
  }
  # stratify across datasets in proportion to cohort F-control share
  ds_vec <- as.character(info_df$dataset[fctrl_idx])
  ds_tab <- table(ds_vec)
  raw    <- as.numeric(ds_tab) * (n_target_F_ctrl / n_fctrl)
  n_each <- floor(raw)
  resid  <- raw - n_each
  short  <- n_target_F_ctrl - sum(n_each)
  if (short > 0) {
    # assign remainder by largest fractional residual, ties broken by RNG
    ord <- order(-resid, runif(length(resid)))
    n_each[ord[seq_len(short)]] <- n_each[ord[seq_len(short)]] + 1
  }
  names(n_each) <- names(ds_tab)
  draw_fctrl <- c()
  for (ds in names(n_each)) {
    pool <- fctrl_idx[ds_vec == ds]
    n_d  <- min(n_each[[ds]], length(pool))
    if (n_d > 0) draw_fctrl <- c(draw_fctrl, sample(pool, n_d, replace = FALSE))
  }
  keep <- c(setdiff(seq_len(nrow(info_df)), fctrl_idx), draw_fctrl)
  sort(keep)
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

base_seed <- 42 + ITER * 1000
drawn_idx <- NULL
max_attempts <- 10

for (attempt in seq_len(max_attempts)) {
  seed <- base_seed + (attempt - 1)
  candidate <- draw_power_matched(info, N_PER_CELL_M, seed)
  err <- validate_draw(info, candidate, male_label, female_label)
  if (is.null(err)) {
    drawn_idx <- candidate
    cat("Power-matched draw succeeded on attempt", attempt, "(seed =", seed, ")\n")
    break
  } else {
    cat("Attempt", attempt, "failed:", err, "\n")
  }
}

if (is.null(drawn_idx)) {
  err_file <- file.path(OUT_DIR, sprintf("ERROR_i%s.txt", iter_tag))
  writeLines(paste("All", max_attempts, "redraw attempts failed"), err_file)
  stop("All ", max_attempts, " power-matched draw attempts failed.")
}

info_sub <- info[drawn_idx, , drop = FALSE]
dge_sub  <- dge_mega[, drawn_idx]

cat("Power-matched group x sex distribution:\n")
pm_tab <- table(info_sub$group_binary, info_sub$inferred_sex)
print(pm_tab)

# ============================================================================
# Interaction dream model on power-matched subsample
# ============================================================================
cat("\n===== Running interaction dream model =====\n")
form_int <- ~ group_binary * inferred_sex + (1|dataset)

res_int_dt <- tryCatch({
  v_int <- suppressWarnings(voomWithDreamWeights(dge_sub, form_int, info_sub, BPPARAM = param))
  fit_int <- suppressWarnings(dream(v_int, form_int, info_sub, BPPARAM = param))
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
  res_dt
}, error = function(e) {
  warning("Interaction dream failed: ", conditionMessage(e))
  return(NULL)
})

# ============================================================================
# Per-sex dream models (same matched-dataset logic as Script 26)
# ============================================================================
run_stratum_dream <- function(sex_label, dge_full, info_full, param) {
  cat("\n--- Stratum:", sex_label, "---\n")
  info_s <- info_full[info_full$inferred_sex == sex_label, , drop = FALSE]
  dge_s <- dge_full[, colnames(dge_full) %in% rownames(info_s)]
  n_ctrl <- sum(info_s$group_binary == "Control")
  n_dis  <- sum(info_s$group_binary == "Disease")
  if (n_ctrl == 0 || n_dis == 0) {
    warning("  Stratum ", sex_label, " missing controls or disease — skipping.")
    return(NULL)
  }
  ds_counts <- table(info_s$dataset, info_s$group_binary)
  keep_ds <- rownames(ds_counts)[apply(ds_counts, 1, function(x) all(x > 0))]
  if (length(keep_ds) < nrow(ds_counts)) {
    info_s <- info_s[info_s$dataset %in% keep_ds, , drop = FALSE]
    info_s$dataset <- droplevels(info_s$dataset)
    dge_s  <- dge_s[, colnames(dge_s) %in% rownames(info_s)]
  }
  if (nlevels(info_s$dataset) < 2) {
    warning("  <2 datasets in stratum ", sex_label, " — skipping.")
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
    res_s_dt
  }, error = function(e) {
    warning("  Dream failed for stratum ", sex_label, ": ", conditionMessage(e))
    return(NULL)
  })
}

res_male   <- run_stratum_dream(male_label,   dge_sub, info_sub, param)
res_female <- run_stratum_dream(female_label, dge_sub, info_sub, param)

# ============================================================================
# Interaction-based classification (same v2 logic as Script 26)
# ============================================================================
cat("\n===== Classification =====\n")

INT_PADJ_THRESH   <- 0.05
DIVERGENT_MIN_LFC <- 0.1

cls_dt <- NULL
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
  cat("Power-matched sex_class distribution:\n")
  print(cls_dt[, .N, by = sex_class][order(-N)])
}

# ============================================================================
# Save outputs
# ============================================================================
cat("\n===== Saving =====\n")

if (!is.null(cls_dt)) {
  out_cls <- cls_dt[, .(gene, sex_class, sex_dimorphic, interaction_logFC, interaction_padj,
                        logFC_M, padj_M, logFC_F, padj_F, lfc_diff)]
  fwrite(out_cls, file.path(OUT_DIR, sprintf("pm_class_i%s.csv", iter_tag)))
  cat("Saved: pm_class_i", iter_tag, ".csv\n", sep = "")
}

# per-iteration metrics row
metrics <- data.table(
  iter        = ITER,
  n_total     = nrow(info_sub),
  n_F_ctrl    = sum(info_sub$inferred_sex == female_label & info_sub$group_binary == "Control"),
  n_M_ctrl    = sum(info_sub$inferred_sex == male_label   & info_sub$group_binary == "Control"),
  n_F_dis     = sum(info_sub$inferred_sex == female_label & info_sub$group_binary == "Disease"),
  n_M_dis     = sum(info_sub$inferred_sex == male_label   & info_sub$group_binary == "Disease"),
  n_dimorphic = if (!is.null(cls_dt)) sum(cls_dt$sex_dimorphic, na.rm = TRUE) else NA_integer_,
  n_female_biased = if (!is.null(cls_dt)) cls_dt[sex_class == "Female_biased", .N] else NA_integer_,
  n_male_biased   = if (!is.null(cls_dt)) cls_dt[sex_class == "Male_biased", .N]   else NA_integer_,
  n_divergent     = if (!is.null(cls_dt)) cls_dt[sex_class == "Divergent", .N]     else NA_integer_,
  n_concordant    = if (!is.null(cls_dt)) cls_dt[sex_class == "Concordant", .N]    else NA_integer_,
  elapsed_min     = as.numeric((proc.time() - t0)[3]) / 60
)
fwrite(metrics, file.path(OUT_DIR, sprintf("pm_metrics_i%s.csv", iter_tag)))
cat("Saved: pm_metrics_i", iter_tag, ".csv\n", sep = "")

cat(sprintf("\nDone. Elapsed: %.1f min (iter=%s)\n",
            (proc.time() - t0)[3] / 60, iter_tag))
