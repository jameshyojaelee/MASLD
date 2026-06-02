#!/usr/bin/env Rscript
# 60_one_vs_rest_stage_dream.R
# One-vs-rest dream mega-analyses for stage-specific gene signatures
#
# For each NAS level (0-7+) and fibrosis stage (F0-F4), fits a dream model
# comparing that stage vs all other stages to identify stage-UNIQUE genes.
#
# Outputs:
#   - one_vs_rest_nas_dream.csv: One-vs-rest DEG results for each NAS level
#   - one_vs_rest_fibrosis_dream.csv: One-vs-rest DEG results for each fibrosis stage
#   - stage_unique_signatures.csv: Genes uniquely significant in one stage only
#   - stage_specificity_scores.csv: Per-gene specificity metrics across stages
#
# Usage: Rscript 60_one_vs_rest_stage_dream.R
# SLURM: bigmem, 16 CPUs, 200GB RAM, 48h

suppressPackageStartupMessages({
  library(reformulas)
  library(lme4)
  library(data.table)
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

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results")
OUTDIR <- file.path(RDIR, "staging_classifier")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== 60: One-vs-Rest Stage Dream Mega-Analysis ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# --- Setup parallel backend ---
ncpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using", ncpus, "CPU cores\n")
param <- if (ncpus > 1) MulticoreParam(ncpus) else SerialParam()

# --- Load data ---
cat("Loading merged counts...\n")
dge <- readRDS(file.path(RDIR, "integration/merged_dge.rds"))

cat("Loading matched metadata...\n")
meta_matched_file <- file.path(RDIR, "integration/meta_matched.rds")
if (!file.exists(meta_matched_file)) {
  stop("meta_matched.rds not found at: ", meta_matched_file,
       "\nRun integration scripts (03_integrate_counts.R) first.")
}
meta_matched <- readRDS(meta_matched_file)

meta_unified <- fread(file.path(INT, "metadata/unified_metadata.csv"))
cat("Unified metadata:", nrow(meta_unified), "samples\n")

qc <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_samples <- qc[pass_technical == TRUE, sample_id]
cat("QC-passing samples:", length(pass_samples), "\n\n")

# ============================================================
# Helper: Run one-vs-rest dream for a single stage level
# ============================================================
run_ovr_dream <- function(dge_subset, meta_subset, meta_matched_full,
                          stage_col, stage_value, stage_label, param) {
  # Create binary indicator: this stage vs all others
  meta_subset[, ovr := fifelse(get(stage_col) == stage_value, "target", "rest")]

  # Check sample sizes
  n_target <- sum(meta_subset$ovr == "target")
  n_rest   <- sum(meta_subset$ovr == "rest")
  cat("  ", stage_label, ": n_target =", n_target, ", n_rest =", n_rest, "\n")

  if (n_target < 5) {
    cat("  SKIPPING: too few target samples (<5)\n")
    return(NULL)
  }

  # Subset DGE to these samples
  sample_ids <- meta_subset$sample_id
  idx <- colnames(dge_subset) %in% sample_ids
  dge_ovr <- dge_subset[, idx]

  # Reorder metadata to match DGE columns
  meta_ovr <- meta_subset[match(colnames(dge_ovr), meta_subset$sample_id)]

  # Get inferred sex from meta_matched
  sex_vals <- meta_matched_full$inferred_sex[match(colnames(dge_ovr), meta_matched_full$sample_id)]

  # Build info dataframe
  info <- data.frame(
    ovr          = factor(meta_ovr$ovr, levels = c("rest", "target")),
    dataset      = factor(meta_ovr$dataset),
    inferred_sex = factor(sex_vals),
    row.names    = colnames(dge_ovr),
    stringsAsFactors = FALSE
  )

  # Remove samples with missing sex
  valid <- !is.na(info$inferred_sex)
  if (sum(valid) < nrow(info)) {
    cat("  Removing", sum(!valid), "samples with missing sex\n")
    dge_ovr <- dge_ovr[, valid]
    info <- info[valid, , drop = FALSE]
  }

  # Drop unused dataset levels (some datasets may have been removed)
  info$dataset <- droplevels(info$dataset)

  # Check if we still have enough data
  if (sum(info$ovr == "target") < 5) {
    cat("  SKIPPING after sex filter: too few target samples\n")
    return(NULL)
  }

  # Need >=2 datasets for random effect
  n_datasets <- nlevels(info$dataset)
  if (n_datasets < 2) {
    cat("  SKIPPING: only", n_datasets, "dataset(s) — need >=2 for random effect\n")
    return(NULL)
  }

  # Filter by expression
  keep <- filterByExpr(dge_ovr, group = info$ovr)
  dge_ovr <- dge_ovr[keep, , keep.lib.sizes = FALSE]
  dge_ovr <- calcNormFactors(dge_ovr, method = "TMM")

  # Dream formula: one-vs-rest indicator + sex + random intercept per dataset
  form <- ~ ovr + inferred_sex + (1 | dataset)

  # Fit dream
  v <- voomWithDreamWeights(dge_ovr, form, info, BPPARAM = param)
  fit <- dream(v, form, info, BPPARAM = param)

  # Extract results for the "target" coefficient
  coef_name <- "ovrtarget"
  if (!(coef_name %in% colnames(coef(fit)))) {
    cat("  WARNING: coefficient", coef_name, "not found\n")
    return(NULL)
  }

  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$stage_value <- stage_value
  tt$stage_label <- stage_label

  n_sig <- sum(tt$adj.P.Val < 0.1, na.rm = TRUE)
  n_sig_strict <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, na.rm = TRUE)
  cat("  DEGs: padj<0.1 =", n_sig, ", padj<0.05 & |LFC|>0.5 =", n_sig_strict, "\n")

  return(as.data.table(tt))
}

# ============================================================
# ANALYSIS 1: One-vs-Rest NAS Score
# ============================================================
cat("=== ANALYSIS 1: One-vs-Rest NAS Score ===\n\n")

NAS_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066")

# Check if OVR results already exist (skip expensive dream re-runs)
nas_ovr_file <- file.path(OUTDIR, "one_vs_rest_nas_dream.csv")
if (file.exists(nas_ovr_file) && file.size(nas_ovr_file) > 1000) {
  cat("NAS OVR results already exist, loading from cache:", nas_ovr_file, "\n")
  nas_ovr_results <- fread(nas_ovr_file)
  cat("Loaded:", nrow(nas_ovr_results), "rows\n\n")
} else {
  meta_nas <- meta_unified[
    dataset %in% NAS_DATASETS &
    !is.na(nas_score) &
    sample_id %in% pass_samples
  ]

  # Collapse NAS >= 7 into "7+"
  meta_nas[, nas_group := fifelse(nas_score >= 7, 7L, as.integer(nas_score))]

  cat("NAS-annotated QC-passing samples:", nrow(meta_nas), "\n")
  nas_dist <- meta_nas[, .N, by = nas_group][order(nas_group)]
  cat("\nSample distribution by NAS group:\n")
  print(nas_dist)
  cat("\n")

  # Subset DGE to NAS samples
  nas_samples <- meta_nas$sample_id
  dge_nas <- dge[, colnames(dge) %in% nas_samples]
  cat("DGE subset:", ncol(dge_nas), "samples\n\n")

  # Run one-vs-rest for each NAS level
  nas_levels <- sort(unique(meta_nas$nas_group))
  all_nas_ovr <- list()

  for (lvl in nas_levels) {
    label <- paste0("NAS", lvl, "_vs_rest")
    cat("Running:", label, "\n")

    # Fresh copy of metadata for each iteration
    meta_iter <- copy(meta_nas)

    result <- tryCatch(
      run_ovr_dream(dge_nas, meta_iter, meta_matched,
                    "nas_group", lvl, label, param),
      error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n")
        return(NULL)
      }
    )

    if (!is.null(result)) {
      all_nas_ovr[[as.character(lvl)]] <- result
    }
    cat("\n")
  }

  nas_ovr_results <- rbindlist(all_nas_ovr, fill = TRUE)
  setnames(nas_ovr_results, "adj.P.Val", "padj", skip_absent = TRUE)
  fwrite(nas_ovr_results, nas_ovr_file)
  cat("NAS one-vs-rest results saved:", nrow(nas_ovr_results), "rows across",
      length(all_nas_ovr), "contrasts\n\n")
}

# ============================================================
# ANALYSIS 2: One-vs-Rest Fibrosis Stage
# ============================================================
cat("=== ANALYSIS 2: One-vs-Rest Fibrosis Stage ===\n\n")

FIB_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478",
                   "GSE193066", "GSE240729")

fib_ovr_file <- file.path(OUTDIR, "one_vs_rest_fibrosis_dream.csv")
if (file.exists(fib_ovr_file) && file.size(fib_ovr_file) > 1000) {
  cat("Fibrosis OVR results already exist, loading from cache:", fib_ovr_file, "\n")
  fib_ovr_results <- fread(fib_ovr_file)
  cat("Loaded:", nrow(fib_ovr_results), "rows\n\n")
} else {
  meta_fib <- meta_unified[
    dataset %in% FIB_DATASETS &
    !is.na(fibrosis_stage) &
    sample_id %in% pass_samples
  ]
  meta_fib[, fib_stage := as.integer(fibrosis_stage)]
  meta_fib <- meta_fib[fib_stage %in% 0:4]

  cat("Fibrosis-annotated QC-passing samples:", nrow(meta_fib), "\n")
  fib_dist <- meta_fib[, .N, by = fib_stage][order(fib_stage)]
  cat("\nSample distribution by fibrosis stage:\n")
  print(fib_dist)
  cat("\n")

  # Subset DGE
  fib_samples <- meta_fib$sample_id
  dge_fib <- dge[, colnames(dge) %in% fib_samples]
  cat("DGE subset:", ncol(dge_fib), "samples\n\n")

  # Run one-vs-rest for each fibrosis stage
  fib_levels <- sort(unique(meta_fib$fib_stage))
  all_fib_ovr <- list()

  for (lvl in fib_levels) {
    label <- paste0("F", lvl, "_vs_rest")
    cat("Running:", label, "\n")

    meta_iter <- copy(meta_fib)

    result <- tryCatch(
      run_ovr_dream(dge_fib, meta_iter, meta_matched,
                    "fib_stage", lvl, label, param),
      error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n")
        return(NULL)
      }
    )

    if (!is.null(result)) {
      all_fib_ovr[[as.character(lvl)]] <- result
    }
    cat("\n")
  }

  fib_ovr_results <- rbindlist(all_fib_ovr, fill = TRUE)
  setnames(fib_ovr_results, "adj.P.Val", "padj", skip_absent = TRUE)
  fwrite(fib_ovr_results, fib_ovr_file)
  cat("Fibrosis one-vs-rest results saved:", nrow(fib_ovr_results), "rows across",
      length(all_fib_ovr), "contrasts\n\n")
}

# ============================================================
# ANALYSIS 3: Stage Uniqueness & Specificity Scores
# ============================================================
cat("=== ANALYSIS 3: Stage Uniqueness & Specificity Scores ===\n\n")

compute_specificity <- function(ovr_results, stage_col, padj_unique = 0.1, padj_absent = 0.3) {
  # For each gene, determine which stages it is significant in (one-vs-rest)
  # A gene is "stage-unique" if significant in exactly one stage's OVR contrast
  # AND not significant (padj > padj_absent) in all other stages

  stages <- unique(ovr_results$stage_value)
  genes <- unique(ovr_results$gene)

  # Build significance matrix: gene x stage
  # Use safe aggregation functions that handle empty/NA vectors
  safe_min <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    min(x)
  }
  safe_lfc <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    x[which.min(abs(x))]
  }
  sig_matrix <- dcast(ovr_results, gene ~ stage_label, value.var = "padj", fun.aggregate = safe_min)
  lfc_matrix <- dcast(ovr_results, gene ~ stage_label, value.var = "logFC", fun.aggregate = safe_lfc)

  # For each gene, compute:
  # 1. n_stages_sig: number of stages where padj < padj_unique
  # 2. specificity_score: max(-log10(padj)) - mean(-log10(padj)) across stages
  # 3. tau (adapted from cell-type specificity): 1 - mean(x_i / max(x_i))

  stage_cols <- setdiff(names(sig_matrix), "gene")

  specificity <- data.table(gene = sig_matrix$gene)

  # Count stages where significant
  padj_mat <- as.matrix(sig_matrix[, ..stage_cols])
  lfc_mat <- as.matrix(lfc_matrix[, ..stage_cols])

  specificity[, n_stages_sig := rowSums(padj_mat < padj_unique, na.rm = TRUE)]
  specificity[, n_stages_absent := rowSums(padj_mat > padj_absent, na.rm = TRUE)]

  # Identify unique genes: significant in exactly 1, absent in all others
  n_total_stages <- length(stage_cols)
  specificity[, is_stage_unique := (n_stages_sig == 1) & (n_stages_absent == n_total_stages - 1)]

  # Which stage is it unique to?
  specificity[, unique_stage := NA_character_]
  for (i in seq_along(stage_cols)) {
    idx <- specificity$is_stage_unique & padj_mat[, i] < padj_unique
    idx[is.na(idx)] <- FALSE
    specificity[idx, unique_stage := stage_cols[i]]
  }

  # Tau specificity index (adapted from Yanai et al. 2005)
  # Using -log10(padj) as the "expression" analog
  nlog_padj <- -log10(pmax(padj_mat, 1e-300))
  nlog_padj[is.na(nlog_padj)] <- 0
  max_nlog <- apply(nlog_padj, 1, max)
  max_nlog[max_nlog == 0] <- 1  # avoid division by zero

  tau_parts <- sweep(nlog_padj, 1, max_nlog, "/")
  specificity[, tau := 1 - rowMeans(tau_parts, na.rm = TRUE) * n_total_stages / (n_total_stages - 1)]
  specificity[tau < 0, tau := 0]

  # Most specific stage (highest -log10 padj)
  specificity[, most_specific_stage := stage_cols[apply(nlog_padj, 1, which.max)]]
  specificity[, most_specific_padj := apply(padj_mat, 1, min, na.rm = TRUE)]

  # Max absolute LFC across stages
  specificity[, max_abs_lfc := apply(abs(lfc_mat), 1, max, na.rm = TRUE)]

  # Direction in most specific stage
  for (i in seq_along(stage_cols)) {
    idx <- specificity$most_specific_stage == stage_cols[i]
    idx[is.na(idx)] <- FALSE
    specificity[idx, direction := fifelse(lfc_mat[idx, i] > 0, "up", "down")]
  }

  return(specificity)
}

# Compute for NAS
if (nrow(nas_ovr_results) > 0) {
  cat("Computing NAS specificity scores...\n")
  nas_specificity <- compute_specificity(nas_ovr_results, "nas_group")
  nas_specificity[, staging_axis := "NAS"]

  n_unique <- sum(nas_specificity$is_stage_unique)
  cat("NAS stage-unique genes:", n_unique, "\n")
  cat("Unique genes per NAS level:\n")
  print(nas_specificity[is_stage_unique == TRUE, .N, by = unique_stage])
  cat("\n")
}

# Compute for fibrosis
if (nrow(fib_ovr_results) > 0) {
  cat("Computing fibrosis specificity scores...\n")
  fib_specificity <- compute_specificity(fib_ovr_results, "fib_stage")
  fib_specificity[, staging_axis := "Fibrosis"]

  n_unique <- sum(fib_specificity$is_stage_unique)
  cat("Fibrosis stage-unique genes:", n_unique, "\n")
  cat("Unique genes per fibrosis stage:\n")
  print(fib_specificity[is_stage_unique == TRUE, .N, by = unique_stage])
  cat("\n")
}

# Combine specificity results
all_specificity <- rbindlist(list(nas_specificity, fib_specificity), fill = TRUE)
fwrite(all_specificity, file.path(OUTDIR, "stage_specificity_scores.csv"))
cat("Specificity scores saved:", nrow(all_specificity), "gene-axis pairs\n")

# Extract stage-unique signature genes
unique_sigs <- all_specificity[is_stage_unique == TRUE]
fwrite(unique_sigs, file.path(OUTDIR, "stage_unique_signatures.csv"))
cat("Stage-unique signatures saved:", nrow(unique_sigs), "genes\n")

# ============================================================
# Summary statistics
# ============================================================
cat("\n=== Summary ===\n")
cat("NAS contrasts:", length(unique(nas_ovr_results$stage_label)), "\n")
cat("Fibrosis contrasts:", length(unique(fib_ovr_results$stage_label)), "\n")
cat("Total NAS OVR results:", nrow(nas_ovr_results), "\n")
cat("Total Fibrosis OVR results:", nrow(fib_ovr_results), "\n")
cat("NAS stage-unique genes:", sum(nas_specificity$is_stage_unique), "\n")
cat("Fibrosis stage-unique genes:", sum(fib_specificity$is_stage_unique), "\n")
cat("Mean NAS tau:", round(mean(nas_specificity$tau, na.rm = TRUE), 3), "\n")
cat("Mean fibrosis tau:", round(mean(fib_specificity$tau, na.rm = TRUE), 3), "\n")

cat("\n=== 60_one_vs_rest_stage_dream.R completed:", as.character(Sys.time()), "===\n")
