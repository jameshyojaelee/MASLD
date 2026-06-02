#!/usr/bin/env Rscript
# run_bayesprism_celltype_de.R
# ---------------------------------------------------------------------------
# Per-cell-type differential expression on the four CONSECUTIVE F-transitions
# (F1_vs_F0, F2_vs_F1, F3_vs_F2, F4_vs_F3) using BayesPrism per-cell-type
# posterior expression matrices (1,444 bulk donors x 13 cell types).
#
# This is the bulk-side counterpart to scRNA pseudobulk DE: it tells us which
# cell types own each bulk F-transition without requiring F-stage annotation
# on scRNA donors.
#
# Inputs
#   bayesprism_<celltype>.csv.gz    -- samples x genes posterior expression
#                                       (rows = samples, cols = HGNC genes)
#   unified_metadata.csv             -- sample_id, dataset, fibrosis_stage,
#                                       inferred_sex (via per_sample_dysregulation.csv)
#   sample_qc_report.csv             -- pass_technical filter
#
# Method
#   BayesPrism posteriors are real-valued and on a count-like scale (already
#   apportioned by the cell-type proportion). They are NOT raw counts and
#   floor at 0 with a heavy tail, so we:
#     1. log2(x + 1) transform
#     2. run limma with a fixed-effect design ~ stage_factor + inferred_sex + dataset
#     3. test each consecutive contrast via topTable on the relevant coefficient
#
#   Genes are kept if at least 50% of samples in the contrasted strata have
#   nonzero posterior expression (avoids fitting on rare cell types where
#   most posterior mass is zero).
#
# Outputs (Analysis/SingleCell/results_gpu_v2/pseudobulk_de/)
#   <celltype>_<contrast>_bayesprism_de.csv -- gene, logFC, padj
#   bayesprism_celltype_de_log.txt          -- run-level summary
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

EXPR_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results",
  "progression/cibersortx_celltype_expression")
META_PATH <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
QC_PATH <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")
DYSREG_PATH <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results",
  "disease_signatures/per_sample_dysregulation.csv")

OUT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_TYPES <- c("Hepatocyte", "Macrophage", "Stellate", "Endothelial",
                "Cholangiocyte", "T_cell", "B_cell", "Plasma_cell",
                "NK_cell", "Monocyte", "DC", "Neutrophil", "Other_immune")

CONTRASTS <- list(
  F1_vs_F0 = c(0L, 1L),
  F2_vs_F1 = c(1L, 2L),
  F3_vs_F2 = c(2L, 3L),
  F4_vs_F3 = c(3L, 4L)
)

cat("=== run_bayesprism_celltype_de.R ===\n")
cat(sprintf("Started: %s\n", Sys.time()))

# ---------------------------------------------------------------------------
# 1. Build sample annotation: sample_id, dataset, fibrosis_stage, inferred_sex
# ---------------------------------------------------------------------------
meta <- fread(META_PATH)
qc   <- fread(QC_PATH)
dys  <- fread(DYSREG_PATH, select = c("sample_id", "inferred_sex"))

meta <- merge(meta, qc[, .(sample_id, pass_technical)], by = "sample_id", all.x = TRUE)
meta <- merge(meta, dys, by = "sample_id", all.x = TRUE)

meta[, fib_stage := suppressWarnings(as.integer(fibrosis_stage))]
meta[, inferred_sex := ifelse(is.na(inferred_sex) | inferred_sex == "",
                              ifelse(is.na(sex) | sex == "", NA_character_, sex),
                              inferred_sex)]

cat(sprintf("Total samples in metadata: %d\n", nrow(meta)));
cat(sprintf("pass_technical = TRUE: %d\n", sum(meta$pass_technical, na.rm = TRUE)))
cat(sprintf("With fib_stage in 0..4: %d\n", sum(!is.na(meta$fib_stage) & meta$fib_stage %in% 0:4)))
cat(sprintf("With inferred_sex non-NA: %d\n", sum(!is.na(meta$inferred_sex))))

usable <- meta[
  isTRUE(pass_technical) | pass_technical == TRUE
][!is.na(fib_stage) & fib_stage %in% 0:4 & !is.na(inferred_sex) & inferred_sex != ""]
cat(sprintf("Usable samples (pass_technical & fib_stage & inferred_sex): %d\n",
            nrow(usable)))
cat("Stage distribution among usable samples:\n")
print(usable[, .N, by = fib_stage][order(fib_stage)])

# ---------------------------------------------------------------------------
# 2. Per-cell-type, per-contrast DE
# ---------------------------------------------------------------------------
log_lines <- character()
add_log <- function(...) {
  msg <- paste0(...)
  message(msg)
  log_lines <<- c(log_lines, msg)
}

run_one_de <- function(expr_log2, sample_meta, stage_lo, stage_hi, contrast_name,
                       cell_type) {
  # Subset to samples in the two adjacent stages
  keep_samples <- sample_meta[fib_stage %in% c(stage_lo, stage_hi)]
  if (nrow(keep_samples) < 6) {
    add_log(sprintf("    %s/%s: too few samples (%d), skipping",
                    cell_type, contrast_name, nrow(keep_samples)))
    return(NULL)
  }
  expr_sub <- expr_log2[, keep_samples$sample_id, drop = FALSE]

  # Drop genes that are essentially zero in this contrast (>=50% nonzero
  # in either arm to keep)
  nonzero_lo <- rowSums(expr_sub[, keep_samples[fib_stage == stage_lo]$sample_id,
                                 drop = FALSE] > 0)
  nonzero_hi <- rowSums(expr_sub[, keep_samples[fib_stage == stage_hi]$sample_id,
                                 drop = FALSE] > 0)
  n_lo <- sum(keep_samples$fib_stage == stage_lo)
  n_hi <- sum(keep_samples$fib_stage == stage_hi)
  keep_genes <- (nonzero_lo / max(n_lo, 1) >= 0.5) | (nonzero_hi / max(n_hi, 1) >= 0.5)
  expr_sub <- expr_sub[keep_genes, , drop = FALSE]

  if (nrow(expr_sub) < 100) {
    add_log(sprintf("    %s/%s: too few genes after filter (%d), skipping",
                    cell_type, contrast_name, nrow(expr_sub)))
    return(NULL)
  }

  # Design: stage as factor (hi vs lo), inferred_sex, dataset
  keep_samples[, stage_f := factor(ifelse(fib_stage == stage_hi, "hi", "lo"),
                                   levels = c("lo", "hi"))]
  keep_samples[, inferred_sex := factor(inferred_sex)]
  keep_samples[, dataset := factor(dataset)]

  # Drop singletons in dataset (limma requires >=2 levels with >=1 obs each;
  # but we also want each dataset to have both stages or it absorbs the contrast)
  ds_counts <- keep_samples[, .(n_lo = sum(stage_f == "lo"),
                                n_hi = sum(stage_f == "hi")), by = dataset]
  bad_ds <- ds_counts[n_lo == 0 | n_hi == 0]$dataset
  if (length(bad_ds) > 0) {
    keep_samples <- keep_samples[!dataset %in% bad_ds]
    expr_sub <- expr_sub[, keep_samples$sample_id, drop = FALSE]
    keep_samples[, dataset := droplevels(dataset)]
  }
  if (nrow(keep_samples) < 6 || sum(keep_samples$stage_f == "lo") < 3 ||
      sum(keep_samples$stage_f == "hi") < 3) {
    add_log(sprintf("    %s/%s: too few samples after dataset filter (lo=%d, hi=%d), skipping",
                    cell_type, contrast_name,
                    sum(keep_samples$stage_f == "lo"),
                    sum(keep_samples$stage_f == "hi")))
    return(NULL)
  }

  # Drop singleton sex level if present
  if (length(unique(keep_samples$inferred_sex)) < 2) {
    design_form <- ~ stage_f + dataset
  } else {
    keep_samples[, inferred_sex := droplevels(inferred_sex)]
    design_form <- ~ stage_f + inferred_sex + dataset
  }
  design <- model.matrix(design_form, data = keep_samples)

  # Drop rank-deficient columns (e.g. dataset levels with all one stage after filtering)
  qr_d <- qr(design)
  if (qr_d$rank < ncol(design)) {
    keep_cols <- qr_d$pivot[seq_len(qr_d$rank)]
    design <- design[, keep_cols, drop = FALSE]
  }

  if (!"stage_fhi" %in% colnames(design)) {
    add_log(sprintf("    %s/%s: stage_fhi coefficient missing after rank filter, skipping",
                    cell_type, contrast_name))
    return(NULL)
  }

  fit <- tryCatch(
    {
      f <- lmFit(expr_sub, design)
      eBayes(f, trend = TRUE, robust = TRUE)
    },
    error = function(e) {
      add_log(sprintf("    %s/%s: lmFit/eBayes ERROR: %s",
                      cell_type, contrast_name, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  tt <- topTable(fit, coef = "stage_fhi", number = Inf, sort.by = "none")
  out_df <- data.table(
    gene  = rownames(tt),
    logFC = tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    P.Value = tt$P.Value,
    padj = tt$adj.P.Val
  )

  add_log(sprintf("    %s/%s: n_lo=%d n_hi=%d genes=%d sig(padj<0.05,|lfc|>0.5)=%d",
                  cell_type, contrast_name,
                  sum(keep_samples$stage_f == "lo"),
                  sum(keep_samples$stage_f == "hi"),
                  nrow(out_df),
                  sum(out_df$padj < 0.05 & abs(out_df$logFC) > 0.5, na.rm = TRUE)))
  out_df
}

results_summary <- list()
for (ct in CELL_TYPES) {
  expr_path <- file.path(EXPR_DIR, sprintf("bayesprism_%s.csv.gz", ct))
  if (!file.exists(expr_path)) {
    add_log(sprintf("MISSING: %s", expr_path))
    next
  }
  add_log(sprintf("\n[%s] loading %s", ct, basename(expr_path)))
  expr_dt <- fread(expr_path)
  # samples x genes; first column = sample_id
  sample_ids <- expr_dt[[1]]
  expr_mat <- as.matrix(expr_dt[, -1])
  rownames(expr_mat) <- sample_ids
  expr_mat <- t(expr_mat) # genes x samples

  # log2(x + 1)
  expr_log2 <- log2(expr_mat + 1)

  # Restrict to usable samples present in this matrix
  ct_samples <- usable[sample_id %in% colnames(expr_log2)]
  if (nrow(ct_samples) == 0) {
    add_log(sprintf("  %s: no usable samples in matrix, skipping", ct))
    next
  }
  expr_log2 <- expr_log2[, ct_samples$sample_id, drop = FALSE]
  add_log(sprintf("  %s: %d genes x %d samples (after intersect with usable)",
                  ct, nrow(expr_log2), ncol(expr_log2)))

  for (contrast_name in names(CONTRASTS)) {
    pair <- CONTRASTS[[contrast_name]]
    res <- run_one_de(expr_log2, copy(ct_samples), pair[1], pair[2],
                      contrast_name, ct)
    if (is.null(res)) next
    out_path <- file.path(OUT_DIR,
                          sprintf("%s_%s_bayesprism_de.csv", ct, contrast_name))
    fwrite(res[, .(gene, logFC, AveExpr, t, P.Value, padj)], out_path)
    sig_n <- sum(res$padj < 0.05 & abs(res$logFC) > 0.5, na.rm = TRUE)
    results_summary[[length(results_summary) + 1L]] <- data.table(
      cell_type = ct,
      contrast  = contrast_name,
      n_genes   = nrow(res),
      n_sig     = sig_n,
      out_path  = out_path
    )
  }
}

# ---------------------------------------------------------------------------
# 3. Persist run log + per-(cell_type, contrast) significant counts
# ---------------------------------------------------------------------------
log_path <- file.path(OUT_DIR, "bayesprism_celltype_de_log.txt")
writeLines(log_lines, log_path)
cat(sprintf("\nWrote log to: %s\n", log_path))

if (length(results_summary) > 0) {
  summary_dt <- rbindlist(results_summary)
  fwrite(summary_dt,
         file.path(OUT_DIR, "bayesprism_celltype_de_runinfo.csv"))
  cat("Per-(cell_type, contrast) result counts:\n")
  print(summary_dt[, .(cell_type, contrast, n_genes, n_sig)])
}

cat(sprintf("\nDone: %s\n", Sys.time()))
