#!/usr/bin/env Rscript
# 352_dream_per_celltype_sex.R — Dream sex interaction model per cell type.
#
# For one cell type (passed via CELLTYPE env / argv), fits two models on
# donor pseudobulk:
#   PART 1 — interaction:  ~ condition_binary * inferred_sex + (1|dataset)
#   PART 2 — stratified male / female (matched datasets)
# Then classifies each gene as Concordant / Female_biased / Male_biased / Divergent
# using the same logic as RNA-seq Script 26.
#
# Inputs:
#   Analysis/SingleCell/results_gpu_v2/sex_celltype/<ct>_counts.csv.gz
#   Analysis/SingleCell/results_gpu_v2/sex_celltype/<ct>_meta.csv
#
# Outputs:
#   Analysis/SingleCell/results_gpu_v2/sex_celltype/<ct>_sex_interaction_dream.csv
#   Analysis/SingleCell/results_gpu_v2/sex_celltype/<ct>_sex_dream_results.csv  (final per-gene class)
#   Analysis/SingleCell/results_gpu_v2/sex_celltype/<ct>_sex_dream_summary.txt

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(variancePartition)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
ct <- if (length(args) >= 1) args[[1]] else Sys.getenv("CELLTYPE", "hepatocytes")
proj_root <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
setwd(proj_root)

ddir <- "Analysis/SingleCell/results_gpu_v2/sex_celltype"
cfile <- file.path(ddir, paste0(ct, "_counts.csv.gz"))
mfile <- file.path(ddir, paste0(ct, "_meta.csv"))
stopifnot(file.exists(cfile), file.exists(mfile))

cat("=== Cell type:", ct, "===\n")
cat("counts:", cfile, "\nmeta:", mfile, "\n")

ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
param <- SnowParam(workers = ncpu, type = "SOCK", progressbar = TRUE,
                  RNGseed = 42L)
cat("Parallel workers:", ncpu, "\n")

# --- load counts (genes x donors, possibly sparse but small enough) ---
counts <- as.data.frame(fread(cfile, header = TRUE))
rownames(counts) <- counts[[1]]
counts[[1]] <- NULL
counts <- as.matrix(counts)
storage.mode(counts) <- "double"
cat("Counts:", nrow(counts), "genes x", ncol(counts), "donors\n")

meta <- fread(mfile)
stopifnot(all(meta$sample == colnames(counts)))

# Standardize types
meta[, condition_binary := factor(condition_binary, levels = c("Control", "Disease"))]
meta[, inferred_sex := factor(inferred_sex, levels = c("Female", "Male"))]
meta[, dataset := factor(dataset)]

cat("\nSample counts by condition x sex:\n")
print(table(meta$condition_binary, meta$inferred_sex, useNA = "ifany"))

# Drop donors with NA sex (k-means always assigns, but defensive)
keep <- !is.na(meta$inferred_sex) & !is.na(meta$condition_binary)
if (any(!keep)) {
  cat("Dropping", sum(!keep), "donors with NA sex/condition\n")
  meta <- meta[keep]
  counts <- counts[, meta$sample, drop = FALSE]
}

# Need both sexes and both conditions present
sex_x_cond <- table(meta$condition_binary, meta$inferred_sex)
if (any(sex_x_cond < 2)) {
  cat("Insufficient donors for sex x condition contrast in", ct,
      "— at least one cell is < 2. Skipping interaction model.\n")
  print(sex_x_cond)
  # write empty placeholders so the chain doesn't fail
  fwrite(data.table(),
         file.path(ddir, paste0(ct, "_sex_interaction_dream.csv")))
  fwrite(data.table(),
         file.path(ddir, paste0(ct, "_sex_dream_results.csv")))
  quit(status = 0)
}

# --- filter low expression: at least 5 counts in >=20% donors ---
keep_g <- rowSums(counts >= 5) >= max(3, floor(0.20 * ncol(counts)))
cat("Genes passing filter:", sum(keep_g), "/", length(keep_g), "\n")
counts <- counts[keep_g, , drop = FALSE]

# --- DGEList + voomWithDreamWeights ---
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")

# Random-effect dataset only if more than one level
multi_ds <- nlevels(droplevels(meta$dataset)) > 1
form_int <- if (multi_ds) {
  ~ condition_binary * inferred_sex + (1 | dataset)
} else {
  ~ condition_binary * inferred_sex
}
cat("Interaction formula:", deparse(form_int), "\n")

cat("Running voomWithDreamWeights...\n")
v_int <- variancePartition::voomWithDreamWeights(dge, form_int, data = meta,
                                                 BPPARAM = param)
cat("Running dream interaction...\n")
fit_int <- variancePartition::dream(v_int, form_int, data = meta, BPPARAM = param)
fit_int <- variancePartition::eBayes(fit_int)
coefs <- colnames(fit_int$coefficients)
cat("Coefficients:", paste(coefs, collapse = ", "), "\n")
int_coef <- grep("condition_binary.*inferred_sex|inferred_sex.*condition_binary",
                 coefs, value = TRUE)
if (length(int_coef) == 0) {
  stop("No interaction coefficient found")
}
cat("Using interaction coef:", int_coef, "\n")

tt_int <- topTable(fit_int, coef = int_coef, number = Inf, sort.by = "none")
res_int <- data.table(gene = rownames(tt_int),
                      logFC = tt_int$logFC,
                      AveExpr = tt_int$AveExpr,
                      t = tt_int$t,
                      pvalue = tt_int$P.Value,
                      padj = tt_int$adj.P.Val)
fwrite(res_int, file.path(ddir, paste0(ct, "_sex_interaction_dream.csv")))
cat("Wrote interaction dream results.\n")

# --- Part 2: stratified dream (male / female separately) ---
run_stratum <- function(sex_label) {
  m_s <- meta[inferred_sex == sex_label]
  n_per_cond <- table(m_s$condition_binary)
  if (any(n_per_cond < 2)) {
    cat("Stratum", sex_label, "has <2 donors per condition; skipping.\n")
    return(NULL)
  }
  c_s <- counts[, m_s$sample, drop = FALSE]
  # Per-stratum re-filter: drop genes that are all-zero or zero-variance in this
  # sex stratum (e.g. Y-chr genes when stratum=Female). Without this, voom's
  # lowess hits 'delta must be finite and > 0' on stratified subsets.
  keep_strat <- rowSums(c_s >= 1) >= 3 & matrixStats::rowSds(c_s) > 0
  if (sum(keep_strat) < nrow(c_s)) {
    cat("Stratum", sex_label, "drop", sum(!keep_strat),
        "zero-var / low-count genes (",
        sum(keep_strat), "remaining)\n")
    c_s <- c_s[keep_strat, , drop = FALSE]
  }
  dge_s <- DGEList(counts = c_s)
  dge_s <- calcNormFactors(dge_s, method = "TMM")
  multi_ds_s <- nlevels(droplevels(m_s$dataset)) > 1
  form_s <- if (multi_ds_s) {
    ~ condition_binary + (1 | dataset)
  } else {
    ~ condition_binary
  }
  cat("Stratum", sex_label, "formula:", deparse(form_s), "n=", nrow(m_s), "\n")
  v_s <- variancePartition::voomWithDreamWeights(dge_s, form_s, data = m_s,
                                                 BPPARAM = param)
  fit_s <- variancePartition::dream(v_s, form_s, data = m_s, BPPARAM = param)
  fit_s <- variancePartition::eBayes(fit_s)
  cs <- colnames(fit_s$coefficients)
  disease_coef <- grep("condition_binaryDisease|condition_binary.Disease",
                       cs, value = TRUE)
  if (length(disease_coef) == 0) {
    cat("No disease coefficient in stratum", sex_label, "\n")
    return(NULL)
  }
  tt <- topTable(fit_s, coef = disease_coef[1], number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, padj = tt$adj.P.Val)
}

cat("\n--- Stratified dream (Female) ---\n")
res_f <- run_stratum("Female")
cat("\n--- Stratified dream (Male) ---\n")
res_m <- run_stratum("Male")

# --- Part 3: sex_class assignment ---
cls <- copy(res_int)
setnames(cls, c("logFC", "padj"), c("interaction_logFC", "interaction_padj"))
if (!is.null(res_f)) {
  cls <- merge(cls,
               res_f[, .(gene, logFC_F = logFC, padj_F = padj)],
               by = "gene", all.x = TRUE)
} else {
  cls[, `:=`(logFC_F = NA_real_, padj_F = NA_real_)]
}
if (!is.null(res_m)) {
  cls <- merge(cls,
               res_m[, .(gene, logFC_M = logFC, padj_M = padj)],
               by = "gene", all.x = TRUE)
} else {
  cls[, `:=`(logFC_M = NA_real_, padj_M = NA_real_)]
}

INT_PADJ <- 0.05
DIV_MIN_LFC <- 0.1
cls[, sex_dimorphic := !is.na(interaction_padj) & interaction_padj < INT_PADJ]
cls[, sex_class := fcase(
  !sex_dimorphic, "Concordant",
  sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
    sign(logFC_F) != sign(logFC_M) &
    abs(logFC_F) > DIV_MIN_LFC & abs(logFC_M) > DIV_MIN_LFC,
  "Divergent",
  sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
    abs(logFC_F) >= abs(logFC_M), "Female_biased",
  sex_dimorphic & !is.na(logFC_F) & !is.na(logFC_M) &
    abs(logFC_M) > abs(logFC_F),  "Male_biased",
  default = "Unclassified"
)]

cat("\nsex_class breakdown:\n")
print(cls[, .N, by = sex_class][order(-N)])

fwrite(cls, file.path(ddir, paste0(ct, "_sex_dream_results.csv")))

# Summary
sink(file.path(ddir, paste0(ct, "_sex_dream_summary.txt")))
cat("Cell type:", ct, "\n")
cat("N donors:", nrow(meta), "\n")
cat("Sex x Condition table:\n"); print(table(meta$condition_binary, meta$inferred_sex))
cat("\nN genes tested:", nrow(cls), "\n")
cat("Sex-dimorphic (interaction padj <", INT_PADJ, "):", sum(cls$sex_dimorphic), "\n")
cat("Class breakdown:\n"); print(cls[, .N, by = sex_class][order(-N)])
sink()
cat("Done", ct, "\n")
