#!/usr/bin/env Rscript
# 17c_edger_celltypes.R
# ============================================================================
# Donor-level edgeR-QLF DA for the three non-hepatocyte Fig S4g compartments
# (Stellate, Macrophage/Kupffer, Cholangiocyte). Same method as 17_pseudobulk_
# de_retest.R arm (a) — the CANONICAL donor-level test (Squair et al. 2021:
# pseudobulk + edgeR empirical Bayes; effective n = donors, 13 MASLD vs 5
# NORMAL) — just looped over the count matrices 07c_corrected_da_celltypes.py
# exports. Hepatocytes already have scatac_da_corrected_hep_edger.csv (from 17);
# this completes the 4-compartment matched set the figure reads.
#
# Env: rnaseq (edgeR + limma). Reads {label}_pseudobulk_{counts.tsv.gz,
# coldata.tsv}; writes scatac_da_corrected_{label}_edger.csv (schema
# feature,logFC,PValue,FDR — identical to the hepatocyte file).
# ============================================================================

suppressPackageStartupMessages({
  library(edgeR)
})

PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
SNAP_DIR <- file.path(PROJECT_ROOT, "Analysis", "ATAC", "Human_Multiome",
                      "results", "snapatac2")

LABELS <- c("stellate", "macrophage", "cholangiocyte")

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# FAST wide-matrix reader: read.delim takes ~30-40 min on these ~200K-COLUMN
# pseudobulk TSVs (per-column type inference); the matrix is only ~18 rows, so
# readLines+strsplit is ~1000x faster. Output is identical to
# read.delim(gzfile(path), row.names = 1, check.names = FALSE).
fast_read_counts <- function(path) {
  con <- gzfile(path, "r"); L <- readLines(con); close(con)
  peaks <- strsplit(L[1], "\t", fixed = TRUE)[[1]][-1]
  body  <- L[-1]; nd <- length(body)
  donors <- character(nd); mat <- matrix(0, nrow = nd, ncol = length(peaks))
  for (i in seq_len(nd)) {
    f <- strsplit(body[i], "\t", fixed = TRUE)[[1]]
    donors[i] <- f[1]; mat[i, ] <- as.numeric(f[-1])
  }
  rownames(mat) <- donors; colnames(mat) <- peaks
  as.data.frame(mat, check.names = FALSE)
}

run_edger <- function(label) {
  counts_path  <- file.path(SNAP_DIR, sprintf("%s_pseudobulk_counts.tsv.gz", label))
  coldata_path <- file.path(SNAP_DIR, sprintf("%s_pseudobulk_coldata.tsv", label))
  out_path     <- file.path(SNAP_DIR, sprintf("scatac_da_corrected_%s_edger.csv", label))

  if (!file.exists(counts_path) || !file.exists(coldata_path)) {
    msg("  [SKIP] %s: pseudobulk inputs missing (run 07c first)", label)
    return(invisible(NULL))
  }

  counts_in <- fast_read_counts(counts_path)
  coldata   <- read.delim(coldata_path, stringsAsFactors = FALSE)
  rownames(coldata) <- as.character(coldata$donor_id)
  coldata <- coldata[rownames(counts_in), , drop = FALSE]
  condition <- as.integer(coldata$condition)  # 1 = MASLD, 0 = NORMAL
  stopifnot(!any(is.na(condition)))
  msg("  [%s] %d donors (%d MASLD, %d NORMAL) x %d peaks",
      label, nrow(counts_in), sum(condition == 1L), sum(condition == 0L),
      ncol(counts_in))

  mat <- t(as.matrix(counts_in))              # peaks x donors
  storage.mode(mat) <- "double"

  dge <- DGEList(counts = mat)
  dge <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~ condition)         # coef 2 = MASLD effect
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)
  tt  <- topTags(qlf, n = Inf, sort.by = "none")$table

  res <- data.frame(
    feature = rownames(mat),
    logFC   = tt$logFC,
    PValue  = tt$PValue,
    FDR     = tt$FDR,
    stringsAsFactors = FALSE
  )
  res <- res[order(res$FDR, res$PValue), , drop = FALSE]
  write.csv(res, out_path, row.names = FALSE)

  n_sig <- sum(res$FDR < 0.05, na.rm = TRUE)
  n_up  <- sum(res$FDR < 0.05 & res$logFC > 0, na.rm = TRUE)
  n_dn  <- sum(res$FDR < 0.05 & res$logFC < 0, na.rm = TRUE)
  msg("  [%s] wrote %s", label, out_path)
  msg("  [%s] n sig FDR<0.05: %d (%d open / %d close) of %d peaks",
      label, n_sig, n_up, n_dn, nrow(res))
  invisible(res)
}

msg("============================================================")
msg("17c_edger_celltypes.R — donor-level edgeR-QLF (3 compartments)")
msg("  effective n = DONORS (13 MASLD vs 5 NORMAL)")
msg("============================================================")
for (lab in LABELS) run_edger(lab)
msg("")
msg("Done.")
