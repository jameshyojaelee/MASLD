#!/usr/bin/env Rscript
# 17_pseudobulk_de_retest.R
# ============================================================================
# Powerful-and-valid DONOR-level differential test for the MASLD scATAC pipeline
# (Squair et al. 2021, Nat Commun: pseudobulk + edgeR/limma empirical Bayes).
#
# WHY -----------------------------------------------------------------------
# The corrected scATAC differential tests collapse cells -> donors (the donor is
# the true replicate; 18 donors = 13 disease [MASL+MASH] vs 5 control [NORMAL]),
# which fixes pseudoreplication but leaves us with the OPPOSITE problem: a plain
# Wilcoxon-on-donor-means / per-feature OLS over n=18 has almost no power and
# returns ~0 significant features (an underpowered floor, NOT a biological null).
#
# Squair et al. 2021 show that the correct, more-powerful donor-level test is
# pseudobulk + edgeR/limma with empirical-Bayes moderation: information is shared
# across features to stabilise per-feature dispersion/variance, recovering power
# that the per-feature Wilcoxon/OLS throws away, WITHOUT reintroducing
# pseudoreplication (the model still sees n = n_donors).
#
# This script re-tests the three differential layers and writes *_edger / *_limma
# outputs NEXT TO the conservative originals, printing "n sig FDR<0.05" for each
# so the gain over the conservative ~0 floor is explicit. It changes NONE of the
# original outputs.
#
# ARMS ----------------------------------------------------------------------
#   (a) DA peaks            edgeR-QLF on donor x peak pseudobulk COUNTS
#   (b) SCENIC regulons     limma-eBayes on donor x regulon CONTINUOUS means
#   (c) chromVAR per CT     limma-eBayes on donor x TF CONTINUOUS means
#
# Effective n is the DONOR count throughout. With only 5 controls this is still
# the correct, maximally-powered valid test; whatever it finds is reported
# honestly (eBayes does not manufacture significance from nothing).
#
# Env: rnaseq (edgeR 4.4.2 + limma 3.62.2).
# Run: /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript \
#        Analysis/ATAC/Human_Multiome/scripts/17_pseudobulk_de_retest.R
# ============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
ATAC_DIR  <- file.path(PROJECT_ROOT, "Analysis", "ATAC", "Human_Multiome")
RES_DIR   <- file.path(ATAC_DIR, "results")
SNAP_DIR  <- file.path(RES_DIR, "snapatac2")
SCENIC_DIR <- file.path(ATAC_DIR, "scenic_plus")
CHROMVAR_DIR <- file.path(RES_DIR, "chromvar_v2")
META_PATH <- file.path(ATAC_DIR, "metadata", "donor_metadata_curated.tsv")

# Arm (a) inputs (written by 07b_corrected_da_hepatocytes.py)
COUNTS_PATH  <- file.path(SNAP_DIR, "hep_pseudobulk_counts.tsv.gz")
COLDATA_PATH <- file.path(SNAP_DIR, "hep_pseudobulk_coldata.tsv")
# Arm (b) inputs
REG_ACT_PATH <- file.path(SCENIC_DIR, "regulon_activity_scores.csv")
CELLMAP_PATH <- file.path(SNAP_DIR, "cell_donor_condition_map.tsv.gz")  # from 07b
# Arm (c) input
CHROMVAR_PER_DONOR <- file.path(CHROMVAR_DIR, "chromvar_per_donor.tsv.gz")

# Outputs
OUT_A <- file.path(SNAP_DIR, "scatac_da_corrected_hep_edger.csv")
OUT_B <- file.path(SCENIC_DIR, "disease_regulons_limma.csv")
OUT_C <- file.path(CHROMVAR_DIR, "chromvar_limma_per_ct.csv")

# condition -> MASLD(1)/NORMAL(0). MASLD = {MASL, MASH}; control = {NORMAL}.
MASLD_SET  <- c("MASLD", "MASH", "MASL", "NASH", "NAFLD", "NAFL")
NORMAL_SET <- c("NORMAL", "HEALTHY", "CONTROL")

cond_to_masld <- function(x) {
  xu <- toupper(trimws(as.character(x)))
  out <- rep(NA_real_, length(xu))
  out[xu %in% MASLD_SET]  <- 1
  out[xu %in% NORMAL_SET] <- 0
  out
}

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# Load donor -> condition(0/1) map from curated metadata (donor_id = D01..D18).
load_donor_condition <- function() {
  meta <- read.delim(META_PATH, stringsAsFactors = FALSE)
  stopifnot(all(c("donor_id", "condition") %in% colnames(meta)))
  data.frame(
    donor_id  = as.character(meta$donor_id),
    condition = cond_to_masld(meta$condition),
    stringsAsFactors = FALSE
  )
}

# ===========================================================================
# ARM (a): DA peaks — edgeR-QLF on donor x peak pseudobulk COUNTS
# ===========================================================================
run_arm_a <- function() {
  msg("")
  msg("=== ARM (a): DA peaks — edgeR-QLF on pseudobulk COUNTS ===")
  if (!file.exists(COUNTS_PATH) || !file.exists(COLDATA_PATH)) {
    msg("  [SKIP] pseudobulk counts/coldata not found:")
    msg("         %s", COUNTS_PATH)
    msg("         %s", COLDATA_PATH)
    msg("         -> run 07b_corrected_da_hepatocytes.py first (it exports them).")
    return(invisible(NULL))
  }

  # counts: rows = donors, cols = peaks (donor_id index col). edgeR wants
  # features (peaks) as ROWS and samples (donors) as COLUMNS -> transpose.
  counts_in <- read.delim(COUNTS_PATH, row.names = 1, check.names = FALSE)
  coldata   <- read.delim(COLDATA_PATH, stringsAsFactors = FALSE)
  msg("  Loaded counts: %d donors x %d peaks", nrow(counts_in), ncol(counts_in))

  # Align coldata to the count-matrix donor order.
  rownames(coldata) <- as.character(coldata$donor_id)
  coldata <- coldata[rownames(counts_in), , drop = FALSE]
  condition <- as.integer(coldata$condition)  # 1 = MASLD, 0 = NORMAL
  stopifnot(!any(is.na(condition)))
  msg("  Donors: %d MASLD, %d NORMAL",
      sum(condition == 1L), sum(condition == 0L))

  mat <- t(as.matrix(counts_in))             # peaks x donors
  storage.mode(mat) <- "double"

  # The peaks were already filtered upstream in 07b (>= min_total_reads total AND
  # detected in >= min_donor_detect donors). Keep that filtered set as-is.
  dge <- DGEList(counts = mat)
  dge <- calcNormFactors(dge, method = "TMM")

  design <- model.matrix(~ condition)        # coef 2 = MASLD effect
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
  write.csv(res, OUT_A, row.names = FALSE)

  n_sig <- sum(res$FDR < 0.05, na.rm = TRUE)
  msg("  Wrote %s", OUT_A)
  msg("  ARM (a) n sig FDR<0.05: %d  (of %d peaks tested)", n_sig, nrow(res))
  invisible(res)
}

# ===========================================================================
# ARM (b): SCENIC regulon activity — limma-eBayes on CONTINUOUS donor means
# ===========================================================================
run_arm_b <- function() {
  msg("")
  msg("=== ARM (b): SCENIC regulons — limma-eBayes on donor-mean activity ===")
  if (!file.exists(REG_ACT_PATH)) {
    msg("  [SKIP] regulon activity matrix not found: %s", REG_ACT_PATH)
    return(invisible(NULL))
  }
  # The regulon activity CSV is keyed by bare cell barcode with NO donor column.
  # We need a barcode -> donor map. 07b exports one (cell_donor_condition_map.tsv.gz)
  # because recovering donor labels otherwise requires reading the 3.5 GB
  # label-transferred h5ad obs in Python (anndata is not usable cleanly in R).
  if (!file.exists(CELLMAP_PATH)) {
    msg("  [BLOCKED] cell->donor map not found: %s", CELLMAP_PATH)
    msg("            This is emitted by 07b_corrected_da_hepatocytes.py. The")
    msg("            regulon_activity_scores.csv has bare barcodes only; the")
    msg("            barcode->donor->condition mapping lives in the h5ad obs.")
    msg("            Re-run 07b (it now writes the map), then re-run arm (b).")
    msg("  ARM (b) n sig adj.P.Val<0.05: NA (stubbed — donor mapping missing)")
    return(invisible(NULL))
  }

  # cells x regulons (first col = barcode).
  act <- read.csv(REG_ACT_PATH, row.names = 1, check.names = FALSE)
  msg("  Loaded regulon activity: %d cells x %d regulons", nrow(act), ncol(act))

  cmap <- read.delim(CELLMAP_PATH, stringsAsFactors = FALSE)
  stopifnot(all(c("barcode", "donor_id") %in% colnames(cmap)))
  rownames(cmap) <- cmap$barcode

  common <- intersect(rownames(act), rownames(cmap))
  msg("  Barcodes matched to donor map: %d / %d", length(common), nrow(act))
  if (length(common) == 0L) {
    msg("  [BLOCKED] zero barcode overlap between regulon CSV and cell->donor map.")
    msg("  ARM (b) n sig adj.P.Val<0.05: NA (stubbed — no barcode overlap)")
    return(invisible(NULL))
  }
  act  <- act[common, , drop = FALSE]
  donor_of_cell <- cmap[common, "donor_id"]

  # Aggregate cells -> donor MEAN activity (the donor is the unit).
  donors <- sort(unique(donor_of_cell))
  donor_x_reg <- t(vapply(donors, function(d) {
    colMeans(act[donor_of_cell == d, , drop = FALSE])
  }, FUN.VALUE = numeric(ncol(act))))
  rownames(donor_x_reg) <- donors                  # donors x regulons
  colnames(donor_x_reg) <- colnames(act)

  # Donor -> condition(0/1).
  dc <- load_donor_condition()
  rownames(dc) <- dc$donor_id
  condition <- dc[rownames(donor_x_reg), "condition"]
  keep_d <- !is.na(condition)
  donor_x_reg <- donor_x_reg[keep_d, , drop = FALSE]
  condition   <- condition[keep_d]
  msg("  Donors used: %d MASLD, %d NORMAL",
      sum(condition == 1), sum(condition == 0))
  if (sum(condition == 1) < 2L || sum(condition == 0) < 2L) {
    msg("  [BLOCKED] need >=2 donors per group.")
    msg("  ARM (b) n sig adj.P.Val<0.05: NA (stubbed — too few donors/group)")
    return(invisible(NULL))
  }

  # limma on the CONTINUOUS donor x regulon matrix. lmFit wants features (regulons)
  # as ROWS, samples (donors) as COLUMNS -> transpose. eBayes moderates the
  # per-regulon residual variance across regulons (the power recovery).
  design <- model.matrix(~ condition)              # coef 2 = MASLD effect
  fit <- lmFit(t(donor_x_reg), design)
  fit <- eBayes(fit)
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  # logFC here = MASLD - NORMAL mean difference in activity units (continuous).
  res <- data.frame(
    tf_name    = rownames(tt),
    logFC      = tt$logFC,
    P.Value    = tt$P.Value,
    adj.P.Val  = tt$adj.P.Val,
    n_donors   = nrow(donor_x_reg),
    stringsAsFactors = FALSE
  )
  res <- res[order(res$adj.P.Val, res$P.Value), , drop = FALSE]
  write.csv(res, OUT_B, row.names = FALSE)

  n_sig <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
  msg("  Wrote %s", OUT_B)
  msg("  ARM (b) n sig adj.P.Val<0.05: %d  (of %d regulons tested)",
      n_sig, nrow(res))
  invisible(res)
}

# ===========================================================================
# ARM (c): chromVAR — limma-eBayes per cell type on donor-mean deviations
# ===========================================================================
run_arm_c <- function() {
  msg("")
  msg("=== ARM (c): chromVAR — limma-eBayes per cell type ===")
  if (!file.exists(CHROMVAR_PER_DONOR)) {
    msg("  [SKIP] chromVAR per-donor file not found: %s", CHROMVAR_PER_DONOR)
    return(invisible(NULL))
  }

  # Long format: donor_id_atac | cell_type | TF | TF_symbol | mean_deviation | ...
  # NB: despite its name, donor_id_atac holds D## values matching metadata donor_id.
  long <- read.delim(gzfile(CHROMVAR_PER_DONOR), stringsAsFactors = FALSE)
  stopifnot(all(c("donor_id_atac", "cell_type", "TF", "mean_deviation")
                %in% colnames(long)))
  msg("  Loaded chromVAR long table: %d rows, %d cell types, %d TFs",
      nrow(long), length(unique(long$cell_type)), length(unique(long$TF)))

  dc <- load_donor_condition()
  rownames(dc) <- dc$donor_id

  all_res <- list()
  for (ct in sort(unique(long$cell_type))) {
    sub <- long[long$cell_type == ct, , drop = FALSE]

    # Pivot to donor x TF (mean_deviation). Use unique TF column id.
    tfs    <- sort(unique(sub$TF))
    donors <- sort(unique(sub$donor_id_atac))
    M <- matrix(NA_real_, nrow = length(donors), ncol = length(tfs),
                dimnames = list(donors, tfs))
    M[cbind(match(sub$donor_id_atac, donors),
            match(sub$TF, tfs))] <- sub$mean_deviation

    # Donor -> condition(0/1); drop donors with unknown condition or all-NA rows.
    condition <- dc[rownames(M), "condition"]
    keep_d <- !is.na(condition)
    M <- M[keep_d, , drop = FALSE]
    condition <- condition[keep_d]
    n_masld  <- sum(condition == 1)
    n_normal <- sum(condition == 0)
    if (n_masld < 2L || n_normal < 2L || nrow(M) < 3L) {
      msg("  [%s] skipped: MASLD=%d NORMAL=%d (need >=2/group)",
          ct, n_masld, n_normal)
      next
    }

    # Drop TFs that are constant / all-NA across the retained donors (limma
    # cannot fit a moderated variance for a zero-variance feature).
    Y <- t(M)                                  # TFs x donors
    finite_per_tf <- rowSums(is.finite(Y))
    var_per_tf <- apply(Y, 1, function(v) {
      v <- v[is.finite(v)]
      if (length(v) < 2L) return(0)
      stats::var(v)
    })
    keep_tf <- (finite_per_tf == ncol(Y)) & (var_per_tf > 0)
    if (sum(keep_tf) < 1L) {
      msg("  [%s] skipped: no testable TFs after NA/zero-variance filter", ct)
      next
    }
    Y <- Y[keep_tf, , drop = FALSE]

    design <- model.matrix(~ condition)        # coef 2 = MASLD effect
    fit <- lmFit(Y, design)
    fit <- eBayes(fit)
    tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

    res_ct <- data.frame(
      cell_type = ct,
      TF        = rownames(tt),
      logFC     = tt$logFC,             # MASLD - NORMAL mean deviation difference
      P.Value   = tt$P.Value,
      adj.P.Val = tt$adj.P.Val,
      n_donors  = nrow(M),
      stringsAsFactors = FALSE
    )
    all_res[[ct]] <- res_ct
    n_sig_ct <- sum(res_ct$adj.P.Val < 0.05, na.rm = TRUE)
    msg("  [%s] %d TFs tested, %d donors (MASLD=%d NORMAL=%d) -> n sig FDR<0.05: %d",
        ct, nrow(res_ct), nrow(M), n_masld, n_normal, n_sig_ct)
  }

  if (length(all_res) == 0L) {
    msg("  [SKIP] no cell type had enough donors per group to test.")
    return(invisible(NULL))
  }
  out <- do.call(rbind, all_res)
  out <- out[order(out$cell_type, out$adj.P.Val, out$P.Value), , drop = FALSE]
  write.csv(out, OUT_C, row.names = FALSE)
  msg("  Wrote %s", OUT_C)
  msg("  ARM (c) total n sig adj.P.Val<0.05 (all cell types): %d",
      sum(out$adj.P.Val < 0.05, na.rm = TRUE))
  invisible(out)
}

# ===========================================================================
# Main
# ===========================================================================
main <- function() {
  msg("============================================================")
  msg("17_pseudobulk_de_retest.R — donor-level eBayes re-test")
  msg("  Squair et al. 2021: pseudobulk + edgeR/limma empirical Bayes")
  msg("  Effective n = DONORS (13 disease vs 5 control)")
  msg("============================================================")
  run_arm_a()
  run_arm_b()
  run_arm_c()
  msg("")
  msg("Done.")
}

if (sys.nframe() == 0L) {
  main()
}
