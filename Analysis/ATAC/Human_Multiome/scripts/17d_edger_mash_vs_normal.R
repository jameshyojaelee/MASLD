#!/usr/bin/env Rscript
# 17d_edger_mash_vs_normal.R
# ============================================================================
# Documented ENDPOINT contrast for Fig S4g: MASH vs NORMAL, donor-level edgeR-QLF.
#
# WHY (2026-07-09) -----------------------------------------------------------
# A literal F0-vs-F4 contrast is NOT available in GSE244832: all 18 donors have
# fibrosis_stage_documented = NA; the only F-stage labels are scVI-INFERRED
# (F_stage_source == scvi_predicted), which the project rule forbids as a
# cross-cohort analysis axis (validated only on the 58 documented Andrews
# donors; this is the Kim/Rosenthal cohort). The defensible extreme-severity
# endpoint that IS documented is the source paper's own primary comparison:
#   MASH (n=9, documented CRN fibrosis stage 2-4) vs NORMAL (n=5, CRN<3 ~ F0),
# dropping the intermediate MASL steatosis donors. No inference, no rule break.
#
# The donor x peak pseudobulk COUNT matrices are contrast-agnostic (07c/07b
# already wrote them for all 18 donors), so this only re-groups donors and
# re-runs edgeR — no h5ad re-aggregation. Same method as 17/17c arm (a).
#
# Env: rnaseq (edgeR). Writes scatac_da_mashvsnormal_{label}_edger.csv
# (schema feature,logFC,PValue,FDR — matches the MASLD-vs-control files).
# ============================================================================

suppressPackageStartupMessages({ library(edgeR) })

PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
ATAC_DIR <- file.path(PROJECT_ROOT, "Analysis", "ATAC", "Human_Multiome")
SNAP_DIR <- file.path(ATAC_DIR, "results", "snapatac2")
META_PATH <- file.path(ATAC_DIR, "metadata", "donor_metadata_curated.tsv")

# label -> (pseudobulk-count basename, output label). Hepatocyte counts use the
# 'hep_pseudobulk' basename (07b); the other three use '{label}_pseudobulk' (07c).
COMPARTMENTS <- list(
  list(label = "hep",           counts = "hep_pseudobulk"),
  list(label = "stellate",      counts = "stellate_pseudobulk"),
  list(label = "macrophage",    counts = "macrophage_pseudobulk"),
  list(label = "cholangiocyte", counts = "cholangiocyte_pseudobulk")
)

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# donor_id -> documented condition (MASH / MASL / NORMAL).
donor_condition <- local({
  meta <- read.delim(META_PATH, stringsAsFactors = FALSE)
  setNames(toupper(trimws(as.character(meta$condition))),
           as.character(meta$donor_id))
})

run_edger_mvn <- function(comp) {
  counts_path <- file.path(SNAP_DIR, sprintf("%s_counts.tsv.gz", comp$counts))
  out_path    <- file.path(SNAP_DIR, sprintf("scatac_da_mashvsnormal_%s_edger.csv", comp$label))
  if (!file.exists(counts_path)) {
    msg("  [SKIP] %s: counts missing (%s)", comp$label, counts_path)
    return(invisible(NULL))
  }

  counts_in <- read.delim(gzfile(counts_path), row.names = 1, check.names = FALSE)
  donors <- rownames(counts_in)
  cond   <- donor_condition[donors]

  # MASH (1) vs NORMAL (0); drop MASL + any donor without a documented condition.
  grp <- rep(NA_integer_, length(donors))
  grp[cond == "MASH"]   <- 1L
  grp[cond == "NORMAL"] <- 0L
  keep <- !is.na(grp)
  n_mash <- sum(grp[keep] == 1L); n_norm <- sum(grp[keep] == 0L)
  msg("  [%s] %d MASH vs %d NORMAL (dropped %d MASL/other) x %d peaks",
      comp$label, n_mash, n_norm, sum(!keep), ncol(counts_in))
  if (n_mash < 2L || n_norm < 2L) {
    msg("  [%s] [SKIP] need >=2 donors/group", comp$label)
    return(invisible(NULL))
  }

  mat <- t(as.matrix(counts_in[keep, , drop = FALSE]))  # peaks x donors
  storage.mode(mat) <- "double"
  condition <- grp[keep]

  dge <- DGEList(counts = mat)
  dge <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~ condition)     # coef 2 = MASH effect
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)
  tt  <- topTags(qlf, n = Inf, sort.by = "none")$table

  res <- data.frame(
    feature = rownames(mat), logFC = tt$logFC,
    PValue = tt$PValue, FDR = tt$FDR, stringsAsFactors = FALSE
  )
  res <- res[order(res$FDR, res$PValue), , drop = FALSE]
  write.csv(res, out_path, row.names = FALSE)

  n_sig <- sum(res$FDR < 0.05, na.rm = TRUE)
  n_up  <- sum(res$FDR < 0.05 & res$logFC > 0, na.rm = TRUE)
  n_dn  <- sum(res$FDR < 0.05 & res$logFC < 0, na.rm = TRUE)
  msg("  [%s] wrote %s", comp$label, out_path)
  msg("  [%s] n sig FDR<0.05: %d (%d open / %d close) of %d peaks",
      comp$label, n_sig, n_up, n_dn, nrow(res))
  invisible(res)
}

msg("============================================================")
msg("17d_edger_mash_vs_normal.R — MASH vs NORMAL donor-level edgeR-QLF")
msg("  documented endpoint (MASH fibrosis 2-4 vs NORMAL ~F0); no F-stage inference")
msg("============================================================")
for (comp in COMPARTMENTS) run_edger_mvn(comp)
msg("")
msg("Done.")
