#!/usr/bin/env Rscript
# 80_celltype_intrinsic_attribution.R
#
# Analysis A1 (v1) — Per-cell-type intrinsic attribution of bulk DEGs.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
# Plan: /gpfs/commons/home/jameslee/.claude/plans/resilient-purring-koala.md
#
# Strategy (v1, uses existing data):
#   1. Load bulk dream DE (dream_results_ashr.csv)
#   2. Load per-cell-type pseudobulk DE for 11 cell types.
#   3. Load MuSiC attribution (Script 25; hepatocyte-specific).
#   4. For each bulk DEG, build per-cell-type attribution:
#        - signed_score = sign(bulk_lfc) * scRNA_t_stat  (positive=concordant)
#        - sig_concordant = pvalue < 0.05 & concordant direction
#        - Primary celltype = cell type with max positive signed_score
#   5. Tiered classification:
#        - "<CT>_intrinsic_strong":  sig_concordant in best CT AND unique best
#        - "<CT>_intrinsic_likely":  concordant in best CT (not sig, unique best)
#        - "multi_celltype":         sig_concordant in >=2 CTs
#        - "bulk_only":              no CT concordant
#   6. Cross with MuSiC hepatocyte class.
#
# Note: scRNA per-cell-type pseudobulk DE has limited power (fewer donors per CT),
# so nominal p-value is used for significance. v2 (CARseq + bMIND) will provide
# bulk-level CT-specific DE with full bulk power.
#
# Env: rnaseq
# Outputs: RNA-seq/results/celltype_attribution/

suppressPackageStartupMessages({
  library(data.table)
})

BASE    <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_RES <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
CAUSAL  <- file.path(BASE, "RNA-seq/results/causal_inference")
SC_DE   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
OUTDIR  <- file.path(BASE, "RNA-seq/results/celltype_attribution")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Thresholds
PVAL_CT     <- 0.05   # nominal pvalue for per-cell-type concordance
PADJ_BULK   <- 0.05   # bulk dream DEG cutoff
LFC_MIN     <- 0.3    # bulk |logFC|

message("[1] Loading bulk dream DE...")
bulk <- fread(file.path(INT_RES, "dream_results_ashr.csv"),
              select = c("gene", "symbol", "logFC", "t", "padj"))
setnames(bulk, c("logFC","t","padj"), c("bulk_lfc","bulk_t","bulk_padj"))
bulk[, ensg_base := sub("\\.\\d+$", "", gene)]
message(sprintf("  %d genes in bulk dream", nrow(bulk)))

message("[2] Loading per-cell-type pseudobulk DE (11 cell types)...")
sc_de_files <- list.files(SC_DE, pattern = "_de\\.csv$", full.names = TRUE)
sc_de_files <- sc_de_files[!grepl("allcell_", basename(sc_de_files))]
celltype_from_file <- function(f) sub("_de\\.csv$", "", basename(f))

sc_list <- lapply(sc_de_files, function(f) {
  ct <- celltype_from_file(f)
  dt <- fread(f)
  # Handle both limma (t, P.Value) and legacy (t_stat, pvalue) column names
  if ("t_stat" %in% names(dt)) setnames(dt, "t_stat", "t", skip_absent = TRUE)
  if ("pvalue" %in% names(dt)) setnames(dt, "pvalue", "P.Value", skip_absent = TRUE)
  if ("adj.P.Val" %in% names(dt) && !"padj" %in% names(dt)) setnames(dt, "adj.P.Val", "padj")
  dt <- dt[, .SD, .SDcols = intersect(c("gene","logFC","t","P.Value","padj"), names(dt))]
  setnames(dt, c("logFC","t","P.Value","padj"),
           paste0(c("lfc_","t_","pval_","padj_"), ct), skip_absent = TRUE)
  dt[, ensg_base := sub("\\.\\d+$", "", gene)]
  dt[, gene := NULL]
  dt
})
names(sc_list) <- sapply(sc_de_files, celltype_from_file)
celltypes <- names(sc_list)
message(sprintf("  Cell types (%d): %s", length(celltypes), paste(celltypes, collapse = ", ")))

sc_wide <- Reduce(function(a,b) merge(a, b, by = "ensg_base", all = TRUE), sc_list)

message("[3] Loading MuSiC hepatocyte attribution (Script 25)...")
music <- fread(file.path(CAUSAL, "deconv_attribution_scores.csv"))
music[, ensg_base := sub("\\.\\d+$", "", gene)]
music_sub <- music[, .(ensg_base,
                       music_category = category,
                       music_score    = attribution_raw)]

message("[4] Merging bulk + scRNA + MuSiC...")
merged <- merge(bulk, sc_wide, by = "ensg_base", all.x = TRUE)
merged <- merge(merged, music_sub, by = "ensg_base", all.x = TRUE)
message(sprintf("  Joined: %d genes", nrow(merged)))

message("[5] Computing signed attribution score per cell type...")
# signed_score = sign(bulk_lfc) * scRNA_t  (positive => concordant direction)
for (ct in celltypes) {
  t_col    <- paste0("t_",    ct)
  pval_col <- paste0("pval_", ct)
  lfc_col  <- paste0("lfc_",  ct)

  score_col <- paste0("score_", ct)
  sig_col   <- paste0("sig_concordant_", ct)

  merged[, (score_col) := sign(bulk_lfc) * get(t_col)]
  merged[, (sig_col)   := !is.na(get(pval_col)) & get(pval_col) < PVAL_CT &
                           !is.na(get(lfc_col)) & !is.na(bulk_lfc) &
                           sign(get(lfc_col)) == sign(bulk_lfc)]
}

message("[6] Primary attribution + tiered classification...")
score_cols <- paste0("score_",           celltypes)
sig_cols   <- paste0("sig_concordant_",  celltypes)

classify_row <- function(scores, sigs) {
  # scores: numeric per celltype; sigs: logical per celltype
  scores[is.na(scores)] <- -Inf
  order_idx <- order(scores, decreasing = TRUE)
  best <- order_idx[1]
  best_score <- scores[best]
  best_ct <- celltypes[best]

  n_sig_concordant <- sum(sigs, na.rm = TRUE)

  class <- if (best_score <= 0 || is.infinite(best_score)) {
    "bulk_only"
  } else if (n_sig_concordant >= 2) {
    "multi_celltype"
  } else if (isTRUE(sigs[best])) {
    paste0(best_ct, "_intrinsic_strong")
  } else {
    paste0(best_ct, "_intrinsic_likely")
  }
  c(best_ct, as.character(best_score), as.character(n_sig_concordant), class)
}

score_mat <- as.matrix(merged[, ..score_cols])
sig_mat   <- as.matrix(merged[, ..sig_cols])

class_out <- vapply(seq_len(nrow(merged)),
                    function(i) classify_row(score_mat[i,], sig_mat[i,]),
                    character(4))
merged[, primary_celltype       := class_out[1, ]]
merged[, primary_score          := as.numeric(class_out[2, ])]
merged[, n_sig_concordant_ct    := as.integer(class_out[3, ])]
merged[, attribution_class      := class_out[4, ]]

# Restrict attribution_class label to bulk DEGs; otherwise NS
merged[!(bulk_padj < PADJ_BULK & abs(bulk_lfc) > LFC_MIN),
       attribution_class := "NS_bulk"]
merged[primary_score <= 0, primary_celltype := NA]

message("[7] Writing outputs...")
fwrite(merged, file.path(OUTDIR, "celltype_attribution_matrix.csv"))
message(sprintf("  Wrote celltype_attribution_matrix.csv (%d rows)", nrow(merged)))

primary <- merged[, .(ensg_base, symbol,
                      bulk_lfc, bulk_padj,
                      primary_celltype, primary_score,
                      n_sig_concordant_ct,
                      attribution_class,
                      music_category, music_score)]
fwrite(primary, file.path(OUTDIR, "celltype_primary_attribution.csv"))

# Summary
deg <- merged[bulk_padj < PADJ_BULK & abs(bulk_lfc) > LFC_MIN]
summ <- deg[, .N, by = attribution_class][order(-N)]
multict <- deg[attribution_class == "multi_celltype", .N]
strong_by_ct <- deg[grepl("_intrinsic_strong$", attribution_class), .N, by = attribution_class][order(-N)]
likely_by_ct <- deg[grepl("_intrinsic_likely$", attribution_class), .N, by = attribution_class][order(-N)]
hep_agree <- deg[!is.na(music_category) & !is.na(primary_celltype),
                 .N, by = .(music_category, primary_celltype)][order(-N)][1:20]

summary_lines <- c(
  sprintf("Bulk DEGs (padj<%.2f, |LFC|>%.1f): %d",
          PADJ_BULK, LFC_MIN, nrow(deg)),
  "",
  "Per-class attribution:",
  capture.output(print(summ, nrows = 40)),
  "",
  sprintf("Multi-celltype DEGs (sig+concordant in >=2 CTs): %d", multict),
  "",
  "Strong intrinsic (sig+concordant, unique best) by cell type:",
  capture.output(print(strong_by_ct)),
  "",
  "Likely intrinsic (concordant, not sig) by cell type:",
  capture.output(print(likely_by_ct)),
  "",
  "MuSiC hepatocyte class vs primary celltype agreement (top 20):",
  capture.output(print(hep_agree)),
  "",
  sprintf("Thresholds: PVAL_CT=%.2f (nominal), PADJ_BULK=%.2f, LFC_MIN=%.1f",
          PVAL_CT, PADJ_BULK, LFC_MIN),
  "",
  "NOTE: scRNA pseudobulk DE has limited power (fewer donors per CT).",
  "      v2 (CARseq + bMIND) will provide bulk-level CT-specific DE.")
writeLines(summary_lines, file.path(OUTDIR, "celltype_attribution_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
