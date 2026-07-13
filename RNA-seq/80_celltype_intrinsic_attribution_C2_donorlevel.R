#!/usr/bin/env Rscript
# 80_celltype_intrinsic_attribution_C2_donorlevel.R
#
# Donor-collapse re-run of Script 80 (per-cell-type intrinsic attribution of
# bulk DEGs). IDENTICAL logic to 80_celltype_intrinsic_attribution_C2.R; the
# ONLY changes are:
#   (1) SC_DE dir + OUTDIR read from env so run-level vs donor-collapsed
#       pseudobulk DE can be swapped side-by-side WITHOUT touching the canonical
#       c2_treat outputs.
#   (2) The per-cell-type glob is restricted to BASE cell-type files
#       (exclude allcell_, `_vs_` stage-contrast, and `bayesprism` files).
#       The canonical script globs `_de\.csv$` and unintentionally swept in the
#       Jun-8 stage-contrast files (e.g. Hepatocytes_Steatohepatitis_vs_Steatosis)
#       and May-1 bayesprism per-stage files that were deposited in the same dir
#       AFTER the script was written — 56 of its 60 `intrinsic_strong` genes map
#       to those stage-contrast pseudo-celltypes. Those files were NOT part of the
#       run->donor pseudoreplication fix, so an apples-to-apples before/after of
#       the donor-collapse EFFECT must hold the cell-type set fixed to the base
#       per-cell-type DE. An optional CT_ALLOWLIST pins the exact set.
#
# Env: rnaseq
# Env vars:
#   SC_DE_DIR      : dir with {CellType}_de.csv (default = canonical run-level)
#   ATTR_OUTDIR    : output dir (default = .../celltype_attribution/donorlevel/base)
#   CT_ALLOWLIST   : optional comma-sep base cell-type names to restrict to
#                    (basenames minus `_de.csv`). If unset, all base cell types
#                    in SC_DE_DIR are used.

suppressPackageStartupMessages({
  library(data.table)
})

BASE    <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT_RES <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
CAUSAL  <- file.path(BASE, "RNA-seq/results/causal_inference")

SC_DE   <- Sys.getenv("SC_DE_DIR",
                      file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de"))
OUTDIR  <- Sys.getenv("ATTR_OUTDIR",
                      file.path(BASE, "RNA-seq/results/celltype_attribution/donorlevel/base"))
CT_ALLOWLIST <- Sys.getenv("CT_ALLOWLIST", "")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message(sprintf("SC_DE_DIR = %s", SC_DE))
message(sprintf("OUTDIR    = %s", OUTDIR))

# Thresholds
PVAL_CT     <- 0.05   # nominal pvalue for per-cell-type concordance
TREAT_FDR   <- 0.05   # bulk C2 Tier-1 DEG gate (TREAT canonical; the 1,918 set)

message("[1] Loading bulk C2 canonical DE (limma-voom-qw, TREAT Tier-1)...")
bulk <- fread(file.path(INT_RES, "canonical_deg_results.csv"),
              select = c("gene", "symbol", "logFC", "t", "padj", "treat_fdr"))
setnames(bulk, c("logFC","t","padj"), c("bulk_lfc","bulk_t","bulk_padj"))
bulk[, ensg_base := sub("\\.\\d+$", "", gene)]
message(sprintf("  %d genes in bulk C2 canonical (%d Tier-1 TREAT DEGs at treat_fdr<%.2f)",
                nrow(bulk), sum(bulk$treat_fdr < TREAT_FDR, na.rm = TRUE), TREAT_FDR))
stopifnot(sum(bulk$treat_fdr < TREAT_FDR, na.rm = TRUE) == 1918L)

message("[2] Loading per-cell-type pseudobulk DE (BASE cell types only)...")
sc_de_files <- list.files(SC_DE, pattern = "_de\\.csv$", full.names = TRUE)
# Restrict to base per-cell-type DE: drop allcell_, stage-contrast (_vs_),
# and bayesprism per-stage files that share the dir.
sc_de_files <- sc_de_files[!grepl("allcell_|_vs_|bayesprism", basename(sc_de_files))]
celltype_from_file <- function(f) sub("_de\\.csv$", "", basename(f))
if (nzchar(CT_ALLOWLIST)) {
  allow <- trimws(strsplit(CT_ALLOWLIST, ",")[[1]])
  sc_de_files <- sc_de_files[celltype_from_file(sc_de_files) %in% allow]
}

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

message("[3] Loading MuSiC hepatocyte attribution (Script 25, C2 recount)...")
music <- fread(file.path(CAUSAL, "c2_recount", "deconv_attribution_scores.csv"))
music[, ensg_base := sub("\\.\\d+$", "", gene)]
music_sub <- music[, .(ensg_base,
                       music_category = category,
                       music_score    = attribution_raw)]

message("[4] Merging bulk + scRNA + MuSiC...")
merged <- merge(bulk, sc_wide, by = "ensg_base", all.x = TRUE)
merged <- merge(merged, music_sub, by = "ensg_base", all.x = TRUE)
message(sprintf("  Joined: %d genes", nrow(merged)))

message("[5] Computing signed attribution score per cell type...")
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

merged[!(treat_fdr < TREAT_FDR), attribution_class := "NS_bulk"]
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
deg <- merged[treat_fdr < TREAT_FDR]
summ <- deg[, .N, by = attribution_class][order(-N)]
multict <- deg[attribution_class == "multi_celltype", .N]
strong_by_ct <- deg[grepl("_intrinsic_strong$", attribution_class), .N, by = attribution_class][order(-N)]
likely_by_ct <- deg[grepl("_intrinsic_likely$", attribution_class), .N, by = attribution_class][order(-N)]
hep_agree <- deg[!is.na(music_category) & !is.na(primary_celltype),
                 .N, by = .(music_category, primary_celltype)][order(-N)][1:20]

summary_lines <- c(
  sprintf("SC_DE_DIR: %s", SC_DE),
  sprintf("Cell types (%d): %s", length(celltypes), paste(celltypes, collapse = ", ")),
  "",
  sprintf("Bulk Tier-1 TREAT DEGs (treat_fdr<%.2f): %d", TREAT_FDR, nrow(deg)),
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
  sprintf("Thresholds: PVAL_CT=%.2f (nominal), bulk Tier-1 gate = TREAT treat_fdr<%.2f (1,918 set)",
          PVAL_CT, TREAT_FDR))
writeLines(summary_lines, file.path(OUTDIR, "celltype_attribution_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
