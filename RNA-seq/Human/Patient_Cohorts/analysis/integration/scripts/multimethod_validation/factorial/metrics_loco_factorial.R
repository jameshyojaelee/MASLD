#!/usr/bin/env Rscript
# =============================================================================
# metrics_loco_factorial.R
# -----------------------------------------------------------------------------
# CONSUMER scoring for the factorial DE-method benchmark -- the R1 axis of the
# Step-2 selection composite (preregistration.yaml step2_rank.R1_loco_repro-
# ducibility). PURE consumer: it reads the per-fold TRAIN DEG tables already
# written by the degx engine `--loco` mode and compares each to the held-out
# cohort's INDEPENDENT per-study DE. NO model is (re)fit here.
#
# THE METRIC (per the frozen pre-registration, R1):
#   For each cell, for each leave-one-COhort-out fold, the cell was fit on the
#   4 NON-held-out cohorts; we compare that TRAIN result to the held-out
#   cohort's per-study DE (the held-out cohort never trained the cell), then
#   AVERAGE the four sub-metrics across folds. Each cell yields one R1 scalar:
#     loco_repro_scalar = mean( loco_jaccard,
#                               (loco_lfc_spearman + 1)/2,
#                               loco_dir_concord / 100,
#                               loco_auc )
#   Aggregation order (prereg R1 aggregation_across_folds): MEAN each sub-metric
#   across folds first, THEN average the four mean sub-metrics into the scalar.
#
# REUSE (studied, then CALLED -- not re-derived):
#   * de_validation_helpers.R::compute_loocv_metrics (L378-449) -- the canonical
#     per-fold metric kernel used by the production LOOCV harness. We SOURCE the
#     helper module and call it so the four sub-metrics are byte-for-byte the
#     same computation the live pipeline uses.
#       compute_loocv_metrics(method_dt, ho_dt, method_name) returns a 1-row
#       data.table whose columns map to the four prereg sub-metrics as:
#         jaccard_01            -> loco_jaccard        (Jaccard of {padj<0.1})
#         lfc_spearman          -> loco_lfc_spearman   (Spearman train vs HO logFC)
#         direction_concordance -> loco_dir_concord    (% sign match, 0-100)
#         auc_replication       -> loco_auc            (AUC |train logFC| ~ HO padj<0.05)
#       (it also returns fisher_or / pi1_heldout / lfc_spearman_sig etc., which
#        are NOT part of the R1 scalar and are ignored here.)
#     method_dt schema expected by the kernel: .(gene, logFC, padj[, stat]).
#       The degx LOCO fold files carry gene,logFC,SE,pval,padj -> a superset, so
#       they feed the kernel directly (the kernel selects logFC + padj).
#     ho_dt schema expected by the kernel: .(gene, ho_logFC, ho_padj, ho_pval).
#
# HELD-OUT REFERENCE DE (resolved to match loocv_multimethod.R exactly):
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study/
#     {COHORT}_de_results.csv  with columns logFC, adj.P.Val, P.Value, gene.
#   -> ho_dt <- .(gene, ho_logFC = logFC, ho_padj = adj.P.Val, ho_pval = P.Value)
#   This is the SAME reference + the SAME column convention loocv_multimethod.R
#   uses for its held-out target (see that script, section 4).
#
# GENE-ID CONVENTION: both the LOCO fold files and the per-study reference carry
#   VERSIONED ensembl ids (ENSG....N), so the merge in compute_loocv_metrics
#   joins them directly. As a safety net, if the raw versioned join yields no
#   common genes (e.g. an unversioned reference), we retry on the version-
#   stripped base id. The training-side table is never altered beyond id form.
#
# Env: micromamba run -n rnaseq Rscript metrics_loco_factorial.R [args]
# CLI:
#   --loco_dir   dir of per-fold CSVs (default RNA-seq/results/degx_factorial/loco)
#                files: cell_<cell_id>__heldout_<COHORT>.csv
#                optional manifest: loco_manifest_disease_vs_control.csv
#   --out_dir    output dir (default RNA-seq/results/degx_factorial)
#                writes metrics_loco.csv
#
# OUTPUT metrics_loco.csv:
#   cell_id, n_folds, loco_jaccard, loco_lfc_spearman, loco_dir_concord,
#   loco_auc, loco_repro_scalar
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# 0. paths + CLI
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEGX_DIR <- file.path(BASE, "RNA-seq/results/degx_factorial")
INT_DIR  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PERSTUDY_DIR <- file.path(INT_DIR, "results/per_study")
HELPERS  <- file.path(INT_DIR,
  "scripts/multimethod_validation/de_validation_helpers.R")

# Source the helper module to reuse compute_loocv_metrics. The module also loads
# variancePartition/edgeR/etc.; sourcing is heavier than strictly needed for a
# CSV-only scorer but guarantees we call the IDENTICAL metric kernel rather than
# re-deriving it (the whole point of the reuse mandate). If the heavy deps are
# unavailable, fall back to a self-contained kernel transcribed from L378-449.
have_helper_module <- tryCatch({
  suppressWarnings(suppressMessages(source(HELPERS)))
  is.function(get0("compute_loocv_metrics"))
}, error = function(e) {
  cat("[warn] could not source de_validation_helpers.R (", conditionMessage(e),
      "); using transcribed metric kernel.\n", sep = "")
  FALSE
})

# Self-contained fallback: a verbatim transcription of the four R1 sub-metrics
# from compute_loocv_metrics (L395-422). Used only when the helper module's
# heavy Bioconductor deps cannot load on the current node. Returns the same
# four named columns the kernel does (plus n_genes_common), so the downstream
# extraction is identical regardless of which path produced the row.
compute_loco_kernel <- function(method_dt, ho_dt, method_name) {
  md <- method_dt[!is.na(gene), .(gene, train_logFC = logFC, train_padj = padj)]
  common <- merge(md, ho_dt, by = "gene")
  n_common <- nrow(common)
  na_row <- data.table(method = method_name, n_genes_common = n_common,
                       jaccard_01 = NA_real_, lfc_spearman = NA_real_,
                       direction_concordance = NA_real_, auc_replication = NA_real_)
  if (n_common < 10) return(na_row)

  train_sig01 <- common[!is.na(train_padj) & train_padj < 0.1, gene]
  ho_sig01    <- common[!is.na(ho_padj)    & ho_padj    < 0.1, gene]
  un <- length(union(train_sig01, ho_sig01))
  jaccard_01 <- if (un > 0) length(intersect(train_sig01, ho_sig01)) / un else 0

  lfc_spearman <- suppressWarnings(cor(common$train_logFC, common$ho_logFC,
                    method = "spearman", use = "pairwise.complete.obs"))
  tsig <- common[!is.na(train_padj) & train_padj < 0.1]
  dir_conc <- if (nrow(tsig) > 0)
    mean(sign(tsig$train_logFC) == sign(tsig$ho_logFC)) * 100 else NA_real_

  ok_auc  <- !is.na(common$ho_padj) & is.finite(common$train_logFC)
  ho_bin  <- as.integer(common$ho_padj[ok_auc] < 0.05)
  lfc_auc <- abs(common$train_logFC[ok_auc])
  pos <- lfc_auc[ho_bin == 1]; neg <- lfc_auc[ho_bin == 0]
  auc_rep <- if (length(pos) > 0 && length(neg) > 0) {
    wt <- suppressWarnings(wilcox.test(pos, neg))
    as.numeric(wt$statistic) / (length(pos) * length(neg))
  } else NA_real_

  data.table(method = method_name, n_genes_common = n_common,
             jaccard_01 = round(jaccard_01, 4),
             lfc_spearman = round(lfc_spearman, 4),
             direction_concordance = round(dir_conc, 1),
             auc_replication = round(auc_rep, 4))
}

# Single entry point: prefer the sourced kernel, else the transcription.
score_fold <- function(method_dt, ho_dt, cell_id) {
  if (have_helper_module) compute_loocv_metrics(method_dt, ho_dt, cell_id)
  else                    compute_loco_kernel(method_dt, ho_dt, cell_id)
}

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  out <- list(loco_dir = file.path(DEGX_DIR, "loco"),
              out_dir  = DEGX_DIR)
  i <- 1
  while (i <= length(a)) {
    key <- sub("^--", "", a[i])
    if (key %in% names(out)) { out[[key]] <- a[i + 1]; i <- i + 2 }
    else { i <- i + 1 }
  }
  out
}
args <- parse_args()
dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

strip_v <- function(x) sub("[.][0-9]+$", "", x)

cat("=== metrics_loco_factorial.R (CONSUMER, R1 reproducibility) ===\n")
cat("loco_dir :", args$loco_dir, "\n")
cat("out_dir  :", args$out_dir, "\n")
cat("metric kernel:", if (have_helper_module)
      "de_validation_helpers.R::compute_loocv_metrics (sourced)"
    else "transcribed fallback (helper deps unavailable)", "\n\n")

if (!dir.exists(args$loco_dir))
  stop("loco_dir does not exist: ", args$loco_dir,
       "\n(The degx --loco mode has not written fold files yet.)")

# ---------------------------------------------------------------------------
# 1. held-out reference DE (per-study), loaded lazily + cached per cohort
# ---------------------------------------------------------------------------
# Matches loocv_multimethod.R section 4 exactly: ho_dt columns gene, ho_logFC,
# ho_padj, ho_pval from {COHORT}_de_results.csv (logFC, adj.P.Val, P.Value).
.ho_cache <- new.env(parent = emptyenv())
load_heldout <- function(cohort) {
  if (!is.null(.ho_cache[[cohort]])) return(.ho_cache[[cohort]])
  f <- file.path(PERSTUDY_DIR, paste0(cohort, "_de_results.csv"))
  if (!file.exists(f)) {
    .ho_cache[[cohort]] <- NULL
    return(NULL)
  }
  ho <- fread(f)
  req <- c("gene", "logFC", "adj.P.Val", "P.Value")
  if (!all(req %in% names(ho))) {
    cat(sprintf("  [skip] held-out reference %s missing cols: %s\n",
                basename(f), paste(setdiff(req, names(ho)), collapse = ",")))
    .ho_cache[[cohort]] <- NULL
    return(NULL)
  }
  ho_dt <- ho[, .(gene, ho_logFC = logFC, ho_padj = adj.P.Val, ho_pval = P.Value)]
  ho_dt <- ho_dt[!is.na(gene) & gene != ""]
  .ho_cache[[cohort]] <- ho_dt
  ho_dt
}

# ---------------------------------------------------------------------------
# 2. fold discovery (manifest preferred; glob fallback for SMOKE)
#    File convention: cell_<cell_id>__heldout_<COHORT>.csv
# ---------------------------------------------------------------------------
FOLD_RE <- "^cell_(.+)__heldout_([^.]+)\\.csv$"

discover_folds <- function() {
  manifest <- file.path(args$loco_dir, "loco_manifest_disease_vs_control.csv")
  if (file.exists(manifest)) {
    mf <- fread(manifest)
    need <- c("cell_id", "heldout_cohort")
    if (all(need %in% names(mf))) {
      if ("status" %in% names(mf)) mf <- mf[status == "ok"]
      mf <- mf[, .(cell_id = as.character(cell_id),
                   heldout_cohort = as.character(heldout_cohort))]
      mf[, file := file.path(args$loco_dir,
            sprintf("cell_%s__heldout_%s.csv", cell_id, heldout_cohort))]
      mf <- mf[file.exists(file)]
      cat(sprintf("Manifest mode: %d ok folds with existing files\n", nrow(mf)))
      if (nrow(mf) > 0) return(mf[, .(cell_id, heldout_cohort, file)])
      cat("  (manifest had no usable rows -> falling back to glob)\n")
    } else {
      cat("  (manifest lacks cell_id/heldout_cohort -> glob mode)\n")
    }
  }
  fs <- list.files(args$loco_dir, pattern = FOLD_RE, full.names = TRUE)
  fs <- fs[!grepl("manifest", basename(fs), ignore.case = TRUE)]
  if (length(fs) == 0) {
    cat("Glob mode: no cell_*__heldout_*.csv fold files found.\n")
    return(data.table(cell_id = character(), heldout_cohort = character(),
                      file = character()))
  }
  bn  <- basename(fs)
  cid <- sub(FOLD_RE, "\\1", bn)
  hoc <- sub(FOLD_RE, "\\2", bn)
  cat(sprintf("Glob mode (no manifest): %d fold files in %s\n",
              length(fs), args$loco_dir))
  data.table(cell_id = cid, heldout_cohort = hoc, file = fs)
}

folds <- discover_folds()
if (nrow(folds) == 0)
  stop("No LOCO fold files to score in ", args$loco_dir,
       " (expected cell_<id>__heldout_<COHORT>.csv).")

setorder(folds, cell_id, heldout_cohort)
cat(sprintf("Discovered %d fold(s) across %d cell(s).\n\n",
            nrow(folds), uniqueN(folds$cell_id)))

# ---------------------------------------------------------------------------
# 3. read a TRAIN fold table -> .(gene, logFC, padj) feed for the kernel
# ---------------------------------------------------------------------------
read_train_fold <- function(f) {
  dt <- fread(f)
  req <- c("gene", "logFC", "padj")          # SE / pval present but unused here
  miss <- setdiff(req, names(dt))
  if (length(miss))
    stop(sprintf("fold %s missing columns: %s", basename(f),
                 paste(miss, collapse = ",")))
  dt[!is.na(gene) & gene != "", .(gene, logFC, padj)]
}

# Run the metric kernel for one fold; retry on version-stripped ids if the raw
# versioned join produced < 10 common genes (defensive: unversioned reference).
score_one <- function(train_dt, ho_dt, cell_id) {
  res <- score_fold(train_dt, ho_dt, cell_id)
  if (is.na(res$n_genes_common) || res$n_genes_common >= 10) return(res)
  td2 <- copy(train_dt)[, gene := strip_v(gene)][!duplicated(gene)]
  ho2 <- copy(ho_dt)[,  gene := strip_v(gene)][!duplicated(gene)]
  res2 <- score_fold(td2, ho2, cell_id)
  if (!is.na(res2$n_genes_common) && res2$n_genes_common >= 10) {
    cat(sprintf("    [%s] versioned join sparse (n=%s); used version-stripped join (n=%d).\n",
                cell_id, res$n_genes_common, res2$n_genes_common))
    return(res2)
  }
  res
}

# ---------------------------------------------------------------------------
# 4. per-fold metrics -> average across folds per cell -> R1 scalar
# ---------------------------------------------------------------------------
perfold_rows <- list()
for (k in seq_len(nrow(folds))) {
  cid <- folds$cell_id[k]; hoc <- folds$heldout_cohort[k]; f <- folds$file[k]

  ho_dt <- load_heldout(hoc)
  if (is.null(ho_dt)) {
    cat(sprintf("  [skip fold] cell=%s heldout=%s : no held-out reference DE.\n", cid, hoc))
    next
  }
  train_dt <- tryCatch(read_train_fold(f), error = function(e) {
    cat(sprintf("  [skip fold] cell=%s heldout=%s : %s\n", cid, hoc, conditionMessage(e)))
    NULL
  })
  if (is.null(train_dt) || nrow(train_dt) == 0) {
    cat(sprintf("  [skip fold] cell=%s heldout=%s : empty/failed train table.\n", cid, hoc))
    next
  }

  m <- tryCatch(score_one(train_dt, ho_dt, cid), error = function(e) {
    cat(sprintf("  [skip fold] cell=%s heldout=%s : metric error %s\n",
                cid, hoc, conditionMessage(e)))
    NULL
  })
  if (is.null(m)) next

  perfold_rows[[length(perfold_rows) + 1L]] <- data.table(
    cell_id           = cid,
    heldout_cohort    = hoc,
    n_genes_common    = as.numeric(m$n_genes_common),
    loco_jaccard      = as.numeric(m$jaccard_01),
    loco_lfc_spearman = as.numeric(m$lfc_spearman),
    loco_dir_concord  = as.numeric(m$direction_concordance),
    loco_auc          = as.numeric(m$auc_replication))
}

if (length(perfold_rows) == 0)
  stop("No fold produced metrics (all skipped/failed).")

perfold <- rbindlist(perfold_rows)

# Average each sub-metric across folds (prereg: mean across folds FIRST), counting
# only folds where the sub-metric is defined (n_folds = folds that scored at all).
agg <- perfold[, .(
  n_folds           = .N,
  loco_jaccard      = mean(loco_jaccard,      na.rm = TRUE),
  loco_lfc_spearman = mean(loco_lfc_spearman, na.rm = TRUE),
  loco_dir_concord  = mean(loco_dir_concord,  na.rm = TRUE),
  loco_auc          = mean(loco_auc,          na.rm = TRUE)
), by = cell_id]

# NaN (a sub-metric NA in every fold) -> NA for honest reporting.
for (col in c("loco_jaccard", "loco_lfc_spearman", "loco_dir_concord", "loco_auc"))
  agg[is.nan(get(col)), (col) := NA_real_]

# R1 scalar: average the four bounded constituents AFTER folding.
#   jaccard in [0,1]; (spearman+1)/2 in [0,1]; dir_concord/100 in [0,1]; auc in [0,1].
agg[, loco_repro_scalar := rowMeans(
  cbind(loco_jaccard,
        (loco_lfc_spearman + 1) / 2,
        loco_dir_concord / 100,
        loco_auc),
  na.rm = TRUE)]
agg[is.nan(loco_repro_scalar), loco_repro_scalar := NA_real_]

# round for readability (keep enough precision for ranking)
for (col in c("loco_jaccard", "loco_lfc_spearman", "loco_auc", "loco_repro_scalar"))
  agg[, (col) := round(get(col), 4)]
agg[, loco_dir_concord := round(loco_dir_concord, 1)]

setcolorder(agg, c("cell_id", "n_folds", "loco_jaccard", "loco_lfc_spearman",
                   "loco_dir_concord", "loco_auc", "loco_repro_scalar"))
setorder(agg, cell_id)

# ---------------------------------------------------------------------------
# 5. write outputs
# ---------------------------------------------------------------------------
out_path     <- file.path(args$out_dir, "metrics_loco.csv")
perfold_path <- file.path(args$out_dir, "metrics_loco_perfold.csv")
fwrite(agg,     out_path)
fwrite(perfold, perfold_path)

cat("\n=== wrote ===\n")
cat(" R1 LOCO reproducibility :", out_path, sprintf("(%d cells)\n", nrow(agg)))
cat(" per-fold detail         :", perfold_path, sprintf("(%d folds)\n", nrow(perfold)))

cat("\n--- metrics_loco.csv ---\n")
print(agg)

# sanity: every reported value bounded as the prereg requires.
chk <- agg[!is.na(loco_repro_scalar)]
ok_bounds <- all(
  chk$loco_jaccard      >= 0 & chk$loco_jaccard      <= 1,
  chk$loco_lfc_spearman >= -1 & chk$loco_lfc_spearman <= 1,
  chk$loco_dir_concord  >= 0 & chk$loco_dir_concord  <= 100,
  chk$loco_auc          >= 0 & chk$loco_auc          <= 1,
  chk$loco_repro_scalar >= 0 & chk$loco_repro_scalar <= 1,
  na.rm = TRUE)
cat(sprintf("\nSanity: all reported sub-metrics within bounds: %s\n",
            if (isTRUE(ok_bounds)) "YES" else "NO -- INVESTIGATE"))
cat("Done.\n")
