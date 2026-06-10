#!/usr/bin/env Rscript
# ===========================================================================
# 14b_nas_fibrosis_stage_lvqw.R
# ---------------------------------------------------------------------------
# HARMONIZED replacement for the GRANULAR staging analyses
# 14b_nas_fibrosis_stage_dream.R + 14d_consecutive_stage_dream.R, swapped from
# dream() to the canonical limma-voom-quality-weighted (LVQW) engine
# (de_engine_lvqw.R). This was MISSED in the 2026-06-08 WS2 harmonization; a PI
# spot-check caught it.
#
# What this produces (OVERWRITES the dream-named files in-place after backup):
#   nas_score_dream.csv            — per-NAS-score signatures   NAS{k}_vs_NAS0
#   fibrosis_stage_dream.csv       — per-fibrosis-stage sigs     F{k}_vs_F0
#   nas_consecutive_dream.csv      — consecutive NAS transitions NAS{k}_vs_NAS{k-1}
#   fibrosis_consecutive_dream.csv — consecutive F transitions   F{k}_vs_F{k-1}
#
# Engine swap (dream -> LVQW)
# ---------------------------
#   * (1|dataset) random intercept -> FIXED `dataset` effect (method consistency
#     with the C2 canonical). inferred_sex stays a nuisance covariate.
#   * voomWithDreamWeights+dream -> voomWithQualityWeights -> lmFit -> eBayes,
#     ashr-shrunk (mixcompdist="normal"), all via de_engine_lvqw.R.
#
# Staging definitions PRESERVED VERBATIM from 14b/14d
# ---------------------------------------------------
#   * NAS universe: GSE130970/135251/162694/174478/193066; nas_score collapsed
#     7+ -> 7; reference = NAS0.
#   * Fibrosis universe: GSE130970/135251/162694/174478/193066/240729 (the
#     original FIB_DATASETS — GSE213621 excluded as it has grouped, not
#     individual, fibrosis stages); fibrosis_stage in {0..4}; reference = F0.
#   * Same QC gate: pass_technical == TRUE.
#
# COHORT-EXPANDED (per task + 14_lvqw principle)
# ----------------------------------------------
#   These are disease-internal staging contrasts. For the cumulative-vs-reference
#   stage signatures (NAS{k}_vs_NAS0, F{k}_vs_F0) we fit ONE multi-level factor
#   design over the full annotated universe (mirrors 05c_stage_contrasts_lvqw.R:
#   one design, one coef per non-reference level). For consecutive transitions we
#   pool, per adjacent pair, EVERY cohort containing BOTH stages (mirrors
#   05b/14_lvqw: binary low/high, build_design_guarded drops degenerate dataset
#   terms). NO per-contrast re-TMM / re-filterByExpr — merged_dge is already
#   normalized/filtered (05h/14_lvqw input pattern).
#
# fibrosis_consecutive OVERLAP with fibrosis_pairwise.csv
# -------------------------------------------------------
#   The consecutive fibrosis transitions (F0->F1 ... F3->F4) are IDENTICAL biology
#   to the already-harmonized fibrosis_pairwise.csv (method=limma_voom_qw_C2,
#   2026-06-08; same binary low/high LVQW fits). To avoid a divergent second
#   analysis we DERIVE fibrosis_consecutive_dream.csv directly from
#   fibrosis_pairwise.csv (rename transition F{lo}_to_F{hi} -> contrast
#   F{hi}_vs_F{lo}; carry gene/logFC/SE/padj/shrunk_logFC/lfsr). nas_consecutive
#   has no existing LVQW equivalent -> fit fresh here.
#
# STAGING ONLY for review; OVERWRITES the dream-named granular files in place
# (the 5 fig2 panels read those exact paths). Backs up each dream CSV ->
# *_premigration_backup.csv first.
# ===========================================================================

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(ashr); library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR  <- file.path(INT, "results/integration")
SIGS  <- file.path(INT, "results/disease_signatures")
SCRDIR<- file.path(INT, "scripts")
stopifnot(dir.exists(SIGS), dir.exists(RDIR))

source(file.path(SCRDIR, "de_engine_lvqw.R"))

METHOD_TAG <- "limma_voom_qw"   # task contract: method="limma_voom_qw"
cat("=== 14b/14d LVQW: GRANULAR staging harmonization (dream -> LVQW) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# Input loading — 05h / 14_lvqw pattern. Subset only; NO re-TMM / re-filterByExpr.
# Staging fields (nas_score, fibrosis_stage) from unified_metadata.csv (the
# original 14b/14d source); inferred_sex + dataset from meta_matched.rds.
# ---------------------------------------------------------------------------
dge_all <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta_m  <- as.data.table(readRDS(file.path(RDIR, "meta_matched.rds")))
ucols   <- fread(file.path(INT, "metadata/unified_metadata.csv"),
                 select = c("sample_id", "dataset", "fibrosis_stage", "nas_score"))
qc      <- fread(file.path(INT, "qc/sample_qc_report.csv"))
pass_ids <- qc[pass_technical == TRUE, sample_id]
cat(sprintf("[in] %d genes x %d samples in merged_dge; %d pass_technical\n",
            nrow(dge_all), ncol(dge_all), length(pass_ids)))

# unified staging table restricted to QC-passing samples present in the DGE
stage_meta <- merge(
  ucols,
  meta_m[, .(sample_id, inferred_sex)],
  by = "sample_id", all.x = TRUE)
stage_meta <- stage_meta[sample_id %in% pass_ids & sample_id %in% colnames(dge_all)]

# helper: subset the (already normalized) merged DGE; attach sample covariates.
subset_dge <- function(sample_ids, cov_dt) {
  idx <- colnames(dge_all) %in% sample_ids
  d   <- dge_all[, idx]
  cov <- cov_dt[match(colnames(d), cov_dt$sample_id)]
  d$samples$dataset      <- factor(cov$dataset)
  d$samples$inferred_sex <- factor(cov$inferred_sex)
  d$samples$stage_val    <- cov$stage_val
  d
}

spearman <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
}

sanity_rows <- list()
add_sanity <- function(contrast, new_dt, old_path, old_contrast_filter = NULL) {
  # new_dt: data.table with gene, logFC, padj. old_path: dream CSV (pre-backup).
  if (!file.exists(old_path)) {
    sanity_rows[[length(sanity_rows) + 1]] <<- data.table(
      contrast = contrast,
      new_n_deg = sum(new_dt$padj < 0.05, na.rm = TRUE),
      old_n_deg = NA_integer_, n_shared = 0L,
      dir_concord = NA_real_, logfc_spearman = NA_real_,
      note = "dream reference absent")
    return(invisible())
  }
  od <- fread(old_path)
  cc <- if ("contrast" %in% names(od)) "contrast" else NA
  if (!is.na(cc) && !is.null(old_contrast_filter))
    od <- od[get(cc) == old_contrast_filter]
  setnames(od, "padj", "old_padj", skip_absent = TRUE)
  setnames(od, "logFC", "old_logFC", skip_absent = TRUE)
  mg <- merge(new_dt[, .(gene, new_logFC = logFC, new_padj = padj)],
              od[, .(gene, old_logFC, old_padj)], by = "gene")
  dir_conc <- mean(sign(mg$new_logFC) == sign(mg$old_logFC), na.rm = TRUE)
  rho <- spearman(mg$new_logFC, mg$old_logFC)
  sanity_rows[[length(sanity_rows) + 1]] <<- data.table(
    contrast       = contrast,
    new_n_deg      = sum(new_dt$padj < 0.05, na.rm = TRUE),
    old_n_deg      = sum(od$old_padj < 0.05, na.rm = TRUE),
    n_shared       = nrow(mg),
    dir_concord    = round(dir_conc, 4),
    logfc_spearman = round(rho, 4),
    note           = "LVQW vs dream (cohort-expanded; dream z.std format)")
  cat(sprintf("  [sanity %s] new=%d old(dream)=%d shared=%d dir=%.3f rho=%.3f\n",
              contrast, sum(new_dt$padj < 0.05, na.rm = TRUE),
              sum(od$old_padj < 0.05, na.rm = TRUE), nrow(mg), dir_conc, rho))
}

# ===========================================================================
# ANALYSIS 1: Per-NAS-Score signatures (cumulative vs NAS0) — nas_score_dream.csv
#   Single multi-level factor design over the full NAS-annotated universe;
#   one coef per non-reference NAS level (cohort-expanded; dataset FIXED).
# ===========================================================================
cat("\n===== (1) Per-NAS-score signatures (LVQW, NAS{k}_vs_NAS0) =====\n")
NAS_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066")

m_nas <- stage_meta[dataset %in% NAS_DATASETS & !is.na(nas_score)]
m_nas[, nas_group := fifelse(nas_score >= 7, 7L, as.integer(nas_score))]
cat("NAS-annotated QC-passing samples:", nrow(m_nas), "\n")
print(m_nas[, .N, by = nas_group][order(nas_group)])

m_nas[, stage_val := nas_group]
dge_nas <- subset_dge(m_nas$sample_id, m_nas)
nas_lvls <- sort(unique(dge_nas$samples$stage_val))
info_nas <- data.frame(
  dataset      = factor(dge_nas$samples$dataset),
  inferred_sex = factor(dge_nas$samples$inferred_sex),
  nas_group    = factor(dge_nas$samples$stage_val, levels = nas_lvls),
  row.names    = colnames(dge_nas))
info_nas$nas_group <- relevel(info_nas$nas_group, ref = "0")

gd_nas <- build_design_guarded(info_nas, c("dataset", "inferred_sex", "nas_group"))
cat("NAS stage design:", gd_nas$formula_used, "\n")
if (length(gd_nas$dropped)) cat("  dropped:", paste(gd_nas$dropped, collapse = ", "), "\n")

nas_out <- list()
for (lvl in setdiff(levels(info_nas$nas_group), "0")) {
  cf <- paste0("nas_group", lvl)
  if (!cf %in% colnames(gd_nas$design)) { cat("  [skip] coef", cf, "absent\n"); next }
  res <- fit_lvqw(dge_nas, gd_nas$design, coef = cf, do_ashr = TRUE)
  res[, nas_level := as.integer(lvl)]
  res[, contrast  := paste0("NAS", lvl, "_vs_NAS0")]
  res[, method    := METHOD_TAG]
  nas_out[[lvl]] <- res
  cat(sprintf("  NAS%s_vs_NAS0: %d DEGs (padj<0.05) | %d (padj<0.05,|LFC|>0.5)\n",
              lvl, sum(res$padj < 0.05, na.rm = TRUE),
              sum(res$padj < 0.05 & abs(res$logFC) > 0.5, na.rm = TRUE)))
}
nas_sig <- rbindlist(nas_out, use.names = TRUE, fill = TRUE)
setcolorder(nas_sig, c("gene", "logFC", "SE", "t", "P.Value", "padj",
                       "shrunk_logFC", "lfsr", "AveExpr",
                       "nas_level", "contrast", "method"))

# ===========================================================================
# ANALYSIS 2: Per-Fibrosis-Stage signatures (cumulative vs F0) — fibrosis_stage_dream.csv
# ===========================================================================
cat("\n===== (2) Per-fibrosis-stage signatures (LVQW, F{k}_vs_F0) =====\n")
FIB_DATASETS <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478",
                  "GSE193066", "GSE240729")

m_fib <- stage_meta[dataset %in% FIB_DATASETS & !is.na(fibrosis_stage)]
m_fib[, fib_int := suppressWarnings(as.integer(as.character(fibrosis_stage)))]
m_fib <- m_fib[fib_int %in% 0:4]
cat("Fibrosis-annotated QC-passing samples:", nrow(m_fib), "\n")
print(m_fib[, .N, by = fib_int][order(fib_int)])

m_fib[, stage_val := fib_int]
dge_fib <- subset_dge(m_fib$sample_id, m_fib)
fib_lvls <- sort(unique(dge_fib$samples$stage_val))
info_fib <- data.frame(
  dataset      = factor(dge_fib$samples$dataset),
  inferred_sex = factor(dge_fib$samples$inferred_sex),
  fib_stage    = factor(dge_fib$samples$stage_val, levels = fib_lvls),
  row.names    = colnames(dge_fib))
info_fib$fib_stage <- relevel(info_fib$fib_stage, ref = "0")

gd_fib <- build_design_guarded(info_fib, c("dataset", "inferred_sex", "fib_stage"))
cat("Fibrosis stage design:", gd_fib$formula_used, "\n")
if (length(gd_fib$dropped)) cat("  dropped:", paste(gd_fib$dropped, collapse = ", "), "\n")

fib_out <- list()
for (lvl in setdiff(levels(info_fib$fib_stage), "0")) {
  cf <- paste0("fib_stage", lvl)
  if (!cf %in% colnames(gd_fib$design)) { cat("  [skip] coef", cf, "absent\n"); next }
  res <- fit_lvqw(dge_fib, gd_fib$design, coef = cf, do_ashr = TRUE)
  res[, fib_stage := as.integer(lvl)]
  res[, contrast  := paste0("F", lvl, "_vs_F0")]
  res[, method    := METHOD_TAG]
  fib_out[[lvl]] <- res
  cat(sprintf("  F%s_vs_F0: %d DEGs (padj<0.05) | %d (padj<0.05,|LFC|>0.5)\n",
              lvl, sum(res$padj < 0.05, na.rm = TRUE),
              sum(res$padj < 0.05 & abs(res$logFC) > 0.5, na.rm = TRUE)))
}
fib_sig <- rbindlist(fib_out, use.names = TRUE, fill = TRUE)
setcolorder(fib_sig, c("gene", "logFC", "SE", "t", "P.Value", "padj",
                       "shrunk_logFC", "lfsr", "AveExpr",
                       "fib_stage", "contrast", "method"))

# ===========================================================================
# ANALYSIS 3: NAS consecutive transitions (NAS{k}_vs_NAS{k-1}) — nas_consecutive_dream.csv
#   Binary low/high per adjacent pair, cohort-expanded (mirror 05b/14_lvqw).
# ===========================================================================
cat("\n===== (3) NAS consecutive transitions (LVQW, NAS{k}_vs_NAS{k-1}) =====\n")
nas_pairs <- lapply(1:7, function(k) c(lo = k - 1L, hi = k))
nas_consec_out <- list()
for (pr in nas_pairs) {
  lo <- pr[["lo"]]; hi <- pr[["hi"]]
  cname <- sprintf("NAS%d_vs_NAS%d", hi, lo)
  both <- m_nas[nas_group %in% c(lo, hi),
                .(ns = uniqueN(nas_group)), by = dataset][ns == 2L, dataset]
  mp <- m_nas[dataset %in% both & nas_group %in% c(lo, hi)]
  if (nrow(mp) < 6L || uniqueN(mp$nas_group) < 2L) {
    cat(sprintf("  [skip] %s (n=%d, cohorts={%s})\n",
                cname, nrow(mp), paste(sort(both), collapse = ","))); next
  }
  mp[, stage_val := nas_group]
  dge_p <- subset_dge(mp$sample_id, mp)
  grp <- factor(ifelse(dge_p$samples$stage_val == hi, "high", "low"),
                levels = c("low", "high"))
  info_p <- data.frame(
    dataset      = factor(dge_p$samples$dataset),
    inferred_sex = factor(dge_p$samples$inferred_sex),
    stage_grp    = grp, row.names = colnames(dge_p))
  gd_p <- build_design_guarded(info_p, c("dataset", "inferred_sex", "stage_grp"))
  if (!"stage_grphigh" %in% colnames(gd_p$design)) {
    cat(sprintf("  [skip] %s coef absent (%s)\n", cname, gd_p$formula_used)); next }
  res <- fit_lvqw(dge_p, gd_p$design, coef = "stage_grphigh", do_ashr = TRUE)
  res[, contrast  := cname]
  res[, n_cohorts := length(both)]
  res[, n_samples := nrow(mp)]
  res[, method    := METHOD_TAG]
  nas_consec_out[[cname]] <- res
  cat(sprintf("  %s: %d cohorts {%s}, n=%d -> %d DEGs (padj<0.05)\n",
              cname, length(both), paste(sort(both), collapse = ","),
              nrow(mp), sum(res$padj < 0.05, na.rm = TRUE)))
}
nas_consec <- rbindlist(nas_consec_out, use.names = TRUE, fill = TRUE)
setcolorder(nas_consec, c("gene", "logFC", "SE", "t", "P.Value", "padj",
                          "shrunk_logFC", "lfsr", "AveExpr",
                          "contrast", "method", "n_cohorts", "n_samples"))

# ===========================================================================
# ANALYSIS 4: Fibrosis consecutive transitions (F{k}_vs_F{k-1}) — fibrosis_consecutive_dream.csv
#   DERIVED from the already-harmonized fibrosis_pairwise.csv (same LVQW binary
#   low/high fits; do NOT re-fit). Rename transition F{lo}_to_F{hi} ->
#   contrast F{hi}_vs_F{lo}; carry gene/logFC/SE/padj/shrunk_logFC/lfsr.
# ===========================================================================
cat("\n===== (4) Fibrosis consecutive transitions (DERIVED from fibrosis_pairwise.csv) =====\n")
pw_path <- file.path(SIGS, "fibrosis_pairwise.csv")
stopifnot(file.exists(pw_path))
pw <- fread(pw_path)
# F{lo}_to_F{hi}  ->  F{hi}_vs_F{lo}  (panel contract: F1_vs_F0, F2_vs_F1, ...)
pw[, contrast := gsub("F([0-9])_to_F([0-9])", "F\\2_vs_F\\1", transition)]
fib_consec <- pw[, .(gene, logFC, SE, t, P.Value, padj, shrunk_logFC, lfsr,
                     AveExpr, contrast, method = METHOD_TAG,
                     n_cohorts, n_samples)]
for (cn in unique(fib_consec$contrast)) {
  s <- fib_consec[contrast == cn]
  cat(sprintf("  %s: %d DEGs (padj<0.05) | %d (padj<0.05,|LFC|>0.25 cascade-cut)\n",
              cn, sum(s$padj < 0.05, na.rm = TRUE),
              sum(s$padj < 0.05 & abs(s$logFC) > 0.25, na.rm = TRUE)))
}

# ===========================================================================
# SANITY GATE — dream vs LVQW per contrast (BEFORE overwriting)
# ===========================================================================
cat("\n===== SANITY GATE (dream vs LVQW, before overwrite) =====\n")
nas_dream_path <- file.path(SIGS, "nas_score_dream.csv")
fib_dream_path <- file.path(SIGS, "fibrosis_stage_dream.csv")
nasc_dream_path<- file.path(SIGS, "nas_consecutive_dream.csv")
fibc_dream_path<- file.path(SIGS, "fibrosis_consecutive_dream.csv")

for (lvl in setdiff(levels(info_nas$nas_group), "0")) {
  cn <- paste0("NAS", lvl, "_vs_NAS0")
  if (cn %in% nas_sig$contrast)
    add_sanity(cn, nas_sig[contrast == cn], nas_dream_path, cn)
}
for (lvl in setdiff(levels(info_fib$fib_stage), "0")) {
  cn <- paste0("F", lvl, "_vs_F0")
  if (cn %in% fib_sig$contrast)
    add_sanity(cn, fib_sig[contrast == cn], fib_dream_path, cn)
}
for (cn in unique(nas_consec$contrast))
  add_sanity(cn, nas_consec[contrast == cn], nasc_dream_path, cn)
for (cn in unique(fib_consec$contrast))
  add_sanity(cn, fib_consec[contrast == cn], fibc_dream_path, cn)

sanity <- rbindlist(sanity_rows, use.names = TRUE, fill = TRUE)
fwrite(sanity, file.path(SIGS, "granular_staging_lvqw_sanity.csv"))
cat("Saved sanity: granular_staging_lvqw_sanity.csv\n")
print(sanity, row.names = FALSE)

# ===========================================================================
# BACKUP dream files, then OVERWRITE the dream-named files with LVQW content
# ===========================================================================
cat("\n===== BACKUP + OVERWRITE =====\n")
overwrite_map <- list(
  list(out = nas_score_dream <- nas_dream_path,  dt = nas_sig),
  list(out = fibrosis_stage_dream <- fib_dream_path, dt = fib_sig),
  list(out = nas_consecutive_dream <- nasc_dream_path, dt = nas_consec),
  list(out = fibrosis_consecutive_dream <- fibc_dream_path, dt = fib_consec))

for (it in overwrite_map) {
  out <- it$out
  bak <- sub("\\.csv$", "_premigration_backup.csv", out)
  if (file.exists(out) && !file.exists(bak)) {
    file.copy(out, bak, overwrite = FALSE)
    cat(sprintf("  backed up %s -> %s\n", basename(out), basename(bak)))
  } else if (file.exists(bak)) {
    cat(sprintf("  backup exists, leaving: %s\n", basename(bak)))
  } else {
    cat(sprintf("  [warn] no dream file to back up at %s\n", basename(out)))
  }
  fwrite(it$dt, out)
  cat(sprintf("  OVERWROTE %s  (%d rows, %d contrasts)\n",
              basename(out), nrow(it$dt), uniqueN(it$dt$contrast)))
}

cat("\n=== 14b/14d LVQW complete:", as.character(Sys.time()), "===\n")
