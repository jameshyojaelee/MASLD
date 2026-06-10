#!/usr/bin/env Rscript
# ===========================================================================
# 14b_stage_vs_healthy_lvqw.R
# ---------------------------------------------------------------------------
# Migration of 14b_stage_vs_healthy_dream.R to the canonical limma-voom-
# quality-weighted (LVQW) engine (de_engine_lvqw.R). This per-stage-VS-CONTROL
# bulk DEG analysis was MISSED in the 2026-06-08 WS2 dream->LVQW harmonization;
# a user spot-check caught it. This script closes that gap.
#
# What it computes (mirrors 14b dream, engine swapped + cohort-expanded):
#   ANALYSIS 1 — PER-FIBROSIS-STAGE vs strictly healthy controls:
#       F1_vs_Ctrl, F2_vs_Ctrl, F3_vs_Ctrl, F4_vs_Ctrl
#   ANALYSIS 2 — PER-NAS-SCORE vs strictly healthy controls:
#       NAS1_vs_Ctrl ... NAS7_vs_Ctrl  (NAS>=7 collapsed into NAS7)
#
# STAGING DEFINITIONS — PRESERVED VERBATIM from 14b_stage_vs_healthy_dream.R:
#   * Control ("Ctrl") = condition == "Control" (strictly healthy; NOT all
#     F0/NAS0 samples — F0/NAS0 disease patients are EXCLUDED from the contrast).
#   * Fibrosis disease groups = samples with fibrosis_stage in 1:4, condition!=Control.
#   * NAS disease groups = samples with nas_score >= 1, condition != Control;
#     NAS score >=7 collapsed into the NAS7 bin.
#   Samples failing pass_technical QC are excluded.
#
# What CHANGED vs 14b dream
# -------------------------
#   * ENGINE: dream() (mixed model, ~ grp + inferred_sex + (1|dataset), all
#     stage levels in ONE multi-level factor fit) -> LVQW (voomWithQualityWeights
#     -> lmFit -> eBayes -> ashr) via fit_lvqw(). Each stage bin is a SEPARATE
#     binary (bin vs Ctrl) fit, matching the WS2 LVQW transition pattern
#     (14_fibrosis_progression_lvqw.R).
#   * DESIGN per contrast: ~ dataset + inferred_sex + group  (group = this stage
#     bin vs control, binary; coef = grouphigh). dataset is a FIXED effect
#     (replaces dream's (1|dataset) random intercept), for method consistency
#     with the C2 canonical.
#   * COHORT-EXPANDED: dream 14b hard-coded FIB_DATASETS / NAS_DATASETS lists
#     (a degenerate subset). Here, EACH contrast pools EVERY cohort that has
#     BOTH >=2 samples in the stage bin AND >=2 controls (within-cohort contrast
#     info, so the `dataset` fixed effect is identifiable). build_design_guarded()
#     drops `dataset` automatically if only ONE cohort survives, and drops
#     `inferred_sex` if it collapses to a single level.
#   * INPUT: subset the canonical merged_dge.rds (already TMM-normalized +
#     filterByExpr'd) -> NO re-TMM, NO re-filterByExpr (05h / WS2 pattern). Stage
#     fields joined from meta_matched.rds by sample_id.
#
# OUTPUT (figure-compatible) — OVERWRITES the two canonical dream files after
# backing them up to *_premigration_backup.csv (copy):
#   results/disease_signatures/nas_stage_vs_ctrl_dream.csv
#   results/disease_signatures/fibrosis_stage_vs_ctrl_dream.csv
#   results/disease_signatures/nas_stage_vs_ctrl_sample_sizes.csv      (regenerated)
#   results/disease_signatures/fibrosis_stage_vs_ctrl_sample_sizes.csv (regenerated)
# The four fig2 panels read columns gene, logFC, padj, contrast — ALL preserved;
# this script ADDS SE, shrunk_logFC, lfsr, method="limma_voom_qw".
# A sanity gate vs the pre-migration dream files is written to:
#   results/disease_signatures/stage_vs_ctrl_lvqw_sanity.csv
#
# Usage (SLURM, never login node):
#   env -u SLURM_JOB_ID sbatch --partition=cpu --qos=interactive --mem=64G \
#     --cpus-per-task=4 --time=12:00:00 --job-name=limma \
#     --output=<scripts>/logs/14b_lvqw_%j.log \
#     --wrap="<rnaseq>/bin/Rscript 14b_stage_vs_healthy_lvqw.R"
# ===========================================================================
set.seed(42)
suppressMessages({
  library(edgeR); library(limma); library(ashr); library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/integration")
SIG  <- file.path(INT, "results/disease_signatures")
SCR  <- file.path(INT, "scripts")
stopifnot(dir.exists(RDIR), dir.exists(SIG))

source(file.path(SCR, "de_engine_lvqw.R"))

METHOD_TAG <- "limma_voom_qw"   # value the figure-data CSV carries in `method`

cat("=== 14b LVQW: per-stage vs strictly-healthy-control DEGs (cohort-expanded) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ---------------------------------------------------------------------------
# 1. Input loading — subset the canonical merged_dge.rds (NO re-TMM / re-filter).
#    Stage fields joined from meta_matched.rds; QC gate = pass_technical.
# ---------------------------------------------------------------------------
dge_all <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta    <- as.data.table(readRDS(file.path(RDIR, "meta_matched.rds")))
qc      <- fread(file.path(INT, "qc/sample_qc_report.csv"))
keep_ids <- qc[pass_technical == TRUE, sample_id]
meta <- meta[sample_id %in% keep_ids & sample_id %in% colnames(dge_all)]
cat(sprintf("[in] %d genes x %d samples in merged_dge; %d pass_technical & in DGE\n",
            nrow(dge_all), ncol(dge_all), nrow(meta)))

# Symbol lookup (drop-in parity with WS2 LVQW stage files; harmless extra col,
# panels ignore it).
sym_gm <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
sym_gm[, eb := sub("[.][0-9]+$", "", gene_id)]

# helper: subset the (already normalized) merged DGE to a sample set and attach
# the per-sample covariates needed for the design. NO calcNormFactors / filterByExpr.
subset_dge <- function(sample_ids, cov_dt) {
  idx <- colnames(dge_all) %in% sample_ids
  d   <- dge_all[, idx]
  cov <- cov_dt[match(colnames(d), cov_dt$sample_id)]
  stopifnot(identical(cov$sample_id, colnames(d)))
  d$samples$dataset      <- factor(as.character(cov$dataset))
  d$samples$inferred_sex <- factor(as.character(cov$inferred_sex))
  d$samples$grp          <- factor(cov$grp, levels = c("Ctrl", "high"))
  d
}

# ---------------------------------------------------------------------------
# 2. Generic per-bin LVQW runner.
#    For one stage bin: pool cohorts with >=2 in the bin AND >=2 controls,
#    fit ~ dataset + inferred_sex + group (group = bin vs Ctrl, binary).
#    Returns the fixed two-tier LVQW schema + contrast/method/provenance cols,
#    or NULL if no cohort has both arms (>=2 each).
# ---------------------------------------------------------------------------
MIN_PER_ARM <- 2L

run_bin_vs_ctrl <- function(meta_ctrl, meta_bin, contrast_name) {
  # cohorts that have BOTH >=MIN_PER_ARM bin samples AND >=MIN_PER_ARM controls
  n_bin  <- meta_bin [, .N, by = dataset]
  n_ctrl <- meta_ctrl[, .N, by = dataset]
  both   <- merge(n_bin, n_ctrl, by = "dataset", suffixes = c("_bin", "_ctrl"))
  both   <- both[N_bin >= MIN_PER_ARM & N_ctrl >= MIN_PER_ARM]
  cohorts <- sort(both$dataset)

  if (length(cohorts) == 0L) {
    cat(sprintf("  %s: SKIPPED (no cohort with >=%d in bin AND >=%d controls)\n",
                contrast_name, MIN_PER_ARM, MIN_PER_ARM))
    return(NULL)
  }

  mb <- meta_bin [dataset %in% cohorts];  mb[, grp := "high"]
  mc <- meta_ctrl[dataset %in% cohorts];  mc[, grp := "Ctrl"]
  m_pair <- rbind(mc, mb, fill = TRUE)

  n_hi <- nrow(mb); n_ct <- nrow(mc)
  cat(sprintf("  %s: %d cohorts {%s} | bin=%d ctrl=%d (n=%d)\n",
              contrast_name, length(cohorts), paste(cohorts, collapse = ","),
              n_hi, n_ct, nrow(m_pair)))

  dge_p <- subset_dge(m_pair$sample_id, m_pair)
  info <- data.frame(
    dataset      = factor(dge_p$samples$dataset),
    inferred_sex = factor(dge_p$samples$inferred_sex),
    group        = factor(dge_p$samples$grp, levels = c("Ctrl", "high")),
    row.names    = colnames(dge_p),
    check.names  = FALSE)

  # build_design_guarded drops `dataset` if only 1 cohort survives and
  # `inferred_sex` if it collapses to a single level. `group` is the term of
  # interest and is placed LAST so it is never dropped; verify its coef.
  gd <- build_design_guarded(info, rhs_terms = c("dataset", "inferred_sex", "group"))
  coef_name <- "grouphigh"
  if (length(gd$dropped))
    cat(sprintf("    design: %s ; dropped: %s\n", gd$formula_used,
                paste(gd$dropped, collapse = ", ")))
  else
    cat(sprintf("    design: %s\n", gd$formula_used))
  stopifnot(coef_name %in% colnames(gd$design))

  res <- fit_lvqw(dge_p, gd$design, coef = coef_name, do_ashr = TRUE)
  res[, contrast  := contrast_name]
  res[, method    := METHOD_TAG]
  res[, n_cohorts := length(cohorts)]
  res[, n_bin     := n_hi]
  res[, n_ctrl    := n_ct]
  res[, symbol    := sym_gm[match(sub("[.][0-9]+$", "", gene), eb), gene_name]]

  n_sig <- sum(res$padj < 0.05 & abs(res$logFC) > 0.5, na.rm = TRUE)
  cat(sprintf("    DEGs (padj<0.05, |LFC|>0.5): %d\n", n_sig))
  res
}

# Final column order: figure-read cols first (gene, logFC, padj, contrast all
# present), ashr + provenance after. Panels only require gene/logFC/padj/contrast.
SCHEMA <- c("gene", "symbol", "logFC", "SE", "t", "P.Value", "padj",
            "shrunk_logFC", "lfsr", "AveExpr", "contrast", "method",
            "n_cohorts", "n_bin", "n_ctrl")

# ===========================================================================
# ANALYSIS 1: Per-Fibrosis-Stage vs Healthy Controls (F1..F4)
# ===========================================================================
cat("=== ANALYSIS 1: Fibrosis stages vs healthy controls ===\n")
meta_ctrl_fib <- meta[condition == "Control"]
meta_fib_dis  <- meta[!is.na(fibrosis_stage) &
                      as.integer(fibrosis_stage) %in% 1:4 &
                      condition != "Control"]
meta_fib_dis[, fib_bin := as.integer(fibrosis_stage)]

fib_results <- list()
for (st in 1:4) {
  cn  <- sprintf("F%d_vs_Ctrl", st)
  mb  <- meta_fib_dis[fib_bin == st]
  res <- run_bin_vs_ctrl(meta_ctrl_fib, mb, cn)
  if (!is.null(res)) fib_results[[cn]] <- res
}
fib_out <- rbindlist(fib_results, use.names = TRUE, fill = TRUE)
setcolorder(fib_out, SCHEMA)

# Sample-sizes table (figure subtitle reads grp=="Ctrl" -> N). Ctrl count is the
# union of controls across all cohorts that contributed to >=1 contrast; per-bin
# counts are the within-contrast pooled sizes.
fib_sizes <- rbindlist(c(
  list(data.table(grp = "Ctrl",
                  N = uniqueN(meta_ctrl_fib[dataset %in% unique(unlist(
                    lapply(names(fib_results), function(cn) {
                      st <- as.integer(sub("F(\\d)_vs_Ctrl", "\\1", cn))
                      mb <- meta_fib_dis[fib_bin == st]
                      nb <- mb[, .N, by = dataset]; nc <- meta_ctrl_fib[, .N, by = dataset]
                      both <- merge(nb, nc, by = "dataset")
                      both[N.x >= MIN_PER_ARM & N.y >= MIN_PER_ARM, dataset]
                    })))]$sample_id))),
  lapply(names(fib_results), function(cn) {
    data.table(grp = sub("_vs_Ctrl", "", cn), N = fib_results[[cn]]$n_bin[1])
  })
))
setorder(fib_sizes, grp)
fwrite(fib_sizes, file.path(SIG, "fibrosis_stage_vs_ctrl_sample_sizes.csv"))
cat("\nFibrosis sample sizes:\n"); print(fib_sizes)

# ===========================================================================
# ANALYSIS 2: Per-NAS-Score vs Healthy Controls (NAS1..NAS7; >=7 -> NAS7)
# ===========================================================================
cat("\n=== ANALYSIS 2: NAS stages vs healthy controls ===\n")
meta_ctrl_nas <- meta[condition == "Control"]
meta_nas_dis  <- meta[!is.na(nas_score) & nas_score >= 1 & condition != "Control"]
meta_nas_dis[, nas_bin := fifelse(nas_score >= 7, 7L, as.integer(nas_score))]

nas_results <- list()
for (st in 1:7) {
  cn  <- sprintf("NAS%d_vs_Ctrl", st)
  mb  <- meta_nas_dis[nas_bin == st]
  res <- run_bin_vs_ctrl(meta_ctrl_nas, mb, cn)
  if (!is.null(res)) nas_results[[cn]] <- res
}
nas_out <- rbindlist(nas_results, use.names = TRUE, fill = TRUE)
setcolorder(nas_out, SCHEMA)

nas_sizes <- rbindlist(c(
  list(data.table(grp = "Ctrl",
                  N = uniqueN(meta_ctrl_nas[dataset %in% unique(unlist(
                    lapply(names(nas_results), function(cn) {
                      st <- as.integer(sub("NAS(\\d)_vs_Ctrl", "\\1", cn))
                      mb <- meta_nas_dis[nas_bin == st]
                      nb <- mb[, .N, by = dataset]; nc <- meta_ctrl_nas[, .N, by = dataset]
                      both <- merge(nb, nc, by = "dataset")
                      both[N.x >= MIN_PER_ARM & N.y >= MIN_PER_ARM, dataset]
                    })))]$sample_id))),
  lapply(names(nas_results), function(cn) {
    data.table(grp = sub("_vs_Ctrl", "", cn), N = nas_results[[cn]]$n_bin[1])
  })
))
setorder(nas_sizes, grp)
fwrite(nas_sizes, file.path(SIG, "nas_stage_vs_ctrl_sample_sizes.csv"))
cat("\nNAS sample sizes:\n"); print(nas_sizes)

# ===========================================================================
# 3. BACKUP the canonical dream files (copy), THEN overwrite with LVQW content.
# ===========================================================================
cat("\n=== BACKUP + OVERWRITE canonical stage-vs-ctrl files ===\n")

backup_and_write <- function(new_dt, fname) {
  src <- file.path(SIG, fname)
  if (file.exists(src)) {
    dst <- file.path(SIG, sub("\\.csv$", "_premigration_backup.csv", fname))
    ok  <- file.copy(src, dst, overwrite = FALSE)  # do not clobber a prior backup
    if (file.exists(dst))
      cat(sprintf("[backup] %s -> %s (%s)\n", basename(src), basename(dst),
                  if (ok) "copied" else "already existed; kept"))
  } else {
    cat(sprintf("[backup] WARN: %s missing; no backup\n", fname))
  }
  fwrite(new_dt, src)
  cat(sprintf("[write] overwrote %s (%d rows, LVQW)\n", fname, nrow(new_dt)))
}
backup_and_write(fib_out, "fibrosis_stage_vs_ctrl_dream.csv")
backup_and_write(nas_out, "nas_stage_vs_ctrl_dream.csv")

# ===========================================================================
# 4. SANITY GATE vs the pre-migration dream files.
#    Per contrast: dream vs LVQW DEG count (padj<0.05), direction concordance,
#    logFC Spearman rho on shared genes. Also the figure DEG def (|LFC|>0.5).
# ===========================================================================
cat("\n=== SANITY GATE (dream vs LVQW) ===\n")
spearman <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
}
sanity_one <- function(new_dt, dream_bk_path) {
  dr <- fread(dream_bk_path)
  drp <- if ("padj" %in% names(dr)) "padj" else "adj.P.Val"
  setnames(dr, drp, "dr_padj")
  rows <- list()
  for (cn in unique(new_dt$contrast)) {
    lv <- new_dt[contrast == cn]
    dd <- dr[contrast == cn]
    mm <- merge(lv[, .(gene, lv_lfc = logFC, lv_padj = padj)],
                dd[, .(gene, dr_lfc = logFC, dr_padj)], by = "gene")
    rows[[cn]] <- data.table(
      contrast              = cn,
      nDEG_dream_p05        = sum(dd$dr_padj < 0.05, na.rm = TRUE),
      nDEG_lvqw_p05         = sum(lv$padj   < 0.05, na.rm = TRUE),
      nDEG_dream_p05_lfc0.5 = sum(dd$dr_padj < 0.05 & abs(dd$logFC) > 0.5, na.rm = TRUE),
      nDEG_lvqw_p05_lfc0.5  = sum(lv$padj   < 0.05 & abs(lv$logFC) > 0.5, na.rm = TRUE),
      n_shared              = nrow(mm),
      direction_concordance = round(mean(sign(mm$lv_lfc) == sign(mm$dr_lfc), na.rm = TRUE), 4),
      logFC_spearman        = round(spearman(mm$lv_lfc, mm$dr_lfc), 4),
      n_cohorts_lvqw        = lv$n_cohorts[1],
      n_bin_lvqw            = lv$n_bin[1],
      n_ctrl_lvqw           = lv$n_ctrl[1])
    cat(sprintf("[%s] dream=%d lvqw=%d (p05) | dream=%d lvqw=%d (p05,|LFC|>.5) | shared=%d dir=%.3f rho=%.3f | %d cohorts\n",
                cn, rows[[cn]]$nDEG_dream_p05, rows[[cn]]$nDEG_lvqw_p05,
                rows[[cn]]$nDEG_dream_p05_lfc0.5, rows[[cn]]$nDEG_lvqw_p05_lfc0.5,
                nrow(mm), rows[[cn]]$direction_concordance, rows[[cn]]$logFC_spearman,
                rows[[cn]]$n_cohorts_lvqw))
  }
  rbindlist(rows)
}
san_fib <- sanity_one(fib_out, file.path(SIG, "fibrosis_stage_vs_ctrl_dream_premigration_backup.csv"))
san_fib[, analysis := "fibrosis"]
san_nas <- sanity_one(nas_out, file.path(SIG, "nas_stage_vs_ctrl_dream_premigration_backup.csv"))
san_nas[, analysis := "nas"]
sanity <- rbindlist(list(san_fib, san_nas), use.names = TRUE)
setcolorder(sanity, c("analysis", "contrast"))
fwrite(sanity, file.path(SIG, "stage_vs_ctrl_lvqw_sanity.csv"))
cat("\nSaved sanity table:", file.path(SIG, "stage_vs_ctrl_lvqw_sanity.csv"), "\n")
print(sanity, row.names = FALSE)

cat("\n=== 14b LVQW complete:", as.character(Sys.time()), "===\n")
