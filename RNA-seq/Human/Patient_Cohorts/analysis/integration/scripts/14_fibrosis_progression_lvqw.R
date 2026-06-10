#!/usr/bin/env Rscript
# ===========================================================================
# 14_fibrosis_progression_lvqw.R
# ---------------------------------------------------------------------------
# Migration of 14_fibrosis_progression_de.R to the canonical limma-voom-
# quality-weighted (LVQW) engine (de_engine_lvqw.R). STAGING ONLY — writes to
# results/disease_signatures/fibrosis_*_lvqw.csv and a sanity-gate CSV; does
# NOT overwrite any canonical file.
#
# What it computes (mirrors the original 14, engine swapped):
#   (1) FIBROSIS ORDINAL SLOPE — fibrosis_stage as a NUMERIC covariate, pooled
#       across cohorts with `dataset` + `inferred_sex` as fixed-effect nuisance.
#       Design ~ dataset + inferred_sex + fibrosis_stage ; coef = fibrosis_stage.
#   (2) PAIRWISE CONSECUTIVE TRANSITIONS — F0->F1, F1->F2, F2->F3, F3->F4, each a
#       binary low/high contrast pooled across the cohorts containing BOTH stages.
#       Design ~ dataset + inferred_sex + fib_group ; coef = fib_grouphigh.
#
# COHORT-SET PRINCIPLE (contrast-specific, NOT include_in_mega):
#   These are disease-internal progression contrasts. We pool EVERY cohort that
#   contains the relevant arms WITHIN-cohort so the `dataset` fixed effect is
#   identifiable:
#     - ORDINAL: every cohort with >=2 distinct fibrosis stages.
#     - EACH TRANSITION: every cohort containing BOTH stages of that pair.
#   build_design_guarded() drops `dataset` automatically if only one cohort
#   qualifies (no within-cohort contrast info), and drops `inferred_sex` if it
#   collapses to a single level. The original 14 hard-coded only 3 cohorts
#   (GSE130970/GSE135251/GSE213621) — a DEGENERATE subset; meta_matched.rds in
#   fact carries fibrosis_stage for 7 cohorts. We use the full identifiable set
#   here and flag the difference in the sanity gate.
#
# Input loading mirrors 05h (merged_dge.rds + meta_matched.rds; subset only,
# NO re-TMM / NO re-filterByExpr — the merged DGE is already normalized/filtered).
# ===========================================================================
suppressMessages({
  library(edgeR)
  library(limma)
  library(ashr)
  library(data.table)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
INT  <- file.path(BASE, "analysis/integration")
SDIR <- file.path(INT, "scripts")
RDIR <- file.path(INT, "results/integration")
ODIR <- file.path(INT, "results/disease_signatures")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(SDIR, "de_engine_lvqw.R"))

# ---------------------------------------------------------------------------
# 0) Input loading (05h template, lines 19-33) — subset only, no re-norm/filter
# ---------------------------------------------------------------------------
dge_all  <- readRDS(file.path(RDIR, "merged_dge.rds"))
meta_new <- as.data.table(readRDS(file.path(RDIR, "meta_matched.rds")))
qc       <- fread(file.path(INT, "qc/sample_qc_report.csv"))

# original 14 QC gate: pass_technical
keep_ids <- qc[pass_technical == TRUE, sample_id]
meta_new <- meta_new[sample_id %in% keep_ids & sample_id %in% colnames(dge_all)]

# fibrosis-annotated universe
meta_fib <- meta_new[!is.na(fibrosis_stage)]
cat(sprintf("[in] %d genes x %d samples in merged_dge; %d fibrosis-annotated (pass_technical)\n",
            nrow(dge_all), ncol(dge_all), nrow(meta_fib)))
cat("[in] per cohort x stage:\n")
print(meta_fib[, .N, by = .(dataset, fibrosis_stage)][order(dataset, fibrosis_stage)])

# helper: subset the (already normalized) merged DGE to a sample set, attach
# sample covariates from meta. NO calcNormFactors / NO filterByExpr re-run.
subset_dge <- function(sample_ids, cov_dt) {
  idx <- colnames(dge_all) %in% sample_ids
  d   <- dge_all[, idx]
  cov <- cov_dt[match(colnames(d), cov_dt$sample_id)]
  d$samples$dataset      <- factor(cov$dataset)
  d$samples$inferred_sex <- factor(cov$inferred_sex)
  d$samples$fibrosis_stage <- cov$fibrosis_stage
  d
}

# ===========================================================================
# (1) FIBROSIS ORDINAL SLOPE — pooled, dataset + sex fixed, stage NUMERIC
# ===========================================================================
cat("\n===== (1) FIBROSIS ORDINAL SLOPE (LVQW) =====\n")

# cohorts with >=2 distinct fibrosis stages (within-cohort slope info)
ord_cohorts <- meta_fib[, .(n_stages = uniqueN(fibrosis_stage)), by = dataset][n_stages >= 2, dataset]
m_ord <- meta_fib[dataset %in% ord_cohorts]
cat(sprintf("Ordinal cohort set (>=2 stages): {%s}\n", paste(sort(unique(m_ord$dataset)), collapse = ", ")))
cat(sprintf("Ordinal samples: %d across %d cohorts\n", nrow(m_ord), uniqueN(m_ord$dataset)))

dge_ord <- subset_dge(m_ord$sample_id, m_ord)
info_ord <- data.frame(
  dataset        = factor(dge_ord$samples$dataset),
  inferred_sex   = factor(dge_ord$samples$inferred_sex),
  fibrosis_stage = as.numeric(dge_ord$samples$fibrosis_stage)
)
rownames(info_ord) <- colnames(dge_ord)

# build_design_guarded: fibrosis_stage is NUMERIC -> always survives; dataset /
# inferred_sex are guarded against single-level collapse.
gd_ord <- build_design_guarded(info_ord, c("dataset", "inferred_sex", "fibrosis_stage"))
cat(sprintf("Ordinal design: %s\n", gd_ord$formula_used))
if (length(gd_ord$dropped)) cat(sprintf("  dropped: %s\n", paste(gd_ord$dropped, collapse = ", ")))
stopifnot("fibrosis_stage" %in% colnames(gd_ord$design))

ord_res <- fit_lvqw(dge_ord, gd_ord$design, coef = "fibrosis_stage", do_ashr = TRUE)
ord_res[, method := "limma_voom_qw_C2"]
fwrite(ord_res, file.path(ODIR, "fibrosis_ordinal_lvqw.csv"))
n_ord_sig <- sum(ord_res$padj < 0.05, na.rm = TRUE)
cat(sprintf("Ordinal DEGs (padj<0.05): %d (%d up / %d down)\n",
            n_ord_sig,
            sum(ord_res$padj < 0.05 & ord_res$logFC > 0, na.rm = TRUE),
            sum(ord_res$padj < 0.05 & ord_res$logFC < 0, na.rm = TRUE)))
cat("Saved: fibrosis_ordinal_lvqw.csv\n")

# ===========================================================================
# (2) PAIRWISE CONSECUTIVE TRANSITIONS — binary low/high, dataset + sex fixed
# ===========================================================================
cat("\n===== (2) PAIRWISE CONSECUTIVE TRANSITIONS (LVQW) =====\n")

stage_pairs <- list(c(0L, 1L), c(1L, 2L), c(2L, 3L), c(3L, 4L))
pairwise_out <- list()

for (pr in stage_pairs) {
  f_low  <- pr[1]; f_high <- pr[2]
  tname  <- sprintf("F%d_to_F%d", f_low, f_high)

  # cohorts containing BOTH stages of this pair (within-cohort contrast info)
  both <- meta_fib[fibrosis_stage %in% c(f_low, f_high),
                   .(n_stage = uniqueN(fibrosis_stage)), by = dataset][n_stage == 2L, dataset]
  m_pair <- meta_fib[dataset %in% both & fibrosis_stage %in% c(f_low, f_high)]

  if (nrow(m_pair) < 6L || uniqueN(m_pair$fibrosis_stage) < 2L) {
    cat(sprintf("%s: SKIPPED (n=%d, cohorts={%s})\n", tname, nrow(m_pair), paste(both, collapse = ",")))
    next
  }

  cat(sprintf("%s: %d cohorts {%s}, %d samples (low=%d, high=%d)\n",
              tname, length(both), paste(sort(both), collapse = ","), nrow(m_pair),
              sum(m_pair$fibrosis_stage == f_low), sum(m_pair$fibrosis_stage == f_high)))

  dge_p <- subset_dge(m_pair$sample_id, m_pair)
  fib_group <- factor(ifelse(dge_p$samples$fibrosis_stage == f_high, "high", "low"),
                      levels = c("low", "high"))
  info_p <- data.frame(
    dataset      = factor(dge_p$samples$dataset),
    inferred_sex = factor(dge_p$samples$inferred_sex),
    fib_group    = fib_group
  )
  rownames(info_p) <- colnames(dge_p)

  # build_design_guarded drops `dataset` if only 1 cohort has both stages, and
  # drops `inferred_sex` if it collapses to a single level.
  gd_p <- build_design_guarded(info_p, c("dataset", "inferred_sex", "fib_group"))
  coef_p <- "fib_grouphigh"
  if (length(gd_p$dropped)) cat(sprintf("  design: %s ; dropped: %s\n", gd_p$formula_used, paste(gd_p$dropped, collapse = ", ")))
  else cat(sprintf("  design: %s\n", gd_p$formula_used))
  stopifnot(coef_p %in% colnames(gd_p$design))

  res_p <- fit_lvqw(dge_p, gd_p$design, coef = coef_p, do_ashr = TRUE)
  res_p[, transition := tname]
  res_p[, n_cohorts  := length(both)]
  res_p[, n_samples  := nrow(m_pair)]
  res_p[, method     := "limma_voom_qw_C2"]
  pairwise_out[[tname]] <- res_p

  cat(sprintf("  DEGs (padj<0.05, no LFC filter): %d\n", sum(res_p$padj < 0.05, na.rm = TRUE)))
}

all_pairwise <- rbindlist(pairwise_out, use.names = TRUE, fill = TRUE)
# column order: tidy schema first, then transition/provenance
setcolorder(all_pairwise,
  c("gene", "logFC", "SE", "t", "P.Value", "padj", "shrunk_logFC", "lfsr", "AveExpr",
    "method", "transition", "n_cohorts", "n_samples"))
fwrite(all_pairwise, file.path(ODIR, "fibrosis_pairwise_lvqw.csv"))
cat("Saved: fibrosis_pairwise_lvqw.csv\n")

# ===========================================================================
# (3) BACKUP original dream/limma outputs (copy only; do NOT overwrite canon)
# ===========================================================================
cat("\n===== BACKUP pre-migration outputs =====\n")
backup_map <- c(
  "fibrosis_dream.csv"            = "fibrosis_dream_premigration_backup.csv",
  "fibrosis_pairwise.csv"         = "fibrosis_pairwise_premigration_backup.csv",
  "fibrosis_slopes_meta.csv"      = "fibrosis_slopes_meta_premigration_backup.csv",
  "fibrosis_ordinal_per_study.csv"= "fibrosis_ordinal_per_study_premigration_backup.csv"
)
for (src in names(backup_map)) {
  s <- file.path(ODIR, src); d <- file.path(ODIR, backup_map[[src]])
  if (file.exists(s)) {
    file.copy(s, d, overwrite = TRUE)
    cat(sprintf("  backed up %s -> %s\n", src, backup_map[[src]]))
  } else {
    cat(sprintf("  [skip] %s not found\n", src))
  }
}

# ===========================================================================
# (4) SANITY GATE vs original dream/limma outputs
#     - ORDINAL: new LVQW pooled slope vs original dream pooled slope
#       (fibrosis_dream.csv = dream ~ fibrosis_stage + sex + (1|dataset)).
#       NOTE cohort-set difference: original dream used 3 cohorts; LVQW uses the
#       full identifiable set (7).
#     - PAIRWISE: per transition, new LVQW pooled vs the original PER-STUDY
#       limma pairwise (fibrosis_pairwise.csv) collapsed to a per-gene
#       inverse-variance fixed-effect summary across its (<=3) cohorts. This is
#       the closest apples-to-apples summary; the original never pooled.
# ===========================================================================
cat("\n===== SANITY GATE =====\n")
gate_rows <- list()

spearman <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
}

# ---- ordinal gate ----
old_dream_path <- file.path(ODIR, "fibrosis_dream_premigration_backup.csv")
if (file.exists(old_dream_path)) {
  od <- fread(old_dream_path)  # logFC, AveExpr, t, P.Value, adj.P.Val, z.std, gene
  setnames(od, "adj.P.Val", "old_padj", skip_absent = TRUE)
  mg <- merge(ord_res[, .(gene, new_logFC = logFC, new_padj = padj)],
              od[, .(gene, old_logFC = logFC, old_padj)], by = "gene")
  shared <- nrow(mg)
  dir_conc <- mean(sign(mg$new_logFC) == sign(mg$old_logFC), na.rm = TRUE)
  rho <- spearman(mg$new_logFC, mg$old_logFC)
  gate_rows[[length(gate_rows) + 1]] <- data.table(
    analysis     = "ordinal_slope",
    new_n_deg    = sum(ord_res$padj < 0.05, na.rm = TRUE),
    old_n_deg    = sum(od$old_padj < 0.05, na.rm = TRUE),
    n_shared     = shared,
    dir_concord  = round(dir_conc, 4),
    logfc_spearman = round(rho, 4),
    note = "LVQW pooled (7 cohorts) vs original dream pooled (3 cohorts GSE130970/135251/213621) — cohort set differs"
  )
  cat(sprintf("[ordinal] new DEG=%d old(dream) DEG=%d shared=%d dir=%.3f rho=%.3f\n",
              sum(ord_res$padj < 0.05, na.rm = TRUE), sum(od$old_padj < 0.05, na.rm = TRUE),
              shared, dir_conc, rho))
}

# ---- pairwise gate (per transition) ----
old_pw_path <- file.path(ODIR, "fibrosis_pairwise_premigration_backup.csv")
if (file.exists(old_pw_path)) {
  opw <- fread(old_pw_path)  # logFC, AveExpr, t, P.Value, adj.P.Val, B, gene, dataset, contrast
  # original contrast names are F{high}_vs_F{low}; map to F{low}_to_F{high}
  opw[, transition := gsub("F([0-9])_vs_F([0-9])", "F\\2_to_F\\1", contrast)]
  opw[, old_se := abs(logFC / t)]
  opw[!is.finite(old_se) | old_se == 0, old_se := NA_real_]

  for (tname in names(pairwise_out)) {
    new_t <- pairwise_out[[tname]]
    old_t <- opw[transition == tname]
    if (nrow(old_t) == 0) {
      gate_rows[[length(gate_rows) + 1]] <- data.table(
        analysis = tname,
        new_n_deg = sum(new_t$padj < 0.05, na.rm = TRUE),
        old_n_deg = NA_integer_, n_shared = 0L,
        dir_concord = NA_real_, logfc_spearman = NA_real_,
        note = "no original per-study rows for this transition")
      next
    }
    # collapse original per-study to per-gene fixed-effect IVW summary
    old_sum <- old_t[is.finite(old_se), {
      w <- 1 / (old_se^2)
      lf <- sum(w * logFC) / sum(w)
      list(old_logFC = lf, old_padj_min = min(p.adjust(P.Value, "BH")), n_studies = .N)
    }, by = gene]
    # original per-study significance: a gene is "old-sig" if sig in any contributing study
    old_sig <- old_t[, .(any_sig = any(adj.P.Val < 0.05, na.rm = TRUE)), by = gene]
    n_old_sig <- sum(old_sig$any_sig, na.rm = TRUE)

    mg <- merge(new_t[, .(gene, new_logFC = logFC, new_padj = padj)],
                old_sum[, .(gene, old_logFC)], by = "gene")
    shared <- nrow(mg)
    dir_conc <- mean(sign(mg$new_logFC) == sign(mg$old_logFC), na.rm = TRUE)
    rho <- spearman(mg$new_logFC, mg$old_logFC)
    gate_rows[[length(gate_rows) + 1]] <- data.table(
      analysis     = tname,
      new_n_deg    = sum(new_t$padj < 0.05, na.rm = TRUE),
      old_n_deg    = n_old_sig,
      n_shared     = shared,
      dir_concord  = round(dir_conc, 4),
      logfc_spearman = round(rho, 4),
      note = sprintf("LVQW pooled (%d cohorts) vs original per-study limma (<=3 cohorts, IVW-collapsed; old DEG=any-study-sig)",
                     new_t$n_cohorts[1])
    )
    cat(sprintf("[%s] new DEG=%d old(any-study) DEG=%d shared=%d dir=%.3f rho=%.3f\n",
                tname, sum(new_t$padj < 0.05, na.rm = TRUE), n_old_sig, shared, dir_conc, rho))
  }
}

gate <- rbindlist(gate_rows, use.names = TRUE, fill = TRUE)
fwrite(gate, file.path(ODIR, "fibrosis_lvqw_sanity.csv"))
cat("\nSaved: fibrosis_lvqw_sanity.csv\n")

cat("\n=== Script 14_lvqw complete (STAGING ONLY — canonical NOT overwritten) ===\n")
