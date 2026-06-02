#!/usr/bin/env Rscript
# ============================================================================
# 348_stage_ccc_bulk_concordance.R
#
# Cross-modal concordance check: for each LR pair flagged as stage-progressive
# by Script 346, test whether ligand and/or receptor bulk DEG sign agrees
# with the scRNA LIANA effect direction.
#
# Stage-matched anchors (produced by 05c_dream_stage_contrasts.R):
#   sc term `disease_stage_coarseSteatosis`        -> dream_results_stage_steatosis.csv
#   sc term `disease_stage_coarseSteatohepatitis`  -> dream_results_stage_sh.csv
#   sc term `disease_stage_coarseCirrhosis`        -> dream_results_stage_cirrhosis.csv
#   continuous (macrophage pseudotime, F-stage)    -> dream_results_ashr.csv
#       (overall MASLD-vs-Healthy = appropriate match for continuous severity axis)
#
# If a stage-specific bulk anchor is missing, we silently fall back to the
# all-MASLD anchor for that contrast and tag axis name with `_fallback`.
#
# Outputs:
#   lr_bulk_concordance_steatosis.tsv
#   lr_bulk_concordance_sh.tsv
#   lr_bulk_concordance_cirrhosis.tsv
#   lr_bulk_concordance_continuous.tsv
#   lr_bulk_concordance_fstage.tsv          (F-stage axis vs all-MASLD + SH anchors)
#   lr_bulk_concordance.tsv                 (combined; one row per axis-contrast x LR pair)
#   Columns: ct_pair, lr_pair, ligand, receptor, lig_bulk_lfc, lig_bulk_padj,
#            rec_bulk_lfc, rec_bulk_padj, lig_concordant, rec_concordant,
#            both_concordant, axis, sc_estimate
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")

COARSE_TSV    <- file.path(OUT_DIR, "stage_lr_lmm_coarse.tsv")
CONTINUOUS_TSV<- file.path(OUT_DIR, "stage_lr_lmm_continuous.tsv")
FSTAGE_TSV    <- file.path(OUT_DIR, "stage_lr_lmm_fstage.tsv")
INT_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
BULK_ALL       <- file.path(INT_DIR, "dream_results_ashr.csv")
BULK_STEATOSIS <- file.path(INT_DIR, "dream_results_stage_steatosis.csv")
BULK_SH        <- file.path(INT_DIR, "dream_results_stage_sh.csv")
BULK_CIRRHOSIS <- file.path(INT_DIR, "dream_results_stage_cirrhosis.csv")

stopifnot(file.exists(COARSE_TSV))
stopifnot(file.exists(BULK_ALL))

# Protocol-contamination exclusion (GSE136103, Liver_Atlas) is inherited from
# Script 346 — donors flagged exclude_stage_analysis == TRUE are dropped
# before LMM fitting, so stage_lr_lmm_*.tsv consumed below is already clean.
# No donor-level filter is applied here because 348 only reads LMM output
# TSVs + bulk DEG CSVs (no per-donor data). If 346 is re-run with a stale
# metadata TSV missing the flag, the backwards-compat fallback there treats
# the flag as FALSE for all donors.
cat("[info] protocol-contamination filtering is inherited from 346 outputs.\n")

# ---- Bulk loader: returns lookup table keyed by HGNC symbol ----------------
load_bulk <- function(path) {
  if (!file.exists(path)) return(NULL)
  bulk <- fread(path)
  # Prefer `symbol` over Ensembl-versioned `gene` because LIANA uses HGNC.
  gene_col <- intersect(c("symbol", "gene_symbol", "hgnc_symbol",
                          "gene", "Gene", "ID"), names(bulk))[1]
  lfc_col  <- intersect(c("logFC", "log2FoldChange", "lfc"), names(bulk))[1]
  padj_col <- intersect(c("padj", "adj.P.Val", "BH", "FDR"), names(bulk))[1]
  if (is.null(gene_col) || is.null(lfc_col) || is.null(padj_col)) return(NULL)
  out <- bulk[, .(gene = get(gene_col),
                  bulk_lfc = get(lfc_col),
                  bulk_padj = get(padj_col))]
  out <- out[!is.na(gene) & gene != ""]
  # Collapse duplicate symbols by smallest padj (lncRNAs / readthroughs)
  out <- out[order(bulk_padj)][, .SD[1], by = gene]
  out
}

bulk_all       <- load_bulk(BULK_ALL)
bulk_steatosis <- load_bulk(BULK_STEATOSIS)
bulk_sh        <- load_bulk(BULK_SH)
bulk_cirrhosis <- load_bulk(BULK_CIRRHOSIS)

cat(sprintf("[bulk] all-MASLD     : %s (%s rows)\n", BULK_ALL,
            if (is.null(bulk_all)) "MISSING" else nrow(bulk_all)))
cat(sprintf("[bulk] Steatosis     : %s\n",
            if (is.null(bulk_steatosis)) "MISSING (fallback to all-MASLD)"
            else sprintf("%s (%s rows)", BULK_STEATOSIS, nrow(bulk_steatosis))))
cat(sprintf("[bulk] Steatohepatitis: %s\n",
            if (is.null(bulk_sh)) "MISSING (fallback to all-MASLD)"
            else sprintf("%s (%s rows)", BULK_SH, nrow(bulk_sh))))
cat(sprintf("[bulk] Cirrhosis     : %s\n",
            if (is.null(bulk_cirrhosis)) "MISSING (fallback to all-MASLD)"
            else sprintf("%s (%s rows)", BULK_CIRRHOSIS, nrow(bulk_cirrhosis))))

# ---- Helper: score concordance for an LMM slice against a chosen anchor ----
score_concord <- function(lmm, bulk_lookup) {
  if (nrow(lmm) == 0 || is.null(bulk_lookup)) return(data.table())
  lmm <- copy(lmm)
  lmm[, ligand   := sub("_.*", "", ligand_complex)]
  lmm[, receptor := sub("_.*", "", receptor_complex)]

  m_lig <- merge(lmm, bulk_lookup, by.x = "ligand", by.y = "gene", all.x = TRUE)
  setnames(m_lig, c("bulk_lfc", "bulk_padj"),
                  c("lig_bulk_lfc", "lig_bulk_padj"))
  m <- merge(m_lig, bulk_lookup, by.x = "receptor", by.y = "gene", all.x = TRUE)
  setnames(m, c("bulk_lfc", "bulk_padj"),
              c("rec_bulk_lfc", "rec_bulk_padj"))

  # Concordance: sign(sc estimate) == sign(bulk lfc) AND bulk padj < 0.1
  m[, sc_sign := sign(Estimate)]
  m[, lig_concordant := (sign(lig_bulk_lfc) == sc_sign) & (lig_bulk_padj < 0.1)]
  m[, rec_concordant := (sign(rec_bulk_lfc) == sc_sign) & (rec_bulk_padj < 0.1)]
  # Master review M-P0-9: require both sides to be testable (non-NA).
  m[, both_concordant := lig_concordant & rec_concordant &
                         !is.na(lig_bulk_padj) & !is.na(rec_bulk_padj)]

  # Rank-based signed -log10(padj) per gene (Phase 5 of rigor plan):
  # cross-modal score that does NOT rely on logFC magnitude (avoids Rich
  # 2026 pseudocount critique). Used for the Spearman rank cross-modal
  # summary computed downstream.
  m[, lig_bulk_signed_log10 := sign(lig_bulk_lfc) *
                               (-log10(pmax(lig_bulk_padj, 1e-300)))]
  m[, rec_bulk_signed_log10 := sign(rec_bulk_lfc) *
                               (-log10(pmax(rec_bulk_padj, 1e-300)))]
  m
}

# ---- Phase 5: Spearman rank cross-modal summary per (ct_pair) ---------------
# For each cell-type pair, compute Spearman rho between the scRNA effect
# estimate (sc_sign * |Estimate|) and bulk signed -log10(padj) for the
# ligand and receptor genes. This is the rank-only firewall the rigor plan
# Phase 5 requires. Reported alongside the binary both_concordant counts.
rank_summary <- function(conc, axis_tag) {
  if (nrow(conc) == 0) return(data.table())
  conc[, sc_score := Estimate]  # already signed continuous
  per_pair <- conc[!is.na(lig_bulk_signed_log10),
                   .(n = .N,
                     rho_lig_spearman = suppressWarnings(stats::cor(
                       sc_score, lig_bulk_signed_log10,
                       method = "spearman", use = "pairwise.complete.obs")),
                     rho_rec_spearman = suppressWarnings(stats::cor(
                       sc_score, rec_bulk_signed_log10,
                       method = "spearman", use = "pairwise.complete.obs"))),
                   by = .(source, target)]
  per_pair[, axis := axis_tag]
  per_pair
}

# ---- Score each stage contrast --------------------------------------------
coarse <- fread(COARSE_TSV)

score_and_save <- function(term_label, stage_anchor, axis_name, out_fname) {
  slice <- coarse[term == term_label]
  if (nrow(slice) == 0) {
    cat(sprintf("[skip] no rows for term=%s\n", term_label))
    return(data.table())
  }
  # Choose anchor: stage-matched if available, else fallback to all-MASLD
  if (is.null(stage_anchor)) {
    anchor <- bulk_all
    axis_tag <- paste0(axis_name, "_fallback_all_MASLD")
  } else {
    anchor <- stage_anchor
    axis_tag <- axis_name
  }
  conc <- score_concord(slice, anchor)
  if (nrow(conc) == 0) {
    cat(sprintf("[skip] empty concordance for %s\n", term_label))
    return(data.table())
  }
  conc[, axis := axis_tag]
  fwrite(conc, file.path(OUT_DIR, out_fname), sep = "\t")
  testable <- conc[!is.na(both_concordant)]
  cat(sprintf("[output] %s: %d rows -> %s\n",
              axis_tag, nrow(conc), out_fname))
  cat(sprintf("[summary] %s: both_concordant=%d (%.1f%%) of %d testable; ",
              axis_tag,
              sum(testable$both_concordant, na.rm = TRUE),
              100 * mean(testable$both_concordant, na.rm = TRUE),
              nrow(testable)))
  cat(sprintf("lig_concordant=%.1f%%, rec_concordant=%.1f%%\n",
              100 * mean(testable$lig_concordant, na.rm = TRUE),
              100 * mean(testable$rec_concordant, na.rm = TRUE)))
  conc
}

steat_conc <- score_and_save(
  "disease_stage_coarseSteatosis",       bulk_steatosis,
  "Steatosis_vs_Healthy",                "lr_bulk_concordance_steatosis.tsv")
sh_conc    <- score_and_save(
  "disease_stage_coarseSteatohepatitis", bulk_sh,
  "Steatohepatitis_vs_Healthy",          "lr_bulk_concordance_sh.tsv")
cirr_conc  <- score_and_save(
  "disease_stage_coarseCirrhosis",       bulk_cirrhosis,
  "Cirrhosis_vs_Healthy",                "lr_bulk_concordance_cirrhosis.tsv")

# Continuous axis (macrophage pseudotime) stays on the all-MASLD anchor -----
cont_conc <- data.table()
if (file.exists(CONTINUOUS_TSV)) {
  cont <- fread(CONTINUOUS_TSV)
  cont_conc <- score_concord(cont, bulk_all)
  if (nrow(cont_conc) > 0) {
    cont_conc[, axis := "macrophage_pseudotime_all_MASLD"]
    fwrite(cont_conc, file.path(OUT_DIR, "lr_bulk_concordance_continuous.tsv"),
           sep = "\t")
    testable <- cont_conc[!is.na(both_concordant)]
    cat(sprintf("[output] continuous: %d rows -> lr_bulk_concordance_continuous.tsv\n",
                nrow(cont_conc)))
    cat(sprintf("[summary] continuous: both_concordant=%d (%.1f%%) of %d testable\n",
                sum(testable$both_concordant, na.rm = TRUE),
                100 * mean(testable$both_concordant, na.rm = TRUE),
                nrow(testable)))
  }
}

# ---- F-stage continuous axis ---------------------------------------------
# Scored against TWO bulk anchors:
#   (1) all-MASLD dream  -> appropriate for continuous severity (0-4)
#   (2) SH-vs-Healthy    -> secondary check since SH donors span F1-F3
fstage_conc_all <- data.table()
fstage_conc_sh  <- data.table()
if (file.exists(FSTAGE_TSV)) {
  fstage <- fread(FSTAGE_TSV)
  fstage <- fstage[term == "F_stage_numeric"]
  if (nrow(fstage) > 0) {
    # Anchor (1): all-MASLD
    fstage_conc_all <- score_concord(fstage, bulk_all)
    if (nrow(fstage_conc_all) > 0) {
      fstage_conc_all[, axis := "F_stage_continuous_all_MASLD"]
      testable <- fstage_conc_all[!is.na(both_concordant)]
      cat(sprintf("[output] F-stage (all-MASLD anchor): %d rows; ",
                  nrow(fstage_conc_all)))
      cat(sprintf("both_concordant=%d (%.1f%%) of %d testable; lig=%.1f%%, rec=%.1f%%\n",
                  sum(testable$both_concordant, na.rm = TRUE),
                  100 * mean(testable$both_concordant, na.rm = TRUE),
                  nrow(testable),
                  100 * mean(testable$lig_concordant, na.rm = TRUE),
                  100 * mean(testable$rec_concordant, na.rm = TRUE)))
    }
    # Anchor (2): SH-vs-Healthy
    sh_anchor <- if (!is.null(bulk_sh)) bulk_sh else bulk_all
    fstage_conc_sh <- score_concord(fstage, sh_anchor)
    if (nrow(fstage_conc_sh) > 0) {
      fstage_conc_sh[, axis := if (!is.null(bulk_sh))
        "F_stage_continuous_SH_anchor" else
        "F_stage_continuous_SH_anchor_fallback_all_MASLD"]
      testable <- fstage_conc_sh[!is.na(both_concordant)]
      cat(sprintf("[output] F-stage (SH anchor): %d rows; ",
                  nrow(fstage_conc_sh)))
      cat(sprintf("both_concordant=%d (%.1f%%) of %d testable; lig=%.1f%%, rec=%.1f%%\n",
                  sum(testable$both_concordant, na.rm = TRUE),
                  100 * mean(testable$both_concordant, na.rm = TRUE),
                  nrow(testable),
                  100 * mean(testable$lig_concordant, na.rm = TRUE),
                  100 * mean(testable$rec_concordant, na.rm = TRUE)))
    }
    fstage_both <- rbindlist(list(fstage_conc_all, fstage_conc_sh), fill = TRUE)
    fwrite(fstage_both, file.path(OUT_DIR, "lr_bulk_concordance_fstage.tsv"),
           sep = "\t")
    cat(sprintf("[output] %d total F-stage rows -> lr_bulk_concordance_fstage.tsv\n",
                nrow(fstage_both)))
  }
} else {
  cat("[skip] stage_lr_lmm_fstage.tsv missing; F-stage concordance skipped.\n")
}

# Combined headline summary -------------------------------------------------
all_conc <- rbindlist(list(steat_conc, sh_conc, cirr_conc, cont_conc,
                           fstage_conc_all, fstage_conc_sh),
                      fill = TRUE)
fwrite(all_conc, file.path(OUT_DIR, "lr_bulk_concordance.tsv"), sep = "\t")
cat(sprintf("\n[combined] %d total axis x LR rows -> lr_bulk_concordance.tsv\n",
            nrow(all_conc)))

# Phase 5: per-ct_pair Spearman rank cross-modal summary -------------------
# Computed across each axis-contrast separately so per-cell-type-pair rank
# concordance is auditable. Output: rho_lig_spearman, rho_rec_spearman per
# (axis, source, target).
rank_summaries <- list()
if (nrow(steat_conc) > 0) {
  rank_summaries[["Steatosis_vs_Healthy"]] <- rank_summary(steat_conc, "Steatosis_vs_Healthy")
}
if (nrow(sh_conc) > 0) {
  rank_summaries[["Steatohepatitis_vs_Healthy"]] <- rank_summary(sh_conc, "Steatohepatitis_vs_Healthy")
}
if (nrow(cirr_conc) > 0) {
  rank_summaries[["Cirrhosis_vs_Healthy"]] <- rank_summary(cirr_conc, "Cirrhosis_vs_Healthy")
}
if (nrow(cont_conc) > 0) {
  rank_summaries[["macrophage_pseudotime_all_MASLD"]] <- rank_summary(cont_conc, "macrophage_pseudotime_all_MASLD")
}
if (nrow(fstage_conc_all) > 0) {
  rank_summaries[["F_stage_continuous_all_MASLD"]] <- rank_summary(fstage_conc_all, "F_stage_continuous_all_MASLD")
}
if (nrow(fstage_conc_sh) > 0) {
  rank_summaries[["F_stage_continuous_SH_anchor"]] <- rank_summary(fstage_conc_sh, "F_stage_continuous_SH_anchor")
}
rank_all <- rbindlist(rank_summaries, fill = TRUE)
if (nrow(rank_all) > 0) {
  fwrite(rank_all,
         file.path(OUT_DIR, "lr_bulk_concordance_spearman_rank.tsv"),
         sep = "\t")
  cat(sprintf("\n[rank] per-(axis,source,target) Spearman summary -> lr_bulk_concordance_spearman_rank.tsv (%d rows)\n",
              nrow(rank_all)))
  # Headline: mean per-axis Spearman rho for ligand-side rank concordance
  hdr_rank <- rank_all[!is.na(rho_lig_spearman),
                       .(mean_rho_lig = mean(rho_lig_spearman),
                         mean_rho_rec = mean(rho_rec_spearman, na.rm = TRUE),
                         n_ct_pairs   = .N),
                       by = axis]
  cat("\n===== Per-axis MEAN Spearman rho (rank-only Phase 5) =====\n")
  print(hdr_rank)
}

# Per-axis headline -- exactly what gets reported back to the agent ---------
cat("\n===== HEADLINE both_concordant per axis =====\n")
if (nrow(all_conc) > 0) {
  hdr <- all_conc[!is.na(both_concordant),
    .(n_testable = .N,
      n_both     = sum(both_concordant, na.rm = TRUE),
      pct_both   = 100 * mean(both_concordant, na.rm = TRUE),
      pct_lig    = 100 * mean(lig_concordant, na.rm = TRUE),
      pct_rec    = 100 * mean(rec_concordant, na.rm = TRUE)),
    by = axis]
  print(hdr)
}
