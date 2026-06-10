#!/usr/bin/env Rscript
# 130b_fix_progression_lvqw_sanity.R
# ---------------------------------------------------------------------------
# Recompute ONLY the LVQW-vs-dream sanity gate for the 9 migrated progression
# contrasts. READ-ONLY on the staged {contrast_id}_lvqw.csv files — does NOT
# re-fit anything, does NOT touch the LVQW outputs or any canonical file.
#
# Bug fixed: the original gate joined on the raw `gene` column. The LVQW outputs
# key on VERSIONED ENSG (ENSG00000241860.8) while the dream reference files key
# on UNVERSIONED ENSG (ENSG00000241860). 100% of shared genes differed on the
# raw key -> n_shared=0 -> concordance never computed. Fix: strip the version
# suffix sub("[.][0-9]+$","",gene) on BOTH sides and join on the stripped key.
#
# Recomputes per contrast:
#   n_shared                       genes present in both (stripped-key join)
#   direction_concordance_shared   mean(sign(lvqw_logFC)==sign(dream_logFC)) over shared
#   direction_concordance_bothsig  same, restricted to genes padj<0.05 in BOTH
#   n_both_sig                     count of both-significant genes
#   logFC_spearman                 Spearman rho of logFC over shared genes
#   logFC_pearson                  Pearson r of logFC over shared genes
#
# Degenerate-reference flag: the dream script re-ran filterByExpr() on raw counts
# per contrast, yielding gene universes from 26,638 to 40,626 — several LARGER
# than the curated 27,638-gene merged_dge universe the LVQW fits use. A dream
# reference with dream_n_genes far above the universe tested many genes LVQW
# never includes (low overlap fraction is structural, not a regression). Flagged
# when dream_n_genes > 1.10 * 27,638. c13 is additionally flagged as a
# covariate-adjusted (fibrosis-adjusted) contrast where far fewer DEGs and a
# weaker rho are expected by design.
# ---------------------------------------------------------------------------
suppressMessages(library(data.table))

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
ODIR <- file.path(BASE, "analysis/integration/results/progression")

UNIVERSE_N <- 27638L   # merged_dge fixed 9-cohort universe used by all LVQW fits

strip_ver <- function(x) sub("[.][0-9]+$", "", x)

contrasts <- c(
  "c3_adv_vs_early_fib", "c4_nafl_vs_ctrl", "c5_nas_ge5_vs_lt5",
  "c6_extreme_endpoints", "c8_cirrhosis", "c9_f2_inflection",
  "c11_nash_vs_ctrl", "c12_early_vs_late_nash", "c13_nash_vs_nafl_fib_adj"
)

sanity <- rbindlist(lapply(contrasts, function(cid) {
  lv_file <- file.path(ODIR, sprintf("%s_lvqw.csv", cid))
  dr_file <- file.path(ODIR, sprintf("%s_dream.csv", cid))

  if (!file.exists(lv_file)) {
    return(data.table(contrast_id = cid, status = "LVQW_MISSING"))
  }
  lv <- fread(lv_file)
  lv[, eb := strip_ver(gene)]
  # collapse any (rare) duplicate stripped keys: keep the most significant row
  setorder(lv, padj)
  lv <- unique(lv, by = "eb")
  lvqw_n_genes <- nrow(lv)
  lvqw_n_deg05 <- sum(lv$padj < 0.05, na.rm = TRUE)

  if (!file.exists(dr_file)) {
    return(data.table(contrast_id = cid, status = "NO_DREAM_REF",
      lvqw_n_genes = lvqw_n_genes, lvqw_n_deg05 = lvqw_n_deg05))
  }
  dr <- fread(dr_file)
  dr[, eb := strip_ver(gene)]
  setorder(dr, padj)
  dr <- unique(dr, by = "eb")
  dream_n_genes <- nrow(dr)
  dream_n_deg05 <- sum(dr$padj < 0.05, na.rm = TRUE)

  # stripped-key join
  m <- merge(lv[, .(eb, lv_lfc = logFC, lv_padj = padj)],
             dr[, .(eb, dr_lfc = logFC, dr_padj = padj)], by = "eb")
  n_shared <- nrow(m)

  if (n_shared > 2) {
    rho_s <- cor(m$lv_lfc, m$dr_lfc, method = "spearman", use = "complete.obs")
    r_p   <- cor(m$lv_lfc, m$dr_lfc, method = "pearson",  use = "complete.obs")
    dir_shared <- mean(sign(m$lv_lfc) == sign(m$dr_lfc), na.rm = TRUE)
    bs <- m[lv_padj < 0.05 & dr_padj < 0.05]
    n_both_sig <- nrow(bs)
    dir_bothsig <- if (n_both_sig > 0) mean(sign(bs$lv_lfc) == sign(bs$dr_lfc), na.rm = TRUE) else NA_real_
  } else {
    rho_s <- r_p <- dir_shared <- dir_bothsig <- NA_real_
    n_both_sig <- NA_integer_
  }

  # degenerate-reference flags
  flags <- character(0)
  if (dream_n_genes > 1.10 * UNIVERSE_N) {
    flags <- c(flags, sprintf("dream_overfiltered(n=%d>univ*1.1)", dream_n_genes))
  }
  if (cid == "c13_nash_vs_nafl_fib_adj") {
    flags <- c(flags, "fib_adjusted_low_DEG_expected")
  }

  data.table(
    contrast_id   = cid,
    status        = "OK",
    lvqw_n_genes  = lvqw_n_genes,
    lvqw_n_deg05  = lvqw_n_deg05,
    dream_n_genes = dream_n_genes,
    dream_n_deg05 = dream_n_deg05,
    n_shared      = n_shared,
    shared_frac_of_lvqw = round(n_shared / lvqw_n_genes, 4),
    logFC_spearman              = round(rho_s, 4),
    logFC_pearson               = round(r_p, 4),
    direction_concordance_shared  = round(dir_shared, 4),
    n_both_sig                    = n_both_sig,
    direction_concordance_bothsig = round(dir_bothsig, 4),
    degenerate_flag = if (length(flags)) paste(flags, collapse = "; ") else ""
  )
}), fill = TRUE)

out_file <- file.path(ODIR, "progression_lvqw_sanity.csv")
fwrite(sanity, out_file)
cat("=== CORRECTED SANITY GATE (LVQW vs dream, stripped-ENSG join) ===\n\n")
print(sanity, row.names = FALSE)
cat(sprintf("\nRewrote: %s\n", out_file))
