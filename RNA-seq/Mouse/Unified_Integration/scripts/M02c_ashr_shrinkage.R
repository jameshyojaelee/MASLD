#!/usr/bin/env Rscript
# M02c_ashr_shrinkage.R
# ---------------------------------------------------------------------------
# Adaptive shrinkage (ashr; Stephens 2017) on each per-diet DE result.
#
# Why: the diet groups span very different sample sizes (Western n=42 ... NASH
# n=249). A fixed |logFC| cutoff mixes effect size with how precisely each diet
# could measure it. ashr shrinks each gene's logFC toward zero in proportion to
# its standard error, so a uniform threshold on the SHRUNK logFC is fair across
# sample sizes. Mirrors the human dream pipeline (05b_ashr_shrinkage.R).
#
# Adds columns to each per_diet/{diet}_de_results.csv (backed up first):
#   shrunk_logFC, shrunk_se, lfsr, svalue
#
# UP membership downstream: lfsr < 0.05 & shrunk_logFC > <threshold>.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ashr)
})

PROJECT <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEDIR <- file.path(PROJECT, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
BK    <- file.path(DEDIR, "pre_ashr_backup")
dir.create(BK, recursive = TRUE, showWarnings = FALSE)

# Canonical diet axis = 4 groups (2026-05-29): NASH_diet + Western_diet merged
# into one "Western" metabolic-overload group. NASH_diet/Western_diet still exist
# on disk for provenance but are no longer the canonical replication axes.
DIETS <- c("MCD", "CDAHFD", "Western", "HFD")

run_ashr_one <- function(diet) {
  f <- file.path(DEDIR, paste0(diet, "_de_results.csv"))
  dt <- fread(f)

  # betahat = logFC; sebetahat = unmoderated SE if present, else moderated SE
  # (|logFC / t|). Both are valid SE estimates for ash().
  betahat <- dt$logFC
  if ("SE_unmoderated" %in% names(dt)) {
    se <- dt$SE_unmoderated
    se_src <- "SE_unmoderated"
  } else {
    se <- abs(dt$logFC / dt$t)
    se_src <- "moderated (logFC/t)"
  }
  # Guard non-finite / zero SE: replace with the max finite SE (→ heavy
  # shrinkage for those genes, which are LFC≈0/t≈0 anyway).
  bad <- !is.finite(se) | se <= 0
  if (any(bad)) se[bad] <- max(se[is.finite(se) & se > 0], na.rm = TRUE)

  fit <- ash(betahat, se, mixcompdist = "normal", method = "shrink")

  dt[, shrunk_logFC := get_pm(fit)]      # posterior mean
  dt[, shrunk_se    := get_psd(fit)]     # posterior SD
  dt[, lfsr         := get_lfsr(fit)]
  dt[, svalue       := get_svalue(fit)]

  file.copy(f, file.path(BK, basename(f)), overwrite = TRUE)
  fwrite(dt, f)

  up05 <- dt[lfsr < 0.05 & shrunk_logFC > 0.5, .N]
  cat(sprintf("  %-13s SE=%-18s genes=%d  UP(lfsr<.05 & shrunk>0.5)=%d\n",
              diet, se_src, nrow(dt), up05))
  dt[, .(gene, lfsr, shrunk_logFC)]
}

cat("=== ashr per-diet shrinkage ===\n")
shr <- lapply(DIETS, run_ashr_one)
names(shr) <- DIETS

# Threshold scan: UP counts at lfsr<0.05 across shrunk_logFC cutoffs
cat("\n=== UP counts (lfsr<0.05) by shrunk_logFC threshold ===\n")
thr <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
hdr <- sprintf("%-13s", "Diet"); for (t in thr) hdr <- paste0(hdr, sprintf("%8s", paste0(">", t)))
cat(hdr, "\n")
for (d in DIETS) {
  s <- shr[[d]]
  row <- sprintf("%-13s", d)
  for (t in thr) row <- paste0(row, sprintf("%8d", s[lfsr < 0.05 & shrunk_logFC > t, .N]))
  cat(row, "\n")
}

# Cross-diet replication preview at each threshold (genes UP in >=2 / >=3 diets)
cat("\n=== Cross-diet replication (genes UP in N diets) by threshold ===\n")
for (t in thr) {
  upsets <- lapply(DIETS, function(d) shr[[d]][lfsr < 0.05 & shrunk_logFC > t,
                                               sub("[.][0-9]+$", "", gene)])
  allg <- unique(unlist(upsets))
  ndiets <- sapply(allg, function(g) sum(sapply(upsets, function(u) g %in% u)))
  cat(sprintf("  shrunk>%.1f : union=%d  >=2 diets=%d  >=3 diets=%d  >=4=%d  all5=%d\n",
              t, length(allg), sum(ndiets >= 2), sum(ndiets >= 3),
              sum(ndiets >= 4), sum(ndiets == 5)))
}
cat("\nDone. Augmented per_diet CSVs in place (backups in pre_ashr_backup/).\n")
