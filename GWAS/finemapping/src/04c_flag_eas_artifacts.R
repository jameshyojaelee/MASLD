#!/usr/bin/env Rscript
# 04c_flag_eas_artifacts.R
# Flags EAS SuSiE "ghost credible sets" — single non-significant variants
# assigned PIP>0.5 due to LD mismatch (1000G EAS panel has only 504 samples).
# Compares SuSiE vs CARMA top variants per locus and recommends preferred method.
#
# Problem: SuSiE with small EAS LD panels produces spurious credible sets where
# the top PIP variant has a non-significant GWAS p-value (>0.01). CARMA's
# Bayesian model averaging is more robust to LD misspecification.
#
# Outputs:
#   results/eas_quality_flags.csv — per-locus flags + method recommendations
#
# Usage: Rscript src/04c_flag_eas_artifacts.R
# Environment: finemapping conda env with data.table

suppressPackageStartupMessages({
  library(data.table)
})

# ── Configuration ──────────────────────────────────────────────────────────
BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
setwd(FM_DIR)

RESULTS_DIR <- file.path(FM_DIR, "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

GHOST_PVAL_THRESH <- 0.01   # p-value above which a CS top variant is "ghost"
DISTANCE_THRESH   <- 10000  # 10kb: if SuSiE-CARMA top variants differ by >10kb, flag

cat("============================================================\n")
cat("EAS SuSiE Ghost Credible Set Flagging\n")
cat("Ghost p-value threshold:", GHOST_PVAL_THRESH, "\n")
cat("Distance threshold:", DISTANCE_THRESH / 1000, "kb\n")
cat("============================================================\n\n")

# ── Step 1: Discover all EAS SuSiE result files ──────────────────────────
# Two naming patterns:
#   output/BBJ_*/EAS_0.5Mb/susie/*.tsv          (BBJ liver enzyme GWAS)
#   output/*_EAS/EAS_0.5Mb/susie/*.tsv           (EAS ancestry disease GWAS)
susie_files <- c(
  Sys.glob(file.path(FM_DIR, "output", "BBJ_*", "EAS_0.5Mb", "susie", "*.tsv")),
  Sys.glob(file.path(FM_DIR, "output", "*_EAS", "EAS_0.5Mb", "susie", "*.tsv"))
)
susie_files <- unique(susie_files)

if (length(susie_files) == 0) {
  stop("No EAS SuSiE result files found.")
}
cat("Found", length(susie_files), "EAS SuSiE result files\n\n")

# ── Step 2: Process each locus ────────────────────────────────────────────
results <- list()

for (f in susie_files) {
  # Parse study and locus from path
  parts <- strsplit(f, "/")[[1]]
  susie_idx <- which(parts == "susie")
  study <- parts[susie_idx - 2]
  fname <- parts[susie_idx + 1]

  # Extract locus ID (chr.bp) from filename
  locus_match <- regmatches(fname, regexpr("\\d+\\.\\d+", fname))
  if (length(locus_match) == 0) next
  locus <- locus_match[1]
  locus_chr <- as.integer(sub("\\..*", "", locus))
  locus_bp  <- as.integer(sub(".*\\.", "", locus))

  # ── Read SuSiE results ──
  susie_dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(susie_dt) || nrow(susie_dt) == 0) next

  # Overall SuSiE top variant (highest PIP across all CS, including CS=0)
  susie_top_row <- susie_dt[which.max(PIP)]
  susie_top_variant <- susie_top_row$SNP
  susie_top_pip     <- susie_top_row$PIP
  susie_top_pval    <- susie_top_row$pval
  susie_top_pos     <- susie_top_row$position
  susie_converged   <- ifelse("converged" %in% names(susie_dt),
                              all(susie_dt$converged == TRUE, na.rm = TRUE),
                              !grepl("notconverged", fname))
  susie_n_cs        <- length(unique(susie_dt[CS > 0]$CS))

  # ── Ghost CS detection ──
  # A "ghost CS" is one where the top PIP variant within that CS has a
  # non-significant GWAS p-value. This indicates SuSiE's L2 regularization
  # on a small LD panel created a spurious signal.
  ghost_cs <- FALSE
  ghost_cs_ids <- c()
  cs_vals <- sort(unique(susie_dt[CS > 0]$CS))
  for (cs_val in cs_vals) {
    cs_sub <- susie_dt[CS == cs_val]
    cs_top <- cs_sub[which.max(PIP)]
    if (cs_top$pval > GHOST_PVAL_THRESH) {
      ghost_cs <- TRUE
      ghost_cs_ids <- c(ghost_cs_ids, cs_val)
    }
  }
  n_ghost_cs <- length(ghost_cs_ids)

  # ── Read CARMA results ──
  # CARMA text file: {study}_1KG_EAS_0.5Mb_{locus}.txt.gz
  carma_txt <- file.path(FM_DIR, "output", study, "EAS_0.5Mb", "CARMA",
                         paste0(study, "_1KG_EAS_0.5Mb_", locus, ".txt.gz"))
  carma_rds <- file.path(FM_DIR, "output", study, "EAS_0.5Mb", "CARMA",
                         paste0(study, "_1KG_EAS_0.5Mb_", locus, ".rds"))

  carma_top_variant <- NA_character_
  carma_top_pip     <- NA_real_
  carma_top_pos     <- NA_integer_
  carma_n_cs        <- NA_integer_
  carma_available   <- FALSE

  # Prefer txt.gz (tabular, has SNP IDs); fall back to rds
  if (file.exists(carma_txt) && file.size(carma_txt) > 0) {
    carma_dt <- tryCatch(fread(carma_txt), error = function(e) NULL)
    if (!is.null(carma_dt) && nrow(carma_dt) > 0) {
      carma_available <- TRUE
      carma_top_row <- carma_dt[which.max(PIP)]
      carma_top_variant <- carma_top_row$SNP
      carma_top_pip     <- carma_top_row$PIP
      carma_top_pos     <- carma_top_row$position
      carma_n_cs        <- length(unique(carma_dt[CS > 0]$CS))
    }
  } else if (file.exists(carma_rds)) {
    # Fall back to RDS: extract PIPs (unnamed numeric vector indexed by position)
    carma_obj <- tryCatch(readRDS(carma_rds), error = function(e) NULL)
    if (!is.null(carma_obj) && length(carma_obj) > 0) {
      # CARMA RDS structure: list(list(PIPs=..., "Credible set"=..., ...))
      inner <- if (is.list(carma_obj[[1]]) && "PIPs" %in% names(carma_obj[[1]])) {
        carma_obj[[1]]
      } else {
        carma_obj
      }
      if ("PIPs" %in% names(inner)) {
        carma_available <- TRUE
        pips <- inner[["PIPs"]]
        top_idx <- which.max(pips)
        carma_top_pip <- pips[top_idx]
        # RDS PIPs are unnamed — match back to SuSiE variant order for position
        if (top_idx <= nrow(susie_dt)) {
          carma_top_variant <- susie_dt$SNP[top_idx]
          carma_top_pos     <- susie_dt$position[top_idx]
        }
        cs_list <- inner[["Credible set"]]
        carma_n_cs <- if (is.list(cs_list)) length(cs_list) else 0L
      }
    }
  }

  # ── Distance between SuSiE and CARMA top variants ──
  distance <- NA_integer_
  if (!is.na(susie_top_pos) && !is.na(carma_top_pos)) {
    distance <- abs(as.integer(susie_top_pos) - as.integer(carma_top_pos))
  }

  # ── Method recommendation ──
  # CARMA preferred if ghost CS detected or methods disagree by >10kb
  recommended <- if (ghost_cs) {
    "CARMA"
  } else if (!is.na(distance) && distance > DISTANCE_THRESH) {
    "CARMA"
  } else if (!susie_converged) {
    "CARMA"
  } else {
    "either"
  }

  # ── Store result ──
  results[[length(results) + 1]] <- data.table(
    study            = study,
    locus            = locus,
    locus_chr        = locus_chr,
    locus_bp         = locus_bp,
    susie_converged  = susie_converged,
    susie_n_cs       = susie_n_cs,
    susie_top_variant = susie_top_variant,
    susie_top_pip    = susie_top_pip,
    susie_top_pval   = susie_top_pval,
    carma_available  = carma_available,
    carma_n_cs       = carma_n_cs,
    carma_top_variant = carma_top_variant,
    carma_top_pip    = carma_top_pip,
    distance_bp      = distance,
    ghost_cs         = ghost_cs,
    n_ghost_cs       = n_ghost_cs,
    ghost_cs_ids     = paste(ghost_cs_ids, collapse = ";"),
    recommended_method = recommended
  )
}

# ── Step 3: Combine and write ────────────────────────────────────────────
out <- rbindlist(results)
out <- out[order(study, locus)]

fwrite(out, file.path(RESULTS_DIR, "eas_quality_flags.csv"))
cat("Wrote:", file.path(RESULTS_DIR, "eas_quality_flags.csv"), "\n")
cat("  Rows:", nrow(out), "\n\n")

# ── Step 4: Summary ─────────────────────────────────────────────────────
cat("============================================================\n")
cat("EAS SuSiE Quality Summary\n")
cat("============================================================\n\n")

cat("Total EAS loci:", nrow(out), "\n")
cat("  With credible sets (SuSiE):", sum(out$susie_n_cs > 0), "\n")
cat("  SuSiE converged:", sum(out$susie_converged), "\n")
cat("  CARMA available:", sum(out$carma_available), "\n\n")

n_ghost     <- sum(out$ghost_cs)
n_with_cs   <- sum(out$susie_n_cs > 0)
ghost_frac  <- if (n_with_cs > 0) 100 * n_ghost / n_with_cs else 0
cat("Ghost CS flagged:", n_ghost, "of", n_with_cs, "loci with CS",
    sprintf("(%.1f%%)\n", ghost_frac))

n_distant <- sum(!is.na(out$distance_bp) & out$distance_bp > DISTANCE_THRESH)
cat("SuSiE-CARMA top variant >10kb apart:", n_distant, "loci\n")

n_carma_rec <- sum(out$recommended_method == "CARMA")
n_either    <- sum(out$recommended_method == "either")
cat("\nMethod recommendation:\n")
cat("  CARMA preferred:", n_carma_rec, "loci\n")
cat("  Either method:  ", n_either, "loci\n")

# Per-study breakdown
cat("\nPer-study breakdown:\n")
study_stats <- out[, .(
  n_loci      = .N,
  n_ghost     = sum(ghost_cs),
  n_distant   = sum(!is.na(distance_bp) & distance_bp > DISTANCE_THRESH),
  n_carma_rec = sum(recommended_method == "CARMA"),
  median_dist = median(distance_bp, na.rm = TRUE)
), by = study][order(-n_ghost)]

for (i in 1:nrow(study_stats)) {
  s <- study_stats[i]
  cat(sprintf("  %-35s %2d loci, %2d ghost (%.0f%%), %2d distant, median dist=%s bp\n",
              s$study, s$n_loci, s$n_ghost,
              if (s$n_loci > 0) 100 * s$n_ghost / s$n_loci else 0,
              s$n_distant,
              format(s$median_dist, big.mark = ",")))
}

cat("\n============================================================\n")
cat("Recommendation: Use CARMA PIPs for EAS loci flagged as ghost_cs=TRUE\n")
cat("or where SuSiE-CARMA top variants are >10kb apart.\n")
cat("============================================================\n")
