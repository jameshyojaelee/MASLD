#!/usr/bin/env Rscript
# 03_run_fm_per_locus.R
# Runs SuSiE + CARMA fine-mapping for a single locus using local LD matrices
# Usage: Rscript 03_run_fm_per_locus.R <sumstats_name> <ld_pop> <locus> <N_tot> <N_cases> <window_mb> <ancestry>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript 03_run_fm_per_locus.R <sumstats_name> <ld_pop> <locus> <N_tot> <N_cases> <window_mb> <ancestry>")
}

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
setwd(FM_DIR)
source("src/finemapping_functions.R")

library(susieR)
library(CARMA)

sumstats_name <- as.character(args[1])
ld_pop <- as.character(args[2])
locus <- as.character(args[3])
N_tot <- as.numeric(args[4])
N_cases <- as.numeric(args[5])  # NA for quantitative traits
window_mb <- as.numeric(args[6])
ancestry <- as.character(args[7])  # EUR or EAS

if (N_cases == 0 || is.nan(N_cases)) N_cases <- NA

cat("============================================================\n")
cat("Fine-mapping locus:", locus, "\n")
cat("Study:", sumstats_name, "\n")
cat("LD population:", ld_pop, "(", ancestry, ")\n")
cat("N_tot:", N_tot, " N_cases:", ifelse(is.na(N_cases), "NA (quantitative)", N_cases), "\n")
cat("Window:", window_mb, "Mb\n")
cat("============================================================\n")

# Parse locus into CHR and BP
locus_parts <- strsplit(locus, "\\.")[[1]]
CHR <- as.numeric(locus_parts[1])
BP <- as.numeric(locus_parts[2])
window_bp <- window_mb * 1000000
START <- max(1, BP - window_bp)
END <- BP + window_bp

# Read per-locus summary stats
ss_path <- file.path(FM_DIR, "output", sumstats_name,
                     paste0(ld_pop, "_", window_mb, "Mb"), "ss",
                     paste0(sumstats_name, "_", window_mb, "Mb_", locus, ".txt"))

if (!file.exists(ss_path)) {
  stop(paste("Per-locus summary stats not found:", ss_path))
}

ss <- fread(ss_path)
cat("Loaded", nrow(ss), "variants for locus\n")

# Extract LD matrix from local UKBB/1KG files
cat("\n--- Extracting LD matrix ---\n")
dat <- get_ld_per_locus(ss, locus, CHR, START, END, ancestry = ancestry)

if (is.null(dat)) {
  cat("SKIPPING locus", locus, "- could not extract LD\n")
  # Write skip marker
  writeLines(paste("SKIPPED:", locus, "- LD extraction failed"),
             file.path(FM_DIR, "output", sumstats_name,
                       paste0(ld_pop, "_", window_mb, "Mb"),
                       paste0("SKIPPED_", locus, ".txt")))
  quit(save = "no", status = 0)
}

ss_filtered <- dat[[1]]
ld_filtered <- dat[[2]]

cat("After LD matching:", nrow(ss_filtered), "variants\n")
cat("LD matrix dimensions:", nrow(ld_filtered), "x", ncol(ld_filtered), "\n")

LDpanel <- switch(ancestry,
  EUR = "UKBB",
  EAS = "1KG_EAS",
  AFR = "1KG_AFR",
  AMR = "1KG_AMR",
  SAS = "1KG_SAS",
  stop(paste("Unsupported ancestry for LDpanel:", ancestry))
)

# Run SuSiE
cat("\n--- Running SuSiE ---\n")
set.seed(42)  # C7 reproducibility: deterministic per-locus SuSiE
run_susie(ss_filtered, ld_filtered, N_tot, N_cases,
          sumstats_name, ld_pop, window_mb, locus, LDpanel,
          base_dir = FM_DIR)

# Run CARMA (skippable via SKIP_CARMA=1 — e.g. exploratory non-EUR runs where
# SuSiE + the convergence/purity/lambda_s guards suffice and CARMA's per-locus
# cost isn't justified at scale)
if (Sys.getenv("SKIP_CARMA", "0") != "1") {
  cat("\n--- Running CARMA ---\n")
  set.seed(42)  # C7 reproducibility: deterministic per-locus CARMA
  run_CARMA(ss_filtered, ld_filtered,
            sumstats_name, ld_pop, window_mb, locus, LDpanel,
            base_dir = FM_DIR)
} else {
  cat("\n--- CARMA skipped (SKIP_CARMA=1) ---\n")
}

cat("\n=== Completed fine-mapping for locus", locus, "===\n")
