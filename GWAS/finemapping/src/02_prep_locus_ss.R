#!/usr/bin/env Rscript
# 02_prep_locus_ss.R
# Extracts per-locus summary stats around each lead SNP
# Usage: Rscript 02_prep_locus_ss.R <sumstats_name> <path_to_sumstats> <path_to_leadSNPs> <ld_pop> <window_mb>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript 02_prep_locus_ss.R <sumstats_name> <path_to_sumstats> <path_to_leadSNPs> <ld_pop> <window_mb>")
}

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
setwd(FM_DIR)
source("src/finemapping_functions.R")

sumstats_name <- as.character(args[1])
path_to_sumstats <- as.character(args[2])
path_to_leadSNPs <- as.character(args[3])
ld_pop <- as.character(args[4])
window_mb <- as.numeric(args[5])

cat("=== Preparing per-locus summary stats ===\n")
cat("  Study:", sumstats_name, "\n")
cat("  Summary stats:", path_to_sumstats, "\n")
cat("  Lead SNPs:", path_to_leadSNPs, "\n")
cat("  LD population:", ld_pop, "\n")
cat("  Window:", window_mb, "Mb\n")

# Read lead SNPs
df.leadSNP <- fread(path_to_leadSNPs, colClasses = c("numeric", "numeric", "character"))
# Handle header row (skip if first row is header)
if (df.leadSNP$CHR[1] == "CHR" || is.na(as.numeric(df.leadSNP$CHR[1]))) {
  df.leadSNP <- df.leadSNP[-1, ]
  df.leadSNP$CHR <- as.numeric(df.leadSNP$CHR)
  df.leadSNP$BP <- as.numeric(df.leadSNP$BP)
}
cat("  Found", nrow(df.leadSNP), "lead SNPs\n")

window_bp <- window_mb * 1000000
df.leadSNP <- df.leadSNP %>%
  mutate(LOWER = ifelse((BP - window_bp) < 0, 1, BP - window_bp),
         UPPER = BP + window_bp)

# Create output directories
out_base <- file.path(FM_DIR, "output", sumstats_name, paste0(ld_pop, "_", window_mb, "Mb"))
dir.create(file.path(out_base, "ss"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "susie"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_base, "CARMA"), recursive = TRUE, showWarnings = FALSE)

# Read full summary stats once
cat("  Loading summary stats...\n")
ss <- fread(path_to_sumstats, na.strings = c("", " ", "NA"))
cat("  Loaded", nrow(ss), "variants\n")

# Extract per-locus summary stats
for (i in 1:nrow(df.leadSNP)) {
  LOCUS <- as.character(df.leadSNP$locus[i])
  CHR <- df.leadSNP$CHR[i]
  START <- df.leadSNP$LOWER[i]
  END <- df.leadSNP$UPPER[i]
  cat(paste0("\nLocus ", i, "/", nrow(df.leadSNP), ": ", LOCUS,
             " (chr", CHR, ":", START, "-", END, ")\n"))
  save_ss_per_locus(sumstats_name, ss, LOCUS, CHR, START, END, window_mb, ld_pop, FM_DIR)
}

cat("\n=== Done. Per-locus summary stats saved to", file.path(out_base, "ss"), "===\n")
