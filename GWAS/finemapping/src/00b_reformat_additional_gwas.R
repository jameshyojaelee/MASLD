#!/usr/bin/env Rscript
# 00b_reformat_additional_gwas.R
# Reformats 4 additional GWAS missed in the initial run:
#   1. Pazoki 2022 PDFF (GWAS Catalog format, hg38)
#   2. Anstee 2020 NAFLD (GWAS Catalog format, hg38, uses lnOR not beta)
#   3. Ghouse Cirrhosis (GWAS Catalog format, hg38)
#   4. Ghouse HCC (has pos_hg19, like BBJ)
#
# NOTE: Chen 2023 EXCLUDED — only has direction (beta=±1, se=1), not real effect sizes.

library(data.table)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
OUT_DIR <- file.path(FM_DIR, "data/sumstats")
GWAS_DIR <- file.path(BASE_DIR, "GWAS/MR_Data")

CHAIN_FILE <- file.path(BASE_DIR, "data/broadaway_eqtl/hg38ToHg19.over.chain")
chain <- import.chain(CHAIN_FILE)
cat("Loaded liftover chain:", length(chain), "chains\n")

liftover_to_hg19 <- function(dt) {
  gr <- GRanges(seqnames = paste0("chr", dt$chromosome),
                ranges = IRanges(start = dt$position, width = 1))
  lifted <- liftOver(gr, chain)
  keep <- lengths(lifted) == 1
  cat("  Liftover:", sum(keep), "/", length(keep), "variants mapped (1:1)\n")
  dt_out <- dt[keep, ]
  lifted_pos <- unlist(lifted[keep])
  dt_out$position <- start(lifted_pos)
  return(dt_out)
}

cat("\n============================================================\n")
cat("Reformatting 4 additional GWAS\n")
cat("============================================================\n")

# --- 1. Pazoki 2022 PDFF (GWAS Catalog format) ---
cat("\n=== Pazoki 2022 PDFF ===\n")
dt <- fread(file.path(GWAS_DIR, "GCST90267352_PDFF_Pazoki2022.tsv.gz"))
cat("  Loaded", nrow(dt), "variants\n")
dt_fmt <- dt %>%
  select(chromosome, position = base_pair_location,
         allele1 = effect_allele, allele2 = other_allele,
         beta, se = standard_error, pval = p_value) %>%
  filter(!is.na(beta), !is.na(se), !is.na(pval),
         abs(beta) < Inf, abs(beta) > 0,
         nchar(allele1) >= 1, nchar(allele2) >= 1) %>%
  mutate(chromosome = as.integer(chromosome)) %>%
  filter(!is.na(chromosome), chromosome >= 1, chromosome <= 22)
cat("  After QC:", nrow(dt_fmt), "variants\n")
dt_hg19 <- liftover_to_hg19(dt_fmt)
out_path <- file.path(OUT_DIR, "Pazoki_PDFF_reformatted_hg19.tsv")
fwrite(dt_hg19, out_path, sep = "\t")
cat("  Written:", out_path, "(", nrow(dt_hg19), "variants)\n")

# --- 2. Anstee 2020 NAFLD (uses lnOR, not beta) ---
cat("\n=== Anstee 2020 NAFLD ===\n")
dt <- fread(file.path(GWAS_DIR, "Anstee_2020_MASLD_GWAS_harmonised.tsv.gz"))
cat("  Loaded", nrow(dt), "variants\n")
dt_fmt <- dt %>%
  select(chromosome, position = base_pair_location,
         allele1 = effect_allele, allele2 = other_allele,
         beta = lnOR, se = standard_error, pval = p_value) %>%
  filter(!is.na(beta), !is.na(se), !is.na(pval),
         abs(beta) < Inf, abs(beta) > 0, se > 0,
         nchar(allele1) >= 1, nchar(allele2) >= 1) %>%
  mutate(chromosome = as.integer(chromosome)) %>%
  filter(!is.na(chromosome), chromosome >= 1, chromosome <= 22)
cat("  After QC:", nrow(dt_fmt), "variants\n")
dt_hg19 <- liftover_to_hg19(dt_fmt)
out_path <- file.path(OUT_DIR, "Anstee_NAFLD_reformatted_hg19.tsv")
fwrite(dt_hg19, out_path, sep = "\t")
cat("  Written:", out_path, "(", nrow(dt_hg19), "variants)\n")

# --- 3. Ghouse Cirrhosis (GWAS Catalog format) ---
cat("\n=== Ghouse Cirrhosis ===\n")
dt <- fread(file.path(GWAS_DIR, "Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz"))
cat("  Loaded", nrow(dt), "variants\n")
dt_fmt <- dt %>%
  select(chromosome, position = base_pair_location,
         allele1 = effect_allele, allele2 = other_allele,
         beta, se = standard_error, pval = p_value) %>%
  filter(!is.na(beta), !is.na(se), !is.na(pval),
         abs(beta) < Inf, abs(beta) > 0, se > 0,
         nchar(allele1) >= 1, nchar(allele2) >= 1) %>%
  mutate(chromosome = as.integer(chromosome)) %>%
  filter(!is.na(chromosome), chromosome >= 1, chromosome <= 22)
cat("  After QC:", nrow(dt_fmt), "variants\n")
dt_hg19 <- liftover_to_hg19(dt_fmt)
out_path <- file.path(OUT_DIR, "Ghouse_Cirrhosis_reformatted_hg19.tsv")
fwrite(dt_hg19, out_path, sep = "\t")
cat("  Written:", out_path, "(", nrow(dt_hg19), "variants)\n")

# --- 4. Ghouse HCC (has pos_hg19, like BBJ) ---
cat("\n=== Ghouse HCC ===\n")
dt <- fread(file.path(GWAS_DIR, "Ghouse_HCC/GCST90809296_EUR_HCC_harmonised_hg38.tsv.gz"))
cat("  Loaded", nrow(dt), "variants\n")
dt_fmt <- dt %>%
  select(chromosome, position = pos_hg19,
         allele1 = effect_allele, allele2 = other_allele,
         beta, se = standard_error, pval = p_value) %>%
  filter(!is.na(beta), !is.na(se), !is.na(pval), !is.na(position),
         abs(beta) < Inf, abs(beta) > 0, se > 0,
         nchar(allele1) >= 1, nchar(allele2) >= 1) %>%
  mutate(chromosome = as.integer(chromosome)) %>%
  filter(!is.na(chromosome), chromosome >= 1, chromosome <= 22)
cat("  After QC:", nrow(dt_fmt), "variants\n")
out_path <- file.path(OUT_DIR, "Ghouse_HCC_reformatted_hg19.tsv")
fwrite(dt_fmt, out_path, sep = "\t")
cat("  Written:", out_path, "(", nrow(dt_fmt), "variants)\n")

cat("\n============================================================\n")
cat("All 4 additional GWAS reformatted\n")
cat("NOTE: Chen 2023 excluded — only has direction (beta=±1), not real effect sizes\n")
cat("============================================================\n")
