#!/usr/bin/env Rscript
# format_finngen_sex_stratified.R
# ---------------------------------------------------------------------------
# Format FinnGen R12 sex-stratified summary statistics for the SuSiE-COLOC
# pipeline. Mirrors GWAS/finemapping/src/00_reformat_gwas.R (function
# reformat_finngen) to keep column conventions identical.
#
# STATUS (2026-05-11): script is staged but does nothing useful until
# FinnGen sex-stratified gz files actually land in
#   GWAS/MR_Data/FinnGen/sex_stratified/finngen_R12_{NAFLD,NASH,HCC}_{F,M}.gz
# (FinnGen R12 public release does not include these; see B3_blocked.md).
#
# Inputs : GWAS/MR_Data/FinnGen/sex_stratified/finngen_R12_<LABEL>_<SEX>.gz
#                                 LABEL in {NAFLD, NASH, HCC}
#                                 SEX   in {F, M}
# Outputs: GWAS/finemapping/data/sumstats/FinnGen_<LABEL>_<SEX>_reformatted_hg19.tsv
#         + appended row in GWAS/finemapping/config/gwas_registry.tsv
#
# Usage: Rscript format_finngen_sex_stratified.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(rtracklayer)
  library(GenomicRanges)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR     <- file.path(BASE_DIR, "GWAS/finemapping")
SUMSTAT_DIR <- file.path(FM_DIR, "data/sumstats")
SRC_DIR    <- file.path(BASE_DIR, "GWAS/MR_Data/FinnGen/sex_stratified")
REGISTRY   <- file.path(FM_DIR, "config/gwas_registry.tsv")

CHAIN_FILE <- file.path(BASE_DIR, "data/broadaway_eqtl/hg38ToHg19.over.chain")
if (!file.exists(CHAIN_FILE)) {
  stop("Missing liftover chain: ", CHAIN_FILE)
}
chain <- import.chain(CHAIN_FILE)

# Sample sizes per sex stratum (placeholder until FinnGen publishes counts).
# FinnGen R12 pooled-sex NAFLD: 3,504 cases / 496,844 controls (n=500,348).
# Sex split in FinnGen is roughly F:M ~ 56:44; per-stratum numbers below are
# best-guess estimates derived from the pooled manifest and FinnGen R12 cohort
# composition. REPLACE with the per-stratum N from the sandbox export when
# available.
strata_n <- list(
  NAFLD = list(F = list(N_cases =  2200L, N_controls = 282000L),
               M = list(N_cases =  1300L, N_controls = 215000L)),
  NASH  = list(F = list(N_cases =  1050L, N_controls = 280000L),
               M = list(N_cases =   750L, N_controls = 213000L)),
  HCC   = list(F = list(N_cases =   300L, N_controls = 232000L),
               M = list(N_cases =   650L, N_controls = 168000L))
)

liftover_to_hg19 <- function(dt) {
  gr <- GRanges(seqnames = paste0("chr", dt$chromosome),
                ranges   = IRanges(start = dt$position, width = 1))
  lifted <- liftOver(gr, chain)
  keep   <- lengths(lifted) == 1
  cat("  Liftover:", sum(keep), "/", length(keep), "variants mapped (1:1)\n")
  dt_out <- dt[keep, ]
  dt_out$position <- start(unlist(lifted[keep]))
  dt_out
}

reformat_one <- function(label, sex_short) {
  src <- file.path(SRC_DIR, paste0("finngen_R12_", label, "_", sex_short, ".gz"))
  if (!file.exists(src)) {
    cat("  SKIP — missing:", src, "\n")
    return(invisible(NULL))
  }
  cat("\n=== Reformatting FinnGen ", label, " (", sex_short, ") ===\n", sep = "")
  dt <- fread(src)
  cat("  Loaded", nrow(dt), "variants\n")

  chrom_col <- intersect(c("#chrom", "chrom", "X.chrom"), colnames(dt))[1]
  if (is.na(chrom_col)) stop("Cannot find chrom column")
  setnames(dt, chrom_col, "chromosome")

  dt_fmt <- dt %>%
    select(chromosome, position = pos,
           allele1 = alt, allele2 = ref,  # alt = effect allele (FinnGen)
           beta, se = sebeta, pval) %>%
    filter(!is.na(beta), !is.na(se), !is.na(pval), abs(beta) < Inf,
           nchar(allele1) >= 1, nchar(allele2) >= 1) %>%
    mutate(chromosome = as.integer(chromosome)) %>%
    filter(!is.na(chromosome), chromosome >= 1, chromosome <= 22)
  cat("  After QC:", nrow(dt_fmt), "variants\n")

  dt_hg19 <- liftover_to_hg19(dt_fmt)

  out_name <- paste0("FinnGen_", label, "_", sex_short, "_reformatted_hg19.tsv")
  out_path <- file.path(SUMSTAT_DIR, out_name)
  fwrite(dt_hg19, out_path, sep = "\t")
  cat("  Wrote:", out_path, "(", nrow(dt_hg19), "rows )\n")

  # Append to GWAS registry (idempotent)
  reg <- read.delim(REGISTRY, stringsAsFactors = FALSE)
  study_name <- paste0("FINNGEN_", sex_short, "_", label)
  if (study_name %in% reg$study_name) {
    cat("  Registry row already exists:", study_name, "\n")
    return(invisible(NULL))
  }
  s <- strata_n[[label]][[sex_short]]
  new_row <- data.frame(
    study_name    = study_name,
    sumstats_path = file.path("data/sumstats", out_name),
    leadsnps_path = file.path("data/lead_snps",
                              paste0(study_name, "_leadSNPs.tsv")),
    ancestry      = "EUR",
    trait_type    = "binary",
    N_tot         = s$N_cases + s$N_controls,
    N_cases       = s$N_cases,
    ld_panel      = "polyfun",
    window_mb     = 0.5,
    stringsAsFactors = FALSE
  )
  reg2 <- rbind(reg, new_row)
  write.table(reg2, REGISTRY, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Appended registry row:", study_name, "\n")
}

cat("============================================================\n")
cat("FinnGen R12 sex-stratified → finemapping pipeline format\n")
cat("============================================================\n")

for (label in c("NAFLD", "NASH", "HCC")) {
  for (sex_short in c("F", "M")) {
    reformat_one(label, sex_short)
  }
}

cat("\n============================================================\n")
cat("Done. New registry entries (if files were present):\n")
cat("  FINNGEN_F_NAFLD, FINNGEN_M_NAFLD,\n")
cat("  FINNGEN_F_NASH,  FINNGEN_M_NASH,\n")
cat("  FINNGEN_F_HCC,   FINNGEN_M_HCC\n")
cat("Next: GWAS_NAME=FINNGEN_F_NAFLD sbatch GWAS/finemapping/src/06_susie_coloc_bigmem.sh --array=1-22\n")
cat("============================================================\n")
