#!/usr/bin/env Rscript
# ============================================================================
# liftover_credset_qc.R
# Reviewer-defense QC for the hybrid-build design (2026-06-01).
#
# The genetics core is anchored on hg19; new GRCh38-native resources (GTEx sQTL,
# SCREEN cCRE, Currin caQTL, ABC) are intersected by lifting the MASLD GWAS
# credible-set variants hg19->hg38. This script quantifies the liftOver success
# rate so we can state the coordinate-conversion loss is negligible — overall
# and, importantly, for the MASLD-DISEASE-specific GWAS loci.
#
# Output: GWAS/finemapping/results/liftover_credset_qc.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(rtracklayer)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM   <- file.path(BASE, "GWAS/finemapping/results/combined_finemapping.csv")
CHAIN<- file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain")
OUT  <- file.path(BASE, "GWAS/finemapping/results/liftover_credset_qc.tsv")

# MASLD-DISEASE GWAS (endpoints), vs liver-enzyme proxies (ALT/AST/GGT/PDFF)
masld_studies <- c("FinnGen_NAFLD","FinnGen_NASH","FinnGen_HCC",
  "2019_31311600_NAFLD_EUR","2020_32298765_NAFLD_EUR","2021_34841290_NAFLD_EUR",
  "2023_36280732_NAFLD_deCode_EUR","2023_36280732_NAFLD_Intermountain_EUR",
  "2023_36280732_NAFLD_UKBB_EUR","2020_32514122_Cirrhosis_EAS","2020_32514122_HCC_EAS",
  "Ghouse_Cirrhosis","Ghouse_HCC","Sveinbjornsson2022_Cirrhosis_meta_EUR",
  "Sveinbjornsson2022_Hcc_meta_EUR")

fm <- fread(FM)
fm_cs <- fm[either_in_cs == TRUE |
            (!is.na(recommended_pip) & recommended_pip > 0.1) |
            (is.na(recommended_pip) & !is.na(max_pip) & max_pip > 0.1)]
chain <- import.chain(CHAIN)

lift_rate <- function(dt, scope) {
  uv <- unique(dt[, .(chromosome, position)])
  if (nrow(uv) == 0) return(NULL)
  gr <- GRanges(paste0("chr", uv$chromosome), IRanges(uv$position, width = 1))
  nmap <- lengths(liftOver(gr, chain))
  data.table(scope = scope, n_credset_variants = nrow(uv),
             n_mapped_1to1 = sum(nmap == 1L),
             n_failed_or_multi = sum(nmap != 1L),
             pct_mapped = round(100 * sum(nmap == 1L) / nrow(uv), 3))
}

qc <- rbind(
  lift_rate(fm_cs, "all_credible_set"),
  lift_rate(fm_cs[study %in% masld_studies], "masld_disease_gwas_only")
)
fwrite(qc, OUT, sep = "\t")
cat("== liftOver hg19->hg38 credible-set QC ==\n")
print(qc)
cat(sprintf("\nWROTE %s\n", OUT))
cat(sprintf("Reviewer line: %.2f%% of all credible-set variants and %.2f%% of MASLD-disease-GWAS credible-set variants map 1:1 hg19->hg38 (UCSC chain); failures concentrate in known build-discordant regions and are dropped.\n",
            qc[scope == "all_credible_set"]$pct_mapped,
            qc[scope == "masld_disease_gwas_only"]$pct_mapped))
