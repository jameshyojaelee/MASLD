#!/usr/bin/env Rscript
# build_mvp_registry.R
# After download: parse MVP metadata (case/control N), build the harmonization
# manifest + sample-size table, and append 27 per-ancestry strata to the
# finemapping GWAS registry. META strata are recorded for reference but NOT
# registered for COLOC (mixed-ancestry LD invalid for SuSiE).
#
# Strata: phecodes 571.5/571.51/571 (binary; 571.8/571.81 dropped 2026-07-01)
#         x present ancestries
#         + liver biomarkers ALT/AST/Albumin/Platelet (quantitative) x EUR/AFR/AMR/EAS

suppressPackageStartupMessages(library(data.table))

BASE   <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
MVP    <- file.path(BASE, "GWAS/MR_Data/MVP")
PHEDIR <- file.path(MVP, "raw/phecodes")
LABDIR <- file.path(MVP, "raw/labs")
FM     <- file.path(BASE, "GWAS/finemapping")
REG    <- file.path(FM, "config/gwas_registry.tsv")
SUMDIR_REL <- "data/sumstats"   # relative to FM (registry convention)
LEAD_REL   <- "data/lead_snps"

# friendly phenotype names
# Phe_571_8 (liver abscess/sequelae) + Phe_571_81 (portal hypertension) were
# DROPPED 2026-07-01 as off-target/non-MASLD (registry 55->50 GWAS, 27 MVP
# strata). They are intentionally ABSENT from phe_name so the loop at L43 skips
# them -- do NOT re-add or the purge silently reverts on any registry rebuild.
phe_name <- c("Phe_571_5"="NAFLD", "Phe_571_51"="Cirrhosis", "Phe_571"="ChronLiver")
ld_for <- c(EUR="polyfun", EAS="1kg_eas", AFR="1kg_afr", AMR="1kg_amr")

parse_meta <- function(metafile) {
  if (!file.exists(metafile)) return(list(N=NA_integer_, cases=NA_integer_, controls=NA_integer_))
  txt <- paste(readLines(metafile, warn=FALSE), collapse=" ")
  num <- function(pat) { m <- regmatches(txt, regexpr(pat, txt, perl=TRUE)); if (length(m)) as.integer(gsub("[^0-9]","",m)) else NA_integer_ }
  list(N        = num("Total Sample Size=\\d+"),
       cases    = num("Cases=\\d+"),
       controls = num("Controls=\\d+"))
}

rows <- list()
# ---- phecodes (binary) ----
for (gz in list.files(PHEDIR, pattern="\\.txt\\.gz$", full.names=TRUE)) {
  bn <- basename(gz)
  f  <- strsplit(bn, ".", fixed=TRUE)[[1]]      # MVP_R4 . 1000G_AGR . Phe_571_5 . EUR . GIA . dbGaP . txt . gz
  phe <- f[3]; anc <- f[4]
  if (anc == "META" || !(anc %in% names(ld_for))) next
  if (!(phe %in% names(phe_name))) next
  meta <- parse_meta(sub("\\.txt\\.gz$", ".metadata.txt", gz))
  nm  <- sprintf("MVP_%s_%s", phe_name[[phe]], anc)
  rows[[length(rows)+1]] <- data.table(
    gwas_name=nm, input_gz=gz,
    out_sumstats=file.path(SUMDIR_REL, paste0(nm, "_reformatted_hg19.tsv")),
    leadsnps=file.path(LEAD_REL, paste0(nm, "_leadSNPs.tsv")),
    ancestry=anc, trait_type="binary",
    N_tot=meta$N, N_cases=meta$cases, ld_panel=ld_for[[anc]])
}
# ---- liver biomarkers (quantitative) ----
for (gz in list.files(LABDIR, pattern="\\.txt\\.gz$", full.names=TRUE)) {
  bn <- basename(gz)
  f  <- strsplit(bn, ".", fixed=TRUE)[[1]]      # MVP_R4 . 1000G_AGR . ALT_Mean_INT . EUR . GIA . dbGaP . txt . gz
  analyte <- sub("_Mean_INT$", "", f[3]); anc <- f[4]
  if (anc == "META" || !(anc %in% names(ld_for))) next
  meta <- parse_meta(sub("\\.txt\\.gz$", ".metadata.txt", gz))
  nm  <- sprintf("MVP_%s_%s", analyte, anc)
  rows[[length(rows)+1]] <- data.table(
    gwas_name=nm, input_gz=gz,
    out_sumstats=file.path(SUMDIR_REL, paste0(nm, "_reformatted_hg19.tsv")),
    leadsnps=file.path(LEAD_REL, paste0(nm, "_leadSNPs.tsv")),
    ancestry=anc, trait_type="quantitative",
    N_tot=meta$N, N_cases=0L, ld_panel=ld_for[[anc]])
}
man <- rbindlist(rows)
setorder(man, gwas_name)
cat("MVP strata for COLOC:", nrow(man), "\n")
print(man[, .(gwas_name, ancestry, trait_type, N_tot, N_cases)])

# ---- write manifest + sample-size table ----
fwrite(man, file.path(MVP, "mvp_manifest.tsv"), sep="\t")
ss <- man[, .(study=gwas_name, ancestry, trait_type, N_tot, N_cases,
              N_controls=ifelse(trait_type=="binary", N_tot-N_cases, NA_integer_),
              s=ifelse(trait_type=="binary", round(N_cases/N_tot,5), NA_real_))]
fwrite(ss, file.path(MVP, "mvp_sample_sizes.tsv"), sep="\t")

# ---- append to registry (idempotent) ----
reg <- fread(REG, sep="\t")
reg_new <- man[, .(study_name=gwas_name, sumstats_path=out_sumstats, leadsnps_path=leadsnps,
                   ancestry, trait_type, N_tot, N_cases, ld_panel, window_mb=0.5)]
reg2 <- reg[!grepl("^MVP_", study_name)]          # drop any prior MVP rows
reg2 <- rbind(reg2, reg_new, fill=TRUE)
fwrite(reg2, REG, sep="\t")
cat("Registry rows:", nrow(reg), "->", nrow(reg2), "(+", nrow(reg_new), "MVP)\n")
cat("DONE build_mvp_registry\n")
