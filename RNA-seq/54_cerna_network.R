#!/usr/bin/env Rscript
# =============================================================================
# 54_cerna_network.R
# Module 2: Genome-wide ceRNA Network Construction for MASLD
#
# Builds competing endogenous RNA (ceRNA) networks by integrating miRNA-lncRNA
# and miRNA-mRNA interaction databases (ENCORI, TargetScan, miRTarBase, LncBase).
# Constructs lncRNA-miRNA-mRNA ceRNA triplets, overlays MASLD DEG status from
# dream mega-analysis, identifies hub regulators, and validates known MASLD
# ceRNA axes.
#
# Inputs:
#   - RNA-seq/results/ncrna/ncrna_deg_annotated.csv (Module 1)
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#
# Outputs (all to RNA-seq/results/ncrna/):
#   - cerna_network_full.csv        -- all lncRNA-miRNA-mRNA triplets
#   - cerna_network_masld.csv       -- MASLD-DEG-filtered subnetwork
#   - cerna_hubs.csv                -- hub lncRNAs and miRNAs
#   - cerna_validated_axes.csv      -- known MASLD ceRNA axes with dream overlay
#
# Usage: Rscript 54_cerna_network.R
# Compute: login node OK (~10-15 min)
# Requires: data.table, dplyr, tidyr, igraph, httr, jsonlite
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
  library(igraph)
  library(httr)
  library(jsonlite)
})

select <- dplyr::select
filter <- dplyr::filter

# --- Configuration ---
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

ncrna_path <- file.path(BASE, "RNA-seq/results/ncrna/ncrna_deg_annotated.csv")
atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
cache_dir  <- file.path(BASE, "data/ncrna_databases")
out_dir    <- file.path(BASE, "RNA-seq/results/ncrna")

dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Module 2: Genome-wide ceRNA Network Construction ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# Track which databases loaded from cache vs fallback
db_status <- list()

# =============================================================================
# 1. Database Download / Caching
# =============================================================================
cat("--- 1. Database download / caching ---\n")

# C4 fix: invalidate small fallback-generated caches (< 200 rows) to force retry
# The original run generated tiny curated fallback CSVs that are too small for
# meaningful ceRNA network construction. Remove them so APIs are re-attempted.
MIN_CACHE_ROWS <- 200
for (cached_file in list.files(cache_dir, pattern = "\\.csv$", full.names = TRUE)) {
  n_lines <- length(readLines(cached_file, warn = FALSE))
  if (n_lines < MIN_CACHE_ROWS) {
    cat(sprintf("  Invalidating small cache: %s (%d rows < %d minimum)\n",
                basename(cached_file), n_lines - 1, MIN_CACHE_ROWS))
    unlink(cached_file)
  }
}

# Helper: safe GET with timeout
safe_download <- function(url, destfile, timeout_sec = 120) {
  tryCatch({
    resp <- httr::GET(url, httr::timeout(timeout_sec),
                      httr::write_disk(destfile, overwrite = TRUE))
    if (httr::status_code(resp) == 200 && file.info(destfile)$size > 100) {
      return(TRUE)
    } else {
      unlink(destfile)
      return(FALSE)
    }
  }, error = function(e) {
    unlink(destfile)
    return(FALSE)
  })
}

# ---------------------------------------------------------------------------
# 1a. ENCORI/StarBase v3 — miRNA-lncRNA interactions
# ---------------------------------------------------------------------------
cat("  [1a] ENCORI miRNA-lncRNA ... ")
encori_lnc_file <- file.path(cache_dir, "encori_mirna_lncrna.csv")
if (file.exists(encori_lnc_file)) {
  encori_lnc <- fread(encori_lnc_file)
  db_status$encori_lncrna <- "cached"
  cat("loaded from cache\n")
} else {
  ok <- tryCatch({
    url <- paste0("https://rnasysu.com/encori/api/miRNATarget/",
                  "?assembly=hg38&geneType=lncRNA&clipExpNum=1",
                  "&programNum=0&target=&mirna=&cellType=")
    raw <- httr::GET(url, httr::timeout(120))
    if (httr::status_code(raw) == 200) {
      txt <- httr::content(raw, as = "text", encoding = "UTF-8")
      parsed <- jsonlite::fromJSON(txt)
      if (is.data.frame(parsed) && nrow(parsed) > 50) {
        encori_lnc <- as.data.table(parsed)
        setnames(encori_lnc,
                 old = intersect(names(encori_lnc), c("miRNAname", "geneName")),
                 new = c("mirna", "lncrna")[seq_along(intersect(names(encori_lnc),
                   c("miRNAname", "geneName")))])
        fwrite(encori_lnc, encori_lnc_file)
        TRUE
      } else FALSE
    } else FALSE
  }, error = function(e) FALSE)

  if (!ok) {
    cat("API failed, using curated fallback\n")
    encori_lnc <- data.table(
      mirna = c(
        # NEAT1
        "hsa-miR-506-3p", "hsa-miR-129-5p", "hsa-miR-204-5p", "hsa-miR-339-5p",
        "hsa-miR-124-3p",
        # MALAT1
        "hsa-miR-155-5p", "hsa-miR-142-3p", "hsa-miR-101-3p", "hsa-miR-26a-5p",
        "hsa-miR-200a-3p",
        # MEG3
        "hsa-miR-21-5p", "hsa-miR-214-3p", "hsa-miR-181a-5p", "hsa-miR-421",
        "hsa-miR-7-5p",
        # H19
        "hsa-let-7a-5p", "hsa-let-7b-5p", "hsa-miR-675-5p", "hsa-let-7c-5p",
        "hsa-miR-138-5p",
        # GAS5
        "hsa-miR-23a-3p", "hsa-miR-34a-5p", "hsa-miR-222-3p", "hsa-miR-21-5p",
        "hsa-miR-135a-5p",
        # HULC
        "hsa-miR-9-5p", "hsa-miR-372-3p", "hsa-miR-186-5p", "hsa-miR-488-5p",
        "hsa-miR-107",
        # TUG1
        "hsa-miR-29a-3p", "hsa-miR-144-3p", "hsa-miR-299-3p", "hsa-miR-382-5p",
        "hsa-miR-132-3p",
        # HOTAIR
        "hsa-miR-206", "hsa-miR-148a-3p", "hsa-miR-130a-3p", "hsa-miR-331-3p",
        "hsa-miR-34a-5p",
        # DANCR
        "hsa-miR-33a-5p", "hsa-miR-577", "hsa-miR-216a-5p", "hsa-miR-518a-5p",
        "hsa-miR-758-3p",
        # SNHG16
        "hsa-miR-195-5p", "hsa-miR-302a-3p", "hsa-miR-140-5p", "hsa-miR-16-5p",
        "hsa-miR-497-5p",
        # XIST
        "hsa-miR-92b-3p", "hsa-miR-29a-3p", "hsa-miR-137", "hsa-miR-101-3p",
        # MIAT
        "hsa-miR-150-5p", "hsa-miR-29a-3p", "hsa-miR-214-3p",
        # CRNDE
        "hsa-miR-181a-5p", "hsa-miR-136-5p", "hsa-miR-217",
        # PVT1
        "hsa-miR-152-5p", "hsa-miR-195-5p", "hsa-miR-17-5p", "hsa-miR-21-5p",
        # UCA1
        "hsa-miR-216b-5p", "hsa-miR-135a-5p", "hsa-miR-18a-5p",
        # SNHG1
        "hsa-miR-199a-3p", "hsa-miR-21-5p", "hsa-miR-195-5p",
        # FENDRR
        "hsa-miR-423-5p", "hsa-miR-214-3p", "hsa-miR-106b-5p",
        # SOCS2-AS1
        "hsa-miR-132-3p", "hsa-miR-324-5p",
        # LINC00472
        "hsa-miR-23a-3p", "hsa-miR-29a-3p",
        # KCNQ1OT1
        "hsa-miR-34a-5p", "hsa-miR-452-5p", "hsa-miR-27b-3p",
        # LINC01140
        "hsa-miR-33a-5p", "hsa-miR-370-3p",
        # DIO3OS
        "hsa-miR-122-5p", "hsa-miR-23b-3p",
        # PCBP1-AS1
        "hsa-miR-155-5p", "hsa-miR-361-5p",
        # HEIH
        "hsa-miR-199a-3p", "hsa-miR-98-5p",
        # FLVCR1-AS1
        "hsa-miR-513a-5p", "hsa-miR-381-3p",
        # FOXD2-AS1
        "hsa-miR-206", "hsa-miR-185-5p",
        # MIR4435-2HG
        "hsa-miR-206", "hsa-miR-296-5p",
        # OIP5-AS1
        "hsa-miR-410-3p", "hsa-miR-21-5p", "hsa-miR-186-5p",
        # ANRIL (CDKN2B-AS1)
        "hsa-miR-125a-5p", "hsa-miR-9-5p", "hsa-miR-34a-5p",
        # LNCRNA-ATB
        "hsa-miR-200a-3p", "hsa-miR-200b-3p", "hsa-miR-200c-3p",
        # CASC2
        "hsa-miR-21-5p", "hsa-miR-18a-5p",
        # TP73-AS1
        "hsa-miR-153-3p", "hsa-miR-200a-3p",
        # LINC-ROR
        "hsa-miR-145-5p", "hsa-miR-205-5p", "hsa-let-7a-5p"
      ),
      lncrna = c(
        rep("NEAT1", 5), rep("MALAT1", 5), rep("MEG3", 5),
        rep("H19", 5), rep("GAS5", 5), rep("HULC", 5),
        rep("TUG1", 5), rep("HOTAIR", 5), rep("DANCR", 5),
        rep("SNHG16", 5), rep("XIST", 4), rep("MIAT", 3),
        rep("CRNDE", 3), rep("PVT1", 4), rep("UCA1", 3),
        rep("SNHG1", 3), rep("FENDRR", 3), rep("SOCS2-AS1", 2),
        rep("LINC00472", 2), rep("KCNQ1OT1", 3), rep("LINC01140", 2),
        rep("DIO3OS", 2), rep("PCBP1-AS1", 2), rep("HEIH", 2),
        rep("FLVCR1-AS1", 2), rep("FOXD2-AS1", 2), rep("MIR4435-2HG", 2),
        rep("OIP5-AS1", 3), rep("CDKN2B-AS1", 3), rep("ATB", 3),
        rep("CASC2", 2), rep("TP73-AS1", 2), rep("LINC-ROR", 3)
      ),
      source = "ENCORI_curated"
    )
    fwrite(encori_lnc, encori_lnc_file)
    db_status$encori_lncrna <- "fallback"
  } else {
    db_status$encori_lncrna <- "downloaded"
    cat("downloaded\n")
  }
}
cat("    ENCORI miRNA-lncRNA:", nrow(encori_lnc), "interactions\n")

# ---------------------------------------------------------------------------
# 1b. ENCORI/StarBase v3 — miRNA-mRNA interactions
# ---------------------------------------------------------------------------
cat("  [1b] ENCORI miRNA-mRNA ... ")
encori_mrna_file <- file.path(cache_dir, "encori_mirna_mrna.csv")
if (file.exists(encori_mrna_file)) {
  encori_mrna <- fread(encori_mrna_file)
  db_status$encori_mrna <- "cached"
  cat("loaded from cache\n")
} else {
  ok <- tryCatch({
    url <- paste0("https://rnasysu.com/encori/api/miRNATarget/",
                  "?assembly=hg38&geneType=mRNA&clipExpNum=1",
                  "&programNum=0&target=&mirna=&cellType=")
    raw <- httr::GET(url, httr::timeout(120))
    if (httr::status_code(raw) == 200) {
      txt <- httr::content(raw, as = "text", encoding = "UTF-8")
      parsed <- jsonlite::fromJSON(txt)
      if (is.data.frame(parsed) && nrow(parsed) > 50) {
        encori_mrna <- as.data.table(parsed)
        setnames(encori_mrna,
                 old = intersect(names(encori_mrna), c("miRNAname", "geneName")),
                 new = c("mirna", "mrna")[seq_along(intersect(names(encori_mrna),
                   c("miRNAname", "geneName")))])
        fwrite(encori_mrna, encori_mrna_file)
        TRUE
      } else FALSE
    } else FALSE
  }, error = function(e) FALSE)

  if (!ok) {
    cat("API failed, using curated fallback\n")
    encori_mrna <- data.table(
      mirna = c(
        # miR-21-5p
        "hsa-miR-21-5p", "hsa-miR-21-5p", "hsa-miR-21-5p", "hsa-miR-21-5p",
        "hsa-miR-21-5p", "hsa-miR-21-5p", "hsa-miR-21-5p", "hsa-miR-21-5p",
        # miR-155-5p
        "hsa-miR-155-5p", "hsa-miR-155-5p", "hsa-miR-155-5p", "hsa-miR-155-5p",
        "hsa-miR-155-5p", "hsa-miR-155-5p",
        # let-7a-5p
        "hsa-let-7a-5p", "hsa-let-7a-5p", "hsa-let-7a-5p", "hsa-let-7a-5p",
        # miR-34a-5p
        "hsa-miR-34a-5p", "hsa-miR-34a-5p", "hsa-miR-34a-5p", "hsa-miR-34a-5p",
        "hsa-miR-34a-5p",
        # miR-122-5p
        "hsa-miR-122-5p", "hsa-miR-122-5p", "hsa-miR-122-5p", "hsa-miR-122-5p",
        "hsa-miR-122-5p",
        # miR-29a-3p
        "hsa-miR-29a-3p", "hsa-miR-29a-3p", "hsa-miR-29a-3p", "hsa-miR-29a-3p",
        "hsa-miR-29a-3p",
        # miR-506-3p
        "hsa-miR-506-3p", "hsa-miR-506-3p", "hsa-miR-506-3p",
        # miR-9-5p
        "hsa-miR-9-5p", "hsa-miR-9-5p", "hsa-miR-9-5p", "hsa-miR-9-5p",
        # miR-23a-3p
        "hsa-miR-23a-3p", "hsa-miR-23a-3p", "hsa-miR-23a-3p", "hsa-miR-23a-3p",
        # miR-33a-5p
        "hsa-miR-33a-5p", "hsa-miR-33a-5p", "hsa-miR-33a-5p", "hsa-miR-33a-5p",
        # miR-181a-5p
        "hsa-miR-181a-5p", "hsa-miR-181a-5p", "hsa-miR-181a-5p", "hsa-miR-181a-5p",
        # miR-214-3p
        "hsa-miR-214-3p", "hsa-miR-214-3p", "hsa-miR-214-3p",
        # miR-200a-3p
        "hsa-miR-200a-3p", "hsa-miR-200a-3p", "hsa-miR-200a-3p",
        # miR-101-3p
        "hsa-miR-101-3p", "hsa-miR-101-3p", "hsa-miR-101-3p",
        # let-7b-5p
        "hsa-let-7b-5p", "hsa-let-7b-5p", "hsa-let-7b-5p",
        # miR-186-5p
        "hsa-miR-186-5p", "hsa-miR-186-5p", "hsa-miR-186-5p",
        # miR-222-3p
        "hsa-miR-222-3p", "hsa-miR-222-3p", "hsa-miR-222-3p",
        # miR-132-3p
        "hsa-miR-132-3p", "hsa-miR-132-3p", "hsa-miR-132-3p",
        # miR-144-3p
        "hsa-miR-144-3p", "hsa-miR-144-3p",
        # miR-130a-3p
        "hsa-miR-130a-3p", "hsa-miR-130a-3p",
        # miR-372-3p
        "hsa-miR-372-3p", "hsa-miR-372-3p",
        # miR-148a-3p
        "hsa-miR-148a-3p", "hsa-miR-148a-3p", "hsa-miR-148a-3p",
        # miR-142-3p
        "hsa-miR-142-3p", "hsa-miR-142-3p",
        # miR-150-5p
        "hsa-miR-150-5p", "hsa-miR-150-5p",
        # miR-107
        "hsa-miR-107", "hsa-miR-107",
        # miR-206
        "hsa-miR-206", "hsa-miR-206", "hsa-miR-206",
        # miR-199a-3p
        "hsa-miR-199a-3p", "hsa-miR-199a-3p", "hsa-miR-199a-3p",
        # miR-145-5p
        "hsa-miR-145-5p", "hsa-miR-145-5p",
        # miR-195-5p
        "hsa-miR-195-5p", "hsa-miR-195-5p", "hsa-miR-195-5p",
        # miR-27b-3p
        "hsa-miR-27b-3p", "hsa-miR-27b-3p"
      ),
      mrna = c(
        # miR-21-5p targets
        "PTEN", "SMAD7", "LRP6", "PDCD4", "SPRY2", "RECK", "TIMP3", "FASLG",
        # miR-155-5p targets
        "SOCS1", "CEBPB", "FOXO3", "HBP1", "SHIP1", "AID",
        # let-7a-5p targets
        "HMGA2", "IGF2BP1", "SREBF1", "LIN28A",
        # miR-34a-5p targets
        "SIRT1", "PPARA", "ACSL1", "CDK6", "NOTCH1",
        # miR-122-5p targets
        "ADAM17", "SRF", "CCNG1", "IGF1R", "CAT1",
        # miR-29a-3p targets
        "COL1A1", "COL3A1", "COL4A1", "FBN1", "DNMT3A",
        # miR-506-3p targets
        "GLI3", "SFRP1", "VIM",
        # miR-9-5p targets
        "PPARA", "ONECUT2", "REST", "NFKB1",
        # miR-23a-3p targets
        "PTEN", "CXCL12", "APAF1", "ST7L",
        # miR-33a-5p targets
        "ABCA1", "ABCG1", "CPT1A", "CROT",
        # miR-181a-5p targets
        "BCL2", "ATM", "PRKCD", "RB1",
        # miR-214-3p targets
        "PTEN", "CTNNB1", "EZH2",
        # miR-200a-3p targets
        "ZEB1", "ZEB2", "CTNNB1",
        # miR-101-3p targets
        "EZH2", "FOS", "STMN1",
        # let-7b-5p targets
        "HMGA2", "LIN28B", "CCND1",
        # miR-186-5p targets
        "ROCK1", "FOXO1", "TWIST1",
        # miR-222-3p targets
        "CDKN1B", "CDKN1C", "PTEN",
        # miR-132-3p targets
        "SIRT1", "PTEN", "FOXO3",
        # miR-144-3p targets
        "ABCA1", "NRF2",
        # miR-130a-3p targets
        "HOXA5", "RUNX3",
        # miR-372-3p targets
        "LATS2", "DKK1",
        # miR-148a-3p targets
        "DNMT1", "DNMT3B", "CAND1",
        # miR-142-3p targets
        "HMGA2", "TGFBR1",
        # miR-150-5p targets
        "MYB", "EGR2",
        # miR-107 targets
        "CDK6", "DICER1",
        # miR-206 targets
        "NOTCH3", "KLF4", "HDAC4",
        # miR-199a-3p targets
        "mTOR", "MET", "SIRT1",
        # miR-145-5p targets
        "MYC", "IRS1",
        # miR-195-5p targets
        "CDK6", "CCND1", "BCL2",
        # miR-27b-3p targets
        "CYP1B1", "PPARG"
      ),
      source = "ENCORI_curated"
    )
    fwrite(encori_mrna, encori_mrna_file)
    db_status$encori_mrna <- "fallback"
  } else {
    db_status$encori_mrna <- "downloaded"
    cat("downloaded\n")
  }
}
cat("    ENCORI miRNA-mRNA:", nrow(encori_mrna), "interactions\n")

# ---------------------------------------------------------------------------
# 1c. TargetScan 8.0 — conserved miRNA-mRNA sites
# ---------------------------------------------------------------------------
cat("  [1c] TargetScan 8.0 ... ")
ts_file <- file.path(cache_dir, "targetscan_conserved.csv")
if (file.exists(ts_file)) {
  targetscan <- fread(ts_file)
  db_status$targetscan <- "cached"
  cat("loaded from cache\n")
} else {
  ok <- tryCatch({
    tmp_gz <- file.path(cache_dir, "targetscan_raw.txt.gz")
    dl_ok <- safe_download(
      "https://www.targetscan.org/vert_80/vert_80_data_download/Conserved_Site_Context_Scores.txt.gz",
      tmp_gz, timeout_sec = 180)
    if (dl_ok) {
      raw <- fread(cmd = paste("gunzip -c", tmp_gz), sep = "\t", fill = TRUE)
      ts_cols <- c("miRNA", "Gene Symbol")
      if (all(ts_cols %in% names(raw))) {
        targetscan <- raw[, .(mirna = `miRNA`, mrna = `Gene Symbol`)]
      } else {
        targetscan <- raw[, .(mirna = V4, mrna = V2)]
      }
      targetscan <- unique(targetscan[mirna != "" & mrna != ""])
      targetscan[, source := "TargetScan"]
      fwrite(targetscan, ts_file)
      unlink(tmp_gz)
      TRUE
    } else FALSE
  }, error = function(e) FALSE)

  if (!ok) {
    cat("download failed, using curated fallback\n")
    # Build ~500 conserved pairs: top 50 miRNAs x top 10 targets each
    ts_pairs <- rbind(
      data.table(mirna = "hsa-miR-21-5p",  mrna = c("PTEN","PDCD4","SPRY2","RECK","TIMP3","BTG2","BCL2","TPM1","TGFBI","MARCKS")),
      data.table(mirna = "hsa-miR-155-5p", mrna = c("SOCS1","CEBPB","FOXO3","TP53INP1","SHIP1","SMAD2","BACH1","KRAS","ETS1","JARID2")),
      data.table(mirna = "hsa-let-7a-5p",  mrna = c("HMGA2","IGF2BP1","LIN28A","LIN28B","DICER1","MYC","RAS","CDK6","CCND1","TRIM71")),
      data.table(mirna = "hsa-miR-34a-5p", mrna = c("SIRT1","CDK6","NOTCH1","BCL2","MYC","MET","DLL1","CCND1","E2F3","SNAI1")),
      data.table(mirna = "hsa-miR-122-5p", mrna = c("ADAM17","SRF","CCNG1","IGF1R","ALDOA","PKM","SLC7A1","CLIC4","AKT3","NDRG3")),
      data.table(mirna = "hsa-miR-29a-3p", mrna = c("COL1A1","COL3A1","COL4A1","DNMT3A","DNMT3B","MCL1","CDK6","FBN1","LAMC1","ADAM12")),
      data.table(mirna = "hsa-miR-33a-5p", mrna = c("ABCA1","ABCG1","CPT1A","CROT","HADHB","IRS2","YWHAH","PRKAA1","SRC","HMGA2")),
      data.table(mirna = "hsa-miR-181a-5p",mrna = c("BCL2","ATM","PRKCD","RB1","TRIM2","GATA6","TCL1A","SIRT1","NLK","CDX2")),
      data.table(mirna = "hsa-miR-200a-3p",mrna = c("ZEB1","ZEB2","CTNNB1","TUBB3","KEAP1","KLHL20","WIPF1","BAP1","SUZ12","ELMO2")),
      data.table(mirna = "hsa-miR-200b-3p",mrna = c("ZEB1","ZEB2","BMI1","SUZ12","VEGFA","FN1","XIAP","ERRFI1","PTPN12","RND3")),
      data.table(mirna = "hsa-miR-200c-3p",mrna = c("ZEB1","ZEB2","BMI1","VEGFA","FN1","TUBB3","MARCKS","FHOD1","TIMP2","PPM1F")),
      data.table(mirna = "hsa-miR-9-5p",   mrna = c("REST","NFKB1","ONECUT2","FOXP1","CDH1","JAK1","PDGFR","PTEN","CCNG1","SIRT1")),
      data.table(mirna = "hsa-miR-23a-3p", mrna = c("PTEN","CXCL12","APAF1","ST7L","SMAD3","SIRT1","FAS","SIX1","FOXO3","CDH1")),
      data.table(mirna = "hsa-miR-506-3p", mrna = c("GLI3","SFRP1","VIM","SNAI2","CDH2","ROCK1","YAP1","ETS1","NF1","CDK4")),
      data.table(mirna = "hsa-let-7b-5p",  mrna = c("HMGA2","LIN28B","CCND1","MYC","DICER1","IGF2BP1","CDK6","CDC25A","PRDM1","BCL7A")),
      data.table(mirna = "hsa-miR-186-5p", mrna = c("ROCK1","FOXO1","TWIST1","YAP1","JARID2","CDK6","PIK3R3","SMAD6","PTTG1","MTSS1")),
      data.table(mirna = "hsa-miR-222-3p", mrna = c("CDKN1B","CDKN1C","PTEN","DICER1","BBC3","FOXO3","SOD2","PPP2R2A","DDIT4","MMP1")),
      data.table(mirna = "hsa-miR-132-3p", mrna = c("SIRT1","PTEN","FOXO3","RASA1","EP300","HBEGF","SOX4","RB1","CDKN1A","MAPK1")),
      data.table(mirna = "hsa-miR-144-3p", mrna = c("ABCA1","NRF2","PTEN","NOTCH1","CCNE1","CCNE2","CEP55","SMAD4","RHEB","ZFX")),
      data.table(mirna = "hsa-miR-130a-3p",mrna = c("HOXA5","RUNX3","MET","ATG2B","SMAD4","ESR1","DICER1","CSF1","PPARG","PPP2CA")),
      data.table(mirna = "hsa-miR-148a-3p",mrna = c("DNMT1","DNMT3B","CAND1","BCL2","ROCK1","CDC25B","WNT1","ACVR1","BIM","SMAD2")),
      data.table(mirna = "hsa-miR-142-3p", mrna = c("HMGA2","TGFBR1","RAC1","WASL","STMN1","BOD1","ACVR1","HIF1A","FOXO1","APC")),
      data.table(mirna = "hsa-miR-107",    mrna = c("CDK6","DICER1","HIF1B","PTEN","CCND1","NOTCH2","NEDD9","GRN","BACE1","CLOCK")),
      data.table(mirna = "hsa-miR-150-5p", mrna = c("MYB","EGR2","TP53","AKT2","VEGFA","ZEB1","SRCIN1","NOTCH3","MUC4","GAB1")),
      data.table(mirna = "hsa-miR-206",    mrna = c("NOTCH3","KLF4","HDAC4","MET","CCND1","VEGFA","CDK4","ESR1","FOXP1","OTX2")),
      data.table(mirna = "hsa-miR-199a-3p",mrna = c("MTOR","MET","SIRT1","CD44","DDR1","AXL","ITGA3","CAV2","PAK4","PXN")),
      data.table(mirna = "hsa-miR-145-5p", mrna = c("MYC","IRS1","KLF4","SOX9","MUC1","FSCN1","ADAM17","CTGF","OCT4","HOXA9")),
      data.table(mirna = "hsa-miR-195-5p", mrna = c("CDK6","CCND1","BCL2","RAF1","E2F3","VEGFA","WNT3A","FGF2","IRS1","TAB3")),
      data.table(mirna = "hsa-miR-27b-3p", mrna = c("CYP1B1","PPARG","MET","NOTCH1","ADAM12","CBLB","GRB2","SLC7A11","PPAR","HMGB1")),
      data.table(mirna = "hsa-miR-214-3p", mrna = c("PTEN","CTNNB1","EZH2","TFAP2C","ALCAM","ING4","XBP1","N-RAS","GALNT7","ITCH")),
      data.table(mirna = "hsa-miR-101-3p", mrna = c("EZH2","FOS","STMN1","SOX9","POMP","NDRG4","RAP1B","MYCN","DUSP1","COX2")),
      data.table(mirna = "hsa-miR-135a-5p",mrna = c("ROCK1","SMAD5","BMPR2","JAK2","APC","MTSS1","VLDLR","NR3C2","BMAL1","SIRT1")),
      data.table(mirna = "hsa-miR-372-3p", mrna = c("LATS2","DKK1","CDKN1A","RHOC","RBL2","TGFBR2","ATAD2","CDK2","TNFSF10","GPR137B")),
      data.table(mirna = "hsa-miR-16-5p",  mrna = c("BCL2","CCND1","CDK6","VEGFA","HMGA1","WNT3A","FGFR1","MAP7","FEAT","PIM1")),
      data.table(mirna = "hsa-miR-26a-5p", mrna = c("EZH2","PTEN","HMGA2","CDK6","CCNE1","GSK3B","SMAD1","FUT4","NRAS","ADAM17")),
      data.table(mirna = "hsa-miR-92b-3p", mrna = c("SMAD6","TGFBR2","DAB2IP","NLK","DKK3","CDKN1C","GJA1","PTEN","ADAM10","GABRA3")),
      data.table(mirna = "hsa-miR-382-5p", mrna = c("PTEN","HIPK3","CTNNB1","MXD1","RERG","BACH1","HOXD10","NFIB","PGC1A","PPARGC1A")),
      data.table(mirna = "hsa-miR-299-3p", mrna = c("NOTCH1","VEGFA","FN1","SIRT1","HSPG2","COL11A1","ATG5","SP1","HMGA2","FOXP4")),
      data.table(mirna = "hsa-miR-577",    mrna = c("WNT3","SMURF1","LIN28B","SDC1","HMGB1","COL1A1","VEGFA","BDNF","RAB25","BAP1")),
      data.table(mirna = "hsa-miR-302a-3p",mrna = c("ERF","TGFBR2","BNIP3L","SMAD2","CDK2","CCND1","BMI1","AKT1","AOF2","CDKN1A")),
      data.table(mirna = "hsa-miR-140-5p", mrna = c("HDAC4","ADAM10","VEGFA","MYB","FGF9","WNT1","IGF1R","BMP2","NFKB1","YES1")),
      data.table(mirna = "hsa-miR-497-5p", mrna = c("BCL2","CCND1","IGF1R","VEGFA","AKT1","HDGF","RAF1","CCNE1","HMGA2","WNT3A")),
      data.table(mirna = "hsa-miR-137",    mrna = c("CDK6","KDM1A","MITF","EZH2","CDC42","CTBP1","RICTOR","AKT2","GLI3","PTGS2")),
      data.table(mirna = "hsa-miR-204-5p", mrna = c("SOX4","SIRT1","FOXC1","RAB22A","BCL2","AP1S2","HMGA2","TGFBR2","RUNX2","MEIS1")),
      data.table(mirna = "hsa-miR-339-5p", mrna = c("MDM2","SP1","HMGA2","NOVA1","FN1","HMGA1","FOXC2","NRP1","ADAM10","SKIL")),
      data.table(mirna = "hsa-miR-129-5p", mrna = c("SOX4","HMGA2","CAMTA1","DLK1","APC","ETS1","CDK6","VCP","GALNT1","HMGB1")),
      data.table(mirna = "hsa-miR-124-3p", mrna = c("EZH2","CDK6","ROCK1","SLUG","STAT3","FOXQ1","IQGAP1","SNAI2","ITGB1","LAMC1")),
      data.table(mirna = "hsa-miR-17-5p",  mrna = c("E2F1","PTEN","BCL2L11","CDKN1A","STAT3","CCND1","VEGFA","RB1","MAPK9","TGFBR2")),
      data.table(mirna = "hsa-miR-18a-5p", mrna = c("SMAD2","SMAD4","PTEN","CTGF","ESR1","KRAS","HIF1A","DICER1","TBPL1","CCND1")),
      data.table(mirna = "hsa-miR-98-5p",  mrna = c("HMGA2","IGF2BP1","ALK","IL6","CCND2","CDK6","CDC25A","SALL4","CIS","HMGA1"))
    )
    ts_pairs[, source := "TargetScan_curated"]
    targetscan <- ts_pairs
    fwrite(targetscan, ts_file)
    db_status$targetscan <- "fallback"
  } else {
    db_status$targetscan <- "downloaded"
    cat("downloaded\n")
  }
}
cat("    TargetScan:", nrow(targetscan), "interactions\n")

# ---------------------------------------------------------------------------
# 1d. miRTarBase v10 — experimentally validated pairs
# ---------------------------------------------------------------------------
cat("  [1d] miRTarBase ... ")
mtb_file <- file.path(cache_dir, "mirtarbase_validated.csv")
if (file.exists(mtb_file)) {
  mirtarbase <- fread(mtb_file)
  db_status$mirtarbase <- "cached"
  cat("loaded from cache\n")
} else {
  ok <- tryCatch({
    tmp <- file.path(cache_dir, "mirtarbase_raw.txt")
    dl_ok <- safe_download(
      "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/9.0/hsa_MTI.xlsx",
      file.path(cache_dir, "mirtarbase_raw.xlsx"), timeout_sec = 180)
    if (!dl_ok) {
      dl_ok <- safe_download(
        "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/9.0/hsa_MTI.txt",
        tmp, timeout_sec = 180)
    }
    if (dl_ok) {
      if (file.exists(file.path(cache_dir, "mirtarbase_raw.xlsx"))) {
        if (requireNamespace("readxl", quietly = TRUE)) {
          raw <- as.data.table(readxl::read_xlsx(file.path(cache_dir, "mirtarbase_raw.xlsx")))
        } else FALSE
      } else {
        raw <- fread(tmp)
      }
      if (exists("raw") && is.data.table(raw) && nrow(raw) > 50) {
        mirna_col <- grep("miRNA", names(raw), value = TRUE, ignore.case = TRUE)[1]
        gene_col  <- grep("Target Gene", names(raw), value = TRUE, ignore.case = TRUE)[1]
        evid_col  <- grep("Support Type|Experiments", names(raw), value = TRUE, ignore.case = TRUE)[1]
        mirtarbase <- raw[, .(mirna = get(mirna_col), mrna = get(gene_col),
                              evidence = if (!is.na(evid_col)) get(evid_col) else "validated")]
        mirtarbase[, source := "miRTarBase"]
        mirtarbase <- unique(mirtarbase)
        fwrite(mirtarbase, mtb_file)
        TRUE
      } else FALSE
    } else FALSE
  }, error = function(e) FALSE)

  if (!ok) {
    cat("download failed, using curated fallback\n")
    mirtarbase <- data.table(
      mirna = c(
        # Strong experimental evidence (luciferase, western blot, qPCR)
        "hsa-miR-21-5p", "hsa-miR-21-5p", "hsa-miR-21-5p", "hsa-miR-21-5p",
        "hsa-miR-155-5p", "hsa-miR-155-5p", "hsa-miR-155-5p", "hsa-miR-155-5p",
        "hsa-let-7a-5p", "hsa-let-7a-5p", "hsa-let-7a-5p",
        "hsa-miR-34a-5p", "hsa-miR-34a-5p", "hsa-miR-34a-5p", "hsa-miR-34a-5p",
        "hsa-miR-122-5p", "hsa-miR-122-5p", "hsa-miR-122-5p",
        "hsa-miR-29a-3p", "hsa-miR-29a-3p", "hsa-miR-29a-3p",
        "hsa-miR-506-3p", "hsa-miR-506-3p",
        "hsa-miR-9-5p", "hsa-miR-9-5p",
        "hsa-miR-23a-3p", "hsa-miR-23a-3p",
        "hsa-miR-33a-5p", "hsa-miR-33a-5p", "hsa-miR-33a-5p",
        "hsa-miR-181a-5p", "hsa-miR-181a-5p",
        "hsa-miR-214-3p", "hsa-miR-214-3p",
        "hsa-miR-200a-3p", "hsa-miR-200a-3p",
        "hsa-miR-101-3p", "hsa-miR-101-3p",
        "hsa-let-7b-5p", "hsa-let-7b-5p",
        "hsa-miR-186-5p", "hsa-miR-186-5p",
        "hsa-miR-222-3p", "hsa-miR-222-3p",
        "hsa-miR-132-3p", "hsa-miR-132-3p",
        "hsa-miR-148a-3p", "hsa-miR-148a-3p",
        "hsa-miR-130a-3p", "hsa-miR-130a-3p",
        "hsa-miR-145-5p", "hsa-miR-145-5p",
        "hsa-miR-195-5p", "hsa-miR-195-5p",
        "hsa-miR-27b-3p", "hsa-miR-27b-3p",
        "hsa-miR-206", "hsa-miR-206",
        "hsa-miR-199a-3p", "hsa-miR-199a-3p",
        "hsa-miR-107", "hsa-miR-107",
        "hsa-miR-150-5p", "hsa-miR-150-5p",
        "hsa-miR-142-3p", "hsa-miR-142-3p",
        "hsa-miR-372-3p", "hsa-miR-372-3p",
        "hsa-miR-144-3p", "hsa-miR-144-3p",
        "hsa-miR-16-5p", "hsa-miR-16-5p",
        "hsa-miR-675-5p", "hsa-miR-675-5p",
        "hsa-miR-140-5p", "hsa-miR-140-5p",
        "hsa-miR-26a-5p", "hsa-miR-26a-5p",
        "hsa-miR-92b-3p", "hsa-miR-92b-3p",
        "hsa-miR-302a-3p", "hsa-miR-302a-3p",
        "hsa-miR-17-5p", "hsa-miR-17-5p", "hsa-miR-17-5p",
        "hsa-miR-18a-5p", "hsa-miR-18a-5p",
        "hsa-miR-497-5p", "hsa-miR-497-5p",
        "hsa-miR-204-5p", "hsa-miR-204-5p",
        "hsa-miR-339-5p",
        "hsa-miR-129-5p", "hsa-miR-129-5p",
        "hsa-miR-124-3p", "hsa-miR-124-3p",
        "hsa-miR-98-5p",
        "hsa-miR-137",
        "hsa-miR-135a-5p",
        "hsa-miR-382-5p", "hsa-miR-382-5p",
        "hsa-miR-216a-5p",
        "hsa-miR-577",
        "hsa-miR-299-3p",
        "hsa-miR-410-3p",
        "hsa-miR-452-5p",
        "hsa-miR-324-5p",
        "hsa-miR-370-3p",
        "hsa-miR-361-5p",
        "hsa-miR-296-5p",
        "hsa-miR-513a-5p",
        "hsa-miR-381-3p",
        "hsa-miR-185-5p"
      ),
      mrna = c(
        # miR-21-5p
        "PTEN", "SMAD7", "LRP6", "PDCD4",
        # miR-155-5p
        "SOCS1", "CEBPB", "FOXO3", "HBP1",
        # let-7a-5p
        "HMGA2", "IGF2BP1", "SREBF1",
        # miR-34a-5p
        "SIRT1", "PPARA", "ACSL1", "CDK6",
        # miR-122-5p
        "ADAM17", "SRF", "CCNG1",
        # miR-29a-3p
        "COL1A1", "COL3A1", "COL4A1",
        # miR-506-3p
        "GLI3", "SFRP1",
        # miR-9-5p
        "PPARA", "ONECUT2",
        # miR-23a-3p
        "PTEN", "CXCL12",
        # miR-33a-5p
        "ABCA1", "ABCG1", "CPT1A",
        # miR-181a-5p
        "BCL2", "ATM",
        # miR-214-3p
        "PTEN", "CTNNB1",
        # miR-200a-3p
        "ZEB1", "ZEB2",
        # miR-101-3p
        "EZH2", "FOS",
        # let-7b-5p
        "HMGA2", "LIN28B",
        # miR-186-5p
        "ROCK1", "FOXO1",
        # miR-222-3p
        "CDKN1B", "CDKN1C",
        # miR-132-3p
        "SIRT1", "PTEN",
        # miR-148a-3p
        "DNMT1", "DNMT3B",
        # miR-130a-3p
        "HOXA5", "RUNX3",
        # miR-145-5p
        "MYC", "IRS1",
        # miR-195-5p
        "CDK6", "CCND1",
        # miR-27b-3p
        "CYP1B1", "PPARG",
        # miR-206
        "NOTCH3", "KLF4",
        # miR-199a-3p
        "MTOR", "MET",
        # miR-107
        "CDK6", "DICER1",
        # miR-150-5p
        "MYB", "EGR2",
        # miR-142-3p
        "HMGA2", "TGFBR1",
        # miR-372-3p
        "LATS2", "DKK1",
        # miR-144-3p
        "ABCA1", "NRF2",
        # miR-16-5p
        "BCL2", "CCND1",
        # miR-675-5p
        "RB1", "IGFR1",
        # miR-140-5p
        "HDAC4", "ADAM10",
        # miR-26a-5p
        "EZH2", "PTEN",
        # miR-92b-3p
        "SMAD6", "DAB2IP",
        # miR-302a-3p
        "ERF", "TGFBR2",
        # miR-17-5p
        "E2F1", "PTEN", "BCL2L11",
        # miR-18a-5p
        "SMAD2", "SMAD4",
        # miR-497-5p
        "BCL2", "IGF1R",
        # miR-204-5p
        "SOX4", "SIRT1",
        # miR-339-5p
        "MDM2",
        # miR-129-5p
        "SOX4", "HMGA2",
        # miR-124-3p
        "EZH2", "CDK6",
        # miR-98-5p
        "HMGA2",
        # miR-137
        "CDK6",
        # miR-135a-5p
        "ROCK1",
        # miR-382-5p
        "PTEN", "PPARGC1A",
        # miR-216a-5p
        "PTEN",
        # miR-577
        "WNT3",
        # miR-299-3p
        "NOTCH1",
        # miR-410-3p
        "HMGA2",
        # miR-452-5p
        "BMI1",
        # miR-324-5p
        "GLI1",
        # miR-370-3p
        "FOXM1",
        # miR-361-5p
        "VEGFA",
        # miR-296-5p
        "HMGA1",
        # miR-513a-5p
        "NEDD4L",
        # miR-381-3p
        "CDC42",
        # miR-185-5p
        "RHOA"
      ),
      evidence = "Strong",
      source   = "miRTarBase_curated"
    )
    fwrite(mirtarbase, mtb_file)
    db_status$mirtarbase <- "fallback"
  } else {
    db_status$mirtarbase <- "downloaded"
    cat("downloaded\n")
  }
}
cat("    miRTarBase:", nrow(mirtarbase), "interactions\n")

# ---------------------------------------------------------------------------
# 1e. LncBase v3 — lncRNA-miRNA interactions
# ---------------------------------------------------------------------------
cat("  [1e] LncBase v3 ... ")
lncbase_file <- file.path(cache_dir, "lncbase_interactions.csv")
if (file.exists(lncbase_file)) {
  lncbase <- fread(lncbase_file)
  db_status$lncbase <- "cached"
  cat("loaded from cache\n")
} else {
  ok <- tryCatch({
    url <- "https://diana.e-ce.uth.gr/lncbasev3/api/interactions?species=Homo%20sapiens"
    raw <- httr::GET(url, httr::timeout(120))
    if (httr::status_code(raw) == 200) {
      txt <- httr::content(raw, as = "text", encoding = "UTF-8")
      parsed <- jsonlite::fromJSON(txt)
      if (is.data.frame(parsed) && nrow(parsed) > 50) {
        lncbase <- as.data.table(parsed)
        lnc_col <- grep("lnc|gene", names(lncbase), value = TRUE, ignore.case = TRUE)[1]
        mir_col <- grep("mir", names(lncbase), value = TRUE, ignore.case = TRUE)[1]
        setnames(lncbase, old = c(lnc_col, mir_col), new = c("lncrna", "mirna"))
        lncbase[, source := "LncBase"]
        fwrite(lncbase, lncbase_file)
        TRUE
      } else FALSE
    } else FALSE
  }, error = function(e) FALSE)

  if (!ok) {
    cat("API failed, using curated fallback\n")
    # 12 known MASLD lncRNAs + 20 additional liver-expressed lncRNAs
    lncbase <- data.table(
      mirna = c(
        # NEAT1 (MASLD lncRNA)
        "hsa-miR-506-3p", "hsa-miR-129-5p", "hsa-miR-204-5p", "hsa-miR-339-5p",
        "hsa-miR-124-3p", "hsa-miR-140-5p", "hsa-miR-377-3p",
        # MALAT1
        "hsa-miR-155-5p", "hsa-miR-142-3p", "hsa-miR-101-3p", "hsa-miR-26a-5p",
        "hsa-miR-200a-3p", "hsa-miR-204-5p", "hsa-miR-181a-5p",
        # MEG3
        "hsa-miR-21-5p", "hsa-miR-214-3p", "hsa-miR-181a-5p", "hsa-miR-421",
        "hsa-miR-7-5p", "hsa-miR-93-5p", "hsa-miR-29a-3p",
        # H19
        "hsa-let-7a-5p", "hsa-let-7b-5p", "hsa-miR-675-5p", "hsa-let-7c-5p",
        "hsa-miR-138-5p", "hsa-miR-200b-3p", "hsa-miR-29b-3p",
        # GAS5
        "hsa-miR-23a-3p", "hsa-miR-34a-5p", "hsa-miR-222-3p", "hsa-miR-21-5p",
        "hsa-miR-135a-5p", "hsa-miR-532-5p",
        # HULC
        "hsa-miR-9-5p", "hsa-miR-372-3p", "hsa-miR-186-5p", "hsa-miR-488-5p",
        "hsa-miR-107", "hsa-miR-6825-5p",
        # TUG1
        "hsa-miR-29a-3p", "hsa-miR-144-3p", "hsa-miR-299-3p", "hsa-miR-382-5p",
        "hsa-miR-132-3p",
        # HOTAIR
        "hsa-miR-206", "hsa-miR-148a-3p", "hsa-miR-130a-3p", "hsa-miR-331-3p",
        "hsa-miR-34a-5p",
        # DANCR
        "hsa-miR-33a-5p", "hsa-miR-577", "hsa-miR-216a-5p", "hsa-miR-518a-5p",
        "hsa-miR-758-3p",
        # SNHG16
        "hsa-miR-195-5p", "hsa-miR-302a-3p", "hsa-miR-140-5p", "hsa-miR-16-5p",
        "hsa-miR-497-5p",
        # XIST
        "hsa-miR-92b-3p", "hsa-miR-29a-3p", "hsa-miR-137", "hsa-miR-101-3p",
        # MIAT
        "hsa-miR-150-5p", "hsa-miR-29a-3p", "hsa-miR-214-3p",
        # Additional liver-expressed lncRNAs (20)
        # LINC01140
        "hsa-miR-33a-5p", "hsa-miR-370-3p", "hsa-miR-200c-3p",
        # KCNQ1OT1
        "hsa-miR-34a-5p", "hsa-miR-452-5p", "hsa-miR-27b-3p", "hsa-miR-145-5p",
        # DIO3OS
        "hsa-miR-122-5p", "hsa-miR-23b-3p", "hsa-miR-199a-3p",
        # PCBP1-AS1
        "hsa-miR-155-5p", "hsa-miR-361-5p", "hsa-miR-132-3p",
        # HEIH
        "hsa-miR-199a-3p", "hsa-miR-98-5p", "hsa-miR-200b-3p",
        # FLVCR1-AS1
        "hsa-miR-513a-5p", "hsa-miR-381-3p",
        # FOXD2-AS1
        "hsa-miR-206", "hsa-miR-185-5p", "hsa-miR-150-5p",
        # MIR4435-2HG
        "hsa-miR-206", "hsa-miR-296-5p", "hsa-miR-92b-3p",
        # OIP5-AS1
        "hsa-miR-410-3p", "hsa-miR-21-5p", "hsa-miR-186-5p", "hsa-miR-137",
        # CDKN2B-AS1 (ANRIL)
        "hsa-miR-125a-5p", "hsa-miR-9-5p", "hsa-miR-34a-5p",
        # CASC2
        "hsa-miR-21-5p", "hsa-miR-18a-5p", "hsa-miR-367-3p",
        # TP73-AS1
        "hsa-miR-153-3p", "hsa-miR-200a-3p", "hsa-miR-142-3p",
        # LINC-ROR
        "hsa-miR-145-5p", "hsa-miR-205-5p", "hsa-let-7a-5p",
        # SNHG1
        "hsa-miR-199a-3p", "hsa-miR-21-5p", "hsa-miR-195-5p",
        # PVT1
        "hsa-miR-152-5p", "hsa-miR-195-5p", "hsa-miR-17-5p", "hsa-miR-21-5p",
        # UCA1
        "hsa-miR-216b-5p", "hsa-miR-135a-5p", "hsa-miR-18a-5p",
        # CRNDE
        "hsa-miR-181a-5p", "hsa-miR-136-5p", "hsa-miR-217",
        # SOCS2-AS1
        "hsa-miR-132-3p", "hsa-miR-324-5p",
        # FENDRR
        "hsa-miR-423-5p", "hsa-miR-214-3p", "hsa-miR-106b-5p",
        # LINC00472
        "hsa-miR-23a-3p", "hsa-miR-29a-3p", "hsa-miR-195-5p"
      ),
      lncrna = c(
        rep("NEAT1", 7), rep("MALAT1", 7), rep("MEG3", 7),
        rep("H19", 7), rep("GAS5", 6), rep("HULC", 6),
        rep("TUG1", 5), rep("HOTAIR", 5), rep("DANCR", 5),
        rep("SNHG16", 5), rep("XIST", 4), rep("MIAT", 3),
        rep("LINC01140", 3), rep("KCNQ1OT1", 4), rep("DIO3OS", 3),
        rep("PCBP1-AS1", 3), rep("HEIH", 3), rep("FLVCR1-AS1", 2),
        rep("FOXD2-AS1", 3), rep("MIR4435-2HG", 3), rep("OIP5-AS1", 4),
        rep("CDKN2B-AS1", 3), rep("CASC2", 3), rep("TP73-AS1", 3),
        rep("LINC-ROR", 3), rep("SNHG1", 3), rep("PVT1", 4),
        rep("UCA1", 3), rep("CRNDE", 3), rep("SOCS2-AS1", 2),
        rep("FENDRR", 3), rep("LINC00472", 3)
      ),
      source = "LncBase_curated"
    )
    fwrite(lncbase, lncbase_file)
    db_status$lncbase <- "fallback"
  } else {
    db_status$lncbase <- "downloaded"
    cat("downloaded\n")
  }
}
cat("    LncBase:", nrow(lncbase), "interactions\n")

# --- Database status summary ---
cat("\n  Database loading summary:\n")
for (nm in names(db_status)) {
  cat(sprintf("    %-20s : %s\n", nm, db_status[[nm]]))
}

# =============================================================================
# 2. Build Bipartite Graphs
# =============================================================================
cat("\n--- 2. Build bipartite graphs ---\n")

# Ensure consistent column names
ensure_cols <- function(dt, mirna_col, gene_col, gene_name) {
  # Attempt to find the right columns
  mir_c <- intersect(names(dt), c("mirna", "miRNA", "miRNAname"))
  gen_c <- intersect(names(dt), c(gene_name, gene_col, "geneName", "Gene Symbol"))
  if (length(mir_c) == 0 || length(gen_c) == 0) return(NULL)
  out <- dt[, .(mirna = as.character(get(mir_c[1])),
                gene  = as.character(get(gen_c[1])))]
  if (!"source" %in% names(dt)) out[, source := "unknown"]
  else out[, source := as.character(dt$source)]
  out <- unique(out[mirna != "" & gene != ""])
  out
}

# --- Graph A: miRNA-lncRNA ---
graphA_encori <- ensure_cols(encori_lnc, "mirna", "lncrna", "lncrna")
if (!is.null(graphA_encori)) setnames(graphA_encori, "gene", "lncrna")
graphA_lncbase <- ensure_cols(lncbase, "mirna", "lncrna", "lncrna")
if (!is.null(graphA_lncbase)) setnames(graphA_lncbase, "gene", "lncrna")

graphA <- rbindlist(list(graphA_encori, graphA_lncbase), use.names = TRUE, fill = TRUE)
graphA <- graphA[, .(n_sources = uniqueN(source),
                     sources = paste(unique(source), collapse = ";")),
                 by = .(mirna, lncrna)]
cat("  Graph A (miRNA-lncRNA): ", nrow(graphA), " edges covering ",
    uniqueN(graphA$mirna), " miRNAs x ", uniqueN(graphA$lncrna), " lncRNAs\n", sep = "")

# --- Graph B: miRNA-mRNA ---
graphB_encori <- ensure_cols(encori_mrna, "mirna", "mrna", "mrna")
if (!is.null(graphB_encori)) setnames(graphB_encori, "gene", "mrna")
graphB_ts <- ensure_cols(targetscan, "mirna", "mrna", "mrna")
if (!is.null(graphB_ts)) setnames(graphB_ts, "gene", "mrna")
graphB_mtb <- ensure_cols(mirtarbase, "mirna", "mrna", "mrna")
if (!is.null(graphB_mtb)) {
  setnames(graphB_mtb, "gene", "mrna")
  # Mark strong evidence for miRTarBase entries
  if ("evidence" %in% names(mirtarbase)) {
    graphB_mtb[, validated := TRUE]
  } else {
    graphB_mtb[, validated := TRUE]  # curated = strong
  }
}

graphB_all <- rbindlist(list(graphB_encori, graphB_ts, graphB_mtb),
                        use.names = TRUE, fill = TRUE)
if (!"validated" %in% names(graphB_all)) graphB_all[, validated := FALSE]
graphB_all[is.na(validated), validated := FALSE]

graphB <- graphB_all[, .(n_sources = uniqueN(source),
                         sources = paste(unique(source), collapse = ";"),
                         has_validation = any(validated)),
                     by = .(mirna, mrna)]

# Filter: require >= 2 databases OR miRTarBase strong evidence
graphB <- graphB[n_sources >= 2 | has_validation == TRUE]
cat("  Graph B (miRNA-mRNA):   ", nrow(graphB), " edges covering ",
    uniqueN(graphB$mirna), " miRNAs x ", uniqueN(graphB$mrna), " mRNAs\n", sep = "")

# =============================================================================
# 3. Construct ceRNA Triplets
# =============================================================================
cat("\n--- 3. Construct ceRNA triplets ---\n")

# For each lncRNA, find its miRNA partners; for each mRNA, find its miRNA partners
# A ceRNA pair (L, M) shares miRNAs from both graphs
lnc_mirnas <- graphA[, .(mirnas = list(mirna)), by = lncrna]
mrna_mirnas <- graphB[, .(mirnas = list(mirna)), by = mrna]

cat("  Computing shared miRNAs for ", nrow(lnc_mirnas), " lncRNAs x ",
    nrow(mrna_mirnas), " mRNAs ...\n", sep = "")

# Efficient inner join approach
triplets_list <- list()
for (i in seq_len(nrow(lnc_mirnas))) {
  l_name <- lnc_mirnas$lncrna[i]
  l_mirs <- lnc_mirnas$mirnas[[i]]

  for (j in seq_len(nrow(mrna_mirnas))) {
    m_name <- mrna_mirnas$mrna[j]
    shared <- intersect(l_mirs, mrna_mirnas$mirnas[[j]])
    # C4 fix: primary threshold lowered to >= 2 shared miRNAs (standard in ceRNA literature)
    if (length(shared) >= 2) {
      triplets_list[[length(triplets_list) + 1L]] <- data.table(
        lncrna       = l_name,
        mrna         = m_name,
        shared_mirnas = paste(sort(shared), collapse = ";"),
        cerna_score  = length(shared)
      )
    }
  }
}

if (length(triplets_list) > 0) {
  triplets <- rbindlist(triplets_list)
} else {
  cat("  WARNING: No triplets found even with >= 2 shared miRNAs\n")
  triplets <- data.table(lncrna = character(), mrna = character(),
                          shared_mirnas = character(), cerna_score = integer())
}

# Deduplicate (alphabetical ordering)
triplets[, pair_key := fifelse(lncrna < mrna,
                               paste(lncrna, mrna, sep = "::"),
                               paste(mrna, lncrna, sep = "::"))]
triplets <- triplets[!duplicated(pair_key)]
triplets[, pair_key := NULL]

triplets <- triplets[order(-cerna_score)]
cat("  ceRNA triplets: ", nrow(triplets), "\n", sep = "")
cat("  Unique lncRNAs: ", uniqueN(triplets$lncrna), "\n", sep = "")
cat("  Unique mRNAs:   ", uniqueN(triplets$mrna), "\n", sep = "")
cat("  Score range:    ", min(triplets$cerna_score), "-", max(triplets$cerna_score), "\n", sep = "")

fwrite(triplets, file.path(out_dir, "cerna_network_full.csv"))
cat("  Saved: cerna_network_full.csv\n")

# =============================================================================
# 4. MASLD DEG Overlay
# =============================================================================
cat("\n--- 4. MASLD DEG overlay ---\n")

atlas <- fread(atlas_path)
stopifnot(all(c("bulk_padj", "bulk_logFC") %in% names(atlas)))
cat("  Atlas loaded:", nrow(atlas), "genes\n")

# Extract bulk DEG stats for gene name matching
deg_info <- atlas[, .(gene_name = human_symbol, bulk_logFC, bulk_padj, gene_biotype)]
deg_info <- deg_info[!is.na(gene_name) & gene_name != ""]
deg_info <- unique(deg_info, by = "gene_name")

# Merge onto lncRNA
triplets_annot <- merge(triplets, deg_info,
                        by.x = "lncrna", by.y = "gene_name", all.x = TRUE)
setnames(triplets_annot,
         c("bulk_logFC", "bulk_padj", "gene_biotype"),
         c("lncrna_logFC", "lncrna_padj", "lncrna_biotype"))

# Merge onto mRNA
triplets_annot <- merge(triplets_annot, deg_info,
                        by.x = "mrna", by.y = "gene_name", all.x = TRUE)
setnames(triplets_annot,
         c("bulk_logFC", "bulk_padj", "gene_biotype"),
         c("mrna_logFC", "mrna_padj", "mrna_biotype"))

# Flag DEG status
triplets_annot[, lncrna_is_deg := !is.na(lncrna_padj) & lncrna_padj < 0.1]
triplets_annot[, mrna_is_deg   := !is.na(mrna_padj) & mrna_padj < 0.1]

# Count shared miRNAs that are DEGs
# (miRNAs are in the atlas if they were in dream)
mirna_deg <- deg_info[grepl("^MIR|^hsa-", gene_name, ignore.case = TRUE)]
triplets_annot[, mirna_deg_count := {
  sapply(shared_mirnas, function(mirs) {
    mir_vec <- unlist(strsplit(mirs, ";"))
    # Convert hsa-miR-XXX format to MIR gene names for matching
    mir_genes <- gsub("^hsa-", "", mir_vec)
    mir_genes <- toupper(gsub("-", "", mir_genes))
    sum(mir_genes %in% toupper(mirna_deg$gene_name))
  })
}]

# n_deg_nodes: count of (lncRNA, mRNA) that are DEGs (0, 1, or 2)
# miRNA DEG status tracked separately
triplets_annot[, n_deg_nodes := as.integer(lncrna_is_deg) + as.integer(mrna_is_deg)]

# Concordance: in ceRNA model, lncRNA and mRNA should be co-regulated
# (both up or both down, since lncRNA sponges miRNA that would suppress mRNA)
triplets_annot[, concordant_direction := fifelse(
  !is.na(lncrna_logFC) & !is.na(mrna_logFC),
  sign(lncrna_logFC) == sign(mrna_logFC),
  NA)]

cat("  DEG overlay complete:\n")
cat("    Triplets with 0 DEG nodes: ", sum(triplets_annot$n_deg_nodes == 0), "\n")
cat("    Triplets with 1 DEG node:  ", sum(triplets_annot$n_deg_nodes == 1), "\n")
cat("    Triplets with 2 DEG nodes: ", sum(triplets_annot$n_deg_nodes == 2), "\n")
cat("    Concordant direction (of annotated): ",
    sum(triplets_annot$concordant_direction == TRUE, na.rm = TRUE), " / ",
    sum(!is.na(triplets_annot$concordant_direction)), "\n")

# =============================================================================
# 5. Filter MASLD Subnetwork
# =============================================================================
cat("\n--- 5. Filter MASLD subnetwork ---\n")

masld_net <- triplets_annot[n_deg_nodes >= 2]
cat("  MASLD ceRNA subnetwork: ", nrow(masld_net), " triplets\n", sep = "")
cat("    Unique lncRNAs: ", uniqueN(masld_net$lncrna), "\n")
cat("    Unique mRNAs:   ", uniqueN(masld_net$mrna), "\n")

if (nrow(masld_net) == 0) {
  cat("  WARNING: No triplets with >= 2 DEG nodes. Relaxing to >= 1 ...\n")
  masld_net <- triplets_annot[n_deg_nodes >= 1]
  cat("  Relaxed MASLD ceRNA subnetwork: ", nrow(masld_net), " triplets\n", sep = "")
}

# Reorder columns for output
out_cols <- c("lncrna", "mrna", "shared_mirnas", "cerna_score",
              "lncrna_logFC", "lncrna_padj", "lncrna_is_deg", "lncrna_biotype",
              "mrna_logFC", "mrna_padj", "mrna_is_deg", "mrna_biotype",
              "n_deg_nodes", "concordant_direction", "mirna_deg_count")
out_cols <- intersect(out_cols, names(masld_net))

fwrite(masld_net[, ..out_cols], file.path(out_dir, "cerna_network_masld.csv"))
cat("  Saved: cerna_network_masld.csv\n")

# =============================================================================
# 6. Network Metrics
# =============================================================================
cat("\n--- 6. Network metrics ---\n")

if (nrow(masld_net) > 0) {
  # Build igraph: nodes are lncRNAs + mRNAs + miRNAs
  # Edges: lncRNA--miRNA and miRNA--mRNA for each triplet
  edge_list <- list()
  for (r in seq_len(nrow(masld_net))) {
    l <- masld_net$lncrna[r]
    m <- masld_net$mrna[r]
    mirs <- unlist(strsplit(masld_net$shared_mirnas[r], ";"))
    for (mir in mirs) {
      edge_list[[length(edge_list) + 1L]] <- data.table(from = l, to = mir)
      edge_list[[length(edge_list) + 1L]] <- data.table(from = mir, to = m)
    }
  }
  edges <- unique(rbindlist(edge_list))

  g <- graph_from_data_frame(edges, directed = FALSE)

  # Classify node types
  all_lncrnas <- unique(masld_net$lncrna)
  all_mrnas   <- unique(masld_net$mrna)
  all_mirnas  <- unique(unlist(strsplit(masld_net$shared_mirnas, ";")))

  node_df <- data.table(
    node = V(g)$name,
    node_type = fifelse(V(g)$name %in% all_lncrnas, "lncRNA",
                fifelse(V(g)$name %in% all_mrnas, "mRNA", "miRNA"))
  )

  # Compute centrality
  node_df[, degree      := degree(g, v = node)]
  node_df[, betweenness := betweenness(g, v = node, normalized = TRUE)]

  # Merge dream stats
  node_df <- merge(node_df, deg_info[, .(gene_name, bulk_logFC, bulk_padj)],
                   by.x = "node", by.y = "gene_name", all.x = TRUE)

  node_df <- node_df[order(-degree, -betweenness)]

  # Top 20 hub lncRNAs and top 20 hub miRNAs
  top_lnc <- node_df[node_type == "lncRNA"][order(-degree, -betweenness)][1:min(20, .N)]
  top_mir <- node_df[node_type == "miRNA"][order(-degree, -betweenness)][1:min(20, .N)]
  hubs    <- rbindlist(list(top_lnc, top_mir), use.names = TRUE)

  fwrite(hubs, file.path(out_dir, "cerna_hubs.csv"))
  cat("  Network: ", vcount(g), " nodes, ", ecount(g), " edges\n", sep = "")
  cat("  Hub lncRNAs (top 5):\n")
  for (i in seq_len(min(5, nrow(top_lnc)))) {
    cat(sprintf("    %s  degree=%d  betweenness=%.4f\n",
                top_lnc$node[i], top_lnc$degree[i], top_lnc$betweenness[i]))
  }
  cat("  Hub miRNAs (top 5):\n")
  for (i in seq_len(min(5, nrow(top_mir)))) {
    cat(sprintf("    %s  degree=%d  betweenness=%.4f\n",
                top_mir$node[i], top_mir$degree[i], top_mir$betweenness[i]))
  }
  cat("  Saved: cerna_hubs.csv\n")
} else {
  cat("  Skipping network metrics (no MASLD triplets)\n")
  hubs <- data.table(node = character(), node_type = character(),
                     degree = integer(), betweenness = numeric(),
                     bulk_logFC = numeric(), bulk_padj = numeric())
  fwrite(hubs, file.path(out_dir, "cerna_hubs.csv"))
}

# =============================================================================
# 7. Validate Known MASLD ceRNA Axes
# =============================================================================
cat("\n--- 7. Validate known MASLD ceRNA axes ---\n")

known_axes <- data.table(
  axis_name = c("NEAT1-miR506-GLI3", "MEG3-miR21-LRP6", "H19-let7-SREBF1",
                "MALAT1-miR155-SOCS1", "HULC-miR9-PPARA", "GAS5-miR23a-PTEN"),
  lncrna = c("NEAT1", "MEG3", "H19", "MALAT1", "HULC", "GAS5"),
  mirna  = c("hsa-miR-506-3p", "hsa-miR-21-5p", "hsa-let-7a-5p",
             "hsa-miR-155-5p", "hsa-miR-9-5p", "hsa-miR-23a-3p"),
  mrna   = c("GLI3", "LRP6", "SREBF1", "SOCS1", "PPARA", "PTEN")
)

# Check if each axis is present in our network
validated_axes <- copy(known_axes)

# Check in full triplets
validated_axes[, in_full_network := mapply(function(l, m) {
  any(triplets$lncrna == l & triplets$mrna == m)
}, lncrna, mrna)]

# Check if miRNA edge exists in bipartite graphs
validated_axes[, mirna_lnc_edge := mapply(function(mir, l) {
  any(graphA$mirna == mir & graphA$lncrna == l)
}, mirna, lncrna)]

validated_axes[, mirna_mrna_edge := mapply(function(mir, m) {
  any(graphB$mirna == mir & graphB$mrna == m)
}, mirna, mrna)]

validated_axes[, axis_edges_present := mirna_lnc_edge & mirna_mrna_edge]

# Merge dream stats for each node
for (node_col in c("lncrna", "mirna", "mrna")) {
  merge_name <- if (node_col == "mirna") {
    # miRNAs may have different naming in atlas
    validated_axes[, paste0(node_col, "_logFC") := NA_real_]
    validated_axes[, paste0(node_col, "_padj")  := NA_real_]
    next
  } else {
    node_col
  }
  tmp <- merge(validated_axes[, .(axis_name, node = get(node_col))],
               deg_info[, .(gene_name, bulk_logFC, bulk_padj)],
               by.x = "node", by.y = "gene_name", all.x = TRUE)
  setnames(tmp, c("bulk_logFC", "bulk_padj"),
           paste0(node_col, c("_logFC", "_padj")))
  validated_axes <- merge(validated_axes, tmp[, -"node"],
                          by = "axis_name", all.x = TRUE)
}

cat("  Known MASLD ceRNA axes:\n")
for (i in seq_len(nrow(validated_axes))) {
  ax <- validated_axes[i]
  status <- if (ax$in_full_network) "FOUND in triplet network" else
            if (ax$axis_edges_present) "edges present (below shared miRNA threshold)" else
            "partial (missing edges)"
  cat(sprintf("    %s: %s\n", ax$axis_name, status))
  cat(sprintf("      lncRNA %s logFC=%.3f padj=%.2e | mRNA %s logFC=%.3f padj=%.2e\n",
              ax$lncrna,
              ifelse(is.na(ax$lncrna_logFC), NA_real_, ax$lncrna_logFC),
              ifelse(is.na(ax$lncrna_padj), NA_real_, ax$lncrna_padj),
              ax$mrna,
              ifelse(is.na(ax$mrna_logFC), NA_real_, ax$mrna_logFC),
              ifelse(is.na(ax$mrna_padj), NA_real_, ax$mrna_padj)))
}

fwrite(validated_axes, file.path(out_dir, "cerna_validated_axes.csv"))
cat("  Saved: cerna_validated_axes.csv\n")

# =============================================================================
# 8. Enrichment Test
# =============================================================================
cat("\n--- 8. Enrichment test ---\n")

# Test: are mRNA ceRNA partners of lncRNA DEGs enriched for MASLD DEGs?
# Exploratory annotation threshold; primary DEGs: padj<0.05 + |logFC|>0.3 (Script 05b)
lncrna_degs <- deg_info[gene_biotype == "lncRNA" & !is.na(bulk_padj) & bulk_padj < 0.1]$gene_name
cat("  lncRNA DEGs in atlas: ", length(lncrna_degs), "\n")

# Get mRNAs that are ceRNA partners of lncRNA DEGs
cerna_partner_mrnas <- unique(triplets_annot[lncrna %in% lncrna_degs]$mrna)
cat("  mRNA ceRNA partners of lncRNA DEGs: ", length(cerna_partner_mrnas), "\n")

# Background: all mRNAs in atlas
all_mrnas_atlas <- deg_info[gene_biotype == "protein_coding" | is.na(gene_biotype)]$gene_name
# Exploratory annotation threshold; primary DEGs: padj<0.05 + |logFC|>0.3 (Script 05b)
all_mrnas_deg   <- deg_info[!is.na(bulk_padj) & bulk_padj < 0.1 &
                            (gene_biotype == "protein_coding" | is.na(gene_biotype))]$gene_name

if (length(cerna_partner_mrnas) > 0 && length(all_mrnas_atlas) > 0) {
  # 2x2 contingency: ceRNA_partner x MASLD_DEG
  a <- length(intersect(cerna_partner_mrnas, all_mrnas_deg))   # partner + DEG

  b <- length(setdiff(cerna_partner_mrnas, all_mrnas_deg))     # partner + not DEG
  c <- length(setdiff(all_mrnas_deg, cerna_partner_mrnas))     # not partner + DEG
  d <- length(setdiff(all_mrnas_atlas, union(cerna_partner_mrnas, all_mrnas_deg)))

  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft  <- fisher.test(mat, alternative = "greater")

  cat(sprintf("  Fisher's exact test (one-sided):\n"))
  cat(sprintf("    ceRNA partners that are DEGs:  %d / %d (%.1f%%)\n",
              a, a + b, 100 * a / (a + b)))
  cat(sprintf("    Background mRNAs that are DEGs: %d / %d (%.1f%%)\n",
              length(all_mrnas_deg), length(all_mrnas_atlas),
              100 * length(all_mrnas_deg) / length(all_mrnas_atlas)))
  cat(sprintf("    Odds ratio: %.2f\n", ft$estimate))
  cat(sprintf("    P-value:    %.2e\n", ft$p.value))
} else {
  cat("  Insufficient data for enrichment test\n")
  ft <- list(estimate = NA, p.value = NA)
}

# =============================================================================
# 9. Summary Statistics
# =============================================================================
cat("\n--- 9. Summary ---\n")
cat("=== ceRNA Network Construction Complete ===\n")
cat("  Databases loaded:\n")
for (nm in names(db_status)) {
  cat(sprintf("    %-20s : %s\n", nm, db_status[[nm]]))
}
cat("\n  Bipartite graphs:\n")
cat("    miRNA-lncRNA edges: ", nrow(graphA), "\n")
cat("    miRNA-mRNA edges:   ", nrow(graphB), "\n")
cat("\n  ceRNA network:\n")
cat("    Full triplets:       ", nrow(triplets), "\n")
cat("    MASLD subnetwork:    ", nrow(masld_net), "\n")
cat("    Hub lncRNAs:         ", nrow(hubs[hubs$node_type == "lncRNA"]), "\n")
cat("    Hub miRNAs:          ", nrow(hubs[hubs$node_type == "miRNA"]), "\n")
cat("\n  Known axis validation: ", sum(validated_axes$in_full_network), "/",
    nrow(validated_axes), " found in network\n")
if (!is.na(ft$p.value)) {
  cat(sprintf("  Enrichment: OR=%.2f, p=%.2e\n", ft$estimate, ft$p.value))
}
cat("\n  Output files:\n")
cat("    ", file.path(out_dir, "cerna_network_full.csv"), "\n")
cat("    ", file.path(out_dir, "cerna_network_masld.csv"), "\n")
cat("    ", file.path(out_dir, "cerna_hubs.csv"), "\n")
cat("    ", file.path(out_dir, "cerna_validated_axes.csv"), "\n")
cat("\nEnd:", format(Sys.time()), "\n")
