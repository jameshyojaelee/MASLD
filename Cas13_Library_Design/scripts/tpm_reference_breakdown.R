#!/usr/bin/env Rscript
# tpm_reference_breakdown.R -- enumerate exactly which mouse disease samples feed
# the pooled disease-liver TPM reference used by the library TPM gate.
suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("MASLD_PROJECT_ROOT","/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FC_DIR <- file.path(BASE,"RNA-seq/Mouse/Unified_Integration/counts/featurecounts")
PUB_FC <- file.path(BASE,"RNA-seq/Mouse/Public_Diet_Models/counts/featurecounts/gene_counts.txt")
WD_DIR <- file.path(BASE,"RNA-seq/Mouse/Western_Diet_Datasets")
MAIN_META <- file.path(BASE,"RNA-seq/Mouse/Unified_Integration/metadata/unified_mouse_metadata.csv")
strip_v <- function(x) sub("\\..*","",x)
cols_of <- function(path){ hdr <- fread(path, skip="Geneid", nrows=1)
  sc <- setdiff(names(hdr), c("Geneid","Chr","Start","End","Strand","Length")); basename(dirname(sc)) }

meta <- fread(MAIN_META)
dz <- if("group_binary" %in% names(meta)) meta[group_binary!="Control"] else meta[!grepl("Control|LFD|Chow|Ctrl",condition)]
dz_main <- dz$sample_id

# --- unified-integration sources (in-house + public), labelled by diet_model ---
mcd_cols <- unlist(lapply(file.path(FC_DIR,c("gene_counts_inhouse.txt","gene_counts_gse156918.txt","gene_counts_gse205974.txt")), cols_of))
pub_cols <- cols_of(PUB_FC)
unified_dz <- dz[sample_id %in% c(mcd_cols, pub_cols)]
cat("=== Unified-integration disease samples by diet_model & dataset ===\n")
print(unified_dz[, .N, by=.(diet_model, dataset)][order(diet_model, dataset)])

# --- standalone Western-diet datasets ---
g220 <- merge(fread(file.path(WD_DIR,"GSE220575/metadata/sample_metadata.csv")),
              fread(file.path(WD_DIR,"GSE220575/metadata/gsm_to_srr.tsv")), by="gsm")[condition %in% c("MASH","HCC")]
g246 <- merge(fread(file.path(WD_DIR,"GSE246088/metadata/sample_metadata.csv")),
              fread(file.path(WD_DIR,"GSE246088/metadata/gsm_to_srr.tsv")), by.x="geo_accession", by.y="gsm")[genotype=="Plvap_Control" & diet %in% c("Western_Diet","High_Fat_Diet")]
g305 <- merge(fread(file.path(WD_DIR,"GSE305484/metadata/sample_metadata.csv")),
              fread(file.path(WD_DIR,"GSE305484/metadata/gsm_to_srr.tsv")), by.x="gsm_accession", by.y="gsm")[condition!="Chow"]
cat("\n=== Standalone Western-diet datasets ===\n")
cat(sprintf("GSE220575 (Western, MASH+HCC): %d\n", nrow(g220)))
cat(sprintf("GSE246088 (Western_Diet + High_Fat_Diet, Plvap_Control): %d  [WD=%d, HFD=%d]\n",
            nrow(g246), sum(g246$diet=="Western_Diet"), sum(g246$diet=="High_Fat_Diet")))
cat(sprintf("GSE305484 (Western, non-Chow): %d\n", nrow(g305)))

cat(sprintf("\nGRAND TOTAL pooled disease samples = %d\n",
            nrow(unified_dz) + nrow(g220) + nrow(g246) + nrow(g305)))
cat("\nPer-gene reference value = MEDIAN TPM across ALL of these pooled disease columns.\n")
