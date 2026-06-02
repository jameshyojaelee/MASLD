#!/usr/bin/env Rscript
# =============================================================================
# Build combined tximport object for ALL Western Diet datasets
# =============================================================================
# Aggregates Kallisto transcript-level quant -> gene-level across all WD
# datasets into a single txi_wd.rds for use by run_wd_kallisto_de.R.
#
# Datasets:
#   GSE220575 (DIAMOND)    — paired-end
#   GSE246088 (WD/Plvap)   — paired-end
#   GSE305484 (WD+fructose)— paired-end
#   GSE292565 (FFC diet)   — paired-end, DNBSEQ-T7
#   GSE246328 (GAN diet)   — single-end, NextSeq 500
#
# Usage:
#   Rscript build_combined_txi.R
# =============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(tximport)
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
WDDIR <- file.path(BASE, "RNA-seq/Mouse/Western_Diet_Datasets")
TX2GENE_PATH <- file.path(BASE, "data/reference/kallisto/tx2gene_vM38.tsv")
OUT_PATH <- file.path(WDDIR, "tximport/txi_wd.rds")

cat("=== Building combined WD tximport object ===\n\n")

# ---- tx2gene ----
tx2gene <- fread(TX2GENE_PATH)
tx2gene_df <- data.frame(
  tx_id   = sub("\\.[0-9]+$", "", tx2gene$transcript_id),
  gene_id = sub("\\.[0-9]+$", "", tx2gene$gene_id),
  stringsAsFactors = FALSE
)
cat(sprintf("tx2gene: %d transcripts -> %d genes\n\n",
            nrow(tx2gene_df), length(unique(tx2gene_df$gene_id))))

# ---- Collect Kallisto outputs from all datasets ----
DATASETS <- c("GSE220575", "GSE246088", "GSE305484", "GSE292565", "GSE246328")
all_files <- c()

for (ds in DATASETS) {
  quant_dir <- file.path(WDDIR, ds, "quant/kallisto")
  gsm_srr_path <- file.path(WDDIR, ds, "metadata/gsm_to_srr.tsv")

  if (!dir.exists(quant_dir)) {
    cat(sprintf("  %s: quant dir not found, skipping\n", ds))
    next
  }

  gsm_srr <- fread(gsm_srr_path)

  # Find all SRR directories with abundance files
  srr_dirs <- list.dirs(quant_dir, recursive = FALSE, full.names = FALSE)
  found <- 0
  for (srr in srr_dirs) {
    h5_file <- file.path(quant_dir, srr, "abundance.h5")
    tsv_file <- file.path(quant_dir, srr, "abundance.tsv")

    if (file.exists(h5_file)) {
      target <- h5_file
    } else if (file.exists(tsv_file)) {
      target <- tsv_file
    } else {
      next
    }

    # Name by SRR accession (used as column name)
    all_files[srr] <- target
    found <- found + 1
  }
  cat(sprintf("  %s: %d samples found\n", ds, found))
}

cat(sprintf("\nTotal samples: %d\n\n", length(all_files)))

if (length(all_files) == 0) {
  stop("No Kallisto outputs found!")
}

# ---- tximport ----
cat("Running tximport (countsFromAbundance='no')...\n")
txi <- tximport(all_files,
                type = "kallisto",
                tx2gene = tx2gene_df,
                countsFromAbundance = "no",
                ignoreTxVersion = TRUE)
cat(sprintf("  Gene-level: %d genes x %d samples\n",
            nrow(txi$counts), ncol(txi$counts)))

# ---- Save ----
dir.create(dirname(OUT_PATH), showWarnings = FALSE, recursive = TRUE)
saveRDS(txi, OUT_PATH)
cat(sprintf("\nSaved: %s (%.1f MB)\n", OUT_PATH,
            file.size(OUT_PATH) / 1e6))
cat("=== Done ===\n")
