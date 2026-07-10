#!/usr/bin/env Rscript
# 24_prepare_bulk_tpm.R
# ---------------------------------------------------------------------------
# Extract per-cohort gene-symbol TPM matrices from the kallisto canonical
# all_cohorts_gene_tpm.tsv.gz for Rectangle deconvolution.
#
# Rectangle (Wagner et al. 2024) expects bulk expression in TPM (unlike MuSiC
# / BayesPrism which use counts), so this mirrors 20_prepare_kallisto_bulk_counts.R
# but keeps TPM floats (no rounding) and aggregates duplicate symbols by SUM
# (TPM is additive within a sample, so summing versioned Ensembl IDs that
#  collapse to one symbol is the correct aggregation, unlike a mean).
#
# The per-cohort TSVs match the format of the existing {COHORT}_counts.tsv:
#   - Rows = gene symbols, first column header = "SYMBOL"
#   - Columns = SRR sample IDs (taken from the existing counts.tsv header)
#   - Values = TPM floats
#
# Input:
#   - RNA-seq/results/kallisto_sensitivity/all_cohorts_gene_tpm.tsv.gz
#     (header col "gene_id" = versioned Ensembl, e.g. ENSG00000000003.17)
#   - Analysis/Deconvolution/bulk/{COHORT}/{COHORT}_counts.tsv (SRR column list)
#   - GENCODE v49 GTF (Ensembl -> gene symbol mapping)
#
# Output:
#   - Analysis/Deconvolution/bulk/{COHORT}/{COHORT}_tpm_symbol.tsv (one per cohort)
#
# Usage: Rscript 24_prepare_bulk_tpm.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== 24: Prepare per-cohort kallisto bulk TPM for Rectangle ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
TPM_FILE  <- file.path(BASE, "RNA-seq/results/kallisto_sensitivity/all_cohorts_gene_tpm.tsv.gz")
GTF_FILE  <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
BULK_DIR  <- file.path(BASE, "Analysis/Deconvolution/bulk")

# Mega cohorts (the only cohorts in scope for the Disease-vs-Control atlas)
COHORTS <- c("GSE126848", "GSE135251", "GSE130970", "GSE213621", "GSE162694")

stopifnot(file.exists(TPM_FILE), file.exists(GTF_FILE))

# --- 1. Build Ensembl -> gene symbol map from GENCODE v49 ---
cat("Building Ensembl -> symbol map from GENCODE v49 GTF...\n")
gtf_cmd <- paste("zcat", GTF_FILE, "| grep -P '^chr.*\\tgene\\t'")
gtf_lines <- fread(cmd = gtf_cmd, header = FALSE, sep = "\t")
gene_ids   <- gsub('.*gene_id "([^"]+)".*', "\\1", gtf_lines$V9)
gene_names <- gsub('.*gene_name "([^"]+)".*', "\\1", gtf_lines$V9)
# Strip version from GTF gene IDs
gene_ids_base <- sub("\\..*$", "", gene_ids)
ens2sym <- data.table(ensembl = gene_ids_base, symbol = gene_names)
ens2sym <- unique(ens2sym, by = "ensembl")
cat(sprintf("  GTF mapping: %d unique Ensembl -> symbol entries\n", nrow(ens2sym)))

# --- 2. Load kallisto TPM ---
cat("\nLoading kallisto TPM...\n")
tpm <- fread(TPM_FILE)
gene_col <- tpm$gene_id
tpm[, gene_id := NULL]
cat(sprintf("  TPM matrix: %d genes x %d samples\n", length(gene_col), ncol(tpm)))

# Strip version from kallisto gene IDs
gene_base <- sub("\\..*$", "", gene_col)

# Map to symbols
sym <- ens2sym$symbol[match(gene_base, ens2sym$ensembl)]
n_mapped <- sum(!is.na(sym))
cat(sprintf("  Mapped: %d / %d genes to symbols (%.1f%%)\n",
            n_mapped, length(sym), 100 * n_mapped / length(sym)))

# Keep only mapped genes, aggregate duplicates by symbol via SUM (TPM additive)
tpm_dt <- data.table(SYMBOL = sym, tpm)
tpm_dt <- tpm_dt[!is.na(SYMBOL)]

n_before <- nrow(tpm_dt)
value_cols <- names(tpm_dt)[names(tpm_dt) != "SYMBOL"]
tpm_dt <- tpm_dt[, lapply(.SD, sum), by = SYMBOL, .SDcols = value_cols]
cat(sprintf("  After symbol aggregation (SUM): %d rows (was %d)\n",
            nrow(tpm_dt), n_before))

# Sanity: no duplicated symbols in aggregated output
stopifnot(!any(duplicated(tpm_dt$SYMBOL)))

sample_cols_all <- setdiff(names(tpm_dt), "SYMBOL")

# --- 3. Write per-cohort TPM files ---
cat("\nWriting per-cohort TPM files...\n")

for (ds in COHORTS) {
  counts_file <- file.path(BULK_DIR, ds, paste0(ds, "_counts.tsv"))
  if (!file.exists(counts_file)) {
    stop(sprintf("Missing counts file for cohort %s: %s", ds, counts_file))
  }

  # SRR sample IDs = header of the existing counts.tsv (minus SYMBOL)
  ds_samples <- setdiff(names(fread(counts_file, nrows = 0)), "SYMBOL")

  # Assert every requested SRR is present in the TPM source
  missing <- setdiff(ds_samples, sample_cols_all)
  if (length(missing) > 0) {
    stop(sprintf("%s: %d SRR ids absent from TPM source: %s",
                 ds, length(missing), paste(head(missing, 5), collapse = ", ")))
  }

  out_dt <- tpm_dt[, c("SYMBOL", ds_samples), with = FALSE]

  # Per-sample TPM column sums (should be ~1e6 by TPM property, minus
  # unmapped/dropped genes)
  colsums <- colSums(as.matrix(out_dt[, ds_samples, with = FALSE]))
  cat(sprintf("  %s: %d genes x %d samples | col-sum range [%.0f, %.0f], mean %.0f\n",
              ds, nrow(out_dt), length(ds_samples),
              min(colsums), max(colsums), mean(colsums)))

  out_file <- file.path(BULK_DIR, ds, paste0(ds, "_tpm_symbol.tsv"))
  if (file.exists(out_file)) {
    stop(sprintf("Refusing to overwrite existing file: %s", out_file))
  }
  fwrite(out_dt, out_file, sep = "\t", quote = FALSE)
  cat(sprintf("       wrote -> %s\n", out_file))
}

cat(sprintf("\n=== Done (%s) ===\n", Sys.time()))
