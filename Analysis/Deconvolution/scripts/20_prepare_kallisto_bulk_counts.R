#!/usr/bin/env Rscript
# 20_prepare_kallisto_bulk_counts.R
# ---------------------------------------------------------------------------
# Extract per-cohort gene-symbol count matrices from the kallisto canonical
# all_cohorts_gene_counts.tsv.gz for BayesPrism deconvolution.
#
# The per-cohort TSVs match the format expected by 12_run_bayesprism.R:
#   - Rows = gene symbols, first column header = "SYMBOL"
#   - Columns = SRR sample IDs
#   - Values = integer (rounded tximport) counts
#
# Input:
#   - RNA-seq/results/kallisto/all_cohorts_gene_counts.tsv.gz
#   - RNA-seq/Human/.../results/integration/meta_matched.rds (sample->dataset map)
#   - GENCODE v49 GTF (Ensembl -> gene symbol mapping)
#
# Output:
#   - Analysis/Deconvolution/bulk/{COHORT}/{COHORT}_counts.tsv (one per cohort)
#
# Usage: Rscript 20_prepare_kallisto_bulk_counts.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

cat("=== 20: Prepare per-cohort kallisto bulk counts for BayesPrism ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
KALL_FILE <- file.path(BASE, "RNA-seq/results/kallisto/all_cohorts_gene_counts.tsv.gz")
META_FILE <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/meta_matched.rds")
GTF_FILE  <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
BULK_DIR  <- file.path(BASE, "Analysis/Deconvolution/bulk")

stopifnot(file.exists(KALL_FILE), file.exists(META_FILE), file.exists(GTF_FILE))

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

# --- 2. Load kallisto counts ---
cat("\nLoading kallisto counts...\n")
kc <- fread(KALL_FILE)
gene_col <- kc$gene_id
kc[, gene_id := NULL]
cat(sprintf("  Kallisto matrix: %d genes x %d samples\n", length(gene_col), ncol(kc)))

# Strip version from kallisto gene IDs
gene_base <- sub("\\..*$", "", gene_col)

# Map to symbols
sym <- ens2sym$symbol[match(gene_base, ens2sym$ensembl)]
n_mapped <- sum(!is.na(sym))
cat(sprintf("  Mapped: %d / %d genes to symbols (%.1f%%)\n",
            n_mapped, length(sym), 100 * n_mapped / length(sym)))

# Keep only mapped genes, aggregate duplicates by symbol
kc_dt <- data.table(SYMBOL = sym, kc)
kc_dt <- kc_dt[!is.na(SYMBOL)]

# Aggregate duplicate symbols (sum counts)
n_before <- nrow(kc_dt)
kc_dt <- kc_dt[, lapply(.SD, sum), by = SYMBOL, .SDcols = names(kc_dt)[names(kc_dt) != "SYMBOL"]]
cat(sprintf("  After symbol aggregation: %d rows (was %d)\n", nrow(kc_dt), n_before))

# Round to integers (tximport can produce non-integer estimates)
for (col in setdiff(names(kc_dt), "SYMBOL")) {
  set(kc_dt, j = col, value = as.integer(round(kc_dt[[col]])))
}

# --- 3. Load metadata for sample->dataset mapping ---
cat("\nLoading metadata...\n")
meta <- readRDS(META_FILE)
sample_ds <- meta$dataset
names(sample_ds) <- meta$sample_id

# Map kallisto samples to datasets
sample_cols <- setdiff(names(kc_dt), "SYMBOL")
ds_map <- sample_ds[sample_cols]
cat(sprintf("  Mapped: %d / %d kallisto samples to datasets\n",
            sum(!is.na(ds_map)), length(ds_map)))

# --- 4. Write per-cohort count files ---
cat("\nWriting per-cohort count files...\n")
datasets <- sort(unique(na.omit(ds_map)))

for (ds in datasets) {
  ds_samples <- sample_cols[ds_map == ds & !is.na(ds_map)]
  if (length(ds_samples) == 0) {
    cat(sprintf("  %s: no samples, skipping\n", ds))
    next
  }

  out_dir <- file.path(BULK_DIR, ds)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Back up existing STAR-era counts
  out_file <- file.path(out_dir, paste0(ds, "_counts.tsv"))
  backup_file <- file.path(out_dir, paste0(ds, "_counts_star_backup.tsv"))
  if (file.exists(out_file) && !file.exists(backup_file)) {
    file.copy(out_file, backup_file)
    cat(sprintf("  %s: backed up STAR counts -> %s\n", ds, basename(backup_file)))
  }

  # Write kallisto counts
  out_dt <- kc_dt[, c("SYMBOL", ds_samples), with = FALSE]
  fwrite(out_dt, out_file, sep = "\t", quote = FALSE)
  cat(sprintf("  %s: wrote %d genes x %d samples -> %s\n",
              ds, nrow(out_dt), length(ds_samples), out_file))
}

cat(sprintf("\n=== Done (%s) ===\n", Sys.time()))
