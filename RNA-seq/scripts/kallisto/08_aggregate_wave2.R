#!/usr/bin/env Rscript
# Wave 2 - Aggregate kallisto transcript-level abundances to gene-level
# counts + TPM using tximport for the 4 disease-only staging cohorts.
#
# Outputs per cohort:
#   RNA-seq/results/kallisto/<cohort>/gene_counts.tsv.gz
#   RNA-seq/results/kallisto/<cohort>/gene_tpm.tsv.gz
# Plus updated pooled matrices across ALL 9 cohorts:
#   RNA-seq/results/kallisto/all_cohorts_gene_counts.tsv.gz
#   RNA-seq/results/kallisto/all_cohorts_gene_tpm.tsv.gz

suppressPackageStartupMessages({
  library(tximport)
  library(data.table)
  library(rtracklayer)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
KALL_ROOT   <- file.path(PROJ, "RNA-seq/results/kallisto")
DRIVER      <- file.path(PROJ, "RNA-seq/scripts/kallisto/sample_driver_wave2.tsv")
WT_KALL     <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/RNA-seq/results/kallisto"
TX2GENE_SRC <- file.path(WT_KALL, "tx2gene_gencode_v49.tsv.gz")
GTF_PATH    <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"

# -------------------------------------------------------------------
# tx2gene: reuse existing or build from GTF
# -------------------------------------------------------------------
TX2GENE_DST <- file.path(KALL_ROOT, "tx2gene_gencode_v49.tsv.gz")
if (file.exists(TX2GENE_DST)) {
  message("[", Sys.time(), "] Using cached tx2gene: ", TX2GENE_DST)
} else if (file.exists(TX2GENE_SRC)) {
  file.copy(TX2GENE_SRC, TX2GENE_DST)
  message("[", Sys.time(), "] Copied tx2gene from worktree")
} else {
  message("[", Sys.time(), "] Building tx2gene from GTF (~3 min)")
  gtf <- rtracklayer::import(GTF_PATH)
  tx  <- gtf[gtf$type == "transcript"]
  tx2g <- data.table(
    TXNAME     = tx$transcript_id,
    GENEID     = tx$gene_id,
    GENESYMBOL = tx$gene_name
  )
  fwrite(tx2g, TX2GENE_DST, sep = "\t")
}
tx2g    <- fread(TX2GENE_DST)
tx2gene <- tx2g[, .(TXNAME, GENEID)]
message("tx2gene: ", nrow(tx2gene), " transcripts -> ",
        uniqueN(tx2gene$GENEID), " genes")

# -------------------------------------------------------------------
# Wave 2: aggregate 4 new cohorts
# -------------------------------------------------------------------
driver <- fread(DRIVER)
message("Driver: ", nrow(driver), " samples across ",
        uniqueN(driver$cohort), " cohorts")

wave2_count_list <- list()
wave2_tpm_list   <- list()

for (cohort_id in unique(driver$cohort)) {
  rows <- driver[cohort == cohort_id]
  tsv_files <- file.path(KALL_ROOT, cohort_id, rows$sample_id, "abundance.tsv")
  names(tsv_files) <- rows$sample_id
  ok <- file.exists(tsv_files)

  if (sum(ok) == 0) {
    message("[", Sys.time(), "] [", cohort_id, "] NO abundance files -- SKIP")
    next
  }
  if (any(!ok)) {
    message("[", Sys.time(), "] [", cohort_id, "] WARNING: missing ",
            sum(!ok), "/", length(ok), " samples; aggregating present subset")
    message("  Missing: ", paste(rows$sample_id[!ok][1:min(5, sum(!ok))],
                                 collapse = ", "), if (sum(!ok) > 5) "...")
    tsv_files <- tsv_files[ok]
  }

  message("[", Sys.time(), "] [", cohort_id, "] tximport ",
          length(tsv_files), " samples")
  txi    <- tximport(tsv_files, type = "kallisto", tx2gene = tx2gene,
                     ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)
  counts <- as.data.frame(txi$counts)
  tpm    <- as.data.frame(txi$abundance)
  counts$gene_id <- rownames(counts)
  tpm$gene_id    <- rownames(tpm)
  setcolorder(setDT(counts), c("gene_id", setdiff(names(counts), "gene_id")))
  setcolorder(setDT(tpm),    c("gene_id", setdiff(names(tpm),    "gene_id")))

  outdir <- file.path(KALL_ROOT, cohort_id)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fwrite(counts, file.path(outdir, "gene_counts.tsv.gz"), sep = "\t")
  fwrite(tpm,    file.path(outdir, "gene_tpm.tsv.gz"),    sep = "\t")

  wave2_count_list[[cohort_id]] <- counts
  wave2_tpm_list[[cohort_id]]   <- tpm
  message("[", Sys.time(), "] [", cohort_id, "] wrote gene_counts + gene_tpm (",
          nrow(counts), " genes x ", ncol(counts) - 1, " samples)")
}

# -------------------------------------------------------------------
# Merge: 5 existing cohorts (from worktree) + 4 new cohorts = 9 total
# -------------------------------------------------------------------
message("\n[", Sys.time(), "] === Merging all 9 cohorts ===")

# Load existing 5-cohort pooled counts/tpm from worktree
wave1_counts <- fread(file.path(WT_KALL, "all_cohorts_gene_counts.tsv.gz"))
wave1_tpm    <- fread(file.path(WT_KALL, "all_cohorts_gene_tpm.tsv.gz"))
message("Wave 1 (5 cohorts): ", nrow(wave1_counts), " genes x ",
        ncol(wave1_counts) - 1, " samples")

# Combine wave2 into single tables
if (length(wave2_count_list) > 0) {
  wave2_counts <- Reduce(function(a, b) merge(a, b, by = "gene_id", all = TRUE),
                         wave2_count_list)
  wave2_tpm    <- Reduce(function(a, b) merge(a, b, by = "gene_id", all = TRUE),
                         wave2_tpm_list)
  message("Wave 2 (", length(wave2_count_list), " cohorts): ",
          nrow(wave2_counts), " genes x ", ncol(wave2_counts) - 1, " samples")

  # Merge wave1 + wave2
  all_counts <- merge(wave1_counts, wave2_counts, by = "gene_id", all = TRUE)
  all_tpm    <- merge(wave1_tpm,    wave2_tpm,    by = "gene_id", all = TRUE)
} else {
  message("[WARN] No wave 2 cohorts aggregated -- keeping wave 1 only")
  all_counts <- wave1_counts
  all_tpm    <- wave1_tpm
}

# Replace NAs with 0 (genes present in one wave but not another)
for (col in names(all_counts)[-1]) {
  set(all_counts, which(is.na(all_counts[[col]])), col, 0)
}
for (col in names(all_tpm)[-1]) {
  set(all_tpm, which(is.na(all_tpm[[col]])), col, 0)
}

message("All 9 cohorts merged: ", nrow(all_counts), " genes x ",
        ncol(all_counts) - 1, " samples")

# Write merged outputs
fwrite(all_counts, file.path(KALL_ROOT, "all_cohorts_gene_counts.tsv.gz"),
       sep = "\t")
fwrite(all_tpm,    file.path(KALL_ROOT, "all_cohorts_gene_tpm.tsv.gz"),
       sep = "\t")
message("[", Sys.time(), "] Wrote merged counts + tpm")

# Sanity checks
n_wave1 <- ncol(wave1_counts) - 1
n_wave2 <- if (length(wave2_count_list) > 0) ncol(wave2_counts) - 1 else 0
n_total <- ncol(all_counts) - 1
message("\nSanity: wave1=", n_wave1, " + wave2=", n_wave2,
        " = ", n_wave1 + n_wave2, " (merged=", n_total, ")")
if (n_total != n_wave1 + n_wave2) {
  warning("Sample count mismatch! Check for duplicate sample IDs.")
}

message("[", Sys.time(), "] aggregation done")
