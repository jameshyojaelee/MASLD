#!/usr/bin/env Rscript
# B1 step 3 - Aggregate kallisto transcript-level abundances to gene-level
# counts + TPM using tximport (tximeta unavailable in rnaseq env; tximport is
# functionally equivalent for our STAR-vs-kallisto comparison).
#
# Outputs per cohort:
#   RNA-seq/results/kallisto/<cohort>/gene_counts.tsv.gz
#   RNA-seq/results/kallisto/<cohort>/gene_tpm.tsv.gz
# Plus a pooled matrix across 5 cohorts:
#   RNA-seq/results/kallisto/all_cohorts_gene_counts.tsv.gz
#   RNA-seq/results/kallisto/all_cohorts_gene_tpm.tsv.gz

suppressPackageStartupMessages({
  library(tximport)
  library(data.table)
  library(rtracklayer)
})

WT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation"
KALL_ROOT <- file.path(WT, "RNA-seq/results/kallisto")
DRIVER <- fread(file.path(WT, "RNA-seq/scripts/kallisto/sample_driver.tsv"))

GTF_PATH <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
TX2GENE_CACHE <- file.path(KALL_ROOT, "tx2gene_gencode_v49.tsv.gz")

if (!file.exists(TX2GENE_CACHE)) {
  message("[", Sys.time(), "] Building tx2gene from GTF (one-time, ~3 min)")
  gtf <- rtracklayer::import(GTF_PATH)
  tx <- gtf[gtf$type == "transcript"]
  tx2g <- data.table(
    TXNAME = tx$transcript_id,
    GENEID = tx$gene_id,
    GENESYMBOL = tx$gene_name
  )
  fwrite(tx2g, TX2GENE_CACHE, sep = "\t")
} else {
  message("[", Sys.time(), "] Using cached tx2gene")
}
tx2g <- fread(TX2GENE_CACHE)
tx2gene <- tx2g[, .(TXNAME, GENEID)]
message("tx2gene: ", nrow(tx2gene), " transcripts -> ",
        uniqueN(tx2gene$GENEID), " genes")

cohort_count_list <- list()
cohort_tpm_list <- list()

for (cohort_id in unique(DRIVER$cohort)) {
  rows <- DRIVER[cohort == cohort_id]
  cohort <- cohort_id  # keep the rest of the loop body compatible
  # kallisto 0.51 drops h5 by default; use abundance.tsv (tximport supports
  # both transparently when files end in .tsv via type="kallisto").
  h5_files  <- file.path(KALL_ROOT, cohort, rows$sample_id, "abundance.h5")
  tsv_files <- file.path(KALL_ROOT, cohort, rows$sample_id, "abundance.tsv")
  # Prefer h5 if any sample wrote it; else fall back to tsv.
  if (any(file.exists(h5_files))) {
    files <- h5_files
    msg_kind <- "abundance.h5"
  } else {
    files <- tsv_files
    msg_kind <- "abundance.tsv (no h5; kallisto 0.51 plaintext mode)"
  }
  names(files) <- rows$sample_id
  ok <- file.exists(files)
  if (sum(ok) == 0) {
    message("[", Sys.time(), "] [", cohort, "] no ", msg_kind, " files yet -- skip")
    next
  }
  if (any(!ok)) {
    message("[", Sys.time(), "] [", cohort, "] missing ",
            sum(!ok), "/", length(ok), " files; aggregating present subset")
    files <- files[ok]
  }
  message("[", Sys.time(), "] [", cohort, "] tximport ", length(files),
          " samples (", msg_kind, ")")
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene,
                  ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)
  counts <- as.data.frame(txi$counts)
  tpm    <- as.data.frame(txi$abundance)
  counts$gene_id <- rownames(counts)
  tpm$gene_id    <- rownames(tpm)
  setcolorder(setDT(counts), c("gene_id", setdiff(names(counts), "gene_id")))
  setcolorder(setDT(tpm),    c("gene_id", setdiff(names(tpm),    "gene_id")))
  fwrite(counts, file.path(KALL_ROOT, cohort, "gene_counts.tsv.gz"), sep = "\t")
  fwrite(tpm,    file.path(KALL_ROOT, cohort, "gene_tpm.tsv.gz"),    sep = "\t")
  cohort_count_list[[cohort]] <- counts
  cohort_tpm_list[[cohort]]   <- tpm
  message("[", Sys.time(), "] [", cohort, "] wrote gene_counts + gene_tpm")
}

# Pool across cohorts (outer-join on gene_id; samples become columns)
if (length(cohort_count_list) > 0) {
  pool_counts <- Reduce(function(a, b) merge(a, b, by = "gene_id", all = TRUE),
                        cohort_count_list)
  pool_tpm    <- Reduce(function(a, b) merge(a, b, by = "gene_id", all = TRUE),
                        cohort_tpm_list)
  fwrite(pool_counts, file.path(KALL_ROOT, "all_cohorts_gene_counts.tsv.gz"),
         sep = "\t")
  fwrite(pool_tpm,    file.path(KALL_ROOT, "all_cohorts_gene_tpm.tsv.gz"),
         sep = "\t")
  message("[", Sys.time(), "] Pooled: ", nrow(pool_counts), " genes x ",
          ncol(pool_counts) - 1, " samples")
}

message("[", Sys.time(), "] aggregation done")
