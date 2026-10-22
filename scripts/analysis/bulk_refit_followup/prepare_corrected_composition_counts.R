#!/usr/bin/env Rscript
# Per-cohort gene-symbol count matrices from the 774-unit corrected fragment counts,
# for the unchanged InstaPrism script Analysis/Deconvolution/scripts/12_run_bayesprism.R.
#
# Same transformation as Analysis/Deconvolution/scripts/20_prepare_kallisto_bulk_counts.R:
# Ensembl version stripped at the first "." and mapped to the GENCODE v49 gene_name,
# rows sharing a symbol summed, first column named SYMBOL, columns = sample IDs.
# Differences from that script: the input is the corrected count matrix (already integer),
# only the 774 primary units are written, and nothing under Analysis/Deconvolution/bulk
# is touched (the output root must not exist).
#
# Usage: prepare_corrected_composition_counts.R <output_root>

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: prepare_corrected_composition_counts.R <output_root>", call. = FALSE)
out_root <- args[[1]]
if (file.exists(out_root)) stop("Refusing to overwrite: ", out_root, call. = FALSE)
dir.create(out_root, recursive = TRUE)
set.seed(20260923)

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
COUNTS <- file.path(BASE, "figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds")
COUNTS_SHA <- "f931e340546d679dc69dd0d6f712f4e73b1d805e80c741b22d3ed25498c68675"
META <- file.path(BASE, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds")
GTF <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
stopifnot(file.exists(COUNTS), file.exists(META), file.exists(GTF))
sha256_file <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE)[[1L]], "[[:space:]]+")[[1L]][[1L]]
if (sha256_file(COUNTS) != COUNTS_SHA) stop("Corrected counts sha256 differs from the candidate input", call. = FALSE)

gtf_lines <- fread(cmd = paste("zcat", GTF, "| grep -P '^chr.*\\tgene\\t'"), header = FALSE, sep = "\t")
gene_ids <- gsub('.*gene_id "([^"]+)".*', "\\1", gtf_lines$V9)
gene_names <- gsub('.*gene_name "([^"]+)".*', "\\1", gtf_lines$V9)
ens2sym <- unique(data.table(ensembl = sub("\\..*$", "", gene_ids), symbol = gene_names), by = "ensembl")

counts <- readRDS(COUNTS)
meta <- as.data.table(readRDS(META))
stopifnot(setequal(colnames(counts), meta$sample_id), is.numeric(counts), all(counts >= 0), all(counts == round(counts)))
sym <- ens2sym$symbol[match(sub("\\..*$", "", rownames(counts)), ens2sym$ensembl)]
mapped <- !is.na(sym)
dt <- data.table(SYMBOL = sym[mapped], as.data.table(counts[mapped, , drop = FALSE]))
n_before <- nrow(dt)
dt <- dt[, lapply(.SD, sum), by = SYMBOL, .SDcols = setdiff(names(dt), "SYMBOL")]
for (col in setdiff(names(dt), "SYMBOL")) set(dt, j = col, value = as.integer(round(dt[[col]])))
audit <- data.table(
  n_count_rows = nrow(counts), n_mapped_to_symbol = sum(mapped), n_rows_before_aggregation = n_before,
  n_symbols_after_aggregation = nrow(dt), n_samples = ncol(counts),
  total_fragments_input = sum(counts), total_fragments_mapped = sum(counts[mapped, ]))
fwrite(audit, file.path(out_root, "symbol_aggregation_audit.tsv"), sep = "\t")
print(audit)

for (ds in sort(unique(meta$dataset))) {
  ids <- meta[dataset == ds, sample_id]
  dir.create(file.path(out_root, ds), recursive = TRUE)
  fwrite(dt[, c("SYMBOL", ids), with = FALSE], file.path(out_root, ds, paste0(ds, "_counts.tsv")),
         sep = "\t", quote = FALSE)
  cat(sprintf("%s: %d symbols x %d samples\n", ds, nrow(dt), length(ids)))
}
fwrite(meta[, .(sample_id, dataset, group_binary, fibrosis_stage, nas_score, analysis_unit_id)],
       file.path(out_root, "sample_manifest.tsv"), sep = "\t")
fwrite(data.table(input = c("corrected_counts", "meta_matched", "gencode_gtf"), path = c(COUNTS, META, GTF),
                  sha256 = c(COUNTS_SHA, sha256_file(META), sha256_file(GTF))),
       file.path(out_root, "input_manifest.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
cat("COMPOSITION_COUNTS_COMPLETE", out_root, "\n")
