#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- "
collate_salmon_quant.R --quant-dir <DIR> --metadata <CSV> --output-dir <DIR> [--sample-column <col>] [--threads <int>]

Required arguments:
  --quant-dir        Directory containing per-sample Salmon outputs (each with quant.sf)
  --metadata         Metadata CSV/TSV file used for the study (must contain sample identifiers)
  --output-dir       Destination directory for merged matrices

Optional arguments:
  --sample-column    Metadata column with sample identifiers (default: run)
  --threads          Number of threads for tximport (default: 4)
"

if (length(args) == 0) {
  cat(usage, "\n")
  quit(status = 1)
}

parse_args <- function(args) {
  res <- list(sample_column = "run", threads = 4L)
  i <- 1
  while (i <= length(args)) {
    flag <- args[[i]]
    if (flag %in% c("--quant-dir", "--metadata", "--output-dir", "--sample-column", "--threads")) {
      if (i == length(args)) {
        stop("Missing value for ", flag)
      }
      value <- args[[i + 1]]
      if (flag == "--quant-dir") res$quant_dir <- value
      if (flag == "--metadata") res$metadata <- value
      if (flag == "--output-dir") res$output_dir <- value
      if (flag == "--sample-column") res$sample_column <- tolower(value)
      if (flag == "--threads") res$threads <- as.integer(value)
      i <- i + 2
    } else if (flag %in% c("--help", "-h")) {
      cat(usage, "\n")
      quit(status = 0)
    } else {
      stop("Unknown flag: ", flag)
    }
  }
  required <- c("quant_dir", "metadata", "output_dir")
  missing <- setdiff(required, names(res))
  if (length(missing) > 0) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "))
  }
  res
}

opts <- parse_args(args)

if (!dir.exists(opts$quant_dir)) {
  stop("Quant directory does not exist: ", opts$quant_dir)
}
if (!file.exists(opts$metadata)) {
  stop("Metadata file not found: ", opts$metadata)
}

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading metadata: ", opts$metadata)
metadata <- suppressMessages(readr::read_csv(opts$metadata, show_col_types = FALSE))
colnames(metadata) <- make.names(tolower(colnames(metadata)))

if (!opts$sample_column %in% colnames(metadata)) {
  stop("Sample column '", opts$sample_column, "' not present in metadata")
}

sample_ids <- unique(trimws(metadata[[opts$sample_column]]))
sample_ids <- sample_ids[!is.na(sample_ids) & sample_ids != ""]

if (length(sample_ids) == 0) {
  stop("No sample identifiers found in metadata column '", opts$sample_column, "'")
}

message("Expecting quantifications for ", length(sample_ids), " samples")

quant_paths <- file.path(opts$quant_dir, sample_ids, "quant.sf")

missing <- quant_paths[!file.exists(quant_paths)]
if (length(missing) > 0) {
  stop(length(missing), " quant.sf files are missing. Example missing file: ", missing[[1]])
}

names(quant_paths) <- sample_ids

message("Importing Salmon quantifications via tximport (threads = ", opts$threads, ")")
txi <- tximport(
  files = quant_paths,
  type = "salmon",
  txOut = TRUE,
  countsFromAbundance = "no",
  dropInfReps = TRUE,
  ignoreTxVersion = TRUE,
  ignoreAfterBar = TRUE
)

counts <- txi$counts
abundance <- txi$abundance
lengths <- txi$length

stopifnot(ncol(counts) == length(sample_ids))

counts_file <- file.path(opts$output_dir, "transcript_counts_matrix.tsv")
abundance_file <- file.path(opts$output_dir, "transcript_tpm_matrix.tsv")
length_file <- file.path(opts$output_dir, "transcript_length_matrix.tsv")
txi_file <- file.path(opts$output_dir, "tximport_salmon.rds")

message("Writing counts to ", counts_file)
write_tsv(
  as_tibble(counts, rownames = "GeneID"),
  counts_file
)

message("Writing TPMs to ", abundance_file)
write_tsv(
  as_tibble(abundance, rownames = "GeneID"),
  abundance_file
)

message("Writing effective lengths to ", length_file)
write_tsv(
  as_tibble(lengths, rownames = "GeneID"),
  length_file
)

saveRDS(txi, txi_file)

summary_df <- tibble(
  metric = c("transcripts", "samples"),
  value = c(nrow(counts), ncol(counts))
)
summary_file <- file.path(opts$output_dir, "salmon_quant_summary.tsv")
write_tsv(summary_df, summary_file)

message("Collation complete: ", nrow(counts), " transcripts x ", ncol(counts), " samples")
