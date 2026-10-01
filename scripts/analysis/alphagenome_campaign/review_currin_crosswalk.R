#!/usr/bin/env Rscript
# Independent direct source accession check. Emit aggregate counts only.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, nzchar(Sys.getenv("SLURM_JOB_ID")))
root <- args[[1]]; out <- args[[2]]
stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE, mode = "0700")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1]])
stopifnot(file.copy(script, file.path(out, "executed_review_currin_crosswalk.R")))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
base <- file.path(root, "GWAS/finemapping/results/alphagenome_campaign/week1-20260915")
asset <- file.path(base, "data/currin-profile-21773775/private_source/liver_ATAC_peakCounts_DESeqDatasets_WASP-filtered_consensusPeaks_autosomalOnly_refinedBoundaries.RData")
stopifnot(unname(tools::md5sum(asset)) == "6f3af5bac0509909d015e8fd99cc1d5e")
env <- new.env(parent = baseenv()); load(asset, envir = env)
object <- env[["peaks.dds"]]
stopifnot(is(object, "DESeqDataSet"), ncol(object) == 138L)
column_ids <- colnames(object)
accession <- as.character(colData(object)$GEO_ID)
stopifnot(!anyNA(accession), !anyDuplicated(accession), !anyDuplicated(column_ids), all(grepl("^GSM[0-9]+$", accession)))
mapping <- read.delim(file.path(base, "data/currin-donor-crosswalk-21773997/public_source_accession_crosswalk.tsv"), check.names = FALSE, colClasses = "character", na.strings = "")
stopifnot(nrow(mapping) == 138L,
          identical(mapping$deposited_count_column, column_ids),
          identical(mapping$author_GEO_ID, accession),
          identical(mapping$GEO_ID_GSM_token, accession),
          all(is.na(mapping$column_GSM_token)), all(mapping$exact_single_accession_field == "TRUE"))
memberships <- list()
for (series in c("GSE26105", "GSE264684")) {
  path <- file.path(base, "data/source_metadata", paste0(series, "_series_metadata.txt"))
  text <- readLines(path, warn = FALSE)
  candidates <- text[startsWith(text, "!Series_sample_id = ")]
  ids <- unique(trimws(sub("^[^=]*=", "", candidates)))
  stopifnot(length(ids) > 0L, all(grepl("^GSM[0-9]+$", ids)))
  member <- accession %in% ids
  column <- if (series == "GSE26105") "in_GSE26105_genotype_series" else "in_GSE264684_ATAC_series"
  stopifnot(identical(member, mapping[[column]] == "TRUE"))
  memberships[[series]] <- list(public_series_sample_count = length(ids), exact_source_accession_matches = sum(member), metadata_sha256 = digest(file = path, algo = "sha256"))
}
stopifnot(!any(grepl("GSM[0-9]+", column_ids)))
old <- read_json(file.path(base, "data/currin-donor-crosswalk-21773997/receipt.json"), simplifyVector = TRUE)
stopifnot(old$count_columns == 138L, old$both_present_single_token_columns == 0L,
          is.null(old$column_GSM_and_GEO_ID_agree_where_comparable),
          old$tokens_matching_GSE26105 == memberships$GSE26105$exact_source_accession_matches,
          old$tokens_matching_GSE264684 == memberships$GSE264684$exact_source_accession_matches)
write_json(list(status = "independent_direct_source_accession_reconstruction_passed", job_id = Sys.getenv("SLURM_JOB_ID"),
  unique_count_columns = length(column_ids), unique_author_GEO_ID = length(accession), memberships = memberships,
  count_column_GSM_tokens = 0L, comparable_column_GSM_identifiers = 0L, agreement_when_no_comparables = "unavailable",
  units = "Author-deposited public accession crosswalk only; no genotype calls, phase, independently verified biological identity or participant molecular interaction established",
  participant_clinical_or_PC_values_exported = FALSE, individual_accessions_reexported = FALSE,
  source_asset_md5 = unname(tools::md5sum(asset))), file.path(out, "checks.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
