#!/usr/bin/env Rscript
# Author-deposited public accession linkage only; no clinical/genotype-PC export.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, nzchar(Sys.getenv("SLURM_JOB_ID")))
root <- args[[1]]
out <- args[[2]]
if (dir.exists(out)) stop("Refusing to overwrite named output")
dir.create(out, recursive = TRUE, mode = "0700")
executed_script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1]])
stopifnot(file.copy(executed_script, file.path(out, "executed_data_currin_donor_crosswalk.R"), overwrite = FALSE))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
asset <- file.path(root, "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/currin-profile-21773775/private_source/liver_ATAC_peakCounts_DESeqDatasets_WASP-filtered_consensusPeaks_autosomalOnly_refinedBoundaries.RData")
stopifnot(unname(tools::md5sum(asset)) == "6f3af5bac0509909d015e8fd99cc1d5e")
source_dir <- file.path(root, "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data/source_metadata")
sources <- setNames(file.path(source_dir, paste0(c("GSE26105", "GSE264684"), "_series_metadata.txt")), c("GSE26105", "GSE264684"))
gsm_tokens <- function(x) unique(regmatches(x, gregexpr("(?<![A-Za-z0-9])GSM[0-9]+(?![A-Za-z0-9])", x, perl = TRUE))[[1]])
sets <- lapply(sources, function(p) {
  text <- readLines(p, warn = FALSE)
  text <- text[grepl("^!Series_sample_id = ", text)]
  unique(sub("^!Series_sample_id = ", "", trimws(text)))
})
stopifnot(all(lengths(sets) > 0L))
env <- new.env(parent = baseenv())
load(asset, envir = env)
dds <- env[["peaks.dds"]]
stopifnot(is(dds, "DESeqDataSet"), "GEO_ID" %in% names(colData(dds)),
          ncol(dds) == 138L, !anyDuplicated(colnames(dds)), !anyNA(colnames(dds)))
geo <- as.character(colData(dds)$GEO_ID)
ids <- colnames(dds)
geo_tokens <- lapply(geo, function(x) if (is.na(x)) character() else gsm_tokens(x))
column_tokens <- lapply(ids, gsm_tokens)
comparable <- lengths(geo_tokens) == 1L & lengths(column_tokens) == 1L
agreements <- vapply(which(comparable), function(i) identical(geo_tokens[[i]], column_tokens[[i]]), logical(1))
mapping <- do.call(rbind, lapply(seq_along(ids), function(i) {
  g <- geo_tokens[[i]]
  column_gsm <- column_tokens[[i]]
  if (!length(g)) g <- NA_character_
  do.call(rbind, lapply(g, function(s) data.frame(
    deposited_count_column = ids[[i]], author_GEO_ID = geo[[i]], GEO_ID_GSM_token = s,
    column_GSM_token = paste(column_gsm, collapse = ";"),
    in_GSE26105_genotype_series = !is.na(s) && s %in% sets$GSE26105,
    in_GSE264684_ATAC_series = !is.na(s) && s %in% sets$GSE264684,
    exact_single_accession_field = !is.na(s) && identical(geo[[i]], s),
    stringsAsFactors = FALSE)))
}))
write.table(mapping, file.path(out, "public_source_accession_crosswalk.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
summary <- list(status = "author_deposited_public_accession_crosswalk_inspected", job_id = Sys.getenv("SLURM_JOB_ID"),
  count_columns = length(ids), unique_count_columns = length(unique(ids)), missing_GEO_ID = sum(is.na(geo) | geo == ""),
  rows = nrow(mapping), unique_GEO_ID_GSM_tokens = length(unique(na.omit(mapping$GEO_ID_GSM_token))),
  tokens_matching_GSE26105 = sum(mapping$in_GSE26105_genotype_series),
  tokens_matching_GSE264684 = sum(mapping$in_GSE264684_ATAC_series),
  tokens_matching_neither = sum(!is.na(mapping$GEO_ID_GSM_token) & !mapping$in_GSE26105_genotype_series & !mapping$in_GSE264684_ATAC_series),
  columns_with_no_GSM_token_in_count_column_ID = sum(lengths(column_tokens) == 0L),
  columns_with_no_GSM_token_in_GEO_ID = sum(lengths(geo_tokens) == 0L),
  columns_with_multiple_GSM_tokens_in_count_column_ID = sum(lengths(column_tokens) > 1L),
  columns_with_multiple_GSM_tokens_in_GEO_ID = sum(lengths(geo_tokens) > 1L),
  both_present_single_token_columns = sum(comparable),
  both_present_single_token_agreements = sum(agreements),
  both_present_single_token_conflicts = sum(!agreements),
  column_GSM_and_GEO_ID_agree_where_comparable = if (length(agreements)) all(agreements) else NA,
  source_series_sample_counts = as.list(lengths(sets)),
  evidence_scope = "Exact author-deposited count-column GEO_ID and public GEO series accession membership. These fields are not an independently reverified biological identity or genotype comparison. An ATAC-series accession must not be interpreted as a GSE26105 genotype join.",
  phase_limit = "No genotype calls, retained phase, phase sets, imputation quality or final VCF inspected. Public accession linkage alone cannot establish phased molecular backgrounds.",
  privacy = "Only identifiers present in the public source object/metadata are exported; no clinical covariates, genotype PCs, sequencing, or allele calls exported.",
  source_asset_md5 = unname(tools::md5sum(asset)), source_metadata_sha256 = lapply(sources, function(p) digest(file = p, algo = "sha256")),
  public_accession_mapping_R_serialized_sha256 = digest(mapping, algo = "sha256"))
write_json(summary, file.path(out, "receipt.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
cat("Public source accession crosswalk inspection complete\n")
