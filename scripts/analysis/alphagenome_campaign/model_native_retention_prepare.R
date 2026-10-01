#!/usr/bin/env Rscript
# Internal native regional-activity target. No VST or source-wide fitted offsets.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L, nzchar(Sys.getenv("SLURM_JOB_ID")))
source <- args[[1]]; inspection <- args[[2]]; peaks_file <- args[[3]]
python <- args[[4]]; selector <- args[[5]]; out <- args[[6]]
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
stopifnot(!dir.exists(out)); dir.create(out, recursive = TRUE)
audit <- read_json(inspection, simplifyVector = TRUE)
stopifnot(isTRUE(audit$measurement_identity_candidate))
env <- new.env(parent = baseenv()); load(source, envir = env)
dds <- env[["peaks.dds"]]
stopifnot(is(dds, "DESeqDataSet"), "counts" %in% assayNames(dds), ncol(dds) == 138L, nrow(dds) == 349685L,
          !anyDuplicated(rownames(dds)), !anyDuplicated(colnames(dds)))
stopifnot(digest(rownames(dds), algo = "sha256") == audit$objects[["peaks.dds"]]$row_identity$ordered_R_serialized_sha256,
          digest(colnames(dds), algo = "sha256") == audit$objects[["peaks.dds"]]$donor_column_identity$ordered_R_serialized_sha256,
          unname(tools::md5sum(source)) == "6f3af5bac0509909d015e8fd99cc1d5e")
axis <- data.frame(peak_id = rownames(dds))
axis_file <- file.path(out, "source_peak_axis.tsv")
write.table(axis, axis_file, sep = "\t", quote = FALSE, row.names = FALSE)
# The selector sees only coordinates and count-matrix row identities, never values.
code <- system2(python, c(shQuote(selector), "--peaks", shQuote(peaks_file),
                        "--source-axis", shQuote(axis_file), "--out", shQuote(out)))
stopifnot(code == 0L)
panel <- read.delim(file.path(out, "panel.tsv"), check.names = FALSE)
stopifnot(nrow(panel) == 1024L, !anyDuplicated(panel$peak_id), all(panel$fold != 0L))
counts <- assay(dds, "counts")
stopifnot(identical(dim(counts), dim(dds)))
totals <- colSums(counts)
stopifnot(all(is.finite(totals)), all(totals > 0))
index <- match(panel$peak_id, rownames(dds)); stopifnot(!anyNA(index))
raw <- as.matrix(counts[index, , drop = FALSE])
stopifnot(all(is.finite(raw)), all(raw >= 0), all(raw == floor(raw)))
target <- log2(1 + sweep(raw, 2L, totals, "/") * 1e6)
stopifnot(all(is.finite(target)))
private <- file.path(out, "private_targets"); dir.create(private)
for (name in c("raw_counts", "log2cpm")) {
  values <- if (name == "raw_counts") raw else target
  dest <- gzfile(file.path(private, paste0(name, ".tsv.gz")), "wt")
  write.table(data.frame(peak_id = panel$peak_id, values, check.names = FALSE), dest,
              sep = "\t", quote = FALSE, row.names = FALSE); close(dest)
}
write.table(data.frame(donor_id = colnames(dds), full_peak_library_count_total = totals),
            file.path(private, "donor_library_totals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
write_json(list(status = "native_profile_targets_prepared", donors = ncol(dds), regions = nrow(panel),
  full_peak_library_regions = nrow(dds), units = "donor_mean_of_log2_1plus_CPM_from_raw_WASP_filtered_counts",
  denominator = "sum_every_deposited_peak_per_donor_including_fold0; not_total_sequenced_or_offpeak_inclusive_reads",
  source_VST_or_normalization_factors_used = FALSE, source_GC_correction_recreated = FALSE,
  full_library_normalization = "fixed_measurement_definition; no_library_totals_used_as_predictor_features",
  donor_independent_crosswalk_verified = FALSE, source_family = "same_Currin_development_family",
  source_exposure = "foundation_checkpoint_exposure_unresolved",
  coordinates = "source_peakID_join_to_exact_GRCh38_BED0; inspection_retained",
  source_md5 = unname(tools::md5sum(source)), inspection_sha256 = digest(file = inspection, algo = "sha256"),
  donor_axis_R_sha256 = digest(colnames(dds), algo = "sha256"),
  panel_sha256 = digest(file = file.path(out, "panel.tsv"), algo = "sha256"),
  redistribution = "internal_source_derived_donor_values_only_no_public_rehosting"),
  file.path(out, "targets_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
