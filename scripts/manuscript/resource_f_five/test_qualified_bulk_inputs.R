#!/usr/bin/env Rscript
# Scientific input invariants; no disease model is fitted by these checks.
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
source("scripts/manuscript/resource_f_five/qualified_bulk_inputs.R")
root <- tempfile("qualified-bulk-inputs-")
dir.create(root)
counts_path <- file.path(root, "counts.rds")
roster_path <- file.path(root, "roster.tsv")
certificate_path <- file.path(root, "certificate.json")
annotation_path <- file.path(root, "annotation.json")
counts <- matrix(c(1, 2, 3, 4), 2, dimnames = list(c("gene1", "gene2"), c("sample1", "sample2")))
roster <- data.table(sample_id = c("sample1", "sample2"), analysis_unit_id = c("unit1", "unit2"),
                     unit_evidence = "source_linked")
write_inputs <- function() {
  saveRDS(counts, counts_path)
  fwrite(roster, roster_path, sep = "\t")
  write_json(list(annotation_resolution_complete = TRUE, n_unresolved_annotations = 0L),
             annotation_path, auto_unbox = TRUE)
  digest <- strsplit(system2("sha256sum", shQuote(counts_path), stdout = TRUE)[[1L]], " ")[[1L]][[1L]]
  certificate <<- list(annotation_resolution_complete = TRUE, n_unresolved_annotations = 0L,
                       counts_complete = TRUE, counts_sha256 = digest,
                       annotation_summary = annotation_path,
                       annotation_summary_sha256 = strsplit(system2("sha256sum", shQuote(annotation_path), stdout = TRUE)[[1L]], " ")[[1L]][[1L]],
                       primary_sample_roster_sha256 = strsplit(system2("sha256sum", shQuote(roster_path), stdout = TRUE)[[1L]], " ")[[1L]][[1L]])
  write_json(certificate, certificate_path, auto_unbox = TRUE)
}
read_inputs <- function() read_qualified_bulk_inputs(counts_path, roster_path, certificate_path,
                                                     c("sample1", "sample2"))
expect_failure <- function(expr, message) {
  result <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(result, "error"), grepl(message, conditionMessage(result), fixed = TRUE))
}
write_inputs()
stopifnot(identical(read_inputs()$counts, counts))
roster$analysis_unit_id[1L] <- "changed_unit"
fwrite(roster, roster_path, sep = "\t")
expect_failure(read_inputs(), "selected unit roster")
roster$analysis_unit_id[1L] <- "unit1"
write_inputs()
write_json(list(annotation_resolution_complete = FALSE, n_unresolved_annotations = 1L),
           annotation_path, auto_unbox = TRUE)
expect_failure(read_inputs(), "annotation evidence")
certificate$annotation_summary_sha256 <- strsplit(system2("sha256sum", shQuote(annotation_path), stdout = TRUE)[[1L]], " ")[[1L]][[1L]]
write_json(certificate, certificate_path, auto_unbox = TRUE)
expect_failure(read_inputs(), "Source annotation correspondence is incomplete")
write_inputs()
certificate$n_unresolved_annotations <- 1L
write_json(certificate, certificate_path, auto_unbox = TRUE)
expect_failure(read_inputs(), "incomplete")
write_inputs()
certificate$counts_complete <- FALSE
write_json(certificate, certificate_path, auto_unbox = TRUE)
expect_failure(read_inputs(), "incomplete")
write_inputs()
certificate$counts_sha256 <- paste(rep("0", 64), collapse = "")
write_json(certificate, certificate_path, auto_unbox = TRUE)
expect_failure(read_inputs(), "certificate")
roster$analysis_unit_id <- "same_unit"
write_inputs()
expect_failure(read_inputs(), "unique units")
roster$analysis_unit_id <- c("unit1", "unit2")
counts[1L, 1L] <- 0.5
write_inputs()
expect_failure(read_inputs(), "integer fragment")
counts[1L, 1L] <- 1
roster$sample_id[2L] <- "outside_model"
write_inputs()
expect_failure(read_inputs(), "original model roster")
roster$sample_id[2L] <- "sample2"
write_inputs()
meta <- copy(roster)
stopifnot(isTRUE(assert_bulk_unit_alignment(meta, c("sample2", "sample1"))))
meta$analysis_unit_id <- "same_unit"
expect_failure(assert_bulk_unit_alignment(meta, c("sample1", "sample2")), "unique analysis unit")
cat("QUALIFIED_BULK_INPUT_CHECKS\tPASS\t12\n")
