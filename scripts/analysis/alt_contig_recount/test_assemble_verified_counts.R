#!/usr/bin/env Rscript
# Small synthetic inputs test whether evidence and integer molecules survive assembly.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
root <- tempfile("verified-count-assembly-")
dir.create(root)
annotation_path <- file.path(root, "resolution_summary.json")
writeLines(toJSON(list(annotation_resolution_complete = TRUE, n_unresolved_annotations = 0L), auto_unbox = TRUE), annotation_path)
fwrite(data.table(counting_gene_id = "G1"), file.path(root, "annotation_correspondence.tsv"), sep = "\t")
sha <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), "[[:space:]]+")[[1]][[1]]
count_root <- file.path(root, "samples")
dir.create(count_root)
supported <- c("source_labels_plus_stable_filtered_genetic_identity_cluster_not_verified_person",
               "authoritative_source_participant_pairing_not_screened_in_A1")
roster <- data.table(sample_id = c("S1", "S2"), analysis_unit_id = c("U1", "U2"), unit_evidence = supported)
for (sample in roster$sample_id) {
  path <- file.path(count_root, sample)
  dir.create(path)
  x <- data.table(gene_id = "G1", count = 2L)
  setnames(x, "count", sample)
  fwrite(x, file.path(path, "gene_counts.tsv"), sep = "\t")
  writeLines(toJSON(list(sample_id = sample, annotation_summary_sha256 = sha(annotation_path),
                        counts_tsv_sha256 = sha(file.path(path, "gene_counts.tsv")),
                        maximum_contribution_per_molecule = 1L, fractional_renormalization = FALSE),
                    auto_unbox = TRUE), file.path(path, "count_completion.json"))
}
roster_path <- file.path(root, "roster.tsv")
fwrite(roster, roster_path, sep = "\t")
run <- function(label, expected) {
  output <- file.path(root, label)
  log <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", "scripts/analysis/alt_contig_recount/assemble_verified_counts.R",
          shQuote(c(annotation_path, roster_path, count_root, output))), stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  if (is.null(status)) status <- 0L
  stopifnot((status == 0L) == expected)
  output
}
valid <- run("valid", TRUE)
certificate <- fromJSON(file.path(valid, "count_completion.json"))
stopifnot(setequal(certificate$unit_definitions, supported),
          !grepl("verified biological individual", certificate$unit, fixed = TRUE),
          certificate$primary_sample_roster_sha256 == sha(roster_path),
          identical(unname(readRDS(file.path(valid, "counts_raw.rds"))), matrix(2L, 1L, 2L)))
roster$unit_evidence[1] <- "verified_person_without_support"
fwrite(roster, roster_path, sep = "\t")
run("unsupported", FALSE)
roster$unit_evidence[1] <- supported[1]
roster$analysis_unit_id[2] <- "U1"
fwrite(roster, roster_path, sep = "\t")
run("duplicate_unit", FALSE)
roster$analysis_unit_id[2] <- "U2"
roster$unit_evidence[1] <- ""
fwrite(roster, roster_path, sep = "\t")
run("empty_evidence", FALSE)
cat("PASS: mixed supported units retained; unsupported, duplicate and empty units rejected\n")
