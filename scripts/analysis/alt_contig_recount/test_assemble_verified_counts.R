#!/usr/bin/env Rscript
# Small synthetic inputs test whether evidence and integer molecules survive assembly.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
root <- tempfile("verified-count-assembly-")
dir.create(root)
annotation_path <- file.path(root, "resolution_summary.json")
fwrite(data.table(gene_id = c("G1", "reference_artifact", "H1"), counting_gene_id = c("G1", "", "H1"),
                 disposition = c("primary_gene", "excluded_reference", "held_untestable")),
       file.path(root, "annotation_correspondence.tsv"), sep = "\t")
sha <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), "[[:space:]]+")[[1]][[1]]
grouped <- file.path(root, "fixture.gtf.gz")
connection <- gzfile(grouped, "wt")
writeLines('chr1\tfixture\texon\t1\t20\t.\t+\t.\tgene_id "G1";', connection)
close(connection)
grouped_sha <- sha(grouped)
correspondence_sha <- sha(file.path(root, "annotation_correspondence.tsv"))
writeLines(toJSON(list(annotation_resolution_complete = TRUE, n_unresolved_annotations = 0L,
  n_held_untestable_annotations = 1L,
  grouped_counting_annotation = grouped, grouped_counting_annotation_sha256 = grouped_sha,
  annotation_correspondence_sha256 = correspondence_sha), auto_unbox = TRUE), annotation_path)
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
  x[[sample]] <- 3L
  fwrite(x, file.path(path, "gene_counts_single_candidate_credit.tsv"), sep = "\t")
  synthetic_source <- list(path = file.path(path, "fixture_source.bam"), size_bytes = 10, mtime_ns = "123")
  execution <- list(status = "complete", sample_id = sample, layout = "SINGLE",
    featurecounts_returncode = 0L, samtools_sort_returncode = 0L,
    annotation_summary_sha256 = sha(annotation_path), grouped_counting_annotation_sha256 = grouped_sha,
    annotation_correspondence_sha256 = correspondence_sha,
    source_bam_before = synthetic_source, source_bam_after = synthetic_source,
    featurecounts_argv = c("/fixture/featureCounts", "-T", "1", "-s", "2", "-M", "-R", "BAM",
      "-a", grouped, "-o", file.path(path, "fixture_counts.tsv"), synthetic_source$path),
    versions = list(python = "synthetic fixture", featureCounts = "synthetic fixture", samtools = "synthetic fixture"))
  execution_path <- file.path(path, "assignment_execution.json")
  writeLines(toJSON(execution, auto_unbox = TRUE), execution_path)
  writeLines(toJSON(list(sample_id = sample, annotation_summary_sha256 = sha(annotation_path),
                        counts_tsv_sha256 = sha(file.path(path, "gene_counts.tsv")),
                        single_candidate_credit_counts_tsv_sha256 = sha(file.path(path, "gene_counts_single_candidate_credit.tsv")),
                        maximum_contribution_per_molecule = 1L, fractional_renormalization = FALSE,
                        grouped_counting_annotation_sha256 = grouped_sha, annotation_correspondence_sha256 = correspondence_sha,
                        assignment_execution_sha256 = sha(execution_path), source_bam = synthetic_source,
                        layout = "SINGLE", count_reader_source_sha256 = sha("scripts/analysis/alt_contig_recount/count_verified_fragments.py")),
                    auto_unbox = TRUE), file.path(path, "count_completion.json"))
}
roster_path <- file.path(root, "roster.tsv")
fwrite(roster, roster_path, sep = "\t")
run <- function(label, expected, arm = character()) {
  output <- file.path(root, label)
  log <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", "scripts/analysis/alt_contig_recount/assemble_verified_counts.R",
          shQuote(c(annotation_path, roster_path, count_root, output, arm))), stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  if (is.null(status)) status <- 0L
  if ((status == 0L) != expected) cat(log, sep = "\n")
  stopifnot((status == 0L) == expected)
  output
}
valid <- run("valid", TRUE)
certificate <- fromJSON(file.path(valid, "count_completion.json"))
stopifnot(setequal(certificate$unit_definitions, supported),
          !grepl("verified biological individual", certificate$unit, fixed = TRUE),
          certificate$primary_sample_roster_sha256 == sha(roster_path),
          identical(unname(readRDS(file.path(valid, "counts_raw.rds"))), matrix(2L, 1L, 2L)),
          identical(rownames(readRDS(file.path(valid, "counts_raw.rds"))), "G1"),
          certificate$n_held_untestable_annotations == 1L, certificate$count_arm == "primary")
credited <- run("single_candidate_credit", TRUE, "single_candidate_credit")
stopifnot(identical(unname(readRDS(file.path(credited, "counts_raw.rds"))), matrix(3L, 1L, 2L)),
          fromJSON(file.path(credited, "count_completion.json"))$count_arm == "single_candidate_credit")
run("unknown_arm", FALSE, "name_fold")
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
cat("PASS: excluded reference and held rows omitted; both count arms assembled; mixed supported units retained; unsupported, duplicate and empty units rejected\n")
