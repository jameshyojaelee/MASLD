#!/usr/bin/env Rscript
# Compare native window models and assembled edges across CPU counts.
suppressPackageStartupMessages({
  library(cicero); library(monocle); library(Matrix); library(data.table); library(jsonlite)
})
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Run within a compute allocation")
arguments <- commandArgs(trailingOnly = TRUE)
arg <- function(name) arguments[[match(name, arguments) + 1L]]
source_path <- normalizePath(arg("--source"))
real_input <- normalizePath(arg("--real-input"))
output <- arg("--output")
if (file.exists(output)) stop("Output already exists")
dir.create(output, recursive = TRUE)
for (expression in parse(source_path)) {
  if (is.call(expression) && identical(expression[[1L]], quote(`<-`)) &&
      length(expression) == 3L && is.symbol(expression[[2L]]) && is.call(expression[[3L]]) &&
      identical(expression[[3L]][[1L]], quote(`function`))) eval(expression)
}
native <- getFromNamespace("generate_cicero_models", "cicero")
source_hash <- unname(tools::md5sum(source_path))
native_hash <- unname(tools::md5sum(system.file("R", "cicero.rdb", package = "cicero")))
make_cds <- function(counts, peaks) {
  rownames(counts) <- paste(peaks$chrom, peaks$start0, peaks$end, sep = "_")
  colnames(counts) <- paste0("metacell_", seq_len(ncol(counts)))
  cds <- newCellDataSet(counts,
    phenoData = new("AnnotatedDataFrame", data = data.frame(row.names = colnames(counts))),
    featureData = new("AnnotatedDataFrame", data = data.frame(
      gene_short_name = rownames(counts), chr = sub("^chr", "", peaks$chrom),
      bp1 = peaks$start0, bp2 = peaks$end, row.names = rownames(counts))),
    expressionFamily = VGAM::negbinomial.size())
  normalize_cicero_cds(cds)
}
compare <- function(cds, coordinates, weights, label, maximum_windows = 24L) {
  compat <- new.env(parent = asNamespace("cicero"))
  compat$t <- window_transpose
  compat$cov <- function(x, ...) donor_covariance(x, weights)
  eligible <- native_model_windows(cds, 500000, coordinates)
  eligible <- eligible[seq_len(min(length(eligible), maximum_windows))]
  stopifnot(length(eligible) >= 3L)
  compat$generate_windows <- function(...) eligible
  serial <- native
  environment(serial) <- compat
  parallel <- native_model_workers(serial, 3L)
  serial_start <- proc.time()[["elapsed"]]
  models_serial <- serial(cds, distance_parameter = 0.1, window = 500000,
                          genomic_coords = coordinates)
  serial_elapsed <- proc.time()[["elapsed"]] - serial_start
  parallel_start <- proc.time()[["elapsed"]]
  models_parallel <- parallel(cds, distance_parameter = 0.1, window = 500000,
                              genomic_coords = coordinates)
  parallel_elapsed <- proc.time()[["elapsed"]] - parallel_start
  stopifnot(identical(models_serial, models_parallel))
  edges_serial <- assemble_connections(models_serial)
  edges_parallel <- assemble_connections(models_parallel)
  stopifnot(identical(edges_serial, edges_parallel))
  fwrite(edges_serial, file.path(output, paste0(label, "_edges.tsv")), sep = "\t")
  list(label = label, windows = length(eligible), models_identical = TRUE,
       edges_identical = TRUE, edge_rows = nrow(edges_serial),
       serial_elapsed_seconds = serial_elapsed, parallel_elapsed_seconds = parallel_elapsed,
       metacells = ncol(cds), peaks = nrow(cds), donor_weight_range = range(weights))
}
set.seed(42)
fixture_counts <- Matrix(matrix(rpois(800L * 30L, 5), nrow = 800L), sparse = TRUE)
fixture_peaks <- data.frame(chrom = "chr1", start0 = seq(0, by = 10000, length.out = 800L))
fixture_peaks$end <- fixture_peaks$start0 + 500L
fixture <- make_cds(fixture_counts, fixture_peaks)
fixture_donors <- rep(c("a", "b", "c", "d", "e"), c(3, 4, 5, 8, 10))
fixture_weights <- setNames(1 / as.numeric(table(fixture_donors)[fixture_donors]), colnames(fixture))
fixture_result <- compare(fixture, data.frame(chr = "chr1", length = 8500000), fixture_weights, "fixture")
stopifnot(inherits(try(native_model_workers(native, 0), silent = TRUE), "try-error"))
stopifnot(inherits(try(native_model_workers(native, 1.5), silent = TRUE), "try-error"))
counts <- as(readMM(file.path(real_input, "counts.mtx")), "dgCMatrix")
peaks <- fread(file.path(real_input, "peaks.tsv"))
donors <- fread(file.path(real_input, "donors.tsv"), colClasses = c(metacell_id = "character", donor_id = "character"))
stopifnot(nrow(counts) == nrow(peaks), ncol(counts) == nrow(donors))
real <- make_cds(counts, peaks)
weights <- setNames(1 / as.numeric(table(donors$donor_id)[donors$donor_id]), colnames(real))
chromosome <- peaks$chrom[[1L]]
real <- real[peaks$chrom == chromosome, ]
coordinates <- data.frame(chr = chromosome, length = max(peaks$end[peaks$chrom == chromosome]) + 500000)
real_result <- compare(real, coordinates, weights, "actual_training")
stopifnot(identical(source_hash, unname(tools::md5sum(source_path))))
stopifnot(identical(native_hash, unname(tools::md5sum(system.file("R", "cicero.rdb", package = "cicero")))))
write_json(list(status = "passed", tests = list(fixture_result, real_result),
  native_package_version = as.character(packageVersion("cicero")),
  source = source_path, source_md5 = source_hash, installed_package_unchanged = TRUE,
  only_native_worker_literal_changed = TRUE, distance_calibration_remains_serial = TRUE,
  inference = "computation_equivalence_only_no_new_biological_claim"),
  file.path(output, "verification.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
cat("PASS: native models and assembled edges are identical at one and three workers\n")
