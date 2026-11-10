#!/usr/bin/env Rscript
# Native Cicero on training-only donor-confined aggregates; no reaggregation.
suppressPackageStartupMessages({
  library(cicero)
  library(monocle)
  library(Matrix)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
donor_covariance <- function(x, weights) {
  if (is.null(rownames(x)) || !identical(rownames(x), names(weights))) {
    stop("Covariance observations do not match the recorded metacell weight order")
  }
  if (any(!is.finite(weights)) || any(weights <= 0)) stop("Invalid donor covariance weights")
  if (length(unique(weights)) == 1L) return(stats::cov(x))
  total <- sum(weights)
  denominator <- total - sum(weights^2) / total
  if (!is.finite(denominator) || denominator <= 0) stop("Insufficient covariance observations")
  mean <- colSums(x * weights) / total
  centered <- sweep(x, 2L, mean, "-") * sqrt(weights)
  crossprod(centered) / denominator
}
window_transpose <- function(x) {
  if (is(x, "Matrix")) {
    if (nrow(x) > 200L) stop("Refusing dense transpose outside a Cicero window")
    return(base::t(as.matrix(x)))
  }
  base::t(x)
}
normalize_cicero_cds <- function(cds) {
  cds <- getFromNamespace("estimateSizeFactorsSimp", "cicero")(cds)
  size_factors <- Biobase::pData(cds)$Size_Factor
  if (any(!is.finite(size_factors)) || any(size_factors <= 0)) stop("Invalid native Cicero size factors")
  normalized <- Biobase::exprs(cds) %*% Diagonal(x = 1 / size_factors)
  dimnames(normalized) <- dimnames(Biobase::exprs(cds))
  cds <- Biobase::assayDataElementReplace(cds, "exprs", normalized, validate = TRUE)
  cds
}
map_native_edges <- function(connections, mapping) {
  # data.table::melt leaves the second endpoint as a factor. Named-vector
  # indexing must use its labels rather than the underlying integer codes.
  connections$peak1 <- unname(mapping[as.character(connections$peak1)])
  connections$peak2 <- unname(mapping[as.character(connections$peak2)])
  if (anyNA(connections$peak1) || anyNA(connections$peak2)) stop("Unmapped native Cicero endpoint")
  connections
}
check_edge_geometry <- function(connections, peaks, window) {
  first <- match(connections$peak1, peaks$peak_id)
  second <- match(connections$peak2, peaks$peak_id)
  if (anyNA(first) || anyNA(second)) stop("Unknown peak in Cicero edge geometry")
  maximum_width <- max(peaks$end - peaks$start0)
  if (any(peaks$chrom[first] != peaks$chrom[second]) ||
      any(abs(peaks$start0[first] - peaks$start0[second]) > window + maximum_width)) {
    stop("Cicero edge violates chromosome/native-window geometry")
  }
  invisible(TRUE)
}
native_model_windows <- function(cds, window, genomic_coords, max_elements = 200L) {
  grs <- getFromNamespace("generate_windows", "cicero")(window, genomic_coords)
  fd <- Biobase::fData(cds)
  chromosome <- gsub("chr", "", as.character(fd$chr))
  bp1 <- as.numeric(as.character(fd$bp1))
  bp2 <- as.numeric(as.character(fd$bp2))
  win_chromosome <- gsub("chr", "", as.character(GenomicRanges::seqnames(grs)))
  left <- GenomicRanges::start(grs)
  right <- GenomicRanges::end(grs)
  counts <- integer(length(grs))
  for (chr in unique(win_chromosome)) {
    peaks <- which(chromosome == chr)
    if (!length(peaks)) next
    windows <- which(win_chromosome == chr)
    # Match get_genomic_range's strict endpoint-inside test exactly, including
    # its behavior at window boundaries. Do not change native genomic tiling.
    counts[windows] <- vapply(windows, function(i) {
      sum((bp1[peaks] < right[i] & bp1[peaks] > left[i]) |
          (bp2[peaks] < right[i] & bp2[peaks] > left[i]))
    }, integer(1))
  }
  cat(format(Sys.time(), tz = "UTC"), "native_windows", length(grs),
      "eligible_model_windows", sum(counts > 1L & counts <= max_elements), "\n")
  grs[counts > 1L & counts <= max_elements]
}
native_model_workers <- function(fn, workers = 1L) {
  if (length(workers) != 1L || !is.finite(workers) || workers < 1L || workers != as.integer(workers)) {
    stop("Cicero workers must be one positive integer")
  }
  replacements <- 0L
  replace_workers <- function(node) {
    if (!is.call(node)) return(node)
    if (identical(node[[1L]], quote(parallel::mclapply))) {
      if (is.null(node[["mc.cores"]]) || !identical(as.numeric(node[["mc.cores"]]), 1)) {
        stop("Native Cicero window worker setting differs from the inspected source")
      }
      node[["mc.cores"]] <- as.integer(workers)
      replacements <<- replacements + 1L
      return(node)
    }
    for (i in seq_along(node)[-1L]) node[[i]] <- replace_workers(node[[i]])
    node
  }
  body(fn) <- replace_workers(body(fn))
  if (replacements != 1L) stop("Expected exactly one native Cicero window mclapply")
  fn
}
if ("--test-transpose" %in% args) {
  if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Tests require a compute allocation")
  set.seed(42)
  numeric_donors <- fread(text = "metacell_id\tdonor_id\nm1\t331\nm2\t331\nm3\t334",
                         colClasses = c(metacell_id = "character", donor_id = "character"))
  numeric_totals <- table(numeric_donors$donor_id)
  stopifnot(identical(as.numeric(numeric_totals[numeric_donors$donor_id]), c(2, 2, 1)))
  unequal <- matrix(c(0, 0, 2, 4, 4, 2), ncol = 2L, byrow = TRUE,
                    dimnames = list(c("m1", "m2", "m3"), c("x", "y")))
  weighted <- donor_covariance(unequal, c(m1 = .5, m2 = .5, m3 = 1))
  stopifnot(isTRUE(all.equal(unname(weighted), matrix(c(4.4, 1.6, 1.6, 3.2), 2L), tolerance = 1e-14)))
  stopifnot(identical(donor_covariance(unequal, c(m1 = .5, m2 = .5, m3 = .5)), stats::cov(unequal)))
  stopifnot(inherits(try(donor_covariance(unequal, c(m2 = .5, m1 = .5, m3 = 1)), silent = TRUE), "try-error"))
  edge_fixture <- data.frame(peak1 = "site_b", peak2 = factor("site_c"), coaccess = .5)
  mapped_fixture <- map_native_edges(edge_fixture, c(site_a = "a", site_b = "b", site_c = "c"))
  stopifnot(identical(mapped_fixture$peak1, "b"), identical(mapped_fixture$peak2, "c"))
  geometry_fixture <- data.frame(peak_id = c("b", "c"), chrom = c("chr1", "chr1"), start0 = c(0, 10000000), end = c(500, 10000500))
  stopifnot(inherits(try(check_edge_geometry(mapped_fixture, geometry_fixture, 500000), silent = TRUE), "try-error"))
  fixture <- Matrix(matrix(rpois(200L * 20L, 2), nrow = 200L), sparse = TRUE)
  stopifnot(identical(window_transpose(fixture), base::t(as.matrix(fixture))))
  stopifnot(identical(cov(window_transpose(fixture)), cov(base::t(as.matrix(fixture)))))
  stopifnot(inherits(try(window_transpose(Matrix(0, nrow = 201L, ncol = 2L, sparse = TRUE)), silent = TRUE), "try-error"))
  rownames(fixture) <- paste0("p", seq_len(nrow(fixture)))
  colnames(fixture) <- paste0("m", seq_len(ncol(fixture)))
  fixture_cds <- newCellDataSet(fixture,
    phenoData = new("AnnotatedDataFrame", data = data.frame(row.names = colnames(fixture))),
    featureData = new("AnnotatedDataFrame", data = data.frame(
      gene_short_name = rownames(fixture), chr = rep(c("1", "2"), each = 100L),
      bp1 = seq(0, by = 500, length.out = 200),
      bp2 = seq(500, by = 500, length.out = 200), row.names = rownames(fixture))),
    expressionFamily = VGAM::negbinomial.size())
  sized <- getFromNamespace("estimateSizeFactorsSimp", "cicero")(fixture_cds)
  dense_reference <- base::t(base::t(as.matrix(fixture)) / Biobase::pData(sized)$Size_Factor)
  normalized_fixture <- normalize_cicero_cds(fixture_cds)
  stopifnot(isTRUE(all.equal(as.matrix(Biobase::exprs(normalized_fixture)), dense_reference, tolerance = 1e-14)))
  coords <- data.frame(chr = c("chr1", "chr2"), length = c(125000, 125000))
  all_windows <- getFromNamespace("generate_windows", "cicero")(50000, coords)
  counts <- vapply(seq_along(all_windows), function(i) {
    nrow(Biobase::exprs(getFromNamespace("get_genomic_range", "cicero")(all_windows, fixture_cds, i)))
  }, integer(1))
  expected <- all_windows[counts > 1L & counts <= 200L]
  actual <- native_model_windows(fixture_cds, 50000, coords)
  stopifnot(identical(actual, expected))
  cat("PASS: sparse-window covariance, native normalization and eligible windows equal native references; dense limit enforced\n")
  quit(status = 0L)
}
get_arg <- function(key, default = NULL) {
  at <- match(key, args)
  if (is.na(at)) return(default)
  if (at == length(args)) stop("Missing argument after ", key)
  args[[at + 1L]]
}
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Cicero requires a compute allocation")
output <- get_arg("--output")
if (is.null(output) || file.exists(output)) stop("Output missing or already exists")
set.seed(as.integer(as.double(get_arg("--seed", "42")) %% .Machine$integer.max))
counts <- as(readMM(get_arg("--counts")), "dgCMatrix")
peaks <- fread(get_arg("--peaks"))
if (nrow(counts) != nrow(peaks)) stop("Peak dimensions differ")
site <- paste(peaks$chrom, peaks$start0, peaks$end, sep = "_")
if (anyDuplicated(site)) stop("Duplicate peak coordinates")
rownames(counts) <- site
colnames(counts) <- paste0("metacell_", seq_len(ncol(counts)))
donor_path <- get_arg("--donors")
covariance_weights <- NULL
if (!is.null(donor_path)) {
  donor_table <- fread(donor_path, colClasses = c(metacell_id = "character", donor_id = "character"))
  if (!all(c("metacell_id", "donor_id") %in% names(donor_table)) ||
      anyDuplicated(donor_table$metacell_id) || anyNA(donor_table$donor_id) ||
      any(!nzchar(donor_table$donor_id)) ||
      !identical(as.character(donor_table$metacell_id), colnames(counts))) {
    stop("Donor sidecar must match the exact aggregate column order")
  }
  donor_totals <- table(donor_table$donor_id)
  covariance_weights <- setNames(1 / as.numeric(donor_totals[as.character(donor_table$donor_id)]), colnames(counts))
  donor_table$covariance_weight <- unname(covariance_weights)
  totals <- tapply(donor_table$covariance_weight, donor_table$donor_id, sum)
  if (any(abs(totals - 1) > 1e-12)) stop("Donors do not contribute equal total weight")
  fwrite(donor_table, paste0(output, ".donor_weights.tsv"), sep = "\t")
}
fd <- new("AnnotatedDataFrame", data = data.frame(
  site_name = site, gene_short_name = site,
  chr = sub("^chr", "", peaks$chrom), bp1 = peaks$start0,
  bp2 = peaks$end, row.names = site
))
pd <- new("AnnotatedDataFrame", data = data.frame(cells = colnames(counts), row.names = colnames(counts)))
cds <- newCellDataSet(counts, phenoData = pd, featureData = fd,
                     expressionFamily = VGAM::negbinomial.size())
# Match make_cicero_cds' native aggregate normalization while retaining sparsity.
cds <- normalize_cicero_cds(cds)
window <- as.integer(get_arg("--window", "500000"))
genomic_coords <- data.frame(chr = sort(unique(peaks$chrom)))
genomic_coords$length <- vapply(genomic_coords$chr, function(chrom) {
  max(peaks$end[peaks$chrom == chrom]) + window
}, numeric(1))
# Calibrate native distance penalty on the full training peak universe. Restrict
# inference to requested cis windows only after calibration, avoiding almost
# empty chromosome-window sampling in a small technical pilot.
# Cicero 1.29.2 uses base t() on sparse Matrix objects. Keep native calculations
# and supply a window-only transpose in local function environments, without
# changing the installed package or densifying the full aggregate matrix.
compat <- new.env(parent = asNamespace("cicero"))
compat$t <- window_transpose
if (!is.null(covariance_weights)) {
  compat$cov <- function(x, ...) donor_covariance(x, covariance_weights)
}
for (name in c("estimate_distance_parameter", "find_distance_parameter", "generate_cicero_models")) {
  fn <- getFromNamespace(name, "cicero")
  environment(fn) <- compat
  assign(name, fn, envir = compat)
}
distance_samples <- as.integer(get_arg("--distance-samples", "100"))
workers <- as.integer(get_arg("--workers", "1"))
allocated_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (!is.finite(workers) || workers < 1L || workers > allocated_workers) {
  stop("Cicero workers exceed the current CPU allocation")
}
compat$generate_cicero_models <- native_model_workers(compat$generate_cicero_models, workers)
cat(format(Sys.time(), tz = "UTC"), "distance_calibration_start", distance_samples, "windows\n")
flush.console()
distance_parameters <- compat$estimate_distance_parameter(
  cds, genomic_coords = genomic_coords, window = window, sample_num = distance_samples
)
if (!length(distance_parameters) || any(!is.finite(distance_parameters))) {
  stop("No finite native Cicero distance calibration")
}
scoring <- if ("score_peak" %in% names(peaks)) as.logical(peaks$score_peak) else rep(TRUE, nrow(peaks))
cat(format(Sys.time(), tz = "UTC"), "distance_calibration_complete", length(distance_parameters), "window_models_start", sum(scoring), "scoring_peaks\n")
flush.console()
model_cds <- cds[scoring, ]
compat$generate_windows <- function(window, genomic_coords) {
  native_model_windows(model_cds, window, genomic_coords,
    as.integer(formals(cicero::generate_cicero_models)$max_elements))
}
models <- compat$generate_cicero_models(model_cds,
  distance_parameter = mean(distance_parameters),
  genomic_coords = genomic_coords, window = window
)
connections <- assemble_connections(models)
cat(format(Sys.time(), tz = "UTC"), "window_models_complete", nrow(connections), "links\n")
flush.console()
names(connections)[match(c("Peak1", "Peak2"), names(connections))] <- c("peak1", "peak2")
if (!all(c("peak1", "peak2", "coaccess") %in% names(connections))) stop("Unexpected Cicero output")
mapping <- setNames(peaks$peak_id, site)
connections <- map_native_edges(connections, mapping)
if (any(!connections$peak1 %in% peaks$peak_id[scoring]) ||
    any(!connections$peak2 %in% peaks$peak_id[scoring])) stop("Cicero edge outside requested scoring peaks")
check_edge_geometry(connections, peaks, window)
connections <- connections[!is.na(connections$peak1) & !is.na(connections$peak2), c("peak1", "peak2", "coaccess")]
fwrite(connections, output, sep = "\t")
writeLines(capture.output(sessionInfo()), paste0(output, ".sessionInfo.txt"))
writeLines(c("method=Cicero_with_window_transpose_compatibility_shim",
             if (is.null(covariance_weights)) "covariance=native_equal_aggregate_weights" else "covariance=donor_equal_reliability_weighted",
             if (is.null(covariance_weights)) "method_variant=native_covariance" else "method_variant=donor_balanced_Cicero_adaptation",
             "aggregation=donor_confined_before_Cicero",
             "normalization=native_estimateSizeFactorsSimp_column_division",
             paste0("distance_samples=", distance_samples),
             paste0("window_model_workers=", workers),
             "calibration=full_training_peak_universe",
             "inference=requested_gene_cis_windows"), paste0(output, ".method.txt"))
