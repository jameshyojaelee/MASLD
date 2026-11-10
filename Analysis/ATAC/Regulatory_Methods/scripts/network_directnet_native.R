#!/usr/bin/env Rscript
# Bounded native-method diagnostic. No donor-held-out performance claim.
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(Signac)
  library(DIRECTNET)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(key, default = NULL) {
  at <- match(key, args)
  if (is.na(at)) return(default)
  if (at == length(args)) stop("Missing argument: ", key)
  args[[at + 1L]]
}
normalize_native_links <- function(native) {
  required <- c("gene", "Peak1", "Peak2", "Importance")
  if (is.list(native) && all(required %in% names(native))) {
    lengths <- lengths(native)
    if (length(unique(lengths)) == 1L && all(lengths > 0L)) {
      return(as.data.frame(native, stringsAsFactors = FALSE))
    }
  }
  data.frame(gene = character(), Peak1 = character(), Peak2 = character(), Importance = numeric())
}
if ("--test-native-output" %in% args) {
  options(Seurat.object.assay.version = "v3")
  fixture <- CreateSeuratObject(matrix(1, 2L, 3L,
    dimnames = list(c("a", "b"), c("c1", "c2", "c3"))))
  expected <- data.frame(gene = c("a", "b"), Peak1 = c("p", "p"),
                         Peak2 = c("q", "r"), Importance = c(.2, .3))
  Misc(fixture, slot = "direct.net") <- expected
  stopifnot(identical(normalize_native_links(Misc(fixture, slot = "direct.net")), expected))
  stopifnot(nrow(normalize_native_links(list(Starts = numeric(), Ends = numeric()))) == 0L)
  cat("PASS: native Misc column-list output retains all link rows and empty output stays untestable\n")
  quit(status = 0L)
}
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Compute allocation required")
output <- arg("--output")
if (is.null(output) || dir.exists(output)) stop("Use a new output directory")
dir.create(output, recursive = TRUE)
options(error = function() {
  write_json(list(status = "failed", scope = "native_technical_comparison_only",
    reason = geterrmessage(), scientific_adoption = FALSE),
    file.path(output, "completion.json"), pretty = TRUE, auto_unbox = TRUE)
  quit(status = 1L)
})
seed <- as.integer(as.double(arg("--seed", "42")) %% .Machine$integer.max)
set.seed(seed)
options(Seurat.object.assay.version = "v3")
cells <- fread(arg("--cells"))
peaks <- fread(arg("--peaks"))
genes <- fread(arg("--genes"))
counts <- as(readMM(arg("--counts")), "dgCMatrix") # peaks by nuclei
if (nrow(counts) != nrow(peaks) || ncol(counts) != nrow(cells)) stop("Input axes differ")
if (ncol(counts) > 1000L || ncol(counts) < 100L) stop("Native diagnostic requires 100 to 1000 nuclei")
site <- paste(peaks$chrom, peaks$start0, peaks$end, sep = "-")
if (anyDuplicated(site)) stop("Duplicate native peak coordinates")
fwrite(data.frame(peak_id = peaks$peak_id,
  canonical_peak_id = paste0(peaks$chrom, ":", peaks$start0, "-", peaks$end),
  native_assay_feature_id = site, native_model_feature_id = gsub("-", "_", site)),
  file.path(output, "native_peak_id_mapping.tsv"), sep = "\t")
rownames(counts) <- site
colnames(counts) <- paste0("nucleus_", seq_len(ncol(counts)))
cells$native_cell_id <- colnames(counts)
fwrite(cells, file.path(output, "cells.tsv"), sep = "\t")
# The native source anchors every strand at gene Starts and uses the preceding
# 500 bp. Retain that behavior explicitly for this diagnostic.
genome <- data.frame(genes = genes$gene_id, Chrom = genes$chrom,
                     Starts = genes$start0, Ends = genes$end)
geometry <- data.frame(gene_id = genes$gene_id, strand = genes$strand,
  source_start0 = genes$start0, tss0 = genes$tss0,
  native_promoter_start = genes$start0 - 500, native_promoter_end = genes$start0,
  native_cis_start = genes$start0 - 250000, native_cis_end = genes$start0 + 250000)
fwrite(geometry, file.path(output, "native_geometry.tsv"), sep = "\t")
object <- CreateSeuratObject(counts = counts, assay = "ATAC", min.cells = 0, min.features = 0)
Idents(object) <- factor(rep("technical_lineage", ncol(object)))
object <- RunTFIDF(object)
object <- FindTopFeatures(object, min.cutoff = "q0")
ncomponents <- min(15L, ncol(object) - 2L, nrow(object) - 2L)
object <- RunSVD(object, n = ncomponents, reduction.name = "lsi")
object <- RunUMAP(object, reduction = "lsi", dims = 2:ncomponents,
  seed.use = seed, n.neighbors = 30L, n.epochs = 200L, verbose = FALSE)
fwrite(data.frame(native_cell_id = rownames(Embeddings(object, "umap")),
  Embeddings(object, "umap")), file.path(output, "native_umap.tsv"), sep = "\t")
started <- proc.time()[["elapsed"]]
api_compatibility <- "--xgboost-api-compatibility" %in% args
native_runner <- DIRECTNET::Run_DIRECT_NET
if (api_compatibility) {
  compatibility <- new.env(parent = asNamespace("DIRECTNET"))
  compatibility$model_audits <- list()
  compatibility$tryCatch <- function(expr, error, ...) {
    base::tryCatch(expr, error = function(condition) {
      cat(conditionMessage(condition), "\n", file = file.path(output, "native_swallowed_errors.txt"), append = TRUE)
      error(condition)
    }, ...)
  }
  compatibility$xgboost <- function(params, data, label, nrounds, nthread,
                                    objective, verbose = 0, ...) {
    parameters <- utils::modifyList(params, list(objective = objective, nthread = nthread))
    model <- xgboost::xgb.train(params = parameters,
      data = xgboost::xgb.DMatrix(data = data, label = label),
      nrounds = nrounds, verbose = verbose)
    importance <- xgboost::xgb.importance(model = model)
    compatibility$model_audits[[length(compatibility$model_audits) + 1L]] <- list(
      n_aggregates = nrow(data), n_predictors = ncol(data),
      target_min = min(label), target_max = max(label), target_variance = stats::var(as.numeric(label)),
      n_importance_rows = nrow(importance))
    write_json(compatibility$model_audits, file.path(output, "native_XGBoost_model_audits.json"),
               pretty = TRUE, auto_unbox = TRUE)
    model
  }
  environment(native_runner) <- compatibility
  writeLines(deparse(compatibility$xgboost), file.path(output, "local_xgboost_API_shim.R"))
}
result <- native_runner(object, peakcalling = FALSE, k_neigh = 50,
  atacbinary = TRUE, max_overlap = 0.8, size_factor_normalize = FALSE,
  genome.info = genome, focus_markers = genes$gene_id,
  nthread = 1, early_stop = FALSE, rescued = FALSE, seed = seed)
native <- normalize_native_links(Misc(result, slot = "direct.net"))
scored <- nrow(native) > 0L
fwrite(native, file.path(output, "native_links.tsv"), sep = "\t")
coverage <- data.frame(gene_id = genes$gene_id,
  n_native_links = as.integer(table(factor(native$gene, levels = genes$gene_id))))
coverage$evaluation_status <- ifelse(coverage$n_native_links > 0L,
  "native_accessibility_importance_only", "native_gene_not_scored")
fwrite(coverage, file.path(output, "native_coverage.tsv"), sep = "\t")
aggregate <- Misc(result, slot = "aggregated.data")
membership <- as.matrix(aggregate$cell_sample)
members <- lapply(seq_len(nrow(membership)), function(i) {
  indices <- as.integer(membership[i, ])
  indices <- indices[is.finite(indices) & indices > 0 & indices <= nrow(cells)]
  data.frame(native_aggregate = i, source_index = indices,
             cell_id = cells$cell_id[indices], donor_id = cells$donor_id[indices])
})
members <- do.call(rbind, members)
fwrite(members, file.path(output, "native_aggregate_membership.tsv"), sep = "\t")
mixing <- stats::aggregate(members$donor_id, list(native_aggregate = members$native_aggregate),
                   function(x) length(unique(x)))
names(mixing)[2] <- "n_donors"
fwrite(mixing, file.path(output, "native_donor_mixing.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
writeLines(deparse(DIRECTNET::Run_DIRECT_NET), file.path(output, "installed_Run_DIRECT_NET.R"))
write_json(list(status = if (scored) "completed" else "untestable",
  reason = if (scored) "" else "native_source_produced_no_usable_link_rows",
  scope = "native_technical_comparison_only", xgboost_API_compatibility_shim = api_compatibility,
  scientific_adoption = FALSE, donor_held_out = FALSE, independent_donor_performance = "not_evaluated",
  aggregation = "native_overlapping_UMAP_neighborhoods_may_mix_donors",
  native_promoter = "unstranded_gene_start_minus500_to_gene_start",
  harmonized_promoter = FALSE, early_stop = FALSE,
  optional_internal_random_split = "not_invoked_early_stop_FALSE",
  native_atacbinary = TRUE, size_factor_normalize = FALSE, seed = seed,
  n_nuclei = ncol(counts), n_peaks = nrow(counts), n_requested_genes = nrow(genes),
  n_links = nrow(native), elapsed_model_seconds = proc.time()[["elapsed"]] - started),
  file.path(output, "completion.json"), pretty = TRUE, auto_unbox = TRUE)
