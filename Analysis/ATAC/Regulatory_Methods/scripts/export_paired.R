#!/usr/bin/env Rscript
# Export genuinely paired native ATAC nuclei; raw RNA is joined separately by well/barcode.
suppressPackageStartupMessages({library(Matrix); library(jsonlite)})

parse_args <- function(args) {
  defaults <- list(mode = "pilot", lineage = "hepatocyte")
  if (length(args) %% 2L) stop("Arguments must be --name value pairs")
  for (i in seq.int(1L, length(args), by = 2L)) {
    if (!startsWith(args[[i]], "--")) stop("Invalid argument: ", args[[i]])
    defaults[[substring(args[[i]], 3L)]] <- args[[i + 1L]]
  }
  if (is.null(defaults$config) || is.null(defaults$output)) stop("Require --config and --output")
  if (!defaults$mode %in% c("pilot", "full")) stop("mode must be pilot or full")
  if (!defaults$lineage %in% c("hepatocyte", "stellate", "all")) stop("Unknown lineage")
  defaults
}

require_true <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
sha <- function(path) {
  result <- system2("sha256sum", c("--", shQuote(path)), stdout = TRUE)
  require_true(length(result) == 1L && grepl("^[0-9a-f]{64} ", result), "Cannot hash source")
  substr(result, 1L, 64L)
}
write_table <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write_sparse <- function(x, path) {
  # Matrix::writeMM does not support a compressed connection in every installed version.
  temporary <- paste0(path, ".uncompressed")
  Matrix::writeMM(x, temporary)
  input <- file(temporary, "rb"); output <- gzfile(path, "wb")
  on.exit({close(input); close(output)}, add = TRUE)
  repeat {
    block <- readBin(input, what = "raw", n = 1024L * 1024L)
    if (!length(block)) break
    writeBin(block, output)
  }
  # Keep the temporary only until both streams have been closed.
  close(input); close(output); on.exit(NULL)
  unlink(temporary)
}

select_cells <- function(cells, mode, seed, cap = 1000L) {
  set.seed(seed)
  eligible <- names(which(table(cells$donor_id) >= 100L))
  eligible <- sort(eligible)
  if (mode == "pilot") eligible <- sample(eligible, min(10L, length(eligible)))
  selected <- unlist(lapply(eligible, function(donor) {
    ii <- which(cells$donor_id == donor)
    sample(ii, min(length(ii), if (mode == "pilot") 100L else cap))
  }), use.names = FALSE)
  sort(selected)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  require_true(grepl("^[0-9]+$", Sys.getenv("SLURM_JOB_ID")) && nzchar(Sys.getenv("SLURM_JOB_NODELIST")),
    "Paired matrix export requires a SLURM compute allocation")
  cfg <- read_json(args$config, simplifyVector = TRUE)
  root <- cfg$project_root
  if (is.null(root)) root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
  settings <- cfg$paired
  default_identity <- file.path(root, "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/codex_gse296875_aligned_native_pseudobulk_21999585")
  identity <- if (is.null(settings$identity_root)) default_identity else settings$identity_root
  source <- if (is.null(settings$rds)) file.path(root,
    "Analysis/SingleCell/continual_integration/candidates/external-gse296875-download-v44-2026-08-18T000000Z/GSE296875_processed.rds") else settings$rds
  expected <- if (is.null(settings$rds_sha256)) "6af74b7dbdd3839841f90cb7f3e47a3184b4dd3dd02242f6f3bdbd14533cfcab" else settings$rds_sha256
  require_true(!file.exists(args$output), "Refusing to overwrite output directory")
  cells <- read.delim(file.path(identity, "cell_identity.tsv"), colClasses = "character", check.names = FALSE)
  geometry <- read.delim(file.path(identity, "native_peak_geometry.tsv"), colClasses = "character", check.names = FALSE)
  genes <- read.delim(file.path(identity, "genes.tsv"), colClasses = "character", check.names = FALSE)
  require_true(nrow(cells) == 68398L && !anyDuplicated(cells$cell_id) && length(unique(cells$donor_id)) == 39L,
    "Authoritative paired cell/donor identity differs")
  require_true(!anyDuplicated(paste(cells$well_id, cells$raw_barcode)), "Duplicated well/barcode identities")
  require_true(nrow(geometry) == 306706L && !anyDuplicated(geometry$peak_id), "Native peak identity differs")
  require_true(nrow(genes) == 36601L && !anyDuplicated(genes$ensembl_id), "Native RNA stable-ID identity differs")
  before <- file.info(source)[, c("size", "mtime")]
  require_true(before$size == 8534583098 && identical(sha(source), expected), "Paired author RDS source identity differs")
  object <- readRDS(source)
  # Read serialized slots directly; no Seurat/Signac accessor or installation is needed.
  assays <- attributes(object)[["assays"]]
  require_true("peaks" %in% names(assays), "Author object lacks native peaks")
  storage <- attributes(assays[["peaks"]])
  counts <- if ("counts" %in% names(storage)) storage$counts else storage$layers$counts
  require_true(inherits(counts, "dgCMatrix") && identical(dim(counts), c(306706L, 68398L)), "Native sparse ATAC dimensions differ")
  require_true(identical(colnames(counts), cells$cell_id) && identical(rownames(counts), geometry$peak_id),
    "Native ATAC axes disagree with authoritative pairing identities")
  require_true(all(is.finite(counts@x) & counts@x >= 0 & counts@x == floor(counts@x)), "ATAC entries are invalid counts")
  peaks <- data.frame(peak_id = geometry$peak_id, chrom = geometry$chromosome,
    start0 = as.numeric(geometry$bed_start_0based), end = as.numeric(geometry$bed_end_half_open))
  require_true(all(peaks$start0 >= 0 & peaks$end > peaks$start0), "Invalid native BED coordinates")
  native_parts <- regmatches(peaks$peak_id, regexec("^(chr[^-]+)-([0-9]+)-([0-9]+)$", peaks$peak_id))
  require_true(all(lengths(native_parts) == 4L), "Malformed native GRanges peak identifiers")
  require_true(all(vapply(native_parts, `[[`, "", 2L) == peaks$chrom) &&
    all(as.numeric(vapply(native_parts, `[[`, "", 3L)) - 1 == peaks$start0) &&
    all(as.numeric(vapply(native_parts, `[[`, "", 4L)) == peaks$end),
    "Native GRanges identifiers and BED coordinates differ")
  # Mesenchymal is preserved as source-dependent stellate context, not a pure stellate label.
  labels <- c(hepatocyte = "Hepatocytes", stellate = "Mesenchymal")
  lineages <- if (args$lineage == "all") names(labels) else args$lineage
  require_true(dir.create(args$output, recursive = TRUE), "Cannot create new paired export")
  seed <- if (is.null(cfg$seed)) 42L else as.integer(cfg$seed)
  summaries <- list()
  for (lineage in lineages) {
    pool <- which(cells$source_label == labels[[lineage]])
    local <- select_cells(cells[pool, , drop = FALSE], args$mode, seed)
    selected <- pool[local]
    require_true(length(selected) > 0L, paste("No eligible paired nuclei:", lineage))
    directory <- file.path(args$output, lineage)
    require_true(dir.create(directory), "Cannot create lineage output")
    metadata <- cells[selected, c("donor_id", "cell_id", "well_id", "raw_barcode", "source_label"), drop = FALSE]
    metadata$cohort <- "GSE296875"; metadata$lineage <- lineage
    metadata$condition <- "not_qualified"
    metadata$lineage_state <- if (lineage == "stellate") "source_dependent" else "source_label"
    metadata$source_lineage <- if (lineage == "stellate") "mesenchymal" else "hepatocyte"
    metadata <- metadata[, c("cohort", "donor_id", "lineage", "cell_id", "well_id", "raw_barcode", "source_label", "source_lineage", "lineage_state", "condition")]
    matrix <- Matrix::t(counts[, selected, drop = FALSE])
    require_true(identical(rownames(matrix), metadata$cell_id), "Exported paired cell axis differs")
    write_sparse(matrix, file.path(directory, "ATAC_counts.mtx.gz"))
    write_table(metadata, file.path(directory, "cells.tsv"))
    write_table(peaks, file.path(directory, "peaks.tsv"))
    write_table(data.frame(gene_id = genes$ensembl_id, gene_symbol = genes$native_gene_symbol), file.path(directory, "genes.tsv"))
    info <- list(cohort = "GSE296875", mode = args$mode, lineage = lineage,
      source_label = labels[[lineage]], lineage_state = unique(metadata$lineage_state),
      cells = nrow(matrix), donors = length(unique(metadata$donor_id)), peaks = ncol(matrix),
      seed = seed, minimum_cells_per_donor = 100L, maximum_cells_per_donor = if (args$mode == "pilot") 100L else 1000L,
      ATAC_count_unit = "native author peak-count events; unique-fragment equivalence not assumed",
      ATAC_total_counts = sum(matrix), same_nucleus_RNA_status = "awaiting_raw_H5_exact_well_barcode_export",
      disease_contrast_qualified = FALSE, external_untouched = FALSE,
      source_RDS = source, source_RDS_sha256 = expected, identity_root = identity,
      matrices_orientation = "cells_by_features", coordinate_system = "GRCh38 BED 0-based half-open",
      peak_id_coordinate_system = "GRanges_1based_closed",
      job_id = Sys.getenv("SLURM_JOB_ID"))
    write_json(info, file.path(directory, "paired_export.json"), pretty = TRUE, auto_unbox = TRUE)
    summaries[[lineage]] <- info
    message(lineage, ": ", nrow(matrix), " paired nuclei across ", info$donors, " donors")
  }
  require_true(identical(file.info(source)[, c("size", "mtime")], before), "Author RDS changed during export")
  write_json(summaries, file.path(args$output, "export_summary.json"), pretty = TRUE, auto_unbox = TRUE)
  capture.output(sessionInfo(), file = file.path(args$output, "R_sessionInfo.txt"))
}

if (sys.nframe() == 0L) main()
