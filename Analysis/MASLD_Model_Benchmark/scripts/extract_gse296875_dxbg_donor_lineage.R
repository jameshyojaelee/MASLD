#!/usr/bin/env Rscript
# Extract the authors' decontaminated dXbg assay and aggregate it to the same
# donor-by-lineage units as the raw fixture.
#
# Direction of the question, stated so the result is read correctly: hepatocyte
# ambient RNA would MANUFACTURE a spurious non-hepatocyte signal, not destroy a
# real one.  What decontamination could plausibly change is the opposite case,
# where heavy ambient dilution depresses signal to noise enough to mask a
# genuine lineage-specific effect.  The raw-count ambient audit argues against
# that, since non-hepatocyte units carry roughly a tenth of the hepatocyte
# marker fraction.  This is therefore a bounded confirmatory check: if the
# negative survives decontamination it is substantially stronger, and if a
# signal appears the decontaminated analysis becomes the load-bearing one.
#
# The feature crosswalk is the real risk.  An unmapped dXbg feature is missing,
# never zero.  If the crosswalk cannot be established cleanly this script stops
# and reports rather than approximating it.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
names(args) <- sub("^--", "", args)
get_arg <- function(flag) {
  hit <- which(args == paste0("--", flag))
  if (length(hit) != 1L) stop(sprintf("missing required argument --%s", flag))
  args[[hit + 1L]]
}

rds_path        <- get_arg("rds")
rds_sha256      <- get_arg("rds-sha256")
membership_path <- get_arg("membership")
folds_path      <- get_arg("folds")
features_path   <- get_arg("raw-features")
output_dir      <- get_arg("output")

PRIMARY   <- c("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
SECONDARY <- c("endothelial_cell", "b_cell")
LINEAGES  <- c(PRIMARY, SECONDARY)
MIN_CELLS <- 20L
EXPECTED_NUCLEI <- 68398L
EXPECTED_DONORS <- 39L

message("verifying source identity (8.5 GB, a few minutes)")
observed_sha <- sub(" .*$", "",
                    system2("sha256sum", shQuote(rds_path), stdout = TRUE))
if (!identical(observed_sha, rds_sha256)) {
  stop(sprintf("processed RDS checksum differs: %s", observed_sha))
}

message("reading membership and folds")
membership <- read.delim(gzfile(membership_path), sep = "\t",
                         stringsAsFactors = FALSE, colClasses = "character")
folds <- read.delim(folds_path, sep = "\t", stringsAsFactors = FALSE,
                    colClasses = "character")
if (nrow(membership) != EXPECTED_NUCLEI) stop("frozen nucleus roster size differs")
if (nrow(folds) != EXPECTED_DONORS) stop("frozen donor roster size differs")
if (!setequal(unique(membership$lineage_id), LINEAGES)) stop("frozen lineage set differs")

raw_features <- read.delim(features_path, sep = "\t", stringsAsFactors = FALSE,
                           colClasses = "character")
if (nrow(raw_features) != 36601L) stop("raw RNA feature roster size differs")

message("loading the processed object (this is the memory peak)")
object <- readRDS(rds_path)
if (!("dXbg" %in% Assays(object))) stop("the dXbg assay is absent from the object")
if (ncol(object) != EXPECTED_NUCLEI) {
  stop(sprintf("object cell count differs: %d", ncol(object)))
}

counts <- SeuratObject::LayerData(object, assay = "dXbg", layer = "counts")
if (is.null(counts) || nrow(counts) == 0L) stop("the dXbg counts layer is empty")
cell_ids <- colnames(counts)

# ---- Cell axis check -------------------------------------------------------
# The well-namespaced barcode is the only unique key: 2,784 bare barcodes recur
# across wells.  A positional or bare-barcode join here would silently
# mis-assign 5,625 nuclei.
if (!setequal(cell_ids, membership$cell_id)) {
  overlap <- length(intersect(cell_ids, membership$cell_id))
  stop(sprintf(
    "dXbg cell identifiers do not match the frozen roster (overlap %d of %d)",
    overlap, EXPECTED_NUCLEI))
}

# ---- Feature crosswalk check ----------------------------------------------
# Unmapped is missing, never zero.  Ambiguity is refused rather than resolved
# by picking a winner.
dxbg_features <- rownames(counts)
symbol_counts <- table(raw_features$symbol)
ambiguous_symbols <- names(symbol_counts)[symbol_counts > 1L]

stable_hit <- match(dxbg_features, raw_features$ensembl_stable_id)
symbol_hit <- match(dxbg_features, raw_features$symbol)
is_ambiguous <- dxbg_features %in% ambiguous_symbols

map_index <- ifelse(
  !is.na(stable_hit), stable_hit,
  ifelse(!is_ambiguous, symbol_hit, NA_integer_)
)
map_state <- ifelse(
  !is.na(stable_hit), "exact_stable_id",
  ifelse(is_ambiguous, "ambiguous_symbol_structurally_missing",
         ifelse(!is.na(symbol_hit), "unique_symbol",
                "absent_from_raw_annotation_structurally_missing"))
)

mapped <- !is.na(map_index)
if (anyDuplicated(map_index[mapped]) != 0L) {
  stop("two dXbg features resolve to one raw feature; crosswalk is not clean")
}
mapped_fraction <- sum(mapped) / length(dxbg_features)
message(sprintf("crosswalk: %d of %d dXbg features mapped (%.4f)",
                sum(mapped), length(dxbg_features), mapped_fraction))
if (mapped_fraction < 0.90) {
  stop(sprintf(
    "dXbg feature crosswalk is not clean: only %.4f mapped. Stopping rather than approximating.",
    mapped_fraction))
}

# ---- Aggregation ----------------------------------------------------------
message("aggregating to donor by lineage units")
donors <- folds$donor_id
units <- expand.grid(lineage_id = LINEAGES, donor_id = donors,
                     stringsAsFactors = FALSE)[, c("donor_id", "lineage_id")]
unit_key <- paste(units$donor_id, units$lineage_id, sep = "|")

cell_unit <- paste(membership$donor_id, membership$lineage_id, sep = "|")
names(cell_unit) <- membership$cell_id
assignment <- factor(cell_unit[cell_ids], levels = unit_key)
if (anyNA(assignment)) stop("a nucleus could not be assigned to a unit")

column <- as.integer(assignment)
group <- sparseMatrix(
  i = seq_along(column), j = column, x = 1,
  dims = c(length(column), length(unit_key)),
  dimnames = list(cell_ids, unit_key)
)
if (ncol(group) != length(unit_key)) stop("a unit column was dropped")
aggregated <- as.matrix(counts %*% group)          # dXbg features by units
nuclei <- as.integer(table(assignment))

# Project onto the raw 36,601-feature axis.  Unmapped raw features are NA, a
# typed missing value, and are never filled with zero.
out <- matrix(NA_real_, nrow = length(unit_key), ncol = nrow(raw_features),
              dimnames = list(unit_key, raw_features$ensembl_stable_id))
out[, map_index[mapped]] <- t(aggregated[mapped, , drop = FALSE])

state <- ifelse(nuclei >= MIN_CELLS, "observed", "insufficient_cells")
libraries <- rowSums(out, na.rm = TRUE)

dir.create(output_dir, recursive = TRUE)
write.table(
  data.frame(
    unit_index = seq_len(nrow(units)) - 1L,
    donor_id = units$donor_id,
    lineage_id = units$lineage_id,
    analysis_role = ifelse(units$lineage_id %in% PRIMARY, "primary", "secondary"),
    nuclei = nuclei,
    minimum_nuclei_per_unit = MIN_CELLS,
    rna_state = state,
    library_dxbg = libraries,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(
    gene_index = seq_len(nrow(raw_features)) - 1L,
    ensembl_stable_id = raw_features$ensembl_stable_id,
    symbol = raw_features$symbol,
    dxbg_observed = seq_len(nrow(raw_features)) %in% map_index[mapped],
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(dxbg_feature = dxbg_features, mapping_state = map_state,
             stringsAsFactors = FALSE),
  file.path(output_dir, "feature_crosswalk.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE
)

connection <- gzfile(file.path(output_dir, "dxbg_counts.tsv.gz"), "wt")
write.table(out, connection, sep = "\t", quote = FALSE, col.names = NA)
close(connection)

audit <- list(
  schema_version = "masld-bench-gse296875-dxbg-donor-lineage-v1",
  dataset_id = "gse296875",
  assay = "dXbg",
  decontaminated = TRUE,
  source_rds_sha256 = rds_sha256,
  nuclei = EXPECTED_NUCLEI,
  donors = EXPECTED_DONORS,
  units = nrow(units),
  dxbg_features = length(dxbg_features),
  dxbg_features_mapped = sum(mapped),
  dxbg_mapped_fraction = mapped_fraction,
  raw_features_covered = sum(seq_len(nrow(raw_features)) %in% map_index[mapped]),
  unmapped_are_missing_never_zero = TRUE,
  counts_are_integer = isTRUE(all.equal(out[!is.na(out)], round(out[!is.na(out)]))),
  cell_axis_matched_frozen_roster = TRUE,
  primary_units_observed = sum(state == "observed" &
                                 units$lineage_id %in% PRIMARY),
  outcomes_read = FALSE,
  model_training_activated = FALSE,
  status = "pass_dxbg_donor_lineage"
)
writeLines(jsonlite::toJSON(audit, auto_unbox = TRUE, pretty = TRUE),
           file.path(output_dir, "audit.json"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat(jsonlite::toJSON(audit, auto_unbox = TRUE), "\n")
