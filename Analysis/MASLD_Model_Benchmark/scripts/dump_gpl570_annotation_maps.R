#!/usr/bin/env Rscript

# Dump the exact GPL570 probe -> Entrez and Entrez -> Ensembl relations from the
# frozen Bioconductor annotation packages.
#
# This script is the annotation authority half of the GPL570 GENCODE v49
# crosswalk. It reads only annotation databases and the frozen GEO platform
# axis. It never opens a CEL file, never reads a label, never reads a GEO series
# matrix, and never fits a model. Ambiguity is dumped, not resolved: every
# probe -> Entrez and Entrez -> Ensembl pair is emitted in long form so that the
# one-to-one decision is made once, in the auditable Python state machine.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: dump_gpl570_annotation_maps.R PLATFORM_AXIS_TSV OUTPUT_DIR")
}
platform_axis_path <- args[[1L]]
output_dir <- args[[2L]]
if (!file.exists(platform_axis_path)) stop("frozen GPL570 platform axis is absent")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260825L)
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
})

EXPECTED <- c(
  "hgu133plus2.db" = "3.13.0",
  "org.Hs.eg.db" = "3.22.0",
  "AnnotationDbi" = "1.72.0"
)
observed <- vapply(
  names(EXPECTED), function(package) as.character(packageVersion(package)), character(1)
)
if (!identical(unname(observed), unname(EXPECTED))) {
  stop(sprintf(
    "annotation package versions differ: %s",
    paste(sprintf("%s=%s", names(observed), observed), collapse = ",")
  ))
}

## --------------------------------------------------- frozen GEO platform axis
axis <- utils::read.delim(
  platform_axis_path, sep = "\t", quote = "", check.names = FALSE,
  colClasses = "character", na.strings = character(0)
)
if (!("ID" %in% names(axis))) stop("frozen platform axis lacks an ID column")
platform_features <- axis[["ID"]]
if (length(platform_features) != 54675L || anyDuplicated(platform_features) > 0L) {
  stop("frozen GPL570 platform axis is not 54675 unique features")
}

## ----------------------------------------- probe -> Entrez, ambiguity dumped
annotation_probes <- AnnotationDbi::keys(hgu133plus2.db, keytype = "PROBEID")
present <- platform_features[platform_features %in% annotation_probes]
absent <- platform_features[!(platform_features %in% annotation_probes)]

probe_entrez <- AnnotationDbi::select(
  hgu133plus2.db, keys = present, columns = "ENTREZID", keytype = "PROBEID"
)
probe_entrez <- probe_entrez[!is.na(probe_entrez$ENTREZID), c("PROBEID", "ENTREZID"), drop = FALSE]
probe_entrez <- unique(probe_entrez)
probe_entrez <- probe_entrez[order(probe_entrez$PROBEID, probe_entrez$ENTREZID), , drop = FALSE]

## --------------------------------------- Entrez -> Ensembl, ambiguity dumped
entrez_ids <- sort(unique(probe_entrez$ENTREZID))
org_entrez <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID")
entrez_present <- entrez_ids[entrez_ids %in% org_entrez]
entrez_absent <- entrez_ids[!(entrez_ids %in% org_entrez)]

entrez_ensembl <- AnnotationDbi::select(
  org.Hs.eg.db, keys = entrez_present, columns = "ENSEMBL", keytype = "ENTREZID"
)
entrez_ensembl <- entrez_ensembl[
  !is.na(entrez_ensembl$ENSEMBL), c("ENTREZID", "ENSEMBL"), drop = FALSE
]
entrez_ensembl <- unique(entrez_ensembl)
entrez_ensembl <- entrez_ensembl[
  order(entrez_ensembl$ENTREZID, entrez_ensembl$ENSEMBL), , drop = FALSE
]

## ------------------------------------------------------------------- outputs
write_tsv <- function(frame, path) {
  utils::write.table(
    frame, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
}
write_tsv(
  data.frame(
    platform_feature_id = probe_entrez$PROBEID,
    entrez_id = probe_entrez$ENTREZID,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "gpl570_probe_entrez.tsv")
)
write_tsv(
  data.frame(
    entrez_id = entrez_ensembl$ENTREZID,
    ensembl_gene_id = entrez_ensembl$ENSEMBL,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "entrez_ensembl.tsv")
)
write_tsv(
  data.frame(
    platform_feature_id = absent,
    state = rep("platform_feature_absent_from_annotation_package", length(absent)),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "gpl570_features_absent_from_annotation.tsv")
)
write_tsv(
  data.frame(
    entrez_id = entrez_absent,
    state = rep("entrez_absent_from_org_Hs_eg_db", length(entrez_absent)),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "entrez_absent_from_org_hs_eg_db.tsv")
)
write_tsv(
  data.frame(
    package = names(observed), version = unname(observed), stringsAsFactors = FALSE
  ),
  file.path(output_dir, "annotation_package_versions.tsv")
)
writeLines(
  capture.output(print(utils::sessionInfo())),
  file.path(output_dir, "R_sessionInfo.txt")
)

receipt <- list(
  schema_version = "masld-bench-gpl570-annotation-dump-v1",
  status = "pass_annotation_relations_dumped",
  platform_id = "GPL570",
  platform_features = length(platform_features),
  annotation_package = "hgu133plus2.db",
  annotation_package_version = unname(observed[["hgu133plus2.db"]]),
  org_hs_eg_db_version = unname(observed[["org.Hs.eg.db"]]),
  annotation_dbi_version = unname(observed[["AnnotationDbi"]]),
  annotation_package_probeid_keys = length(annotation_probes),
  platform_features_present_in_annotation = length(present),
  platform_features_absent_from_annotation = length(absent),
  probe_entrez_pairs = nrow(probe_entrez),
  distinct_entrez_ids = length(entrez_ids),
  entrez_present_in_org_hs_eg_db = length(entrez_present),
  entrez_absent_from_org_hs_eg_db = length(entrez_absent),
  entrez_ensembl_pairs = nrow(entrez_ensembl),
  ambiguity_resolved_in_this_script = FALSE,
  CEL_expression_values_read = FALSE,
  GEO_series_matrix_read = FALSE,
  labels_read = FALSE,
  model_training_activated = FALSE,
  sealed_outcomes_read = FALSE
)
writeLines(
  jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE),
  file.path(output_dir, "annotation_dump_receipt.json")
)
cat(jsonlite::toJSON(receipt, auto_unbox = TRUE), "\n")
