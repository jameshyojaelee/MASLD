#!/usr/bin/env Rscript

# Dump the exact GPL16686 transcript cluster -> Entrez and Entrez -> Ensembl
# relations from the frozen Bioconductor annotation packages.
#
# The runtime capability probe reported 0 of 200 queried transcript clusters
# carrying an Entrez ID, which would make the whole GPL16686 gene lane
# unmappable. That probe queried only the first 200 keys, and on a HuGene array
# the low-numbered transcript clusters are largely controls, so the result is
# most likely a sampling artefact rather than a property of the package. This
# script queries ALL keys and records the answer either way, because a genuine
# zero would block the lane and must not be discovered later.
#
# The platform axis question is settled here too. The frozen requirements record
# platform_rows = 53981 from the GEO GPL16686 table, but the core summarizer
# returns 53617 transcript clusters, and the matrix must be keyed on what the
# summarizer produces. Both counts are recorded and the annotation axis is the
# one carried forward.
#
# Mapping authority per the user's explicit decision:
#   hugene20sttranscriptcluster.db -> org.Hs.eg.db -> GENCODE v49
# The Thermo Fisher HuGene-2_0-st-v1 Release 36 asset stays recorded as
# unavailable provenance. Coordinates and GB_ACC are never used to infer a
# mapping.
#
# Reads no CEL, no label, no GEO series matrix; fits no model.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: dump_gpl16686_annotation_maps.R OUTPUT_DIR")
output_dir <- args[[1L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260825L)
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(hugene20sttranscriptcluster.db)
  library(org.Hs.eg.db)
  library(jsonlite)
})

EXPECTED <- c(
  "hugene20sttranscriptcluster.db" = "8.8.0",
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

db <- hugene20sttranscriptcluster.db
available_columns <- AnnotationDbi::columns(db)
available_keytypes <- AnnotationDbi::keytypes(db)
if (!("ENTREZID" %in% available_columns) || !("PROBEID" %in% available_keytypes)) {
  stop("hugene20sttranscriptcluster.db does not expose PROBEID -> ENTREZID")
}

## ------------------------------------------- every transcript cluster, not 200
platform_features <- AnnotationDbi::keys(db, keytype = "PROBEID")
platform_features <- sort(unique(platform_features))
if (length(platform_features) < 1000L) stop("GPL16686 annotation axis is degenerate")

probe_entrez <- AnnotationDbi::select(
  db, keys = platform_features, columns = "ENTREZID", keytype = "PROBEID"
)
probe_entrez <- probe_entrez[!is.na(probe_entrez$ENTREZID), c("PROBEID", "ENTREZID"), drop = FALSE]
probe_entrez <- unique(probe_entrez)
probe_entrez <- probe_entrez[order(probe_entrez$PROBEID, probe_entrez$ENTREZID), , drop = FALSE]

# Reproduce the probe's sampling so the discrepancy is explained, not just fixed.
first_200 <- head(platform_features, 200L)
first_200_with_entrez <- length(intersect(first_200, unique(probe_entrez$PROBEID)))

## --------------------------------------- Entrez -> Ensembl, ambiguity dumped
entrez_ids <- sort(unique(probe_entrez$ENTREZID))
org_entrez <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID")
entrez_present <- entrez_ids[entrez_ids %in% org_entrez]
entrez_absent <- entrez_ids[!(entrez_ids %in% org_entrez)]

entrez_ensembl <- if (length(entrez_present) > 0L) {
  frame <- AnnotationDbi::select(
    org.Hs.eg.db, keys = entrez_present, columns = "ENSEMBL", keytype = "ENTREZID"
  )
  frame <- frame[!is.na(frame$ENSEMBL), c("ENTREZID", "ENSEMBL"), drop = FALSE]
  unique(frame[order(frame$ENTREZID, frame$ENSEMBL), , drop = FALSE])
} else {
  data.frame(ENTREZID = character(0), ENSEMBL = character(0), stringsAsFactors = FALSE)
}

## ------------------------------------------------------------------- outputs
write_tsv <- function(frame, path) {
  utils::write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
write_tsv(
  data.frame(platform_feature_id = platform_features, stringsAsFactors = FALSE),
  file.path(output_dir, "gpl16686_platform_axis.tsv")
)
write_tsv(
  data.frame(
    platform_feature_id = probe_entrez$PROBEID, entrez_id = probe_entrez$ENTREZID,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "gpl16686_probe_entrez.tsv")
)
write_tsv(
  data.frame(
    entrez_id = entrez_ensembl$ENTREZID, ensembl_gene_id = entrez_ensembl$ENSEMBL,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "entrez_ensembl.tsv")
)
write_tsv(
  data.frame(
    platform_feature_id = character(0), state = character(0), stringsAsFactors = FALSE
  ),
  file.path(output_dir, "gpl16686_features_absent_from_annotation.tsv")
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
  data.frame(package = names(observed), version = unname(observed), stringsAsFactors = FALSE),
  file.path(output_dir, "annotation_package_versions.tsv")
)
writeLines(
  capture.output(print(utils::sessionInfo())),
  file.path(output_dir, "R_sessionInfo.txt")
)

receipt <- list(
  schema_version = "masld-bench-gpl16686-annotation-dump-v1",
  status = "pass_annotation_relations_dumped",
  platform_id = "GPL16686",
  series = "GSE83452",
  cohort_family_id = "antwerp_inserm_shared",
  annotation_package = "hugene20sttranscriptcluster.db",
  annotation_package_version = unname(observed[["hugene20sttranscriptcluster.db"]]),
  org_hs_eg_db_version = unname(observed[["org.Hs.eg.db"]]),
  annotation_dbi_version = unname(observed[["AnnotationDbi"]]),
  available_columns = as.list(available_columns),
  available_keytypes = as.list(available_keytypes),
  platform_features = length(platform_features),
  geo_platform_table_rows = 53981L,
  core_summarizer_features = 53617L,
  axis_carried_forward = "hugene20sttranscriptcluster.db_PROBEID_keys",
  probe_entrez_pairs = nrow(probe_entrez),
  features_with_entrez = length(unique(probe_entrez$PROBEID)),
  first_200_keys_with_entrez = first_200_with_entrez,
  runtime_probe_reported_first_200_with_entrez = 0L,
  first_200_sampling_explains_runtime_probe_zero = (first_200_with_entrez == 0L),
  distinct_entrez_ids = length(entrez_ids),
  entrez_present_in_org_hs_eg_db = length(entrez_present),
  entrez_absent_from_org_hs_eg_db = length(entrez_absent),
  entrez_ensembl_pairs = nrow(entrez_ensembl),
  thermo_release_36_asset = "unavailable_permission_login_gated",
  coordinate_or_GB_ACC_mapping_inference_used = FALSE,
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
cat(sprintf(
  "features=%d with_entrez=%d first200_with_entrez=%d entrez=%d ensembl_pairs=%d\n",
  length(platform_features), length(unique(probe_entrez$PROBEID)),
  first_200_with_entrez, length(entrez_ids), nrow(entrez_ensembl)
))
