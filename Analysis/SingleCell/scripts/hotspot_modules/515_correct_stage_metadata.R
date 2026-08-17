#!/usr/bin/env Rscript
# Build a new metadata candidate that restores source-diagnosis stage labels.
# Existing metadata and program releases are read-only.

suppressPackageStartupMessages(library(data.table))

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
INPUT <- file.path(
  BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
)
PAIRING <- file.path(BASE, "data/GSE202379/metadata/donor_pairing.csv")
OUTPUT <- Sys.getenv("HOTSPOT_V2_METADATA_FILE", unset = "")
if (!nzchar(OUTPUT)) stop("HOTSPOT_V2_METADATA_FILE is required", call. = FALSE)
if (file.exists(OUTPUT)) stop("Refusing to overwrite metadata candidate", call. = FALSE)

metadata <- fread(INPUT)
pairing <- fread(PAIRING, colClasses = "character")
stage_map <- c(
  "Healthy control" = "Healthy",
  "NAFLD" = "Steatosis",
  "NASH w/o cirrhosis" = "Steatohepatitis",
  "NASH with cirrhosis" = "Cirrhosis",
  "end stage" = "Cirrhosis"
)
if (!all(pairing$condition %in% names(stage_map))) {
  stop("Unrecognized GSE202379 diagnosis label", call. = FALSE)
}

run_map <- rbindlist(lapply(seq_len(nrow(pairing)), function(i) {
  runs <- trimws(strsplit(pairing$rna_srrs[i], ";", fixed = TRUE)[[1]])
  data.table(
    sample = runs,
    source_diagnosis = pairing$condition[i],
    corrected_stage = unname(stage_map[pairing$condition[i]])
  )
}))
if (anyDuplicated(run_map$sample)) {
  stop("A GSE202379 run maps to multiple diagnoses", call. = FALSE)
}

metadata[, `:=`(
  disease_stage_source = fifelse(
    dataset == "GSE202379", "GSE202379 published diagnosis", "existing harmonization"
  ),
  source_diagnosis = NA_character_
)]
metadata[run_map, on = "sample", `:=`(
  disease_stage_coarse = i.corrected_stage,
  source_diagnosis = i.source_diagnosis
)]
metadata[disease_stage_coarse == "Healthy", disease_stage_numeric := 0]
metadata[disease_stage_coarse == "Steatosis", disease_stage_numeric := 1]
metadata[disease_stage_coarse == "Steatohepatitis", disease_stage_numeric := 2]
metadata[disease_stage_coarse == "Cirrhosis", disease_stage_numeric := 3]

audit <- unique(metadata[dataset == "GSE202379", .(
  sample, source_diagnosis, disease_stage_coarse
)])
expected <- data.table(
  source_diagnosis = names(stage_map),
  disease_stage_coarse = unname(stage_map)
)
bad <- audit[expected, on = "source_diagnosis"][
  is.na(i.disease_stage_coarse) | disease_stage_coarse != i.disease_stage_coarse
]
if (nrow(bad) || anyNA(audit$source_diagnosis)) {
  stop("GSE202379 diagnosis-to-stage audit failed", call. = FALSE)
}

dir.create(dirname(OUTPUT), recursive = TRUE, showWarnings = FALSE)
fwrite(metadata, OUTPUT, sep = "\t", quote = FALSE, na = "NA")
source(file.path(BASE, "Analysis/SingleCell/scripts/lib_donor_collapse.R"))
donor_map <- build_srr_to_donor_map(BASE)
donor_metadata <- copy(metadata)
donor_metadata[, donor := fifelse(
  sample %in% names(donor_map), unname(donor_map[sample]), sample
)]
donor_conflicts <- donor_metadata[, .(
  n_datasets = uniqueN(dataset[!is.na(dataset) & dataset != ""]),
  n_stages = uniqueN(disease_stage_coarse[
    !is.na(disease_stage_coarse) & disease_stage_coarse != ""
  ])
), by = donor][n_datasets > 1L | n_stages > 1L]
if (nrow(donor_conflicts)) {
  stop("Corrected metadata disagrees within biological donor", call. = FALSE)
}
donor_metadata <- donor_metadata[, .(
  sample = donor,
  dataset = dataset[which(!is.na(dataset) & dataset != "")[1]],
  disease_stage_coarse = disease_stage_coarse[
    which(!is.na(disease_stage_coarse) & disease_stage_coarse != "")[1]
  ],
  exclude_stage_analysis = any(exclude_stage_analysis %in% TRUE),
  n_source_records = .N
), by = donor][, donor := NULL]
fwrite(
  donor_metadata,
  sub("\\.tsv$", "_donor.tsv", OUTPUT),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  audit[, .N, by = .(source_diagnosis, disease_stage_coarse)][
    order(disease_stage_coarse, source_diagnosis)
  ],
  sub("\\.tsv$", "_audit.tsv", OUTPUT), sep = "\t", quote = FALSE
)
cat(sprintf("[515] wrote %s (%d rows)\n", OUTPUT, nrow(metadata)))
