#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))

if (!file.exists(contract_file)) fail("Contract must be frozen before partitioning: ", contract_file)
for (path in required_input_files) if (!file.exists(path)) fail("Missing input: ", path)

if (dir.exists(sealed_input_root)) {
  required <- c(nonholdout_dge, nonholdout_meta, holdout_dge, holdout_meta, partition_manifest)
  required <- c(required, nonholdout_composition, holdout_composition)
  if (all(file.exists(required))) {
    message("Sealed partition already exists; refusing to rewrite")
    quit(save = "no", status = 0)
  }
  fail("Incomplete sealed input directory exists: ", sealed_input_root)
}

tmp <- atomic_dir(sealed_input_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

# This step is intentionally blind to expression values. It uses only dge column
# metadata to split the monolithic corrected object and logs no expression summary.
dge <- readRDS(input_dge)
meta <- as.data.table(readRDS(input_meta))
assert_true(ncol(dge) == 1257L, "Corrected DGE must contain 1,257 QC-passing samples")
assert_true(nrow(meta) == 1281L, "Matched metadata must contain 1,281 rows")
assert_true(all(colnames(dge) %in% meta$sample_id), "Every DGE sample must map to metadata")

dge_dataset <- as.character(dge$samples$dataset)
if (length(dge_dataset) != ncol(dge) || anyNA(dge_dataset)) {
  dge_dataset <- meta$dataset[match(colnames(dge), meta$sample_id)]
}
assert_true(!anyNA(dge_dataset), "Every DGE sample must have a dataset")

holdout_idx <- which(dge_dataset == holdout_cohort)
nonholdout_idx <- which(dge_dataset != holdout_cohort)
assert_true(length(holdout_idx) == 160L, "GSE193066 must contain 160 QC-passing samples")
assert_true(length(nonholdout_idx) == 1097L, "Non-holdout partition must contain 1,097 samples")
assert_true(length(intersect(colnames(dge)[holdout_idx], colnames(dge)[nonholdout_idx])) == 0L,
            "Holdout and non-holdout samples overlap")

dge_nonholdout <- dge[, nonholdout_idx, keep.lib.sizes = TRUE]
dge_holdout <- dge[, holdout_idx, keep.lib.sizes = TRUE]
assert_true(identical(dge_nonholdout$samples$norm.factors,
                      dge$samples$norm.factors[nonholdout_idx]),
            "Frozen non-holdout TMM factors drifted during partitioning")
assert_true(identical(dge_holdout$samples$norm.factors,
                      dge$samples$norm.factors[holdout_idx]),
            "Frozen holdout TMM factors drifted during partitioning")
meta_nonholdout <- meta[sample_id %in% colnames(dge_nonholdout)]
meta_holdout <- meta[sample_id %in% colnames(dge_holdout)]
meta_nonholdout <- meta_nonholdout[match(colnames(dge_nonholdout), sample_id)]
meta_holdout <- meta_holdout[match(colnames(dge_holdout), sample_id)]
assert_true(identical(meta_nonholdout$sample_id, colnames(dge_nonholdout)), "Non-holdout metadata order drift")
assert_true(identical(meta_holdout$sample_id, colnames(dge_holdout)), "Holdout metadata order drift")

# Rebuild the deposited MuSiC table against the corrected fragment-count census.
# Phenotype columns from the historical merged table are deliberately discarded.
composition <- fread(composition_candidate, check.names = FALSE)
composition_contract <- fread(composition_testability)
celltypes <- composition_contract$celltype
assert_true(length(celltypes) == 22L && uniqueN(celltypes) == 22L,
            "Composition contract must enumerate 22 deposited lineages")
assert_true(all(c("sample_id", "dataset", celltypes) %in% names(composition)),
            "Deposited composition table lacks contracted columns")
assert_true(!anyDuplicated(composition$sample_id), "Deposited composition sample IDs are duplicated")
composition <- composition[, c("sample_id", "dataset", celltypes), with = FALSE]
composition_all <- merge(
  data.table(sample_id = colnames(dge), dataset = dge_dataset), composition,
  by = c("sample_id", "dataset"), all.x = TRUE, sort = FALSE
)
composition_all <- composition_all[match(colnames(dge), sample_id)]
assert_true(identical(composition_all$sample_id, colnames(dge)),
            "Corrected composition census order drift")
available <- composition_contract[testability == "testable", celltype]
row_available <- rowSums(!is.na(composition_all[, ..available]))
assert_true(all(row_available %in% c(0L, length(available))),
            "Composition assay missingness is not row-consistent")
assert_true(sum(row_available == length(available)) == composition_expected_available,
            "Corrected composition availability census drift")
assert_true(sum(row_available[holdout_idx] == length(available)) == 160L,
            "Every locked holdout sample must have deposited composition")
composition_nonholdout <- composition_all[nonholdout_idx]
composition_holdout <- composition_all[holdout_idx]

saveRDS(dge_nonholdout, file.path(tmp, basename(nonholdout_dge)), compress = "xz")
saveRDS(meta_nonholdout, file.path(tmp, basename(nonholdout_meta)), compress = "xz")
saveRDS(dge_holdout, file.path(tmp, basename(holdout_dge)), compress = "xz")
saveRDS(meta_holdout, file.path(tmp, basename(holdout_meta)), compress = "xz")
fwrite(composition_nonholdout,
       file.path(tmp, basename(nonholdout_composition)), sep = "\t", quote = FALSE, na = "NA")
fwrite(composition_holdout,
       file.path(tmp, basename(holdout_composition)), sep = "\t", quote = FALSE, na = "NA")

manifest <- rbindlist(list(
  data.table(
    partition = "nonholdout",
    dataset = names(table(dge_dataset[nonholdout_idx])),
    n_samples = as.integer(table(dge_dataset[nonholdout_idx]))
  ),
  data.table(
    partition = "holdout",
    dataset = names(table(dge_dataset[holdout_idx])),
    n_samples = as.integer(table(dge_dataset[holdout_idx]))
  ),
  composition_all[, .(
    partition = "composition_assay",
    n_samples = sum(rowSums(!is.na(.SD)) > 0L)
  ), by = dataset, .SDcols = available]
), use.names = TRUE)
write_tsv(manifest, file.path(tmp, basename(partition_manifest)))

artifacts <- c("nonholdout_dge.rds", "nonholdout_meta.rds", "holdout_dge.rds",
               "holdout_meta.rds", "nonholdout_composition.tsv.gz",
               "holdout_composition.tsv.gz", "partition_manifest.tsv")
write_tsv(data.table(
  artifact = artifacts,
  sha256 = vapply(file.path(tmp, artifacts), sha256_file, character(1))
), file.path(tmp, "artifact_manifest.tsv"))

publish_dir(tmp, sealed_input_root)
Sys.chmod(list.files(sealed_input_root, full.names = TRUE), mode = "0440")
message("Published blind partition: ", sealed_input_root)
