#!/usr/bin/env Rscript

# Blind partition for BULK-PROGRAM-MAP-v9. Forked from the v8 script with one
# change: the composition census is built from the sealed acceptance decision
# (BayesPrism) rather than the deposited MuSiC table, which is degenerate on the
# real bulk substrate. This step reads column metadata only and never inspects an
# expression value.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))

if (!file.exists(composition_accepted_seal)) {
  fail("Composition acceptance must be sealed before partitioning: ", composition_accepted_seal)
}
if (!file.exists(accepted_composition)) {
  fail("Accepted composition table is missing; no estimator qualified: ", accepted_composition)
}
for (path in c(input_dge, input_meta, gse193066_crosswalk)) {
  if (!file.exists(path)) fail("Missing input: ", path)
}

if (dir.exists(sealed_input_root)) {
  required <- c(nonholdout_dge, nonholdout_meta, holdout_dge, holdout_meta,
                partition_manifest, nonholdout_composition, holdout_composition)
  if (all(file.exists(required))) {
    message("Sealed partition already exists; refusing to rewrite")
    quit(save = "no", status = 0)
  }
  fail("Incomplete sealed input directory exists: ", sealed_input_root)
}

tmp <- atomic_dir(sealed_input_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)

dge <- readRDS(input_dge)
meta <- as.data.table(readRDS(input_meta))
assert_true(ncol(dge) == 1257L, "Corrected DGE must contain 1,257 QC-passing samples")
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

dge_nonholdout <- dge[, nonholdout_idx, keep.lib.sizes = TRUE]
dge_holdout <- dge[, holdout_idx, keep.lib.sizes = TRUE]
assert_true(identical(dge_nonholdout$samples$norm.factors,
                      dge$samples$norm.factors[nonholdout_idx]),
            "Frozen non-holdout TMM factors drifted during partitioning")
meta_nonholdout <- meta[sample_id %in% colnames(dge_nonholdout)][
  match(colnames(dge_nonholdout), sample_id)]
meta_holdout <- meta[sample_id %in% colnames(dge_holdout)][
  match(colnames(dge_holdout), sample_id)]
assert_true(identical(meta_nonholdout$sample_id, colnames(dge_nonholdout)),
            "Non-holdout metadata order drift")
assert_true(identical(meta_holdout$sample_id, colnames(dge_holdout)),
            "Holdout metadata order drift")

composition <- fread(accepted_composition, check.names = FALSE)
lineages <- setdiff(names(composition)[vapply(composition, is.numeric, logical(1))],
                    c("sample_id", "dataset"))
assert_true(length(lineages) == 16L,
            paste0("Accepted composition must carry 16 lineages; found ", length(lineages)))
composition <- composition[match(colnames(dge), sample_id)]
assert_true(identical(composition$sample_id, colnames(dge)),
            "Accepted composition census order drift")
row_available <- rowSums(!is.na(composition[, ..lineages]))
assert_true(all(row_available %in% c(0L, length(lineages))),
            "Composition missingness is not row-consistent")
assert_true(sum(row_available == length(lineages)) == 1257L,
            "The accepted estimator must cover every corrected census sample")

saveRDS(dge_nonholdout, file.path(tmp, basename(nonholdout_dge)), compress = "xz")
saveRDS(meta_nonholdout, file.path(tmp, basename(nonholdout_meta)), compress = "xz")
saveRDS(dge_holdout, file.path(tmp, basename(holdout_dge)), compress = "xz")
saveRDS(meta_holdout, file.path(tmp, basename(holdout_meta)), compress = "xz")
fwrite(composition[nonholdout_idx], file.path(tmp, basename(nonholdout_composition)),
       sep = "\t", quote = FALSE, na = "NA")
fwrite(composition[holdout_idx], file.path(tmp, basename(holdout_composition)),
       sep = "\t", quote = FALSE, na = "NA")

manifest <- rbindlist(list(
  data.table(partition = "nonholdout", dataset = names(table(dge_dataset[nonholdout_idx])),
             n_samples = as.integer(table(dge_dataset[nonholdout_idx]))),
  data.table(partition = "holdout", dataset = names(table(dge_dataset[holdout_idx])),
             n_samples = as.integer(table(dge_dataset[holdout_idx])))
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
