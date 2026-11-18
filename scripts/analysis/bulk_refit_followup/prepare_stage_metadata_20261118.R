# Build corrected candidate primary_bulk_metadata_20261118.tsv,
# meta_matched_20261118.rds, and stage-corrected-20261118 release wrapper
# from candidate_20261118/unified_metadata.csv (GSE213621 Control -> fibrosis_stage = NA_integer_),
# while keeping the 774-unit binary discovery roster and merged_dge.rds completely unchanged.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CAND_META <- file.path(ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/candidate_20261118")
ORIG_774 <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
ORIG_REPROD <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls")
ORIG_ROSTER <- file.path(ROOT, "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/bulk_execution_inputs/primary_bulk_metadata.tsv")
STAGE_REL <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/stage-corrected-20261118/primary_source_controls")

dir.create(STAGE_REL, recursive = TRUE, showWarnings = FALSE)

um <- fread(file.path(CAND_META, "unified_metadata.csv"))
stopifnot(nrow(um[dataset == "GSE213621" & group_binary == "Control" & !is.na(fibrosis_stage)]) == 0L)
stopifnot(nrow(um[dataset == "GSE213621" & group_binary == "Control"]) > 0L)

# Update meta_matched.rds cleanly (preserving exact row order and columns of the 774 binary roster)
meta_orig <- readRDS(file.path(ORIG_774, "meta_matched.rds"))
meta_dt <- copy(as.data.table(meta_orig))
idx_um <- match(meta_dt$sample_id, um$sample_id)
stopifnot(!anyNA(idx_um))

n_213621_ctrl <- sum(meta_dt$dataset == "GSE213621" & meta_dt$group_binary == "Control")
stopifnot(n_213621_ctrl == 68L)
meta_dt[, fibrosis_stage := as.integer(um$fibrosis_stage[idx_um])]
stopifnot(all(is.na(meta_dt[dataset == "GSE213621" & group_binary == "Control"]$fibrosis_stage)))
stopifnot(identical(meta_dt$group_binary, as.data.table(meta_orig)$group_binary))
stopifnot(identical(meta_dt$sample_id, as.data.table(meta_orig)$sample_id))

meta_df <- as.data.frame(meta_dt)
rownames(meta_df) <- rownames(meta_orig)
attr(meta_df, "row.names") <- attr(meta_orig, "row.names")

meta_out_rds <- file.path(CAND_META, "meta_matched_20261118.rds")
saveRDS(meta_df, meta_out_rds)
saveRDS(meta_df, file.path(STAGE_REL, "meta_matched.rds"))

# Update primary_bulk_metadata_20261118.tsv
pbm <- fread(ORIG_ROSTER)
idx_pbm <- match(pbm$sample_id, um$sample_id)
stopifnot(!anyNA(idx_pbm))
pbm[, fibrosis_stage := as.integer(um$fibrosis_stage[idx_pbm])]
pbm_out <- file.path(CAND_META, "primary_bulk_metadata_20261118.tsv")
fwrite(pbm, pbm_out, sep = "\t", na = "NA")

# Write candidate stage-corrected-20261118 manifest & VALIDATED.json for build_pi_stage_extensions.R
mm <- fread(file.path(ORIG_774, "model_input_manifest.tsv"))
mm[, meta_path := meta_out_rds]
fwrite(mm, file.path(STAGE_REL, "model_input_manifest.tsv"), sep = "\t")
file.copy(file.path(ORIG_774, "expected_counts.json"), file.path(STAGE_REL, "expected_counts.json"), overwrite = TRUE)
val_obj <- list(
  status = "VALIDATED",
  arm = "primary_source_controls_stage_corrected_20261118",
  fit_root = file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z"),
  note = "GSE213621 Control fibrosis_stage set to NA_integer_; 774 binary roster and merged_dge.rds preserved untouched.",
  created_at_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
)
write_json(val_obj, file.path(STAGE_REL, "VALIDATED.json"), auto_unbox = TRUE, pretty = TRUE)

census <- meta_dt[, .(
  n_total = .N,
  n_control = sum(group_binary == "Control"),
  n_disease = sum(group_binary == "Disease"),
  n_fibrosis_non_na = sum(!is.na(fibrosis_stage)),
  n_nas_non_na = sum(!is.na(nas_score))
), by = dataset]
fwrite(census, file.path(CAND_META, "stage_census_20261118.tsv"), sep = "\t")
print(census)
cat(sprintf("Staged GSE162694+GSE213621 units: %d (GSE162694=%d, GSE213621=%d)\n",
            sum(census[dataset %in% c("GSE162694", "GSE213621")]$n_fibrosis_non_na),
            census[dataset == "GSE162694"]$n_fibrosis_non_na,
            census[dataset == "GSE213621"]$n_fibrosis_non_na))
cat(sprintf("All-cohort staged units: %d\n", sum(census$n_fibrosis_non_na)))
