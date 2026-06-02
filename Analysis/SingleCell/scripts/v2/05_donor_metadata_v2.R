#!/usr/bin/env Rscript
# ============================================================================
# 05_donor_metadata_v2.R
#
# Phase 0.5 v2 equivalent of Script 344. Builds the unified donor-level
# metadata table that drives the v2 stage-CCC chain (Scripts 06/07).
#
# Inputs (read-only on v1 paths; writes only under results_gpu_v2_phase05):
#   - V1 donor_metadata.tsv               (sample x core columns; READ-ONLY)
#   - V1 donor_fstage_documented.tsv      (READ-ONLY documented labels)
#   - V2 donor_fstage_scvi_v2_predicted.tsv
#   - V2 donor_fstage_augmented_v2_predicted.tsv
#   - V2 donor_fstage_cellfrac_v2.tsv     (Script 04 percell output)
#
# Output:
#   - Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

# --- v1 read-only inputs ----------------------------------------------------
DONOR_META_V1 <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
)
FSTAGE_DOC_V1 <- file.path(
  BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv"
)

# --- v2 read-write artifacts ------------------------------------------------
V2_STAGE_DIR <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2"
)
FSTAGE_SCVI_V2 <- file.path(V2_STAGE_DIR, "donor_fstage_scvi_v2_predicted.tsv")
FSTAGE_AUG_V2 <- file.path(V2_STAGE_DIR, "donor_fstage_augmented_v2_predicted.tsv")
FSTAGE_FRAC_V2 <- file.path(V2_STAGE_DIR, "donor_fstage_cellfrac_v2.tsv")

OUT_DIR <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_TSV <- file.path(OUT_DIR, "donor_metadata_v2.tsv")

# Optional v1 pseudotime / hep-subtype / bulk overlays (read-only, may be NA-OK)
PSEUDOTIME_V1 <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime/consensus_pseudotime_all.csv"
)
BULK_META <- file.path(
  BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
)

# ----------------------------------------------------------------------------
# 1. Load core donor metadata (v1 roster)
# ----------------------------------------------------------------------------
donor <- fread(DONOR_META_V1)
cat(sprintf("[meta] loaded %d donors from v1 roster\n", nrow(donor)))

# ----------------------------------------------------------------------------
# 2. v1 documented F-stage (read-only, donor-level — safe to use cross-version)
# ----------------------------------------------------------------------------
if (file.exists(FSTAGE_DOC_V1)) {
  doc <- fread(FSTAGE_DOC_V1)
  doc[, F_stage_documented := as.integer(F_stage_documented)]
  keep_doc <- c("sample", "dataset", "F_stage_documented")
  if ("source" %in% names(doc)) {
    keep_doc <- c(keep_doc, "source")
  }
  donor <- merge(donor, doc[, ..keep_doc],
                 by = c("sample", "dataset"), all.x = TRUE)
  if ("source" %in% names(donor)) {
    setnames(donor, "source", "F_stage_doc_source")
  }
  cat(sprintf("[fstage_v1_doc] %d donors with documented F-stage (v1, read-only)\n",
              sum(!is.na(donor$F_stage_documented))))
} else {
  donor[, F_stage_documented := NA_integer_]
  donor[, F_stage_doc_source := NA_character_]
  cat("[fstage_v1_doc] file missing; F_stage_documented = NA\n")
}

# ----------------------------------------------------------------------------
# 3. v2 scVI-predicted F-stage (343b equivalent)
# ----------------------------------------------------------------------------
if (file.exists(FSTAGE_SCVI_V2)) {
  fs <- fread(FSTAGE_SCVI_V2)
  if ("F_stage_predicted_argmax" %in% names(fs)) {
    fs[, F_stage_predicted_argmax_v2 := as.integer(F_stage_predicted_argmax)]
  }
  pcols <- grep("^P_F[0-4]$", names(fs), value = TRUE)
  if (length(pcols) > 0) {
    setnames(fs, pcols, paste0(pcols, "_scvi_v2"))
  }
  keep_fs <- c("sample", "dataset", "F_stage_predicted_argmax_v2",
               paste0("P_F", 0:4, "_scvi_v2"))
  keep_fs <- intersect(keep_fs, names(fs))
  donor <- merge(donor, fs[, ..keep_fs],
                 by = c("sample", "dataset"), all.x = TRUE)
  cat(sprintf("[fstage_v2_scvi] %d donors with scVI-v2 F-stage\n",
              sum(!is.na(donor$F_stage_predicted_argmax_v2))))
} else {
  cat(sprintf("[fstage_v2_scvi] missing %s\n", FSTAGE_SCVI_V2))
}

# ----------------------------------------------------------------------------
# 4. v2 augmented F-stage (343m equivalent) — PRIMARY downstream driver
# ----------------------------------------------------------------------------
if (file.exists(FSTAGE_AUG_V2)) {
  fa <- fread(FSTAGE_AUG_V2)
  ord_col <- "F_stage_pred_augmented_scvi_v2_ordinal"
  prob_cols <- grep("^P_F[0-4]_augmented_scvi_v2_ordinal$",
                    names(fa), value = TRUE)
  keep_fa <- c("sample", "dataset", ord_col, "origin", prob_cols)
  keep_fa <- intersect(keep_fa, names(fa))
  donor <- merge(donor, fa[, ..keep_fa],
                 by = c("sample", "dataset"), all.x = TRUE)
  if (ord_col %in% names(donor)) {
    donor[, F_stage_augmented_v2 := as.integer(get(ord_col))]
  } else {
    donor[, F_stage_augmented_v2 := NA_integer_]
  }
  cat(sprintf("[fstage_v2_aug] %d donors with augmented v2 F-stage\n",
              sum(!is.na(donor$F_stage_augmented_v2))))
} else {
  donor[, F_stage_augmented_v2 := NA_integer_]
  cat(sprintf("[fstage_v2_aug] missing %s\n", FSTAGE_AUG_V2))
}

# ----------------------------------------------------------------------------
# 5. Choose canonical F_stage_inferred_v2 (preferring augmented > scvi > documented)
# ----------------------------------------------------------------------------
donor[, F_stage_inferred_v2 := fcase(
  !is.na(F_stage_documented), F_stage_documented,
  !is.na(F_stage_augmented_v2), F_stage_augmented_v2,
  !is.na(F_stage_predicted_argmax_v2), F_stage_predicted_argmax_v2,
  default = NA_integer_
)]
donor[, F_stage_source_v2 := fcase(
  !is.na(F_stage_documented), "documented_v1",
  !is.na(F_stage_augmented_v2), "augmented_v2",
  !is.na(F_stage_predicted_argmax_v2), "scvi_v2",
  default = NA_character_
)]
cat(sprintf("[F_stage_inferred_v2] %d / %d donors with F-stage label\n",
            sum(!is.na(donor$F_stage_inferred_v2)), nrow(donor)))

# ----------------------------------------------------------------------------
# 6. v2 per-donor cellfrac (entropy, dominant_fstage)
# ----------------------------------------------------------------------------
if (file.exists(FSTAGE_FRAC_V2)) {
  fr <- fread(FSTAGE_FRAC_V2)
  setnames(fr, "n_cells", "n_hepatocytes_v2", skip_absent = TRUE)
  keep_fr <- intersect(
    c("sample", "dataset", "F0_frac", "F1_frac", "F2_frac",
      "F3_frac", "F4_frac", "n_hepatocytes_v2",
      "dominant_fstage", "entropy"),
    names(fr)
  )
  fr <- fr[, ..keep_fr]
  setnames(fr,
           old = c("F0_frac", "F1_frac", "F2_frac", "F3_frac", "F4_frac",
                   "dominant_fstage", "entropy"),
           new = c("F0_frac_v2", "F1_frac_v2", "F2_frac_v2",
                   "F3_frac_v2", "F4_frac_v2",
                   "dominant_fstage_v2", "entropy_v2"),
           skip_absent = TRUE)
  donor <- merge(donor, fr, by = c("sample", "dataset"), all.x = TRUE)
  cat(sprintf("[cellfrac_v2] %d donors with v2 cellfrac\n",
              sum(!is.na(donor$entropy_v2))))
} else {
  cat(sprintf("[cellfrac_v2] missing %s\n", FSTAGE_FRAC_V2))
}

# ----------------------------------------------------------------------------
# 7. v1 pseudotime overlay (Macrophage + Hepatocyte donor mean) — read-only
# ----------------------------------------------------------------------------
if (file.exists(PSEUDOTIME_V1)) {
  pt <- fread(PSEUDOTIME_V1)
  # First column is cell_id; ensure the schema matches
  if (!"cell_type" %in% names(pt) && ncol(pt) >= 4) {
    setnames(pt, c("V1", "cell_type", "consensus_pseudotime",
                   "aligned_palantir", "aligned_dpt",
                   "aligned_monocle3", "n_methods"),
             skip_absent = TRUE)
    if ("V1" %in% names(pt)) setnames(pt, "V1", "cell_id")
  }
  if (all(c("cell_id", "cell_type", "consensus_pseudotime") %in% names(pt))) {
    pt[, sample := sub("_[^_]+$", "", cell_id)]
    mac_pt <- pt[cell_type == "Macrophages",
                 .(macrophage_pseudotime_mean = mean(consensus_pseudotime, na.rm = TRUE),
                   n_macrophages_with_pt = .N), by = sample]
    hep_pt <- pt[cell_type == "Hepatocytes",
                 .(hepatocyte_pseudotime_mean = mean(consensus_pseudotime, na.rm = TRUE),
                   n_hepatocytes_with_pt = .N), by = sample]
    donor <- merge(donor, mac_pt, by = "sample", all.x = TRUE)
    donor <- merge(donor, hep_pt, by = "sample", all.x = TRUE)
    cat(sprintf("[pseudotime_v1] %d donors with Mac pt; %d with Hep pt\n",
                sum(!is.na(donor$macrophage_pseudotime_mean)),
                sum(!is.na(donor$hepatocyte_pseudotime_mean))))
  }
}

# ----------------------------------------------------------------------------
# 8. Bulk overlay for age, sex_numeric (read-only)
# ----------------------------------------------------------------------------
if (file.exists(BULK_META)) {
  bulk <- fread(BULK_META)
  sample_col <- if ("sample_id" %in% names(bulk)) "sample_id" else "sample"
  bulk_keep_cols <- c(sample_col,
                      intersect(c("age", "sex", "nas_score",
                                  "diagnosis_harmonized"), names(bulk)))
  bulk_keep <- bulk[, ..bulk_keep_cols]
  setnames(bulk_keep, sample_col, "sample")
  if ("sex" %in% names(bulk_keep)) {
    bulk_keep[, sex_numeric := fifelse(
      tolower(sex) %in% c("f", "female"), 1L,
      fifelse(tolower(sex) %in% c("m", "male"), 0L, NA_integer_)
    )]
  }
  donor <- merge(donor, bulk_keep, by = "sample", all.x = TRUE)
  cat(sprintf("[bulk_overlay] %d age / %d sex_numeric\n",
              sum(!is.na(donor$age)),
              sum(!is.na(donor$sex_numeric))))
}

# ----------------------------------------------------------------------------
# 9. Backfill disease_stage_coarse / disease_stage_numeric if blank
# ----------------------------------------------------------------------------
if ("condition" %in% names(donor)) {
  donor[disease_stage_coarse == "" | is.na(disease_stage_coarse),
        disease_stage_coarse := fcase(
          grepl("healthy|control", condition, ignore.case = TRUE), "Healthy",
          !is.na(F_stage_documented) & F_stage_documented == 0L, "Steatosis",
          !is.na(F_stage_documented) & F_stage_documented %in% 1:3, "Steatohepatitis",
          !is.na(F_stage_documented) & F_stage_documented == 4L, "Cirrhosis",
          grepl("masld|nafld|nash|steatohepatitis", condition,
                ignore.case = TRUE), "Steatohepatitis",
          default = NA_character_
        )]
}
donor[disease_stage_coarse == "Healthy",          disease_stage_numeric := 0]
donor[disease_stage_coarse == "Steatosis",        disease_stage_numeric := 1]
donor[disease_stage_coarse == "Steatohepatitis",  disease_stage_numeric := 2]
donor[disease_stage_coarse == "Cirrhosis",        disease_stage_numeric := 3]

# ----------------------------------------------------------------------------
# 10. Write final donor_metadata_v2.tsv
# ----------------------------------------------------------------------------
fwrite(donor, OUT_TSV, sep = "\t")
cat(sprintf("\n[output] %s (%d rows x %d cols)\n",
            OUT_TSV, nrow(donor), ncol(donor)))
cat("\nDisease stage distribution:\n")
print(donor[, .N, by = disease_stage_coarse])
cat("\nF_stage_source_v2 distribution:\n")
print(donor[, .N, by = F_stage_source_v2])
