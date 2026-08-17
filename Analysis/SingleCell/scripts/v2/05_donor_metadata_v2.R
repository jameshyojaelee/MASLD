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
          !is.na(F_stage_documented) & F_stage_documented == 4L, "Cirrhosis",
          grepl("nafld|masl", condition, ignore.case = TRUE), "Steatosis",
          grepl("nash|mash|steatohepatitis", condition,
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

# ============================================================================
# 11. DONOR-COLLAPSE (pseudoreplication fix, 2026-07-12)
#
# obs["sample"] in the v2 atlas is a SEQUENCING RUN (SRR / GSM), not a
# biological donor. The run-level donor_metadata_v2.tsv written above lets the
# downstream stage-CCC chain (06/07/ccc_v3_panels) count runs as independent
# donors -- pseudoreplication. We collapse runs -> biological donor here and
# write a SECOND file keyed on the TRUE donor (sample := biological_donor),
# matching the LIANA donor-collapse in 06_per_donor_liana_v2.py and mirroring
# the v1 fix (345/346). The run-level file above is left intact for before/
# after comparison.
# ============================================================================
DONOR_PAIRING <- c(
  GSE244832 = file.path(BASE, "data/GSE244832/metadata/donor_pairing.csv"),
  GSE202379 = file.path(BASE, "data/GSE202379/metadata/donor_pairing.csv"),
  GSE185477 = file.path(BASE, "data/GSE185477/metadata/donor_pairing.csv"),
  GSE136103 = file.path(BASE, "data/GSE136103/metadata/donor_pairing.csv")
)
srr2donor <- new.env(parent = emptyenv())
for (ds in names(DONOR_PAIRING)) {
  fp <- DONOR_PAIRING[[ds]]
  if (!file.exists(fp)) {
    cat(sprintf("[collapse] WARNING donor_pairing.csv missing for %s: %s\n", ds, fp))
    next
  }
  dp <- fread(fp, colClasses = "character")
  for (i in seq_len(nrow(dp))) {
    did <- paste0(ds, "_", dp$donor_id[i])
    runs <- trimws(strsplit(dp$rna_srrs[i], ";", fixed = TRUE)[[1]])
    runs <- runs[nzchar(runs)]
    for (r in runs) assign(r, did, envir = srr2donor)
  }
}
map_donor <- function(s) {
  s <- as.character(s)
  vapply(s, function(x) if (exists(x, envir = srr2donor, inherits = FALSE))
    get(x, envir = srr2donor) else x, character(1))
}
donor[, biological_donor := map_donor(sample)]

# Report within-donor stage conflicts (should be zero; runs of a donor share a stage)
conflicts <- donor[, .(n_stage = uniqueN(disease_stage_coarse)),
                   by = biological_donor][n_stage > 1]
if (nrow(conflicts) > 0) {
  cat(sprintf("[collapse] WARNING %d donors carry >1 disease_stage_coarse across runs:\n",
              nrow(conflicts)))
  print(conflicts)
} else {
  cat("[collapse] no within-donor disease_stage_coarse conflicts\n")
}

first_non_na <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[[1]] else x[NA_integer_] }
num_cols <- setdiff(names(donor)[vapply(donor, is.numeric, logical(1))],
                    "biological_donor")
chr_cols <- setdiff(names(donor), c(num_cols, "biological_donor"))
# numeric covariates -> donor mean (donor-invariant ones are unchanged);
# categorical -> first non-NA (stage/dataset/F-stage source are donor-invariant)
dc_num <- if (length(num_cols))
  donor[, lapply(.SD, function(x) { m <- mean(x, na.rm = TRUE);
        if (is.nan(m)) NA_real_ else m }),
        by = biological_donor, .SDcols = num_cols] else
  unique(donor[, .(biological_donor)])
dc_chr <- donor[, lapply(.SD, first_non_na),
                by = biological_donor, .SDcols = chr_cols]
dc <- merge(dc_num, dc_chr, by = "biological_donor")
dc[, sample := biological_donor]          # key on TRUE donor for downstream merges
setcolorder(dc, c("sample", setdiff(names(dc), "sample")))
OUT_DC <- file.path(OUT_DIR, "donor_metadata_v2_dc.tsv")
fwrite(dc, OUT_DC, sep = "\t")
cat(sprintf("\n[collapse] %d run-level rows -> %d donor-level rows -> %s\n",
            nrow(donor), nrow(dc), OUT_DC))
cat("\nDonor-collapsed disease stage distribution:\n")
print(dc[, .N, by = disease_stage_coarse])
