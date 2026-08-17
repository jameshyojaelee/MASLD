#!/usr/bin/env Rscript
# ============================================================================
# 344_build_donor_metadata_extended.R
#
# Build the unified donor-level metadata table that drives all downstream
# stage-stratified CCC analyses (Scripts 345-349).
#
# Inputs:
#   - donor_metadata.tsv (sample x core columns)
#   - donor_fstage_documented.tsv (Script 343 output)
#   - consensus_pseudotime_all.csv (per-cell -> donor-mean for Macrophages)
#   - hepatocyte_subtype_metadata.csv (per-cell -> Progressor fraction)
#   - program_donor_scores_wide.tsv.gz (donor x 16 cNMF programs)
#   - unified_metadata.csv (bulk RNA-seq -> overlay age, sex_numeric)
#   - lineage_cell_counts.tsv (Script 345 pre-step, optional; if absent,
#     this script can also compute it directly via Python sidecar; see fallback)
#
# Output:
#   - donor_metadata_extended.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

# Master review M-P0-7: pre-registered bootstrap gate enforcement. If 343q
# wrote BOOTSTRAP_FALLBACK_REQUIRED.flag, F_stage_augmented_v2 failed
# stability acceptance and downstream LMM scripts must NOT treat it as the
# primary axis. Warn here so anyone rebuilding extended metadata sees the
# state up-front.
.boot_flag <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/BOOTSTRAP_FALLBACK_REQUIRED.flag")
if (file.exists(.boot_flag)) {
  warning("[bootstrap gate] BOOTSTRAP_FALLBACK_REQUIRED.flag exists: ",
          paste(readLines(.boot_flag), collapse = " | "),
          "\n  F_stage_augmented_v2 will still be written for sensitivity, ",
          "but downstream scripts treat F_stage_inferred_v2 (343b) as primary.")
}

DONOR_META <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv")
FSTAGE_DOC <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv")
FSTAGE_SCVI <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_scvi_predicted.tsv")
PSEUDOTIME <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime/consensus_pseudotime_all.csv")
HEP_SUBTYPE <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv")
PROGRAM_DON <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mcp/integration/program_donor_scores_wide.tsv.gz")
BULK_META   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
LIN_COUNTS  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/lineage_cell_counts.tsv")
FSTAGE_FRAC <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_cellfrac.tsv")

OUT_TSV <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv")

donor <- fread(DONOR_META)
cat(sprintf("[donor_metadata] %d donors\n", nrow(donor)))

# Documented F-stage (Script 343 output) ------------------------------------
if (file.exists(FSTAGE_DOC)) {
  fstage <- fread(FSTAGE_DOC)
  fstage[, F_stage_documented := as.integer(F_stage_documented)]
  donor <- merge(donor, fstage[, .(sample, dataset, F_stage_documented,
                                   F_stage_doc_source = source)],
                 by = c("sample", "dataset"), all.x = TRUE)
  cat(sprintf("[fstage_documented] %d donors with documented F-stage\n",
              sum(!is.na(donor$F_stage_documented))))
} else {
  donor[, F_stage_documented := NA_integer_]
  donor[, F_stage_doc_source := NA_character_]
  cat("[fstage_documented] file missing; F_stage_documented = NA for all donors\n")
}

# scVI-projected F-stage (Script 343b output) -------------------------------
# Trained on 58 Andrews donors with documented F-stage; LOOCV QWK ~0.74.
# Provides F_stage_predicted_argmax + P_F0..P_F4 for ~260 donors with
# hepatocytes in the integrated atlas.
if (file.exists(FSTAGE_SCVI)) {
  fscvi <- fread(FSTAGE_SCVI)
  fscvi[, F_stage_predicted_argmax := as.integer(F_stage_predicted_argmax)]
  scvi_cols <- intersect(c("F_stage_predicted_argmax",
                            "P_F0", "P_F1", "P_F2", "P_F3", "P_F4",
                            "classifier_qwk_loocv", "n_hepatocytes"),
                          names(fscvi))
  donor <- merge(donor,
                 fscvi[, c("sample", "dataset", scvi_cols), with = FALSE],
                 by = c("sample", "dataset"), all.x = TRUE)
  cat(sprintf("[fstage_scvi] %d donors with scVI-projected F-stage\n",
              sum(!is.na(donor$F_stage_predicted_argmax))))
} else {
  donor[, F_stage_predicted_argmax := NA_integer_]
  cat("[fstage_scvi] file missing; F_stage_predicted_argmax = NA for all donors\n")
}

# F_stage_inferred = documented when available, else scVI argmax ------------
# F_stage_source: "documented" | "scvi_predicted" | "scvi_predicted_dropped" | NA
#
# IMPORTANT (added 2026-05-13 after external validation, Script 343d):
# scVI projection fails to transfer to most non-Andrews cohorts -- it
# hallucinates F4 in healthy donors from Liver_Atlas / GSE185477 / GSE189600,
# and Spearman rho with disease_stage_numeric is only 0.04 in GSE136103.
# Only Andrews (GSE202379, training cohort, rho=0.90) and Wang (GSE244832,
# rho=0.60) transfer reliably. We therefore restrict scVI-projected
# F_stage_inferred to these two cohorts; predictions for other datasets are
# dropped to NA and tagged "scvi_predicted_dropped" so downstream scripts
# treat those donors as F-stage-unlabeled.
SCVI_TRANSFERABLE <- c("GSE202379", "GSE244832")
donor[, F_stage_inferred := fifelse(
        !is.na(F_stage_documented), F_stage_documented,
        F_stage_predicted_argmax)]
donor[, F_stage_source := fcase(
        !is.na(F_stage_documented),         "documented",
        !is.na(F_stage_predicted_argmax),   "scvi_predicted",
        default = NA_character_)]
# Gate scVI predictions to transferable cohorts only.
n_dropped <- sum(donor$F_stage_source == "scvi_predicted" &
                 !donor$dataset %in% SCVI_TRANSFERABLE, na.rm = TRUE)
donor[F_stage_source == "scvi_predicted" &
      !dataset %in% SCVI_TRANSFERABLE,
      `:=`(F_stage_inferred = NA_integer_,
           F_stage_source   = "scvi_predicted_dropped")]
cat(sprintf("[fstage_inferred] %d / %d donors with F_stage_inferred (post-gate)\n",
            sum(!is.na(donor$F_stage_inferred)), nrow(donor)))
cat(sprintf("  - documented:              %d\n", sum(donor$F_stage_source == "documented", na.rm = TRUE)))
cat(sprintf("  - scvi_predicted (kept):   %d (only GSE202379 / GSE244832)\n",
            sum(donor$F_stage_source == "scvi_predicted", na.rm = TRUE)))
cat(sprintf("  - scvi_predicted_dropped:  %d (failed external validation)\n", n_dropped))

# Augmented scVI F-stage (Script 343m output) ------------------------------
# Cross-dataset F0/F4 anchors (clean-healthy and Cirrhosis from
# disease_stage_coarse) added to the 58 Andrews documented donors as
# training set, expanding from 58 -> ~127 anchors with cross-cohort
# coverage. External Spearman rho jumped from 0.391 (vanilla scVI) to 0.639;
# Healthy->F4 hallucinations dropped from 27 to 2. This is now the PRIMARY
# F-stage axis driver across the full atlas (no cohort gating needed -
# augmented training itself solves the OOD transfer problem).
FSTAGE_AUG <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_augmented_predicted.tsv")
if (file.exists(FSTAGE_AUG)) {
  faug <- fread(FSTAGE_AUG)
  # Actual column names in donor_fstage_augmented_predicted.tsv (Script 343m):
  #   F_stage_pred_augmented_scvi_ordinal  (argmax from ordinal logistic, headline)
  #   F_stage_pred_augmented_scvi_knn5     (argmax from kNN-5)
  #   P_F0..P_F4_augmented_scvi_ordinal    (posteriors)
  ordinal_col <- "F_stage_pred_augmented_scvi_ordinal"
  knn_col     <- "F_stage_pred_augmented_scvi_knn5"
  posterior_cols <- grep("^P_F[0-4]_augmented_scvi_ordinal$",
                         names(faug), value = TRUE)
  keep_cols <- intersect(c("sample", "dataset",
                           ordinal_col, knn_col, posterior_cols),
                         names(faug))
  donor <- merge(donor, faug[, ..keep_cols],
                 by = c("sample", "dataset"), all.x = TRUE)
  if (ordinal_col %in% names(donor)) {
    donor[, F_stage_augmented := as.integer(get(ordinal_col))]
  } else {
    donor[, F_stage_augmented := NA_integer_]
  }
  cat(sprintf("[fstage_augmented] %d donors with augmented F-stage prediction (ordinal)\n",
              sum(!is.na(donor$F_stage_augmented))))
} else {
  donor[, F_stage_augmented := NA_integer_]
  cat("[fstage_augmented] file missing; F_stage_augmented = NA for all donors\n")
}

# Macrophage donor-mean consensus pseudotime --------------------------------
pt <- fread(PSEUDOTIME)
setnames(pt, c("V1", "cell_type", "consensus_pseudotime",
               "aligned_palantir", "aligned_dpt", "aligned_monocle3", "n_methods"))
setnames(pt, "V1", "cell_id")
# cell_id format: GSM4041150_AGTGTCAAGCAGCCTC-1 -> sample = GSM4041150
pt[, sample := sub("_[^_]+$", "", cell_id)]
mac_pt <- pt[cell_type == "Macrophages",
             .(macrophage_pseudotime_mean = mean(consensus_pseudotime, na.rm = TRUE),
               n_macrophages_with_pt      = .N),
             by = sample]
hep_pt <- pt[cell_type == "Hepatocytes",
             .(hepatocyte_pseudotime_mean = mean(consensus_pseudotime, na.rm = TRUE),
               n_hepatocytes_with_pt      = .N),
             by = sample]
donor <- merge(donor, mac_pt, by = "sample", all.x = TRUE)
donor <- merge(donor, hep_pt, by = "sample", all.x = TRUE)
cat(sprintf("[pseudotime] %d donors with Mac pseudotime; %d with Hep pseudotime\n",
            sum(!is.na(donor$macrophage_pseudotime_mean)),
            sum(!is.na(donor$hepatocyte_pseudotime_mean))))

# Progressor fraction (hepatocyte meta-subtype) -----------------------------
hep_meta <- fread(HEP_SUBTYPE)
if ("hepatocyte_subtype_label" %in% names(hep_meta)) {
  prog <- hep_meta[, .(progressor_frac = mean(grepl("rogressor",
                                                    hepatocyte_subtype_label,
                                                    ignore.case = TRUE)),
                       n_hepatocytes_hep_atlas = .N),
                   by = sample]
} else {
  prog <- data.table(sample = unique(donor$sample),
                     progressor_frac = NA_real_,
                     n_hepatocytes_hep_atlas = NA_integer_)
}
donor <- merge(donor, prog, by = "sample", all.x = TRUE)
cat(sprintf("[progressor_frac] %d donors with hepatocyte subtype info\n",
            sum(!is.na(donor$progressor_frac))))

# cNMF program scores --------------------------------------------------------
prog_scores <- fread(PROGRAM_DON)
keep_progs <- intersect(c("cnmf_global_k16_P1", "cnmf_global_k16_P5",
                          "cnmf_global_k16_P11", "cnmf_global_k16_P13",
                          "cnmf_global_k16_P14"),
                        names(prog_scores))
if (length(keep_progs) > 0) {
  donor <- merge(donor, prog_scores[, c("sample", keep_progs), with = FALSE],
                 by = "sample", all.x = TRUE)
}
cat(sprintf("[cnmf] joined %d program score columns\n", length(keep_progs)))

# Bulk overlay for age, sex, NAS (where sample matches) ---------------------
if (file.exists(BULK_META)) {
  bulk <- fread(BULK_META)
  sample_col <- if ("sample_id" %in% names(bulk)) "sample_id" else "sample"
  bulk_keep <- bulk[, c(sample_col,
                        intersect(c("age", "sex", "nas_score",
                                    "diagnosis_harmonized"),
                                  names(bulk))),
                    with = FALSE]
  setnames(bulk_keep, sample_col, "sample")
  if ("sex" %in% names(bulk_keep)) {
    bulk_keep[, sex_numeric := fifelse(tolower(sex) %in% c("f", "female"), 1L,
                              fifelse(tolower(sex) %in% c("m", "male"), 0L, NA_integer_))]
  }
  donor <- merge(donor, bulk_keep, by = "sample", all.x = TRUE)
  cat(sprintf("[bulk_overlay] joined; %d donors with age, %d with sex\n",
              sum(!is.na(donor$age)), sum(!is.na(donor$sex_numeric))))
}

# Backfill blank disease_stage_coarse from condition + F_stage_documented --
# Some donors (notably the entire GSE202379 cohort and a few Healthy controls
# from Liver_Atlas / GSE136103) have a blank `disease_stage_coarse` field in
# the upstream donor_metadata.tsv. Backfill rules:
#   condition == "Healthy"                       -> "Healthy"
#   source diagnosis NAFLD/MASL                    -> "Steatosis"
#   source diagnosis NASH/MASH                     -> "Steatohepatitis"
#   condition matches MASLD/NAFL  & F_stage == 4 -> "Cirrhosis"
# Generic MASLD is not itself a stage label and is not backfilled from fibrosis.
donor[disease_stage_coarse == "" | is.na(disease_stage_coarse),
      disease_stage_coarse := fcase(
        grepl("healthy|control", condition, ignore.case = TRUE), "Healthy",
        !is.na(F_stage_documented) & F_stage_documented == 4L, "Cirrhosis",
        grepl("nafld|masl", condition, ignore.case = TRUE), "Steatosis",
        grepl("nash|mash|steatohepatitis", condition,
              ignore.case = TRUE), "Steatohepatitis",
        default = NA_character_)]
donor[disease_stage_coarse == "Healthy",       disease_stage_numeric := 0]
donor[disease_stage_coarse == "Steatosis",     disease_stage_numeric := 1]
donor[disease_stage_coarse == "Steatohepatitis", disease_stage_numeric := 2]
donor[disease_stage_coarse == "Cirrhosis",     disease_stage_numeric := 3]
cat(sprintf("[backfill] disease_stage_coarse populated for %d donors (was %d before)\n",
            sum(!is.na(donor$disease_stage_coarse) & donor$disease_stage_coarse != ""),
            sum(!is.na(donor$disease_stage_coarse) & donor$disease_stage_coarse != "" &
                donor$disease_stage_coarse != "Steatohepatitis")))

# 9-lineage cell counts (from optional sidecar) -----------------------------
if (file.exists(LIN_COUNTS)) {
  lin <- fread(LIN_COUNTS)
  donor <- merge(donor, lin, by = "sample", all.x = TRUE)
  cat("[lineage_counts] joined from Script 345 sidecar\n")
} else {
  cat("[lineage_counts] sidecar missing -> Script 345 will compute per-donor counts at runtime\n")
}

# clean_healthy_flag / dubious_healthy_flag --------------------------------
# Derived from Script 343j's audit logic, computed here so a fresh atlas
# build reproduces it without depending on Script 343j running first.
#   clean_healthy   = Healthy AND dominant_fstage <= 1 AND F4_frac < 0.1
#   dubious_healthy = Healthy AND (dominant_fstage >= 3 OR F4_frac > 0.4)
# Convergent F3/F4 inference across scVI, fib-signature, HSC counting and
# scVI-kNN on these "Healthy" donors flagged them as biologically distinct
# from biopsy-grade healthy liver (see Script 343j + figS_dubious_healthy).
donor[, c("clean_healthy_flag", "dubious_healthy_flag") :=
        .(FALSE, FALSE)]
if (file.exists(FSTAGE_FRAC)) {
  fracdt <- fread(FSTAGE_FRAC)
  fracdt[, dominant_fstage := as.integer(dominant_fstage)]
  fracdt[, F4_frac := as.numeric(F4_frac)]
  donor <- merge(donor,
                 fracdt[, .(sample, dataset, dominant_fstage, F4_frac)],
                 by = c("sample", "dataset"), all.x = TRUE)
  donor[, clean_healthy_flag := disease_stage_coarse == "Healthy" &
                                !is.na(dominant_fstage) &
                                dominant_fstage <= 1 &
                                !is.na(F4_frac) &
                                F4_frac < 0.1]
  donor[, dubious_healthy_flag := disease_stage_coarse == "Healthy" &
                                  ((!is.na(dominant_fstage) &
                                    dominant_fstage >= 3) |
                                   (!is.na(F4_frac) & F4_frac > 0.4))]
  cat(sprintf("[clean_healthy_flag] %d clean / %d dubious among %d Healthy donors\n",
              sum(donor$clean_healthy_flag, na.rm = TRUE),
              sum(donor$dubious_healthy_flag, na.rm = TRUE),
              sum(donor$disease_stage_coarse == "Healthy", na.rm = TRUE)))
} else {
  cat("[clean_healthy_flag] donor_fstage_cellfrac.tsv missing -> all FALSE\n")
}

# === Protocol contamination remediation (added 2026-05-22) ===
# Protocol group classification — 7 datasets total, 3 protocols
donor[, protocol_group := fcase(
  dataset %in% c("GSE244832", "GSE185477", "GSE174748"), "whole_tissue",
  dataset %in% c("GSE136103", "Liver_Atlas"),             "npc_enriched",
  dataset == "GSE202379",                                  "paired_sc_sn",
  dataset == "GSE189600",                                  "hep_specific",
  default = "unknown"
)]
stopifnot(sum(is.na(donor$protocol_group)) == 0)

# Exclude all GSE136103 + Liver_Atlas donors from stage-stratified analyses
donor[, exclude_stage_analysis := dataset %in% c("GSE136103", "Liver_Atlas")]

# Clean F-stage: NA for excluded donors (dubious_healthy_flag is all FALSE per pre-flight,
# but we keep the structure for future-proofing)
donor[, F_stage_augmented_clean := fifelse(
  exclude_stage_analysis |
    (!is.na(dubious_healthy_flag) & dubious_healthy_flag == TRUE & F_stage_augmented >= 3),
  NA_real_,
  as.numeric(F_stage_augmented)
)]
cat(sprintf("[protocol_remediation] protocol_group: %s\n",
            paste(sprintf("%s=%d", names(table(donor$protocol_group)),
                          as.integer(table(donor$protocol_group))), collapse = ", ")))
cat(sprintf("[protocol_remediation] exclude_stage_analysis = TRUE for %d donors\n",
            sum(donor$exclude_stage_analysis)))
cat(sprintf("[protocol_remediation] F_stage_augmented_clean non-NA: %d (was %d)\n",
            sum(!is.na(donor$F_stage_augmented_clean)),
            sum(!is.na(donor$F_stage_augmented))))
# === End protocol remediation ===

# Write -------------------------------------------------------------------
fwrite(donor, OUT_TSV, sep = "\t")
cat(sprintf("\n[output] %s (%d rows x %d cols)\n",
            OUT_TSV, nrow(donor), ncol(donor)))
cat("\nDisease stage distribution:\n")
print(donor[, .N, by = disease_stage_coarse])
cat("\nF-stage documented distribution:\n")
print(donor[, .N, by = F_stage_documented])
