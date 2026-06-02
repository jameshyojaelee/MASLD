#!/usr/bin/env Rscript
# sex_v3/01_load_data.R
# ---------------------------------------------------------------------------
# Module 01 — Load + harmonize inputs for the mashr-based Bayesian
# sex-stratified DEG pipeline (Script 26 v3 refactor).
#
# Reads:
#   - merged_dge.rds                                  (canonical DGEList)
#   - unified_metadata.csv                            (age, diagnosis_harmonized)
#   - persample_celltype_proportions.csv              (composition: Hep/Mac/Endo/Chol)
#   - sex_xci/chrx_stratified_classification.csv      (A9 chrX annotations)
#   - sample_qc_report.csv                            (pass_sex)
#   - meta_matched.rds                                (canonical inferred_sex)
#   - config/human_datasets.yaml                      (include_in_mega filter)
#
# Applies:
#   - yaml include_in_mega = 5-cohort canonical (Suppli/Hoang/Govaere/Bril/Chen)
#   - pass_sex == TRUE filter (drop samples failing sex QC)
#   - Restrict to samples with all 4 composition fractions available
#   - Standardize composition covariates: scale() -> mean 0, sd 1
#
# Writes:
#   - intermediates/sex_v3_input.rds: list(dge, meta, composition, chrx_annot,
#                                          gene_names, mega_cohorts, n_samples)
# ---------------------------------------------------------------------------

set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(yaml)
})

# ---------------------------------------------------------------------------
# Paths + contrast routing
# ---------------------------------------------------------------------------
BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT    <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR   <- file.path(INT, "results/integration")

# Shared utilities (atomic writes + sessionInfo dump + contrast routing)
source(file.path(BASE,
                 "RNA-seq/Human/Patient_Cohorts/analysis/integration",
                 "scripts/sex_v3/sex_v3_utils.R"))

# Contrast-aware paths. Default (disease_vs_ctrl) keeps the legacy layout
# `sex_v3/intermediates/`; new contrasts get sibling subdirs
# `sex_v3/contrast_<name>/intermediates/`.
CSPEC  <- contrast_spec()
CPATHS <- contrast_paths()
SEXV3  <- CPATHS$sexv3
IDIR   <- CPATHS$idir
dir.create(IDIR, recursive = TRUE, showWarnings = FALSE)
cat("Contrast routing: CONTRAST_NAME=", CSPEC$name, " (", CSPEC$pretty_name, ")\n",
    "  IDIR  = ", IDIR, "\n",
    "  SEXV3 = ", SEXV3, "\n", sep = "")

UMETA  <- file.path(INT, "metadata/unified_metadata.csv")
PROP_F <- file.path(BASE, "RNA-seq/results/celltype_attribution",
                    "persample_celltype_proportions.csv")
CHRX_F <- file.path(BASE, "RNA-seq/results/audit_sensitivity",
                    "sex_xci/chrx_stratified_classification.csv")

cat("============================================================\n")
cat("Module 01 — sex_v3 input loader\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")

# ---------------------------------------------------------------------------
# 1) Load DGE + cohort filter (yaml include_in_mega)
# ---------------------------------------------------------------------------
cat("\n[1] Loading merged_dge.rds + applying yaml include_in_mega filter...\n")
dge <- readRDS(file.path(RDIR, "merged_dge.rds"))
cat("  Raw DGE samples:", ncol(dge), "  genes:", nrow(dge), "\n")

ycfg <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_cohorts <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg))
all_cohorts  <- names(ycfg)

if (isTRUE(CSPEC$use_yaml_mega)) {
  cohort_set <- mega_cohorts
  cat("  CONTRAST '", CSPEC$name,
      "' uses yaml include_in_mega cohorts (k=", length(cohort_set), "):\n   ",
      paste(cohort_set, collapse = ", "), "\n", sep = "")
} else {
  # MASH-vs-MASL drops the control requirement, so the 4 control-less cohorts
  # also contribute (per the power analysis: 7 usable cohorts at N=687).
  cohort_set <- all_cohorts
  cat("  CONTRAST '", CSPEC$name,
      "' uses ALL yaml cohorts (k=", length(cohort_set), "):\n   ",
      paste(cohort_set, collapse = ", "), "\n", sep = "")
}

keep_cohort <- dge$samples$dataset %in% cohort_set
dge_mega <- dge[, keep_cohort]
cat("  Samples after cohort filter:", ncol(dge_mega), "\n")

# ---------------------------------------------------------------------------
# 2) pass_sex filter (drop samples failing sex check)
# ---------------------------------------------------------------------------
cat("\n[2] Applying pass_sex filter from sample_qc_report.csv...\n")
qc_report <- fread(file.path(INT, "qc/sample_qc_report.csv"))
sex_pass_ids <- qc_report[pass_sex == TRUE, sample_id]
sex_fail <- !colnames(dge_mega) %in% sex_pass_ids
if (any(sex_fail)) {
  cat("  Removing", sum(sex_fail), "samples failing pass_sex check\n")
  dge_mega <- dge_mega[, !sex_fail]
}
cat("  Samples after pass_sex filter:", ncol(dge_mega), "\n")

# ---------------------------------------------------------------------------
# 3) Inferred sex from canonical meta_matched.rds
# ---------------------------------------------------------------------------
cat("\n[3] Loading meta_matched.rds for canonical inferred_sex...\n")
meta_new <- readRDS(file.path(RDIR, "meta_matched.rds"))
matched_sex <- meta_new$inferred_sex[match(colnames(dge_mega), meta_new$sample_id)]
cat("  Samples with non-NA inferred_sex:", sum(!is.na(matched_sex)), "\n")

# ---------------------------------------------------------------------------
# 4) Unified metadata (age + diagnosis_harmonized)
# ---------------------------------------------------------------------------
cat("\n[4] Loading unified_metadata.csv (age, diagnosis_harmonized)...\n")
umeta <- fread(UMETA)
cat("  unified_metadata rows:", nrow(umeta),
    "  with non-NA age:", sum(!is.na(umeta$age)), "\n")

umeta_sub <- umeta[match(colnames(dge_mega), umeta$sample_id),
                   .(sample_id, age, diagnosis_harmonized)]

# ---------------------------------------------------------------------------
# 5) Composition (4 cell-type fractions)
# ---------------------------------------------------------------------------
cat("\n[5] Loading cell-type proportions (Hep/Mac/Endo/Chol)...\n")
prop <- fread(PROP_F)
cat("  Proportions table:", nrow(prop), "samples,", ncol(prop), "cols\n")

prop_use <- prop[, .(
  sample_id   = sample_id,
  Hepatocytes = `Hepatocytes`,
  Macrophages = `Macrophages`,
  Endothelial = `Endothelial cells`,
  Cholangiocytes = `Cholangiocytes`
)]
cat("  Composition NA counts:\n",
    "    Hep:",  sum(is.na(prop_use$Hepatocytes)),
    "  Mac:",   sum(is.na(prop_use$Macrophages)),
    "  Endo:",  sum(is.na(prop_use$Endothelial)),
    "  Chol:",  sum(is.na(prop_use$Cholangiocytes)), "\n")

# ---------------------------------------------------------------------------
# 6) chrX annotation (A9)
# ---------------------------------------------------------------------------
cat("\n[6] Loading chrX stratified classification (A9)...\n")
chrx_annot <- fread(CHRX_F)
cat("  chrX annotated genes:", nrow(chrx_annot),
    "  chr_category breakdown:\n")
print(table(chrx_annot$chr_category, useNA = "ifany"))

# ---------------------------------------------------------------------------
# 7) Build unified metadata table; drop samples missing required covariates
# ---------------------------------------------------------------------------
cat("\n[7] Merging metadata + composition; require all 4 fractions + sex...\n")
# R1 Issue 6 / Major fix: enforce explicit factor levels c("F","M") so that
# downstream string-matching ("inferred_sexM") and Module 05's hardcoded F/M
# column labels can never silently swap if upstream produces "Female"/"Male"
# or 0/1.
sex_raw <- as.character(matched_sex)
# Coerce common alternate encodings to canonical "F"/"M"
sex_raw[sex_raw %in% c("Female", "female", "f", "0")] <- "F"
sex_raw[sex_raw %in% c("Male",   "male",   "m", "1")] <- "M"
non_canonical <- sex_raw[!is.na(sex_raw) & !sex_raw %in% c("F", "M")]
if (length(non_canonical) > 0) {
  stop("Non-canonical inferred_sex values encountered (",
       length(unique(non_canonical)), " unique): ",
       paste(unique(non_canonical), collapse = ","),
       " — refusing to silently coerce; fix upstream meta_matched.rds.")
}
info0 <- data.frame(
  sample_id    = colnames(dge_mega),
  group_binary = factor(dge_mega$samples$group_binary, levels = c("Control","Disease")),
  dataset      = factor(dge_mega$samples$dataset),
  inferred_sex = factor(sex_raw, levels = c("F", "M")),
  stringsAsFactors = FALSE
)
# Assert canonical levels for downstream safety
stopifnot(identical(levels(info0$inferred_sex), c("F", "M")))
info0 <- merge(info0, umeta_sub, by = "sample_id", all.x = TRUE)
# 2026-05-14 fix: coerce empty-string diagnosis to NA before factoring.
# unified_metadata.csv stores "" for 318/645 Disease samples; without
# coercion R makes "" an implicit factor level that is perfectly collinear
# with group_binary=="Disease". Removed from the dream M2 model (see 04),
# but keep this cleanup so any downstream consumer sees NA rather than "".
info0$diagnosis_harmonized[info0$diagnosis_harmonized == ""] <- NA
info0$diagnosis_harmonized <- factor(info0$diagnosis_harmonized,
                                      exclude = c(NA, ""))
info0 <- merge(info0, as.data.frame(prop_use), by = "sample_id", all.x = TRUE)

# Restrict to samples with sex + all 4 fractions present
has_sex   <- !is.na(info0$inferred_sex)
has_comp  <- !is.na(info0$Hepatocytes) & !is.na(info0$Macrophages) &
             !is.na(info0$Endothelial) & !is.na(info0$Cholangiocytes)
keep_idx  <- has_sex & has_comp
n_before  <- nrow(info0)
cat("  Samples with sex available:", sum(has_sex), "\n")
cat("  Samples with all 4 composition fractions:", sum(has_comp), "\n")
cat("  Samples with BOTH (kept):", sum(keep_idx), "/", n_before, "\n")

# Contrast-specific subset: filter by diagnosis_harmonized and rebuild
# group_binary so the dream formula sees Control = arm 0, Disease = arm 1
# regardless of which biological pairing is being tested.
diag_lvl <- as.character(info0$diagnosis_harmonized)
in_contrast <- diag_lvl %in% c(CSPEC$arm_ctrl_levels, CSPEC$arm_dis_levels) |
                (CSPEC$name == "disease_vs_ctrl" & is.na(diag_lvl))
# disease_vs_ctrl: legacy behavior accepted NA diagnosis_harmonized because
# group_binary already reflects Disease/Control regardless of diagnosis labels.
# New contrasts require diag_lvl to be one of the explicit levels.
if (CSPEC$name != "disease_vs_ctrl") {
  keep_idx <- keep_idx & in_contrast
  cat("  Samples in contrast '", CSPEC$name, "' (",
      paste(c(CSPEC$arm_ctrl_levels, CSPEC$arm_dis_levels), collapse = " | "),
      "): ", sum(in_contrast), "; kept after sex/comp: ", sum(keep_idx),
      "\n", sep = "")
}

info <- info0[keep_idx, , drop = FALSE]
rownames(info) <- info$sample_id

# Rebuild group_binary from diagnosis_harmonized per the contrast spec.
# disease_vs_ctrl preserves the existing group_binary (from dge$samples).
if (CSPEC$name != "disease_vs_ctrl") {
  diag_in_info <- as.character(info$diagnosis_harmonized)
  info$group_binary <- factor(
    fifelse(diag_in_info %in% CSPEC$arm_ctrl_levels, "Control", "Disease"),
    levels = c("Control", "Disease"))
  cat("  Rebuilt group_binary for contrast '", CSPEC$name, "':\n",
      "    Control arm = {", paste(CSPEC$arm_ctrl_levels, collapse = ", "), "}\n",
      "    Disease arm = {", paste(CSPEC$arm_dis_levels,  collapse = ", "), "}\n",
      sep = "")
}

# Filter DGE to matching samples and ensure column order matches info
dge_mega <- dge_mega[, colnames(dge_mega) %in% info$sample_id]
info <- info[match(colnames(dge_mega), info$sample_id), , drop = FALSE]
rownames(info) <- info$sample_id
stopifnot(identical(rownames(info), colnames(dge_mega)))
cat("  Final n_samples:", ncol(dge_mega), "  n_genes:", nrow(dge_mega), "\n")

# Drop unused dataset levels
info$dataset <- droplevels(info$dataset)
info$diagnosis_harmonized <- droplevels(info$diagnosis_harmonized)

cat("\n  Group x sex distribution:\n")
print(table(info$group_binary, info$inferred_sex))
cat("\n  Group x dataset distribution:\n")
print(table(info$group_binary, info$dataset))
cat("\n  diagnosis_harmonized distribution:\n")
print(table(info$diagnosis_harmonized, useNA = "ifany"))
cat("\n  Age summary (raw):\n")
print(summary(info$age))
cat("  NA age:", sum(is.na(info$age)), "/", nrow(info), "\n")

# ---------------------------------------------------------------------------
# 8) Standardize composition covariates (zero-mean unit-SD)
# ---------------------------------------------------------------------------
cat("\n[8] Standardizing 4 composition covariates (scale to mean 0 / sd 1)...\n")
composition <- info[, c("Hepatocytes","Macrophages","Endothelial","Cholangiocytes")]
composition_scaled <- as.data.frame(scale(composition))
# Copy back to info so downstream modules can pull either raw or scaled
info$Hepatocytes    <- composition_scaled$Hepatocytes
info$Macrophages    <- composition_scaled$Macrophages
info$Endothelial    <- composition_scaled$Endothelial
info$Cholangiocytes <- composition_scaled$Cholangiocytes
cat("  Scaled covariates: column means ~",
    paste(round(colMeans(composition_scaled), 4), collapse=", "),
    " sd ~", paste(round(apply(composition_scaled, 2, sd), 4), collapse=", "), "\n")

# ---------------------------------------------------------------------------
# 9) Persist input list
# ---------------------------------------------------------------------------
gene_names <- rownames(dge_mega)
out_list <- list(
  dge          = dge_mega,
  meta         = info,
  composition  = as.matrix(composition_scaled),
  chrx_annot   = chrx_annot,
  gene_names   = gene_names,
  mega_cohorts = mega_cohorts,
  n_samples    = ncol(dge_mega),
  n_genes      = nrow(dge_mega)
)

out_path <- file.path(IDIR, "sex_v3_input.rds")
write_atomic_rds(out_list, out_path)
cat("\n[9] Saved input bundle (atomic tmp -> rename):\n  ", out_path, "\n")
cat("      $dge:         DGEList ", ncol(dge_mega), "samples x ", nrow(dge_mega), "genes\n")
cat("      $meta:        data.frame ", nrow(info), "rows x ", ncol(info), "cols\n")
cat("      $composition: matrix    ", nrow(composition_scaled), "x", ncol(composition_scaled),
    "(scaled fractions)\n")
cat("      $chrx_annot:  data.table ", nrow(chrx_annot), "rows\n")
cat("      $gene_names:  vector    ", length(gene_names), "\n")

# R5 Issue 2 fix: sessionInfo dump for per-module reproducibility audit
dump_session_info(IDIR, "01")

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
