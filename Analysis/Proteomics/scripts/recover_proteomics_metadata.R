#!/usr/bin/env Rscript
# recover_proteomics_metadata.R — Recover disease metadata for proteomics datasets
# Change 0 of the 7→4 figure restructuring plan
# Outputs: pxd052937_disease_metadata.csv
#
# NOTE: GSE276114 was removed. GEO confirms it is bulk RNA-seq
# (Expression profiling by high throughput sequencing), not SomaScan proteomics.

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PROTEO_DIR <- file.path(BASE, "Analysis/Proteomics")
RESULTS_DIR <- file.path(PROTEO_DIR, "results")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  Recover Proteomics Metadata\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Project root:", BASE, "\n")
cat("============================================================\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# Part 2: PXD052937 (DIA MassSpec, 3,346 proteins x 72 plasma samples)
# ══════════════════════════════════════════════════════════════════════════════

cat("── Part 2: PXD052937 (DIA plasma proteomics) ──────────────────────────\n\n")

# ┌──────────────────────────────────────────────────────────────────────────┐
# │ !!! TENTATIVE CONDITION MAPPING — NOT VERIFIED FROM SOURCE PUBLICATION !!! │
# │                                                                            │
# │ The A/B/C/D -> Normal/MASL/MASH/Cirrhosis assignment below is INFERRED     │
# │ from per-group sample sizes and typical MASLD cohort designs ONLY. It has  │
# │ NOT been confirmed against the PXD052937 source publication. Any analysis  │
# │ that relies on these labels is provisional until the mapping is verified.  │
# │                                                                            │
# │ Inferred (UNVERIFIED) mapping:                                             │
# │   A (7 samples)  = Normal/Healthy control                                  │
# │   B (17 samples) = MASL/Steatosis (simple steatosis)                       │
# │   C (38 samples) = MASH/Steatohepatitis                                    │
# │   D (10 samples) = Cirrhosis/Advanced fibrosis                             │
# │                                                                            │
# │ ACTION REQUIRED: verify the letter->condition assignment against          │
# │   Niu et al. 2024 (PRIDE PXD052937) — the sample/condition table in the   │
# │   source publication — BEFORE using these labels in any analysis.          │
# │ A runtime warning() is emitted below whenever this mapping is applied.     │
# └──────────────────────────────────────────────────────────────────────────┘

# Read existing metadata
pxd_meta_file <- file.path(RESULTS_DIR, "pxd052937_metadata.csv")
if (!file.exists(pxd_meta_file)) {
  stop("PXD052937 metadata file not found: ", pxd_meta_file)
}

pxd_raw <- fread(pxd_meta_file)
cat("Loaded PXD052937 metadata:", nrow(pxd_raw), "samples\n")
cat("Conditions in file:", paste(sort(unique(pxd_raw$Condition)), collapse = ", "), "\n")
cat("Samples per condition:\n")
print(pxd_raw[, .N, by = Condition][order(Condition)])
cat("\n")

# TENTATIVE condition mapping — MUST BE VERIFIED FROM PUBLICATION
condition_map <- data.table(
  condition_letter = c("A", "B", "C", "D"),
  condition_label  = c("Normal", "MASL", "MASH", "Cirrhosis"),
  condition_detail = c("Healthy control",
                       "Simple steatosis (MASL/NAFL)",
                       "Steatohepatitis (MASH/NASH)",
                       "Advanced fibrosis / Cirrhosis"),
  is_masld         = c(FALSE, TRUE, TRUE, TRUE)
)

cat("TENTATIVE condition mapping (NEEDS VERIFICATION):\n")
print(condition_map)
cat("\n")

# Build disease metadata table
pxd_disease <- data.table(
  sample_id        = pxd_raw$`Run Label`,
  replicate        = pxd_raw$Replicate,
  condition_letter = pxd_raw$Condition,
  file_name        = pxd_raw$`File Name`
)

# Merge in condition labels.
# Emit a LOUD runtime warning whenever this UNVERIFIED mapping is applied so it
# can never be silently trusted downstream (see TENTATIVE block above).
warning(
  "PXD052937: applying TENTATIVE/UNVERIFIED condition mapping ",
  "(A=Normal, B=MASL, C=MASH, D=Cirrhosis). This was inferred from sample ",
  "sizes only and MUST be verified against Niu et al. 2024 (PRIDE PXD052937) ",
  "before these labels are used in any analysis.",
  call. = FALSE, immediate. = TRUE
)
pxd_disease <- merge(pxd_disease, condition_map, by = "condition_letter", all.x = TRUE)

# Reorder columns
setcolorder(pxd_disease, c("sample_id", "condition_letter", "condition_label",
                           "condition_detail", "is_masld", "replicate", "file_name"))

# Sort by condition then replicate
setorder(pxd_disease, condition_letter, replicate)

cat("PXD052937 disease metadata summary:\n")
cat("  Total samples:", nrow(pxd_disease), "\n")
cat("  MASLD samples:", sum(pxd_disease$is_masld), "\n")
cat("  Control samples:", sum(!pxd_disease$is_masld), "\n")
cat("  Samples per condition:\n")
print(pxd_disease[, .N, by = .(condition_letter, condition_label)])
cat("\n")

# Save
pxd_out_file <- file.path(RESULTS_DIR, "pxd052937_disease_metadata.csv")
fwrite(pxd_disease, pxd_out_file)
cat("Saved:", pxd_out_file, "\n")
cat("Dimensions:", nrow(pxd_disease), "x", ncol(pxd_disease), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("============================================================\n")
cat("  Summary\n")
cat("============================================================\n")
cat("\n")

cat("PXD052937 (DIA, plasma):\n")
cat("  Samples with condition labels:", nrow(pxd_disease), "\n")
cat("  STATUS: TENTATIVE mapping applied (A=Normal, B=MASL, C=MASH, D=Cirrhosis).\n")
cat("          MUST verify from Niu et al. 2024 publication.\n")

cat("\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
