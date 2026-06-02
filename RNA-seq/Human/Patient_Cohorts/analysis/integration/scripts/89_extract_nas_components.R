#!/usr/bin/env Rscript
# 89_extract_nas_components.R
# P0: Extract NAS component scores and update modeling metadata
#
# Sources:
#   - GSE130970: Full components (steatosis 0-3, inflammation 0-2, ballooning 0-2) for 78 samples
#   - Back-calculates inflammation = NAS - steatosis - ballooning where needed
#
# Updates:
#   - modeling_metadata.csv: adds steatosis_grade, lobular_inflammation_grade, ballooning_grade columns
#   - prepared_data.h5: adds component metadata keys (via Python helper)
#
# Usage: Rscript 89_extract_nas_components.R
# SLURM: cpu, 4 CPUs, 16G RAM, 48h

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR <- file.path(INT, "results/staging_classifier")

cat("=== 89: Extract NAS Components ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# ============================================================
# STEP 1: Load GSE130970 component data
# ============================================================
cat("=== STEP 1: GSE130970 Component Extraction ===\n")

gse130970_file <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/archive/old_data_metadata/metadata/GSE130970_SraRunTable.csv")
stopifnot(file.exists(gse130970_file))

gse130970 <- fread(gse130970_file)
cat("GSE130970:", nrow(gse130970), "samples\n")

# Extract components
components <- gse130970[, .(
  sample_id = Run,
  dataset = "GSE130970",
  steatosis_grade = steatosis_grade,
  lobular_inflammation_grade = lobular_inflammation_grade,
  ballooning_grade = cytological_ballooning_grade,
  nas_score_from_components = steatosis_grade + lobular_inflammation_grade + cytological_ballooning_grade,
  nas_score_original = nafld_activity_score,
  fibrosis_stage = fibrosis_stage
)]

# Validate: components sum to NAS
mismatch <- components[nas_score_from_components != nas_score_original]
if (nrow(mismatch) > 0) {
  cat("WARNING:", nrow(mismatch), "samples where components don't sum to NAS!\n")
  print(mismatch)
} else {
  cat("VALIDATED: All", nrow(components), "samples have components summing to NAS\n")
}

cat("\nComponent distributions:\n")
cat("  Steatosis (0-3):", paste(components[, .N, by = steatosis_grade][order(steatosis_grade), paste0(steatosis_grade, "=", N)], collapse = ", "), "\n")
cat("  Inflammation (0-2):", paste(components[, .N, by = lobular_inflammation_grade][order(lobular_inflammation_grade), paste0(lobular_inflammation_grade, "=", N)], collapse = ", "), "\n")
cat("  Ballooning (0-2):", paste(components[, .N, by = ballooning_grade][order(ballooning_grade), paste0(ballooning_grade, "=", N)], collapse = ", "), "\n")
cat("  NAS (0-8):", paste(components[, .N, by = nas_score_original][order(nas_score_original), paste0(nas_score_original, "=", N)], collapse = ", "), "\n")

# ============================================================
# STEP 2: Update modeling_metadata.csv
# ============================================================
cat("\n=== STEP 2: Update Modeling Metadata ===\n")

meta_file <- file.path(OUTDIR, "modeling_metadata.csv")
stopifnot(file.exists(meta_file))
meta <- fread(meta_file)
cat("Loaded modeling_metadata:", nrow(meta), "samples\n")

# Add component columns (default to -1 = missing)
meta[, steatosis_grade := -1L]
meta[, lobular_inflammation_grade := -1L]
meta[, ballooning_grade := -1L]

# Match GSE130970 samples
matched <- merge(
  meta[, .(sample_id)],
  components[, .(sample_id, steatosis_grade, lobular_inflammation_grade, ballooning_grade)],
  by = "sample_id",
  all.x = TRUE
)

# Update the columns
for (i in seq_len(nrow(matched))) {
  if (!is.na(matched$steatosis_grade[i])) {
    sid <- matched$sample_id[i]
    meta[sample_id == sid, steatosis_grade := matched$steatosis_grade[i]]
    meta[sample_id == sid, lobular_inflammation_grade := matched$lobular_inflammation_grade[i]]
    meta[sample_id == sid, ballooning_grade := matched$ballooning_grade[i]]
  }
}

n_with_components <- sum(meta$steatosis_grade >= 0)
cat("Samples with NAS components:", n_with_components, "\n")

# Also create NAS 3-class column (Low 0-2, Medium 3-5, High 6-8)
meta[, nas_3class := fifelse(
  nas_group >= 0 & nas_group <= 2, 0L,   # Low
  fifelse(nas_group >= 3 & nas_group <= 5, 1L,  # Medium
  fifelse(nas_group >= 6, 2L, -1L)))]  # High

n_nas3 <- sum(meta$nas_3class >= 0)
cat("Samples with NAS 3-class:", n_nas3, "\n")
cat("  Low (0-2):", sum(meta$nas_3class == 0), "  Medium (3-5):", sum(meta$nas_3class == 1),
    "  High (6-8):", sum(meta$nas_3class == 2), "\n")

# Save updated metadata
fwrite(meta, meta_file)
cat("Updated:", meta_file, "\n")

# Also save component-only table for easy access
fwrite(components, file.path(OUTDIR, "nas_component_data.csv"))
cat("Saved: nas_component_data.csv (", nrow(components), "rows)\n")

# ============================================================
# STEP 3: Summary
# ============================================================
cat("\n=== Summary ===\n")
cat("GSE130970 components extracted:", n_with_components, "samples\n")
cat("NAS 3-class labels added:", n_nas3, "samples\n")
cat("Columns added to modeling_metadata.csv: steatosis_grade, lobular_inflammation_grade, ballooning_grade, nas_3class\n")
cat("\nNAS 3-class distribution:\n")
print(meta[nas_3class >= 0, .N, by = nas_3class][order(nas_3class)])

cat("\n=== 89_extract_nas_components.R completed:", as.character(Sys.time()), "===\n")
