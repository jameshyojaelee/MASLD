#!/usr/bin/env Rscript
# 00_harmonize_metadata.R
# ---------------------------------------------------------------------------
# Harmonize metadata from all human MASLD cohorts into a unified table.
# Dataset registry and excluded samples are read from:
#   config/human_datasets.yaml   (relative to project root)
#
# TO ADD A NEW DATASET:
#   1. Add an entry to config/human_datasets.yaml (sra_table_path, excluded_samples, etc.)
#   2. Add a case to the harmonize() switch below with the condition mapping
#   3. Re-run this script
#
# Output: analysis/integration/metadata/unified_metadata.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

# Resolve project root: env var override or default absolute path
PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
BASE <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts")
# HARMONIZE_OUT_DIR lets a corrected rerun write beside, not over, the live table.
OUT  <- Sys.getenv("HARMONIZE_OUT_DIR", file.path(BASE, "analysis/integration/metadata"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- Load dataset config ---
cfg_path <- file.path(PROJECT_ROOT, "config/human_datasets.yaml")
if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)
cfg <- yaml.load_file(cfg_path)
datasets_cfg <- cfg$datasets

cat("Datasets in config:", paste(names(datasets_cfg), collapse = ", "), "\n\n")

# --- Load SraRunTables from config-specified paths ---
sra <- lapply(names(datasets_cfg), function(ds) {
  path <- file.path(BASE, datasets_cfg[[ds]]$sra_table_path)
  if (!file.exists(path)) {
    warning("SRA table not found for ", ds, ": ", path, " — skipping")
    return(NULL)
  }
  fread(path)
})
names(sra) <- names(datasets_cfg)
sra <- sra[!sapply(sra, is.null)]

# --- BG-002 fail-closed gate (2026-08-07) ---------------------------------
# Datasets flagged include_in_mega define the canonical pooled estimand. A
# missing SraRunTable, an absent harmonize() case, or a harmonize() error used
# to warn and silently drop the cohort, letting a partial substrate through
# with a valid-looking unified_metadata.csv. Non-mega datasets keep the
# permissive behaviour; mega datasets now fail closed.
mega_datasets <- names(datasets_cfg)[vapply(
  names(datasets_cfg),
  function(ds) isTRUE(datasets_cfg[[ds]]$de$include_in_mega),
  logical(1)
)]
if (!length(mega_datasets)) {
  stop("No include_in_mega datasets found in config/human_datasets.yaml")
}
cat("include_in_mega datasets (", length(mega_datasets), "):",
    paste(mega_datasets, collapse = ", "), "\n\n")

missing_sra <- setdiff(mega_datasets, names(sra))
if (length(missing_sra)) {
  stop("Missing SraRunTable for include_in_mega dataset(s): ",
       paste(missing_sra, collapse = ", "),
       "\n  These define the canonical pooled analysis and may not be dropped.",
       "\n  Fix the sra_table_path in config/human_datasets.yaml or restore the file.")
}

# --- Apply excluded samples from config ---
for (ds in names(sra)) {
  excl <- datasets_cfg[[ds]]$excluded_samples
  if (length(excl) > 0) {
    before <- nrow(sra[[ds]])
    sra[[ds]] <- sra[[ds]][!Run %in% excl]
    cat("Excluded", before - nrow(sra[[ds]]), "samples from", ds,
        "(config: excluded_samples)\n")
  }
}

# ============================================================================
# harmonize()
# Maps dataset-specific SRA columns to the unified schema.
#
# Unified columns: sample_id, dataset, condition, group_binary,
#                  sex, age, fibrosis_stage
#
# ADD NEW DATASET: add a case to this switch matching the dataset ID in config.
# ============================================================================
harmonize <- function(dt, dataset) {
  switch(dataset,

    # ---- GSE126848 ----
    "GSE126848" = {
      dt[, .(
        sample_id     = Run,
        dataset       = "GSE126848",
        condition     = fcase(
          disease == "healthy", "Control",
          disease == "obese",   "Control_Obese",
          disease == "NAFLD",   "NAFL",
          disease == "NASH",    "NASH"
        ),
        group_binary  = fifelse(disease %in% c("healthy", "obese"), "Control", "Disease"),
        sex           = fifelse(gender == "female", "F", "M"),
        age           = NA_real_,
        fibrosis_stage = NA_integer_,
        nas_score      = NA_integer_
      )]
    },

    # ---- GSE167523 ----
    "GSE167523" = {
      dt[, .(
        sample_id     = Run,
        dataset       = "GSE167523",
        condition     = fifelse(disease_subtype == "NAFL", "NAFL", "NASH"),
        group_binary  = "Disease",   # No healthy controls
        sex           = fifelse(gender == "female", "F", "M"),
        age           = as.numeric(AGE),
        fibrosis_stage = NA_integer_,
        nas_score      = NA_integer_
      )]
    },

    # ---- GSE135251 ----
    "GSE135251" = {
      dt[, .(
        sample_id     = Run,
        dataset       = "GSE135251",
        condition     = fcase(
          group_in_paper == "control",    "Control",
          group_in_paper == "NAFL",       "NAFL",
          group_in_paper == "NASH_F0-F1", "NASH",
          group_in_paper == "NASH_F2",    "NASH",
          group_in_paper == "NASH_F3",    "NASH_Fibrosis",
          group_in_paper == "NASH_F4",    "NASH_Fibrosis"
        ),
        group_binary  = fifelse(disease == "Control", "Control", "Disease"),
        sex           = NA_character_,
        age           = NA_real_,
        fibrosis_stage = as.integer(Fibrosis_stage),
        nas_score      = as.integer(nas_score)
      )]
    },

    # ---- GSE130970 ----
    # Controls are the six histologically normal biopsies of Hoang et al. 2019,
    # reconstructed from their Supplementary Table 1 by
    # 00a_gse130970_source_controls.py. Fibrosis stage 0 is NOT control status:
    # 16 of the F0 biopsies are steatotic NAFLD (review item 1, 2026-09-23).
    "GSE130970" = {
      src_path <- file.path(BASE, "metadata/source/GSE130970/GSE130970_source_diagnosis.tsv")
      if (!file.exists(src_path)) stop("GSE130970 source diagnosis table missing: ", src_path)
      src <- fread(src_path, colClasses = "character")
      unresolved <- setdiff(dt$Run, src$run)
      if (length(unresolved)) {
        stop("GSE130970 runs without a source control status: ",
             paste(unresolved, collapse = ", "))
      }
      dt <- merge(dt, src[, .(Run = run, source_control_status, strict_control_nas0)],
                  by = "Run", all.x = TRUE, sort = FALSE)
      if (dt[source_control_status == "Control" & steatosis_grade != 0, .N]) {
        stop("A steatotic GSE130970 biopsy is labelled as a source control")
      }
      dt[, .(
        sample_id     = Run,
        dataset       = "GSE130970",
        condition     = fcase(
          source_control_status == "Control", "Control",
          fibrosis_stage == "0", "Fibrosis_F0",
          fibrosis_stage == "1", "Fibrosis_F1",
          fibrosis_stage %in% c("2", "3"), "Fibrosis_F2F3",
          fibrosis_stage == "4", "Fibrosis_F4"
        ),
        group_binary  = fifelse(source_control_status == "Control", "Control", "Disease"),
        sex           = fifelse(sex == "female", "F", "M"),
        age           = as.numeric(age_at_biopsy),
        fibrosis_stage = as.integer(fibrosis_stage),
        nas_score      = as.integer(nafld_activity_score),
        source_control_status = source_control_status,
        strict_control_nas0   = strict_control_nas0,
        steatosis_grade       = as.integer(steatosis_grade),
        lobular_grade         = as.integer(lobular_inflammation_grade),
        ballooning_grade      = as.integer(cytological_ballooning_grade)
      )]
    },

    # ---- GSE213621 ----
    # Fibrosis-based labeling only (no NAFL/NASH diagnosis).
    "GSE213621" = {
      dt[, .(
        sample_id      = Run,
        dataset        = "GSE213621",
        condition      = fcase(
          fibrotic_stage == "Control", "Control",
          fibrotic_stage == "F0F1",    "Fibrosis_F0F1",
          fibrotic_stage == "F2",      "Fibrosis_F2",
          fibrotic_stage == "F3F4",    "Fibrosis_F3F4"
        ),
        group_binary   = fifelse(fibrotic_stage == "Control", "Control", "Disease"),
        sex            = NA_character_,
        age            = NA_real_,
        fibrosis_stage = fcase(
          fibrotic_stage == "Control", 0L,
          fibrotic_stage == "F0F1",    1L,
          fibrotic_stage == "F2",      2L,
          fibrotic_stage == "F3F4",    3L
        ),
        nas_score      = NA_integer_
      )]
    },

    # ---- GSE162694 ----
    # Bril et al. 2021; bariatric surgery NASH cohort; normal histology = controls
    "GSE162694" = {
      dt[, .(
        sample_id      = Run,
        dataset        = "GSE162694",
        condition      = fcase(
          fibrosis_stage == "normal liver histology", "Control",
          fibrosis_stage == "0", "NASH_F0",
          fibrosis_stage == "1", "NASH_F1",
          fibrosis_stage == "2", "NASH_F2",
          fibrosis_stage == "3", "NASH_F3",
          fibrosis_stage == "4", "NASH_F4"
        ),
        group_binary   = fifelse(fibrosis_stage == "normal liver histology", "Control", "Disease"),
        sex            = fifelse(sex == "female", "F", "M"),
        age            = as.numeric(age),
        fibrosis_stage = fcase(
          fibrosis_stage == "normal liver histology", NA_integer_,
          default = as.integer(fibrosis_stage)
        ),
        nas_score      = as.integer(nas_score)
      )]
    },

    # ---- GSE174478 ----
    # Japanese NAFLD cohort; all NAFLD, no healthy controls; fibrosis 0-4
    "GSE174478" = {
      dt[, .(
        sample_id      = Run,
        dataset        = "GSE174478",
        condition      = fcase(
          fibrosis_stage == "0", "Fibrosis_F0",
          fibrosis_stage == "1", "Fibrosis_F1",
          fibrosis_stage == "2", "Fibrosis_F2",
          fibrosis_stage == "3", "Fibrosis_F3",
          fibrosis_stage == "4", "Fibrosis_F4"
        ),
        group_binary   = "Disease",   # No healthy controls
        sex            = fifelse(sex == "female", "F", "M"),
        age            = as.numeric(age),
        fibrosis_stage = as.integer(fibrosis_stage),
        nas_score      = as.integer(nas_score)
      )]
    },

    # ---- GSE193066 ----
    # Sex is in Sex_y column (Sex_x is empty); fibrosis stage is numeric 0-4
    "GSE193066" = {
      dt[, .(
        sample_id      = Run,
        dataset        = "GSE193066",
        condition      = fcase(
          get("fibrosis stage") == "0" | get("fibrosis stage") == 0, "Fibrosis_F0",
          get("fibrosis stage") == "1" | get("fibrosis stage") == 1, "Fibrosis_F1",
          get("fibrosis stage") == "2" | get("fibrosis stage") == 2, "Fibrosis_F2",
          get("fibrosis stage") == "3" | get("fibrosis stage") == 3, "Fibrosis_F3",
          get("fibrosis stage") == "4" | get("fibrosis stage") == 4, "Fibrosis_F4"
        ),
        group_binary   = "Disease",
        sex            = fifelse(Sex_y == "female", "F", "M"),
        age            = as.numeric(age),
        fibrosis_stage = as.integer(get("fibrosis stage")),
        nas_score      = as.integer(get("nafld activity score"))
      )]
    },

    # ---- GSE240729 ----
    # No sex or age metadata; sex will be inferred from XIST/DDX3Y in script 01
    "GSE240729" = {
      dt[, .(
        sample_id      = Run,
        dataset        = "GSE240729",
        condition      = fcase(
          fibrosisscore == "F0", "Fibrosis_F0",
          fibrosisscore == "F1", "Fibrosis_F1",
          fibrosisscore == "F2", "Fibrosis_F2",
          fibrosisscore == "F3", "Fibrosis_F3",
          fibrosisscore == "F4", "Fibrosis_F4"
        ),
        group_binary   = "Disease",
        sex            = NA_character_,
        age            = NA_real_,
        fibrosis_stage = as.integer(gsub("F", "", fibrosisscore)),
        nas_score      = NA_integer_
      )]
    },

    # ---- FALLBACK — warns loudly for any dataset without a harmonize() case ----
    {
      warning("No harmonize() case for dataset: ", dataset,
              "\nAdd a case to 00_harmonize_metadata.R switch() block.")
      NULL
    }
  )
}

# --- Build unified metadata ---
harmonized <- mapply(
  function(dt, ds) {
    if (is.null(dt)) return(NULL)
    result <- tryCatch(harmonize(dt, ds), error = function(e) {
      warning("harmonize() failed for ", ds, ": ", conditionMessage(e))
      NULL
    })
    result
  },
  sra, names(sra),
  SIMPLIFY = FALSE
)
harmonized <- harmonized[!sapply(harmonized, is.null)]

if (length(harmonized) == 0) stop("No datasets harmonized successfully.")

# --- BG-002 fail-closed gate (2026-08-07) ---------------------------------
# A missing switch() case or a harmonize() error above only warns. For
# include_in_mega cohorts that would silently shrink the canonical estimand,
# so refuse rather than emit a partial substrate.
dropped_mega <- setdiff(mega_datasets, names(harmonized))
if (length(dropped_mega)) {
  stop("harmonize() produced no rows for include_in_mega dataset(s): ",
       paste(dropped_mega, collapse = ", "),
       "\n  Check for a missing switch() case in harmonize() or an error warned above.",
       "\n  Canonical pooled cohorts may not be silently dropped.")
}

unified <- rbindlist(harmonized, use.names = TRUE, fill = TRUE)

# --- Compute diagnosis_harmonized column ---
# NAS-based classification (Kleiner et al. 2005):
#   NAS < 3  → NAFL
#   NAS 3-4  → Borderline (grouped with NASH for analysis)
#   NAS >= 5 → NASH
#   NAS = 0 + Control condition → Control
# For datasets without NAS scores, use original condition labels where possible.
unified[, diagnosis_harmonized := fcase(
  # 0. A source-defined histologically normal control stays Control even when
  #    its GEO NAS is 1 (two GSE130970 controls; see 00a_gse130970_source_controls.py).
  !is.na(source_control_status) & source_control_status == "Control", "Control",

  # 1. Datasets with NAS scores: classify by NAS
  !is.na(nas_score) & nas_score == 0 & condition == "Control", "Control",
  !is.na(nas_score) & nas_score < 3,                          "NAFL",
  !is.na(nas_score) & nas_score %in% 3:4,                     "Borderline",
  !is.na(nas_score) & nas_score >= 5,                         "NASH",

  # 2. GSE126848 / GSE167523: have original NAFL/NASH labels (no NAS)
  condition == "NAFL",                                         "NAFL",
  condition == "NASH",                                         "NASH",

  # 3. GSE135251 NASH_F3/F4: NASH_Fibrosis condition → NASH in harmonized view
  condition == "NASH_Fibrosis",                                "NASH",

  # 4. Controls without NAS
  condition == "Control" | condition == "Control_Obese",       "Control",

  # 5. Fibrosis-only datasets (GSE213621, GSE240729): cannot classify
  # 6. GSE162694 NASH_F* samples that lack a NAS score (26 samples; fibrosis dist
  #    F0=9, F1=5, F2=1, F3=3, F4=8) fall through to NA here and are therefore
  #    EXCLUDED from the NAFL-vs-NASH contrast (Script 13) BY DESIGN: NAFL-vs-NASH
  #    classification is NAS-based (Kleiner 2005) and cannot be assigned without NAS.
  #    117/143 GSE162694 samples are NAS-classifiable (matches the documented count in
  #    docs/dataset_labeling_and_harmonization.md); the 26 NA-NAS are listed in the
  #    "Excluded from NAFL/NASH" tally printed below. (Documented per audit C4-1, 2026-06-09.)
  default = NA_character_
)]

cat("\n========== DIAGNOSIS HARMONIZED ==========\n")
cat("NAS-based reclassification:\n")
print(unified[!is.na(nas_score), .N, by = .(dataset, diagnosis_harmonized)][order(dataset, diagnosis_harmonized)])
cat("\nLabel-based (no NAS):\n")
print(unified[is.na(nas_score) & !is.na(diagnosis_harmonized), .N,
              by = .(dataset, diagnosis_harmonized)][order(dataset, diagnosis_harmonized)])
cat("\nExcluded from NAFL/NASH (no NAS, no diagnosis label):\n")
print(unified[is.na(diagnosis_harmonized), .N, by = dataset])
cat("\nOverall diagnosis_harmonized:\n")
print(unified[, .N, by = diagnosis_harmonized][order(-N)])

# --- Summary ---
cat("\n========== UNIFIED METADATA SUMMARY ==========\n")
cat("Total samples:", nrow(unified), "\n\n")
cat("By dataset:\n");       print(unified[, .N, by = dataset][order(-N)])
cat("\nBy condition:\n");    print(unified[, .N, by = condition][order(-N)])
cat("\nBy group_binary:\n"); print(unified[, .N, by = group_binary])
cat("\nSex available:\n")
print(unified[, .(has_sex = sum(!is.na(sex)), total = .N), by = dataset])
cat("\nFibrosis stage available:\n")
print(unified[, .(has_fib = sum(!is.na(fibrosis_stage)), total = .N), by = dataset])
cat("\nNAS score available:\n")
print(unified[, .(has_nas = sum(!is.na(nas_score)), total = .N), by = dataset])

# --- Write output ---
# BG-002 (2026-08-07): every include_in_mega cohort must still be represented
# in the emitted table, and the write is atomic so a failure cannot leave a
# truncated unified_metadata.csv in place of the previous valid one.
absent_mega <- setdiff(mega_datasets, unique(as.character(unified$dataset)))
if (length(absent_mega)) {
  stop("Unified metadata lacks rows for include_in_mega dataset(s): ",
       paste(absent_mega, collapse = ", "))
}

outfile <- file.path(OUT, "unified_metadata.csv")
tmpfile <- file.path(OUT, paste0(".unified_metadata.csv.", Sys.getpid(), ".tmp"))
fwrite(unified, tmpfile)
if (!file.rename(tmpfile, outfile)) {
  unlink(tmpfile)
  stop("Atomic rename of unified metadata failed: ", tmpfile, " -> ", outfile)
}
cat("\nWritten:", outfile, "\n")
cat("Samples:", nrow(unified), " | Columns:", ncol(unified), "\n")
