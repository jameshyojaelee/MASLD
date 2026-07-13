#!/usr/bin/env Rscript
# M00_harmonize_mouse_metadata.R
# ---------------------------------------------------------------------------
# Harmonize metadata from all mouse MASLD datasets into a unified table.
# Sources: InHouse_MCD, Public_MCD (GSE156918, GSE205974), Public_Diet_Models
# Output:  Unified_Integration/metadata/unified_mouse_metadata.csv
# ---------------------------------------------------------------------------

# ---- Seed pinning (T2.4, 2026-04-22) -----
set.seed(42)
Sys.setenv(R_PARALLEL_SEED = "42")

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

# Resolve project root: env var override or default absolute path
PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
MOUSE  <- file.path(PROJECT_ROOT, "RNA-seq/Mouse")
OUTDIR <- file.path(MOUSE, "Unified_Integration/metadata")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --- Load dataset config ---
# File paths for each dataset's samples.tsv come from config/mouse_datasets.yaml.
# TO ADD A NEW DATASET: add an entry there rather than editing this script.
cfg_path <- file.path(PROJECT_ROOT, "config/mouse_datasets.yaml")
if (file.exists(cfg_path)) {
  mouse_cfg <- yaml.load_file(cfg_path)
  excluded_conditions <- mouse_cfg$excluded_conditions
  cat("Loaded mouse config:", cfg_path, "\n")
  cat("Excluded conditions:", paste(excluded_conditions, collapse = ", "), "\n\n")
} else {
  warning("Mouse config not found: ", cfg_path, " — using defaults")
  mouse_cfg <- list(datasets = list())
  excluded_conditions <- c("High_Fat_Diet_Elafibranor", "High_Fat_Diet_Resmetirom",
                            "FLCN_KO_MCD", "FLCN_KO_Control")
}

cat("=== M00: Harmonize Mouse Metadata ===\n\n")

# --- Helper: evaluate YAML filter expression on a data.table ---
# Returns the filtered data.table.  If filter_expr is NULL or empty, returns
# the input unchanged.  The expression is evaluated with the data.table as the
# parent environment so that column names are directly accessible.
apply_yaml_filter <- function(dt, filter_expr, dataset_id) {
  if (is.null(filter_expr) || !nzchar(trimws(filter_expr))) return(dt)
  n_before <- nrow(dt)
  keep_mask <- tryCatch(
    eval(parse(text = filter_expr), envir = dt),
    error = function(e) {
      warning("Filter failed for ", dataset_id, ": ", conditionMessage(e),
              " — keeping all rows")
      rep(TRUE, nrow(dt))
    }
  )
  dt <- dt[keep_mask]
  cat("  Applied YAML filter for ", dataset_id, ": '", filter_expr, "' -> ",
      n_before, " -> ", nrow(dt), " samples\n", sep = "")
  dt
}

# --- Build per-dataset config lookup from YAML ---
ds_configs <- if (!is.null(mouse_cfg$datasets)) mouse_cfg$datasets else list()

# ============================================================
# 1) InHouse MCD — 12 paired-end samples (weeks 1-3)
# ============================================================
inhouse <- fread(file.path(MOUSE, "InHouse_MCD/metadata/samples.tsv"))
cat("InHouse raw samples:", nrow(inhouse), "\n")

# Apply YAML filter if specified for InHouse_MCD
ds_cfg <- ds_configs[["InHouse_MCD"]]
if (!is.null(ds_cfg) && !is.null(ds_cfg$filter)) {
  inhouse <- apply_yaml_filter(inhouse, ds_cfg$filter, "InHouse_MCD")
}

inhouse_meta <- inhouse[, .(
  sample_id = sample_id,
  dataset   = "InHouse_MCD",
  diet_model = "MCD",
  condition  = diet,                  # "Control" or "MCD"
  group_binary = ifelse(diet == "Control", "Control", "Disease"),
  sex        = tolower(sex),
  layout     = "PAIRED",
  batch      = "InHouse_MCD"
)]

# ============================================================
# 2) Public_MCD — GSE156918 (Cre-only, single-end)
# ============================================================
pmcd <- fread(file.path(MOUSE, "Public_MCD/metadata/samples.tsv"))
cat("Public_MCD raw samples:", nrow(pmcd), "\n")

# GSE156918: apply YAML filter (genotype == 'Cre') + original GSM prefix selection
gse156918 <- pmcd[grepl("GSM47481", sample_id)]
ds_cfg <- ds_configs[["GSE156918"]]
if (!is.null(ds_cfg) && !is.null(ds_cfg$filter)) {
  gse156918 <- apply_yaml_filter(gse156918, ds_cfg$filter, "GSE156918")
}
gse156918_meta <- gse156918[, .(
  sample_id  = sample_id,
  dataset    = "GSE156918",
  diet_model = "MCD",
  condition  = diet,
  group_binary = ifelse(diet == "Control", "Control", "Disease"),
  sex        = "unspecified",
  layout     = "SINGLE",
  batch      = "GSE156918"
)]

# GSE205974: apply YAML filter if any
gse205974 <- pmcd[grepl("GSM623625", sample_id)]
ds_cfg <- ds_configs[["GSE205974"]]
if (!is.null(ds_cfg) && !is.null(ds_cfg$filter)) {
  gse205974 <- apply_yaml_filter(gse205974, ds_cfg$filter, "GSE205974")
}
gse205974_meta <- gse205974[, .(
  sample_id  = sample_id,
  dataset    = "GSE205974",
  diet_model = "MCD",
  condition  = diet,
  group_binary = ifelse(diet == "Control", "Control", "Disease"),
  sex        = "unspecified",
  layout     = "PAIRED",
  batch      = "GSE205974"
)]

cat("  GSE156918 Cre-only:", nrow(gse156918_meta), "\n")
cat("  GSE205974:", nrow(gse205974_meta), "\n")

# ============================================================
# 3) Public_Diet_Models — GEO datasets
# ============================================================
diet <- fread(file.path(MOUSE, "Public_Diet_Models/metadata/samples.tsv"))
cat("Diet Models raw samples:", nrow(diet), "\n")

# --- Step 3a: Apply per-dataset YAML filter expressions ---
# This handles drug-arm exclusion (GSE224069 DDD86481, GSE263273 INT/OCA)
# and any future filter needs, all driven by the YAML config.
diet_datasets <- unique(diet$dataset)
for (ds_id in diet_datasets) {
  ds_cfg <- ds_configs[[ds_id]]
  if (!is.null(ds_cfg) && !is.null(ds_cfg$filter)) {
    mask <- diet$dataset == ds_id
    sub_dt <- diet[mask]
    filtered <- apply_yaml_filter(sub_dt, ds_cfg$filter, ds_id)
    # Remove old rows, add back filtered rows
    diet <- rbindlist(list(diet[!mask], filtered), use.names = TRUE)
  }
}

# --- Step 3b: Exclude conditions from the YAML excluded_conditions list ---
# This catches drug/KO conditions that weren't caught by per-dataset filters
if (length(excluded_conditions) > 0) {
  n_before <- nrow(diet)
  diet <- diet[!condition %in% excluded_conditions]
  n_excluded <- n_before - nrow(diet)
  if (n_excluded > 0) {
    cat("  Excluded", n_excluded, "samples via excluded_conditions list\n")
  }
}

# --- Step 3c: Map conditions to diet_model and group_binary ---
# GSE263273 Vehicle maps to AMLN_ob (Lepob/Obese), NOT HFD
diet[, diet_model := fcase(
  condition %in% c("MCD"),                               "MCD",
  condition %in% c("7w_HFD", "52w_HFD",
                    "High_Fat_Diet_Vehicle",
                    "GAN_diet_D09100310_Research_Diets_USA_DIO_NASH_Vehicle"),
                                                          "HFD",
  condition == "Vehicle" & dataset == "GSE263273",        "AMLN_ob",
  condition %in% c("7wks_CDAHFD"),                       "CDAHFD",
  condition %in% c("20wks_FPC"),                         "FPC",
  condition %in% c("LIDPAD"),                            "LIDPAD",
  condition %in% c("Control", "7w_LFD", "52w_LFD",
                    "7wks_LFD_Ctrl", "20wks_LFD_Ctrl",
                    "Chow_Diet_Vehicle"),                "Control",
  default = NA_character_
)]

diet[, group_binary := fcase(
  diet_model == "Control",           "Control",
  !is.na(diet_model),               "Disease",
  default = NA_character_
)]

# Report any samples that still don't map (should be zero after filter + exclusion)
excluded <- diet[is.na(diet_model)]
if (nrow(excluded) > 0) {
  cat("  Excluded (unmapped conditions):", nrow(excluded), "\n")
  cat("    Conditions:", paste(unique(excluded$condition), collapse = ", "), "\n")
}
diet_keep <- diet[!is.na(diet_model)]

diet_meta <- diet_keep[, .(
  sample_id  = sample_id,
  dataset    = dataset,
  diet_model = diet_model,
  condition  = condition,
  group_binary = group_binary,
  sex        = tolower(sex),
  layout     = ifelse(fastq_2 != "" & !is.na(fastq_2), "PAIRED", "SINGLE"),
  batch      = dataset     # Each GEO dataset is a batch
)]

# For controls in multi-comparison datasets, assign to the diet model
# of their paired disease group
# GSE162876: controls are shared between CDAHFD and FPC comparisons
# GSE274914: controls are shared between 7w and 52w HFD comparisons

cat("  Diet samples kept:", nrow(diet_meta), "\n")

# ============================================================
# 4) Combine all metadata
# ============================================================
unified <- rbindlist(list(inhouse_meta, gse156918_meta, gse205974_meta, diet_meta),
                     use.names = TRUE, fill = TRUE)

# Drop datasets excluded from YAML (commented out = no config entry)
active_datasets <- names(mouse_cfg$datasets)
n_before <- nrow(unified)
unified <- unified[dataset %in% active_datasets]
if (nrow(unified) < n_before) {
  cat("  Dropped", n_before - nrow(unified), "samples from excluded datasets\n")
}

# Normalize sex
unified[sex %in% c("unspecified", ""), sex := NA_character_]

cat("\n===== UNIFIED METADATA SUMMARY =====\n")
cat("Total samples:", nrow(unified), "\n")
cat("\nBy dataset:\n")
print(unified[, .N, by = dataset][order(-N)])
cat("\nBy diet_model:\n")
print(unified[, .N, by = diet_model][order(-N)])
cat("\nBy group_binary:\n")
print(unified[, .N, by = group_binary])
cat("\nBy sex:\n")
print(unified[, .N, by = sex])
cat("\nDisease vs Control by dataset:\n")
print(dcast(unified, dataset ~ group_binary, fun.aggregate = length))

# Save
fwrite(unified, file.path(OUTDIR, "unified_mouse_metadata.csv"))
cat("\nSaved:", file.path(OUTDIR, "unified_mouse_metadata.csv"), "\n")
