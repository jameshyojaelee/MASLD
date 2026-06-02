#!/usr/bin/env Rscript
# 189_recover_nas_subscores.R
# Recover NAS sub-component scores from original publication metadata.
# Currently only GSE130970 (N=76) has subscores in modeling_metadata.csv.
# Attempt to extract from GSE135251 Govaere supplementary.
#
# SLURM: Not needed (runs in seconds on login node)
# Env: micromamba activate rnaseq

library(data.table)
library(readxl)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
STAGING_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier")
OUTDIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("=== 189: NAS Subscore Recovery ===\n")

# Load existing modeling metadata
meta <- fread(file.path(STAGING_DIR, "modeling_metadata.csv"))
cat(sprintf("Loaded %d samples from modeling_metadata.csv\n", nrow(meta)))

# --- GSE130970: Already has subscores ---
gse130970_sub <- meta[dataset == "GSE130970" & steatosis_grade != -1]
cat(sprintf("GSE130970: %d samples with existing subscores\n", nrow(gse130970_sub)))

# --- GSE135251: Try Govaere supplementary ---
govaere_supp <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE135251/metadata",
  "aba4448_data_file_s1.xlsx")

recovered <- data.table()
if (file.exists(govaere_supp)) {
  cat("Attempting GSE135251 subscore recovery from Govaere supplementary...\n")
  sheets <- excel_sheets(govaere_supp)
  cat(sprintf("  Available sheets: %s\n", paste(sheets, collapse=", ")))

  for (sh in sheets) {
    tryCatch({
      dt <- as.data.table(read_excel(govaere_supp, sheet = sh))
      cols <- tolower(names(dt))
      has_steat <- any(grepl("steatosis", cols))
      has_inflam <- any(grepl("inflam", cols))
      has_balloon <- any(grepl("balloon", cols))
      has_id <- any(grepl("patient|sample|id|geo", cols))

      if (has_steat && has_inflam && has_id) {
        cat(sprintf("  Sheet '%s': Found subscores (%d rows)\n", sh, nrow(dt)))
        cat(sprintf("  Columns: %s\n", paste(names(dt), collapse=", ")))
        recovered <- dt
      }
    }, error = function(e) NULL)
  }
} else {
  cat("GSE135251 supplementary not found\n")
}

# --- Summary ---
n_steat <- sum(meta$steatosis_grade != -1)
n_inflam <- sum(meta$lobular_inflammation_grade != -1)
n_balloon <- sum(meta$ballooning_grade != -1)

cat(sprintf("\n=== NAS Subscore Availability ===\n"))
cat(sprintf("Steatosis:    %d samples (datasets: %s)\n",
  n_steat,
  paste(unique(meta[steatosis_grade != -1, dataset]), collapse=", ")))
cat(sprintf("Inflammation: %d samples (datasets: %s)\n",
  n_inflam,
  paste(unique(meta[lobular_inflammation_grade != -1, dataset]), collapse=", ")))
cat(sprintf("Ballooning:   %d samples (datasets: %s)\n",
  n_balloon,
  paste(unique(meta[ballooning_grade != -1, dataset]), collapse=", ")))

if (nrow(recovered) > 0) {
  out_path <- file.path(OUTDIR, "govaere_subscores_raw.csv")
  fwrite(recovered, out_path)
  cat(sprintf("\nRecovered Govaere subscores written to: %s\n", out_path))
  cat("MANUAL REVIEW REQUIRED: Match sample IDs to bulk RNA-seq IDs\n")
} else {
  cat("\nNo additional subscores recovered.\n")
  cat("Proceeding with Tier 1 targets only (fibrosis, NAS composite, disease state).\n")
}

summary <- data.table(
  target = c("steatosis_grade", "lobular_inflammation_grade", "ballooning_grade",
             "fibrosis_stage", "nas_score"),
  n_samples = c(n_steat, n_inflam, n_balloon,
                sum(meta$fibrosis_stage != -1 & !is.na(meta$fibrosis_stage)),
                sum(meta$nas_score != -1 & !is.na(meta$nas_score))),
  n_datasets = c(
    length(unique(meta[steatosis_grade != -1, dataset])),
    length(unique(meta[lobular_inflammation_grade != -1, dataset])),
    length(unique(meta[ballooning_grade != -1, dataset])),
    length(unique(meta[fibrosis_stage != -1 & !is.na(fibrosis_stage), dataset])),
    length(unique(meta[nas_score != -1 & !is.na(nas_score), dataset]))
  ),
  loco_viable = c(FALSE, FALSE, FALSE, TRUE, TRUE)
)
fwrite(summary, file.path(OUTDIR, "subscore_availability.csv"))
cat(sprintf("\nSubscore availability written to: %s\n",
  file.path(OUTDIR, "subscore_availability.csv")))
cat("Done.\n")
