#!/usr/bin/env Rscript
# S4-1b: propeller (Phipson 2022) compositional retest.
#
# propeller fits an empirical-Bayes moderated beta-binomial / arcsin-sqrt
# linear model on per-donor proportions, with the donor as the replicate.
# It is the canonical complement to scCODA.
#
# Inputs (project, not worktree):
#   hepatocyte_subtype_metadata.csv  per-cell subtype + sample
#   meta_subtype_mapping.csv         raw -> meta
#   donor_metadata.tsv               sample x condition
#   subtype_disease_enrichment.csv   original Fisher p
#
# Output:
#   results_gpu_v2/proportion_compositional/propeller_results.tsv

suppressPackageStartupMessages({
  if (!requireNamespace("speckle", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("speckle", ask = FALSE, update = FALSE)
  }
  library(speckle)
  library(data.table)
  library(dplyr)
})

PROJECT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HEP   <- file.path(PROJECT, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes")
DONOR <- file.path(PROJECT, "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv")
OUT   <- file.path(PROJECT, "Analysis/SingleCell/results_gpu_v2/proportion_compositional")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

message("[load] per-cell metadata")
cell <- fread(file.path(HEP, "hepatocyte_subtype_metadata.csv"))
setnames(cell, names(cell)[1], "cell_id")

mapping <- fread(file.path(HEP, "meta_subtype_mapping.csv"))
cell[, meta_subtype := mapping$meta_subtype[match(as.character(hepatocyte_subtype), as.character(mapping$subtype))]]
cell[is.na(meta_subtype), meta_subtype := "Other"]

donor <- fread(DONOR)
cond_col <- ifelse("condition_harmonized" %in% names(donor), "condition_harmonized", "condition")
donor[, condition_binary := ifelse(tolower(get(cond_col)) %in% c("healthy","normal","control"),
                                    "Healthy", "MASLD")]
donor_lookup <- setNames(donor$condition_binary, donor$sample)
cell[, condition_binary := donor_lookup[sample]]
cell <- cell[!is.na(condition_binary)]

run_propeller <- function(group_col, claim) {
  message(sprintf("[propeller] %s  (%d cells, %d donors)", claim,
                  nrow(cell), length(unique(cell$sample))))
  clusters <- as.character(cell[[group_col]])
  samples  <- as.character(cell$sample)
  donors_meta <- unique(cell[, .(sample, condition_binary)])
  group <- donors_meta$condition_binary[match(unique(samples), donors_meta$sample)]
  res <- tryCatch(
    propeller(
      clusters = clusters,
      sample   = samples,
      group    = group,
      transform = "asin",
      robust = TRUE,
      trend = FALSE
    ),
    error = function(e) { message("  propeller error: ", e$message); NULL }
  )
  if (is.null(res)) return(NULL)
  res$claim <- claim
  res$cell_type <- rownames(res)
  data.table(res)
}

res_meta <- run_propeller("meta_subtype", "meta_subtype_Healthy_vs_MASLD")
res_raw  <- run_propeller("hepatocyte_subtype", "raw_subtype_Healthy_vs_MASLD")

# Merge with original Fisher p
orig <- fread(file.path(HEP, "subtype_disease_enrichment.csv"))
orig[, subtype := as.character(subtype)]
if (!is.null(res_raw)) {
  res_raw[, p_original_fisher := orig$fisher_pval[match(cell_type, orig$subtype)]]
}

# Donor-level variance estimate from propeller's mean-variance fit
# (extract via per-claim arcsin transform variance across donors)
donor_variance <- function(group_col) {
  prop <- cell[, .N, by = .(sample, get(group_col))]
  setnames(prop, "get", "label")
  totals <- prop[, .(total = sum(N)), by = sample]
  prop <- merge(prop, totals, by = "sample")
  prop[, p := N / total]
  prop[, asin_p := asin(sqrt(pmin(pmax(p, 0), 1)))]
  prop[, .(asin_var = var(asin_p), n_donors = .N), by = label]
}
dv_meta <- donor_variance("meta_subtype")
dv_raw  <- donor_variance("hepatocyte_subtype")

if (!is.null(res_meta)) {
  res_meta[, donor_asin_variance := dv_meta$asin_var[match(cell_type, dv_meta$label)]]
}
if (!is.null(res_raw)) {
  res_raw[, donor_asin_variance := dv_raw$asin_var[match(cell_type, dv_raw$label)]]
}

combined <- rbindlist(list(res_meta, res_raw), use.names = TRUE, fill = TRUE)
fwrite(combined, file.path(OUT, "propeller_results.tsv"), sep = "\t")
message(sprintf("[done] wrote %s (%d rows)", file.path(OUT, "propeller_results.tsv"), nrow(combined)))
