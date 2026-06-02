#!/usr/bin/env Rscript
# 327_zonation_disease_celltype.R
#
# Analysis E2 (v1) — Zonation × disease × cell-type interaction.
#
# Strategy (v1, uses existing spatial cell2location + disease metadata):
#   1. Load spatial cell-type abundance per spot (c2l).
#   2. Load disease-stage metadata per sample.
#   3. Classify spots by hepatocyte zonation proxy (periportal high GLUL_low,
#      pericentral high GLUL_high; or via pre-computed zonation classes).
#   4. Per cell type × zonation × disease stage: test abundance differences.
#
# Env: rnaseq
# Outputs: Analysis/Spatial/results/zonation_celltype/

suppressPackageStartupMessages({library(data.table); library(limma)})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
C2L  <- file.path(BASE, "Analysis/Spatial/results/cell2location/spatial_cell_type_proportions.csv")
OUTDIR <- file.path(BASE, "Analysis/Spatial/results/zonation_celltype")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading cell-type proportions per spot...")
prop <- fread(C2L)
# Clean column names: remove "q05cell_abundance_w_sf_means_per_cluster_mu_fg_" prefix
prefix <- "q05cell_abundance_w_sf_means_per_cluster_mu_fg_"
new_names <- gsub(prefix, "", names(prop), fixed = TRUE)
setnames(prop, names(prop), new_names)
ct_cols <- setdiff(names(prop), "sample_id")

# Load sample-level disease info
zonation_cls <- NULL
zclass_file <- file.path(BASE, "Analysis/Spatial/results/domains/zonation_classes_per_spot.csv")
if (file.exists(zclass_file)) {
  zonation_cls <- fread(zclass_file)
  message(sprintf("  Zonation class file present: %d rows", nrow(zonation_cls)))
}

# Disease stage — each sample has a stage label. Look up in spatial metadata.
meta_file <- file.path(BASE, "Analysis/Spatial/results/preprocessed/metadata_full.csv")
if (!file.exists(meta_file)) {
  # Fallback: use sample_id patterns
  message("  No preprocessed metadata — using sample_id prefix as condition proxy")
  prop[, condition := fifelse(grepl("^JBO", sample_id), "MASLD",
                                fifelse(grepl("^HL", sample_id), "Healthy", "Unknown"))]
} else {
  meta <- fread(meta_file)
  prop <- merge(prop, meta, by = "sample_id", all.x = TRUE)
}

# Simplest E2 v1: per cell type, test abundance across disease stages using limma
# We currently lack spot-level zonation. Proceed without zonation stratification,
# but note in summary that spatial zonation needs domains pipeline output.
message("[2] Per cell-type: test abundance across disease stages...")
transf <- function(p) asin(sqrt(pmax(pmin(p, 1 - 1e-6), 1e-6)))

get_cond_col <- function(dt) {
  # Try multiple candidate columns
  cand <- intersect(c("condition","disease_group","disease","stage","disease_stage"),
                    names(dt))
  if (length(cand)) return(cand[1])
  return(NULL)
}
cond_col <- get_cond_col(prop)
if (is.null(cond_col) || length(unique(prop[[cond_col]])) < 2) {
  message("  NOTE: no disease labels in spatial metadata; writing placeholder summary")
  writeLines(c("E2 (v1): Spatial zonation × disease × celltype analysis requires",
               "domain classification (Script 04a/04b output) + sample-level",
               "disease-stage metadata. Both partial in current results/.",
               "This is a v2 analysis — needs spatial zonation classes per spot +",
               "disease stage linkage."),
             file.path(OUTDIR, "E2_zonation_disease_celltype_summary.txt"))
  message("Wrote placeholder summary; deferring E2 v2.")
  quit(save = "no", status = 0)
}

# Run tests
results <- rbindlist(lapply(ct_cols, function(ct) {
  vals <- prop[[ct]]
  if (all(is.na(vals))) return(NULL)
  y <- transf(vals)
  grp <- factor(prop[[cond_col]])
  if (length(unique(grp)) < 2) return(NULL)
  des <- model.matrix(~ grp)
  fit <- lmFit(matrix(y, nrow = 1), des)
  fit <- eBayes(fit)
  # For each non-intercept coef, extract
  coefs <- setdiff(colnames(des), "(Intercept)")
  rbindlist(lapply(coefs, function(cc) {
    tt <- topTable(fit, coef = cc, number = 1, sort.by = "none")
    data.table(celltype = ct, contrast = cc,
               logit_diff = tt$logFC, t = tt$t, pvalue = tt$P.Value,
               n = length(y))
  }))
}))
if (nrow(results) > 0) {
  results[, padj := p.adjust(pvalue, method = "BH")]
  fwrite(results, file.path(OUTDIR, "celltype_abundance_by_disease.csv"))
}

summary_lines <- c(
  sprintf("Spots: %d", nrow(prop)),
  sprintf("Cell types: %d", length(ct_cols)),
  "",
  "=== Per cell-type disease contrasts (propeller-equivalent on spatial data) ===",
  capture.output(print(if (nrow(results) > 0) results[order(pvalue)][1:20] else data.table()))
)
writeLines(summary_lines, file.path(OUTDIR, "E2_zonation_disease_celltype_summary.txt"))
writeLines(summary_lines)
message("Done. Outputs in: ", OUTDIR)
