#!/usr/bin/env Rscript
# compute_hep_specificity.R
# ---------------------------------------------------------------------------
# Cache a per-gene hepatocyte-substrate specificity table for the Cas13 v5
# library (readout-aware soft tag). The in-vivo screen reads out CELL-AUTONOMOUS
# HEPATOCYTE lipid content, so a gene can only score if its transcript is present
# in hepatocytes (Cas13 substrate) and not merely a contamination passenger from
# another lineage. We derive, per gene:
#   hep_mean_cpm        mean CPM across hepatocyte donor pseudobulk
#   max_other_cpm       max mean-CPM over all OTHER lineages
#   hep_ratio           hep_mean_cpm / max_other_cpm  (specificity)
#   hep_substrate       high            (hep>=1 CPM AND hep-specific ratio>=1)
#                       ambient_suspect (hep>=1 CPM but ratio<1 -> likely ambient)
#                       absent          (hep<1 CPM -> not a usable Cas13 substrate)
#
# Input : Analysis/SingleCell/results_gpu_v2/pseudobulk/*_pseudobulk.csv
#         (gene x donor RAW counts, one file per cell type)
# Output: Cas13_Library_Design/data/hep_specificity.csv  (keyed on human symbol)
#
# Run on a COMPUTE NODE (reads ~200 MB across 17 files), e.g.:
#   srun --partition=io --qos=interactive --mem=16G --cpus-per-task=2 \
#        --time=4:00:00 micromamba run -n rnaseq Rscript compute_hep_specificity.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PBDIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk")
OUT   <- file.path(BASE, "Cas13_Library_Design/data/hep_specificity.csv")

# Canonical lineage files (dedup the space-named duplicates: prefer underscore /
# Mar-17 versions of "T cells" and "Mono+mono derived cells").
files <- c(
  Hepatocytes        = "Hepatocytes_pseudobulk.csv",
  Cholangiocytes     = "Cholangiocytes_pseudobulk.csv",
  Endothelial        = "Endothelial_cells_pseudobulk.csv",
  Fibroblasts        = "Fibroblasts_pseudobulk.csv",
  Macrophages        = "Macrophages_pseudobulk.csv",
  B_cells            = "B_cells_pseudobulk.csv",
  T_cells            = "T_cells_pseudobulk.csv",
  Plasma_cells       = "Plasma_cells_pseudobulk.csv",
  Basophils          = "Basophils_pseudobulk.csv",
  cDC1s              = "cDC1s_pseudobulk.csv",
  cDC2s              = "cDC2s_pseudobulk.csv",
  Mig_cDCs           = "Mig.cDCs_pseudobulk.csv",
  Mono_derived       = "Mono+mono_derived_cells_pseudobulk.csv",
  Neutrophils        = "Neutrophils_pseudobulk.csv",
  pDCs               = "pDCs_pseudobulk.csv",
  Circulating_NK_NKT = "Circulating_NK_NKT_pseudobulk.csv",
  Resident_NK        = "Resident_NK_pseudobulk.csv"
)

# mean CPM per gene for one lineage (process one file at a time -> low memory)
mean_cpm <- function(path) {
  dt <- fread(path)
  setnames(dt, 1, "gene")
  genes <- dt$gene
  m  <- as.matrix(dt[, -1])
  cs <- colSums(m); cs[cs == 0] <- 1
  cpm <- sweep(m, 2, cs, "/") * 1e6
  setNames(rowMeans(cpm), genes)
}

stopifnot(all(file.exists(file.path(PBDIR, files))))

hep   <- mean_cpm(file.path(PBDIR, files["Hepatocytes"]))
genes <- names(hep)
maxother <- setNames(rep(0,  length(genes)), genes)
topother <- setNames(rep(NA_character_, length(genes)), genes)

for (ct in setdiff(names(files), "Hepatocytes")) {
  v  <- mean_cpm(file.path(PBDIR, files[ct]))
  va <- v[genes]; va[is.na(va)] <- 0
  upd <- va > maxother
  topother[upd] <- ct
  maxother <- pmax(maxother, va)
}

ratio <- hep / pmax(maxother, 1e-6)
sub <- ifelse(hep >= 1 & ratio >= 1, "high",
       ifelse(hep >= 1 & ratio <  1, "ambient_suspect", "absent"))

out <- data.table(gene_symbol = genes,
                  hep_mean_cpm = round(hep, 3),
                  max_other_cpm = round(maxother, 3),
                  hep_ratio = round(ratio, 3),
                  top_other_celltype = topother,
                  hep_substrate = sub)
fwrite(out, OUT)
cat("wrote", OUT, "-", nrow(out), "genes\n")
print(table(out$hep_substrate))
