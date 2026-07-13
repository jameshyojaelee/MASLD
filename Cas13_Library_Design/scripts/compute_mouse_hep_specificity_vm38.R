#!/usr/bin/env Rscript
# compute_mouse_hep_specificity_vm38.R
# ---------------------------------------------------------------------------
# vM38 variant of compute_mouse_hep_specificity.R: cache the per-gene mouse-
# hepatocyte substrate-specificity table for the Cas13 library expression GATE,
# computed on the NEW GENCODE vM38 scRNA re-quant.
#
# Same biology as the vM37 substrate (the in-vivo screen reads out CELL-
# AUTONOMOUS HEPATOCYTE lipid in MOUSE, so a protein-coding target can only score
# if its MOUSE transcript is present in MOUSE hepatocytes), and the SAME output
# schema, so this is a drop-in replacement for mouse_hep_specificity.csv.
#
# The one change vs the vM37 script: gene_id comes DIRECTLY from the vM38
# pseudobulk (Cell Ranger features.tsv col1 = native vM38 ENSMUSG). The vM37
# version had to remap mouse SYMBOL -> Ensembl gene_id through the
# refdata-gex-GRCm39-2024-A GTF; here the pseudobulk already carries the native
# vM38 gene_id, so NO GTF symbol->id mapping is needed. (vM38 = the whole point of
# the re-quant: it recovers the ~68 predicted lncRNA loci vM37 was missing.)
#
# Per gene (mirror of the human/vM37 method exactly):
#   mouse_hep_cpm          mean CPM across hepatocyte sample pseudobulk
#   mouse_max_other_cpm    max mean-CPM over all OTHER mouse lineages
#   mouse_hep_ratio        mouse_hep_cpm / mouse_max_other_cpm  (specificity)
#   mouse_top_other_celltype  lineage with the max other-CPM
#   mouse_hep_substrate    high            (hep>=1 CPM AND hep-specific ratio>=1)
#                          ambient_suspect (hep>=1 CPM but ratio<1 -> likely ambient)
#                          absent          (hep<1 CPM -> not a usable Cas13 substrate)
#
# Input : Analysis/SingleCell/results_gpu_v2/mouse_sc/pseudobulk_vm38/*_pseudobulk.csv
#         (gene x sample RAW counts, one file per mouse lineage; emitted by
#          mouse_pseudobulk_vm38.py). Each CSV has a leading `gene_id` column
#          (native vM38 ENSMUSG) and a `gene_symbol` column, then one column per
#          sample.
# Output: Cas13_Library_Design/data/mouse_hep_specificity_vm38.csv
#         columns: gene_id, gene_symbol, mouse_hep_cpm, mouse_max_other_cpm,
#                  mouse_hep_ratio, mouse_top_other_celltype, mouse_hep_substrate
#         (identical schema to mouse_hep_specificity.csv -> drop-in for the gate).
#
# Run on a COMPUTE NODE, e.g.:
#   srun --partition=io --qos=interactive --mem=16G --cpus-per-task=2 \
#        --time=4:00:00 micromamba run -n rnaseq Rscript compute_mouse_hep_specificity_vm38.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PBDIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mouse_sc/pseudobulk_vm38")
OUT   <- file.path(BASE, "Cas13_Library_Design/data/mouse_hep_specificity_vm38.csv")

# Mouse liver lineage pseudobulk files (Hepatocytes = substrate; rest = "other").
files <- c(
  Hepatocytes    = "Hepatocytes_pseudobulk.csv",
  Cholangiocytes = "Cholangiocytes_pseudobulk.csv",
  LSECs          = "LSECs_pseudobulk.csv",
  Stellate_cells = "Stellate_cells_pseudobulk.csv",
  Kupffer_cells  = "Kupffer_cells_pseudobulk.csv",
  Monocytes      = "Monocytes_pseudobulk.csv",
  NK_cells       = "NK_cells_pseudobulk.csv",
  T_cells        = "T_cells_pseudobulk.csv",
  B_cells        = "B_cells_pseudobulk.csv",
  pDCs           = "pDCs_pseudobulk.csv"
)

# Read one lineage CSV -> data.table keyed on the native vM38 gene_id. The first
# two columns are gene_id + gene_symbol (emitted by mouse_pseudobulk_vm38.py); the
# remaining columns are per-sample raw counts.
read_pb <- function(path) {
  dt <- fread(path)
  stopifnot("gene_id"   %in% names(dt),
            "gene_symbol" %in% names(dt))
  dt
}

# mean CPM per gene_id for one lineage (process one file at a time -> low memory).
# Returns a named numeric vector keyed on gene_id.
mean_cpm <- function(dt) {
  ids <- dt$gene_id
  m   <- as.matrix(dt[, setdiff(names(dt), c("gene_id", "gene_symbol")), with = FALSE])
  storage.mode(m) <- "double"
  cs <- colSums(m); cs[cs == 0] <- 1
  cpm <- sweep(m, 2, cs, "/") * 1e6
  setNames(rowMeans(cpm), ids)
}

stopifnot(all(file.exists(file.path(PBDIR, files))))

# Hepatocyte substrate (anchor gene_id + gene_symbol order from this file).
hep_dt <- read_pb(file.path(PBDIR, files["Hepatocytes"]))
hep    <- mean_cpm(hep_dt)
ids    <- names(hep)
# gene_id -> gene_symbol map (native vM38, no GTF remap).
id2sym <- setNames(hep_dt$gene_symbol, hep_dt$gene_id)

maxother <- setNames(rep(0,  length(ids)), ids)
topother <- setNames(rep(NA_character_, length(ids)), ids)

for (ct in setdiff(names(files), "Hepatocytes")) {
  v  <- mean_cpm(read_pb(file.path(PBDIR, files[ct])))
  va <- v[ids]; va[is.na(va)] <- 0          # align to hepatocyte gene_id order
  upd <- va > maxother
  topother[upd] <- ct
  maxother <- pmax(maxother, va)
}

ratio <- hep / pmax(maxother, 1e-6)
sub <- ifelse(hep >= 1 & ratio >= 1, "high",
       ifelse(hep >= 1 & ratio <  1, "ambient_suspect", "absent"))

# Strip the Ensembl version suffix so gene_id matches the vM37 substrate's
# UNVERSIONED key (ENSMUSG.....) -- rebuild_cas13_library.R joins on the base id
# (strip_v(); gene_id_mouse) and META keys on mouse_ensembl_base, both unversioned.
# vM38 cellranger features carry .version (e.g. ENSMUSG00000051951.6); confirmed
# 1:1 base->versioned (no base-id collisions in the 55,402-feature set), so this
# is lossless. Guard anyway in case a future ref introduces a collision.
ids_base <- sub("[.][0-9]+$", "", ids)
if (anyDuplicated(ids_base))
  stop(sprintf("version-stripped gene_id collision (%d dup base ids) -- aborting",
               sum(duplicated(ids_base))))

out <- data.table(gene_id = ids_base,
                  gene_symbol = id2sym[ids],
                  mouse_hep_cpm = round(hep, 3),
                  mouse_max_other_cpm = round(maxother, 3),
                  mouse_hep_ratio = round(ratio, 3),
                  mouse_top_other_celltype = topother,
                  mouse_hep_substrate = sub)
setcolorder(out, c("gene_id", "gene_symbol", "mouse_hep_cpm",
                   "mouse_max_other_cpm", "mouse_hep_ratio",
                   "mouse_top_other_celltype", "mouse_hep_substrate"))

fwrite(out, OUT)
cat("wrote", OUT, "-", nrow(out), "genes (native vM38 gene_id)\n")
print(table(out$mouse_hep_substrate))

# Sanity panel (canonical hepatocyte genes high; Hkdc1 absent; Col1a1 ambient).
chk <- c("Gck", "Alb", "Scd1", "Fasn", "Gpam", "Hsd17b13", "Pnpla3",
         "Hkdc1", "Hk1", "Hk2", "Col1a1")
cat("\nsanity:\n")
print(out[gene_symbol %in% chk][match(chk, gene_symbol), .(gene_symbol,
      mouse_hep_cpm, mouse_max_other_cpm, mouse_hep_ratio,
      mouse_top_other_celltype, mouse_hep_substrate)])
