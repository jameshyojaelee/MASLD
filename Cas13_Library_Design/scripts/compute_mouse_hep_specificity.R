#!/usr/bin/env Rscript
# compute_mouse_hep_specificity.R
# ---------------------------------------------------------------------------
# MOUSE counterpart of compute_hep_specificity.R: cache a per-gene hepatocyte-
# substrate specificity table in MOUSE for the Cas13 v8 library expression GATE.
#
# The in-vivo screen reads out CELL-AUTONOMOUS HEPATOCYTE lipid content in MOUSE,
# so a protein-coding target can only score if its MOUSE transcript is present in
# MOUSE hepatocytes (Cas13 substrate) -- not merely a whole-liver-bulk passenger
# from another lineage. v7 gated on whole-liver bulk MCD TPM, which (a) retains
# non-parenchymal genes that cannot score (e.g. Col1a1, a stellate gene) and
# (b) is a single thin n=6 experiment. v8 gates on this mouse-hepatocyte table.
#
# Per gene (mirror of the human method exactly):
#   mouse_hep_cpm          mean CPM across hepatocyte sample pseudobulk
#   mouse_max_other_cpm    max mean-CPM over all OTHER mouse lineages
#   mouse_hep_ratio        mouse_hep_cpm / mouse_max_other_cpm  (specificity)
#   mouse_hep_substrate    high            (hep>=1 CPM AND hep-specific ratio>=1)
#                          ambient_suspect (hep>=1 CPM but ratio<1 -> likely ambient)
#                          absent          (hep<1 CPM -> not a usable Cas13 substrate)
#
# Input : Analysis/SingleCell/results_gpu_v2/mouse_sc/pseudobulk/*_pseudobulk.csv
#         (gene x sample RAW counts, one file per mouse lineage; ~114 samples
#          pooled across GSE189600 ALIOS + the Liver Cell Atlas mouse). Keyed on
#          MOUSE SYMBOL (column 1 == "gene").
# Output: Cas13_Library_Design/data/mouse_hep_specificity.csv  (carries BOTH
#         gene_id [Ensembl, via GRCm39-2024-A] and gene_symbol; downstream joins
#         on gene_id to survive symbol drift)
#
# Run on a COMPUTE NODE, e.g.:
#   srun --partition=io --qos=interactive --mem=16G --cpus-per-task=2 \
#        --time=4:00:00 micromamba run -n rnaseq Rscript compute_mouse_hep_specificity.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PBDIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mouse_sc/pseudobulk")
OUT   <- file.path(BASE, "Cas13_Library_Design/data/mouse_hep_specificity.csv")

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
                  mouse_hep_cpm = round(hep, 3),
                  mouse_max_other_cpm = round(maxother, 3),
                  mouse_hep_ratio = round(ratio, 3),
                  mouse_top_other_celltype = topother,
                  mouse_hep_substrate = sub)

# ── gene_id annotation (2026-06-24) ─────────────────────────────────────────
# Map mouse SYMBOL -> Ensembl gene_id through the SAME reference the scRNA atlas
# (hence this pseudobulk) was built on: Cell Ranger refdata-gex-GRCm39-2024-A
# (GENCODE vM37). The pseudobulk `gene` column IS this reference's gene_name, so
# the symbol->gene_id map is a clean 1:1. Downstream (rebuild_cas13_library.R)
# then joins the gate on the STABLE gene_id, not the drift-prone symbol -- which
# fixes false-zeros where the library roster uses a newer Ensembl symbol than the
# atlas (e.g. Akr1b1<->Akr1b3, Atp5f1d<->Atp5d, Lars1<->Lars, Fyb1<->Fyb).
REFGTF <- Sys.getenv("MOUSE_REF_GTF",
  "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz")
gl  <- fread(cmd = paste0("zcat -f ", shQuote(REFGTF), " | awk -F\"\\t\" '$3==\"gene\"'"),
             header = FALSE, sep = "\t", quote = "")
ext <- function(col, key) { m <- regmatches(col, regexpr(paste0(key, ' "[^"]+"'), col))
                            sub(paste0(key, ' "'), "", sub('"$', "", m)) }
refmap <- data.table(gene_id = sub("[.][0-9]+$", "", ext(gl$V9, "gene_id")),
                     gene_symbol = ext(gl$V9, "gene_name"))
refmap <- unique(refmap, by = "gene_symbol")     # collapse rare duplicate symbols (keep first)
out[refmap, gene_id := i.gene_id, on = "gene_symbol"]
setcolorder(out, c("gene_id", "gene_symbol"))
cat(sprintf("gene_id mapped for %d / %d substrate genes (GRCm39-2024-A)\n",
            sum(!is.na(out$gene_id)), nrow(out)))

fwrite(out, OUT)
cat("wrote", OUT, "-", nrow(out), "genes\n")
print(table(out$mouse_hep_substrate))

# Sanity panel (canonical hepatocyte genes high; Hkdc1 absent; Col1a1 ambient).
chk <- c("Gck", "Alb", "Scd1", "Fasn", "Gpam", "Hsd17b13", "Pnpla3",
         "Hkdc1", "Hk1", "Hk2", "Col1a1")
cat("\nsanity:\n")
print(out[gene_symbol %in% chk][match(chk, gene_symbol), .(gene_symbol,
      mouse_hep_cpm, mouse_max_other_cpm, mouse_hep_ratio,
      mouse_top_other_celltype, mouse_hep_substrate)])
