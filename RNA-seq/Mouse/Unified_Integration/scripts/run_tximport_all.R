#!/usr/bin/env Rscript
# tximport aggregation for all 9 existing mouse datasets (Kallisto re-quant)
# Produces: txi.rds for M01 consumption

suppressPackageStartupMessages({
  library(tximport)
  library(data.table)
})

PROJECT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
tx2gene <- fread(file.path(PROJECT, "data/reference/kallisto/tx2gene_vM38.tsv"), header = TRUE)
setnames(tx2gene, c("tx_id", "gene_id", "gene_name"))
tx2gene$tx_id   <- sub("[.][0-9]+$", "", tx2gene$tx_id)
tx2gene$gene_id <- sub("[.][0-9]+$", "", tx2gene$gene_id)
cat("tx2gene:", nrow(tx2gene), "transcripts ->", length(unique(tx2gene$gene_id)), "genes\n")

datasets <- list(
  c("InHouse_MCD",  "RNA-seq/Mouse/InHouse_MCD"),
  c("GSE156918",    "RNA-seq/Mouse/Public_MCD/GSE156918"),
  c("GSE205974",    "RNA-seq/Mouse/Public_MCD/GSE205974"),
  c("GSE159911",    "RNA-seq/Mouse/Public_Diet_Models/GSE159911"),
  c("GSE162876",    "RNA-seq/Mouse/Public_Diet_Models/GSE162876"),
  c("GSE224069",    "RNA-seq/Mouse/Public_Diet_Models/GSE224069"),
  c("GSE225616",    "RNA-seq/Mouse/Public_Diet_Models/GSE225616"),
  c("GSE263273",    "RNA-seq/Mouse/Public_Diet_Models/GSE263273"),
  c("GSE274914",    "RNA-seq/Mouse/Public_Diet_Models/GSE274914")
)

all_files <- c()
all_names <- c()
for (ds in datasets) {
  quant_dir <- file.path(PROJECT, ds[2], "quant/kallisto")
  if (!dir.exists(quant_dir)) {
    cat("SKIP:", ds[1], "\n")
    next
  }
  samples <- list.dirs(quant_dir, recursive = FALSE, full.names = FALSE)
  for (s in samples) {
    f <- file.path(quant_dir, s, "abundance.tsv")
    if (file.exists(f)) {
      all_files <- c(all_files, f)
      all_names <- c(all_names, paste0(ds[1], "__", s))
    }
  }
}
names(all_files) <- all_names
cat("Found", length(all_files), "Kallisto quant files\n")

txi <- tximport(all_files, type = "kallisto",
                tx2gene = tx2gene[, .(tx_id, gene_id)],
                countsFromAbundance = "no",
                ignoreTxVersion = TRUE)

cat("tximport done:", nrow(txi$counts), "genes x", ncol(txi$counts), "samples\n")

out_dir <- file.path(PROJECT, "RNA-seq/Mouse/Unified_Integration/counts/tximport")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(txi, file.path(out_dir, "txi.rds"))
cat("Saved:", file.path(out_dir, "txi.rds"), "\n")
cat("Counts dim:", paste(dim(txi$counts), collapse = " x "), "\n")
cat("Length dim:", paste(dim(txi$length), collapse = " x "), "\n")
