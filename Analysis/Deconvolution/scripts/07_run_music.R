args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript 07_run_music.R <sce_rds> <counts_tsv> <metadata_tsv> <output_dir> <dataset_name>")
}

# Prefer conda env libraries over user libs if provided
lib_sites <- c(Sys.getenv("R_LIBS_SITE"), Sys.getenv("R_LIBS"))
lib_sites <- lib_sites[lib_sites != ""]
if (length(lib_sites) > 0) {
  .libPaths(unique(lib_sites))
}

sce_rds <- args[[1]]
counts_tsv <- args[[2]]
metadata_tsv <- args[[3]]
output_dir <- args[[4]]
dataset_name <- args[[5]]

if (!requireNamespace("MuSiC", quietly = TRUE)) {
  stop("MuSiC is required")
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  stop("SingleCellExperiment is required")
}
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("SummarizedExperiment is required")
}

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

sce <- readRDS(sce_rds)
if (!"counts" %in% SummarizedExperiment::assayNames(sce)) {
  if ("X" %in% SummarizedExperiment::assayNames(sce)) {
    SummarizedExperiment::assay(sce, "counts") <- SummarizedExperiment::assay(sce, "X")
  } else {
    stop("Reference SCE is missing a 'counts' assay")
  }
}

counts <- read.delim(counts_tsv, check.names = FALSE)
rownames(counts) <- counts$SYMBOL
counts$SYMBOL <- NULL

meta <- read.delim(metadata_tsv, check.names = FALSE)
# metadata rownames are already set in prepare step
if (!any(rownames(meta) == "")) {
  # no-op
} else {
  stop("Metadata rownames missing")
}

# Align samples
sample_ids <- intersect(colnames(counts), rownames(meta))
counts <- counts[, sample_ids, drop = FALSE]
meta <- meta[sample_ids, , drop = FALSE]

counts_mat <- as.matrix(counts)
storage.mode(counts_mat) <- "numeric"

common_genes <- intersect(rownames(counts_mat), rownames(sce))
if (length(common_genes) == 0) {
  stop("No overlapping genes between bulk counts and sc reference")
}

counts_mat <- counts_mat[common_genes, , drop = FALSE]
sce <- sce[common_genes, , drop = FALSE]

res <- MuSiC::music_prop(
  bulk.mtx = counts_mat,
  sc.sce = sce,
  clusters = "celltype",
  samples = "sampleID",
  verbose = TRUE
)

out_dir <- output_dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write.table(res$Est.prop.weighted, file.path(out_dir, paste0(dataset_name, "_music_prop_weighted.tsv")),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(res$Est.prop.all, file.path(out_dir, paste0(dataset_name, "_music_prop_all.tsv")),
            sep = "\t", quote = FALSE, row.names = TRUE)

saveRDS(res, file.path(out_dir, paste0(dataset_name, "_music_result.rds")))

cat("MuSiC results written to", out_dir, "\n")
