args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript 06_prepare_bulk_counts.R <counts_path> <metadata_path> <species> <output_dir> <dataset_name>")
}

# Prefer conda env libraries over user libs if provided
lib_sites <- c(Sys.getenv("R_LIBS_SITE"), Sys.getenv("R_LIBS"))
lib_sites <- lib_sites[lib_sites != ""]
if (length(lib_sites) > 0) {
  .libPaths(unique(lib_sites))
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
species <- args[[3]]
output_dir <- args[[4]]
dataset_name <- args[[5]]

if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  stop("AnnotationDbi is required")
}

if (species == "mouse") {
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop("org.Mm.eg.db is required for mouse mapping")
  }
  orgdb <- org.Mm.eg.db::org.Mm.eg.db
  keytype <- "ENSEMBL"
} else if (species == "human") {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db is required for human mapping")
  }
  orgdb <- org.Hs.eg.db::org.Hs.eg.db
  keytype <- "ENSEMBL"
} else {
  stop("species must be mouse or human")
}

# Load counts
counts_raw <- read.delim(counts_path, comment.char = "#", check.names = FALSE)

# If featureCounts format, drop metadata columns
if (all(c("Geneid", "Chr", "Start", "End", "Strand", "Length") %in% colnames(counts_raw))) {
  gene_ids <- counts_raw$Geneid
  counts_mat <- counts_raw[, !(colnames(counts_raw) %in% c("Geneid", "Chr", "Start", "End", "Strand", "Length")), drop = FALSE]

  # Clean sample names from bam paths
  cleaned <- sub(".*/", "", colnames(counts_mat))
  cleaned <- sub("\\.Aligned.*", "", cleaned)
  colnames(counts_mat) <- cleaned
} else if ("GeneID" %in% colnames(counts_raw)) {
  gene_ids <- counts_raw$GeneID
  counts_mat <- counts_raw[, setdiff(colnames(counts_raw), "GeneID"), drop = FALSE]
} else {
  stop("Unknown counts format: missing Geneid or GeneID")
}

# Strip version suffix
gene_ids_clean <- sub("\\..*$", "", gene_ids)

# Map to symbols
mapping <- AnnotationDbi::select(
  orgdb,
  keys = unique(gene_ids_clean),
  keytype = keytype,
  columns = "SYMBOL"
)

mapping <- mapping[!is.na(mapping$SYMBOL), ]

# Merge counts with mapping
counts_df <- data.frame(ENSEMBL = gene_ids_clean, counts_mat, check.names = FALSE)
merged <- merge(mapping, counts_df, by = "ENSEMBL")

# Aggregate by SYMBOL
agg <- aggregate(. ~ SYMBOL, data = merged[, c("SYMBOL", colnames(counts_mat))], FUN = sum)

# Write counts
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
counts_out <- file.path(output_dir, paste0(dataset_name, "_counts.tsv"))
write.table(agg, counts_out, sep = "\t", quote = FALSE, row.names = FALSE)

# Prepare metadata
meta <- read.delim(metadata_path, check.names = FALSE)
if ("sample_id" %in% colnames(meta)) {
  meta$sample_id <- as.character(meta$sample_id)
  rownames(meta) <- meta$sample_id
} else if ("sample" %in% colnames(meta)) {
  meta$sample <- as.character(meta$sample)
  rownames(meta) <- meta$sample
} else {
  stop("Metadata must include sample_id or sample column")
}

# Keep only samples present in counts
sample_cols <- setdiff(colnames(agg), "SYMBOL")
meta <- meta[intersect(rownames(meta), sample_cols), , drop = FALSE]

metadata_out <- file.path(output_dir, paste0(dataset_name, "_metadata.tsv"))
write.table(meta, metadata_out, sep = "\t", quote = FALSE, row.names = TRUE)

cat("Wrote", counts_out, "and", metadata_out, "\n")
