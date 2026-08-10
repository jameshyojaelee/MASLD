#!/usr/bin/env Rscript

# Generic DESeq2 Analysis Script for Snakemake Pipeline
# Updated to handle dynamic design formulas and standard inputs

library(optparse)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# Argument Parsing
option_list <- list(
  make_option(c("-c", "--counts"), type="character", help="Path to featureCounts output (gene_counts.txt)"),
  make_option(c("-m", "--metadata"), type="character", help="Path to sample sheet (metadata)"),
  make_option(c("-o", "--output"), type="character", help="Output directory"),
  make_option(c("-d", "--design"), type="character", default="~ condition", help="Design formula (e.g. '~ condition')"),
  # BG-005 (2026-08-09): a contrast is now MANDATORY. Previously the script called
  # results(dds) with no argument, which silently returns DESeq2's default -- the LAST
  # coefficient in resultsNames(dds). For any multi-level condition that is an
  # arbitrary pairwise comparison, not the registry contrast, and it looks entirely
  # plausible in the output table. Supply factor,numerator,denominator.
  make_option(c("-k", "--contrast"), type="character", default=NULL,
              help="REQUIRED. Contrast as 'factor,numerator,denominator' (e.g. 'condition,NASH,NAFL')")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$counts) || is.null(opt$metadata) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Missing required arguments")
}
if (is.null(opt$contrast)) {
  print_help(opt_parser)
  stop("--contrast is required (BG-005). Refusing to fall back to DESeq2's default\n",
       "  last-coefficient result, which silently reports an unintended comparison.\n",
       "  Note: the canonical per-study human DE is limma-voom via\n",
       "  analysis/integration/scripts/02_per_study_de.R, driven by the de: block in\n",
       "  config/human_datasets.yaml. Use that unless you specifically need DESeq2.",
       call. = FALSE)
}
contrast_parts <- trimws(strsplit(opt$contrast, ",", fixed = TRUE)[[1]])
if (length(contrast_parts) != 3L || any(!nzchar(contrast_parts))) {
  stop("--contrast must have exactly 3 comma-separated parts ",
       "'factor,numerator,denominator'; got: ", opt$contrast, call. = FALSE)
}

# Create output dir
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

# 1. Load Counts
# featureCounts output often has a header line and column names. 
# Rows: Genes. Cols: Samples (paths).
cat("Loading counts from:", opt$counts, "\n")
raw_counts <- read.table(opt$counts, header=TRUE, row.names=1, sep="\t", comment.char="#")

# featureCounts standard output has columns: Chr, Start, End, Strand, Length, [SamplePaths...]
# We need to remove the annotation columns (first 5 usually)
# We can identify sample columns by checking overlap with metadata or just taking columns 6+
cat("Processing count matrix...\n")
# Assuming standard featureCounts output where first 5 cols are metadata
# But safer to match with metadata samples
# Let's load metadata first
meta <- read.table(opt$metadata, header=TRUE, sep="\t", stringsAsFactors=FALSE)

cat("Metadata loaded. Rows:", nrow(meta), "\n")
# Expected sample ID column "sample_id" based on generate_samplesheet.py
if (!"sample_id" %in% colnames(meta)) {
    stop("Metadata must contain 'sample_id' column")
}

# Clean filenames in count matrix columns (often they are full paths)
# featureCounts often implies the filename is the column name
# We need to match `sample_id` to these columns.
cnt_cols <- colnames(raw_counts)
# Helper to simplify path to basename without extension
clean_cols <- basename(cnt_cols)
clean_cols <- gsub("_Aligned.sortedByCoord.out.bam", "", clean_cols)
clean_cols <- gsub(".bam", "", clean_cols)
# Also handle R1/R2 suffixes if present, though featureCounts usually takes the bam name

colnames(raw_counts) <- clean_cols

# Intersect samples
common_samples <- intersect(meta$sample_id, colnames(raw_counts))
if (length(common_samples) == 0) {
    stop("No common samples found between metadata sample_id and count matrix columns!")
}

cat("Matched samples:", length(common_samples), "\n")

count_matrix <- raw_counts[, common_samples, drop=FALSE]
meta_clean <- meta[meta$sample_id %in% common_samples, ]
# Ensure order matches
meta_clean <- meta_clean[match(colnames(count_matrix), meta_clean$sample_id), ]

# 2. Prepare DESeq2 Dataset
cat("Design formula:", opt$design, "\n")
design_formula <- as.formula(opt$design)

# Check if design variables exist in metadata
design_vars <- all.vars(design_formula)
for (v in design_vars) {
    if (!v %in% colnames(meta_clean)) {
        stop(paste("Design variable", v, "not found in metadata columns:", paste(colnames(meta_clean), collapse=", ")))
    }
    # Convert to factor
    meta_clean[[v]] <- factor(meta_clean[[v]])
    # Remove NA
    if (any(is.na(meta_clean[[v]]))) {
        warning(paste("NA values found in design variable", v, "- dropping samples"))
        valid <- !is.na(meta_clean[[v]])
        meta_clean <- meta_clean[valid,]
        count_matrix <- count_matrix[, valid]
    }
}

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta_clean,
                              design = design_formula)

# 3. Run DESeq
cat("Running DESeq...\n")
dds <- DESeq(dds)

# 4. Save basic results
saveRDS(dds, file=file.path(opt$output, "deseq2_results.rds"))

# VST transformation for PCA
vsd <- vst(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup=design_vars, returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() +
  theme_bw()

ggsave(file.path(opt$output, "pca_plot.pdf"), plot=p)

# Write results for the EXPLICIT contrast only (BG-005). Validate the factor and both
# levels against the fitted object before asking DESeq2 for the comparison, so a typo
# fails here rather than silently resolving to something else.
cf <- contrast_parts[1]; num <- contrast_parts[2]; den <- contrast_parts[3]
if (!cf %in% colnames(colData(dds))) {
  stop("Contrast factor '", cf, "' is not a column of colData(dds). Available: ",
       paste(colnames(colData(dds)), collapse = ", "), call. = FALSE)
}
lv <- levels(factor(colData(dds)[[cf]]))
missing_lv <- setdiff(c(num, den), lv)
if (length(missing_lv)) {
  stop("Contrast level(s) not present in '", cf, "': ", paste(missing_lv, collapse = ", "),
       "\n  Observed levels: ", paste(lv, collapse = ", "), call. = FALSE)
}
cat("Contrast:", cf, num, "vs", den, "\n")
res <- results(dds, contrast = c(cf, num, den))
write.csv(as.data.frame(res), file.path(opt$output, "deseq2_results_table.csv"))

cat("Analysis complete. Results saved to", opt$output, "\n")