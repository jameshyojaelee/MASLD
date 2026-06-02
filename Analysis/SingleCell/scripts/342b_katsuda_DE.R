#!/usr/bin/env Rscript
# 342b_katsuda_DE.R - Katsuda 2019 rat 2c/4c/8c hepatocyte Agilent microarray DE via limma.
# Series matrix has normalized expression; treat as standard one-color array.
# Test (8c+4c) vs 2c (polyploid extreme contrast).

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
})
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
setwd(ROOT)

MAT  <- "data/external/katsuda2019/GSE132409_series_matrix.txt.gz"
ORTH <- "data/external/orthologs/rat_human_orthologs.tsv.gz"
OUT_DIR <- "data/ploidy_signatures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time()), "] Parse series matrix\n", sep="")
# Read the entire file, find the data section
lines <- readLines(gzfile(MAT))
# Find the line that starts with "!series_matrix_table_begin"
begin_idx <- grep("series_matrix_table_begin", lines, ignore.case = TRUE)
end_idx   <- grep("series_matrix_table_end",   lines, ignore.case = TRUE)
cat("  data block lines:", begin_idx, "to", end_idx, "\n")

# Parse the sample metadata section
meta_lines <- grep("^!Sample_", lines, value = TRUE)
# Pull title + ploidy
get_meta_row <- function(prefix) {
  ln <- meta_lines[grep(prefix, meta_lines, fixed = TRUE)][1]
  if (is.na(ln)) return(NULL)
  vals <- strsplit(ln, "\t")[[1]][-1]
  gsub('"', '', vals)
}
titles <- get_meta_row("!Sample_title")
gsms   <- get_meta_row("!Sample_geo_accession")
ploidy_lines <- meta_lines[grep("ploidy:", meta_lines)]
ploidy_vals <- gsub('.*ploidy: ([0-9]+c)".*', '\\1', ploidy_lines[1])
# fallback: parse all
ploidy_vec <- strsplit(ploidy_lines[1], "\t")[[1]][-1]
ploidy_vec <- gsub('"|ploidy: ', '', ploidy_vec)

cat("  titles:", titles, "\n")
cat("  GSM ids:", gsms, "\n")
cat("  ploidy values:", ploidy_vec, "\n")

# Parse the expression data block
data_text <- lines[(begin_idx+1):(end_idx-1)]
data_dt <- fread(text = paste(data_text, collapse = "\n"), header = TRUE, quote = "\"")
# First column = probe ID
setnames(data_dt, 1, "PROBE_ID")
cat("  expression matrix:", nrow(data_dt), "probes x", ncol(data_dt) - 1, "samples\n")

expr_mat <- as.matrix(data_dt[, -1])
rownames(expr_mat) <- data_dt$PROBE_ID
colnames(expr_mat) <- gsms

# Map probes -> gene symbol. GPL22740 is Agilent rat array. We can try ANNOTATION via the GPL file.
# For now, query Agilent's gene name from the SOFT annotation file embedded in feature extraction output
# OR fetch GPL via GEOquery. Simplest: assume probe IDs are Agilent probe accessions; map via a
# project-side annotation if available. If not, do DE on probes and map later.

# Sample data frame
sample_dt <- data.table(GSM = gsms, title = titles, ploidy = ploidy_vec)
sample_dt[, replicate := tstrsplit(title, "_", keep = 2)]
sample_dt[, ploidy_factor := factor(ploidy, levels = c("2c", "4c", "8c"))]
sample_dt[, polyploid_binary := factor(ifelse(ploidy == "2c", "2c", "polyploid"), levels = c("2c", "polyploid"))]
cat("\nSample design:\n"); print(sample_dt)

# Make sure expr_mat has no NA rows
ok <- rowSums(is.na(expr_mat)) == 0
expr_mat <- expr_mat[ok, , drop = FALSE]
cat("  probes after NA filter:", nrow(expr_mat), "\n")
cat("  expr range:", round(min(expr_mat), 2), "to", round(max(expr_mat), 2), "\n")
# Series matrix is already normalized (75th pct + median-center per gene per Katsuda). Log2 if not.
if (max(expr_mat) > 50) {
  cat("  Log2-transforming\n")
  expr_mat <- log2(expr_mat + 1)
}

# limma DE: polyploid (4c+8c) vs 2c, blocked on replicate
design <- model.matrix(~ replicate + polyploid_binary, data = sample_dt)
fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "polyploid_binarypolyploid", number = Inf, sort.by = "none")
res$PROBE_ID <- rownames(res)
cat("\nDE summary (4c+8c vs 2c):\n")
cat("  probes with adj.P.Val<0.05:", sum(res$adj.P.Val < 0.05, na.rm=TRUE), "\n")
cat("  probes with |logFC|>0.5 & adj.P.Val<0.05:", sum(res$adj.P.Val < 0.05 & abs(res$logFC) > 0.5, na.rm=TRUE), "\n")
cat("  Top probes:\n")
print(head(res[order(res$adj.P.Val), c("PROBE_ID","logFC","P.Value","adj.P.Val","t")], 10))

# Also continuous ploidy as predictor
sample_dt[, ploidy_numeric := as.integer(gsub("c", "", ploidy))]
design2 <- model.matrix(~ replicate + ploidy_numeric, data = sample_dt)
fit2 <- lmFit(expr_mat, design2)
fit2 <- eBayes(fit2)
res2 <- topTable(fit2, coef = "ploidy_numeric", number = Inf, sort.by = "none")
res2$PROBE_ID <- rownames(res2)
cat("\nContinuous ploidy trend:\n")
cat("  probes with adj.P.Val<0.05:", sum(res2$adj.P.Val < 0.05, na.rm=TRUE), "\n")

# Save raw DE outputs
res$contrast <- "polyploid_vs_2c"
res2$contrast <- "continuous_ploidy"
all_de <- rbind(res, res2)
fwrite(all_de, file.path(OUT_DIR, "katsuda2019_de_full.tsv"), sep = "\t")
cat("\nWrote katsuda2019_de_full.tsv (", nrow(all_de), "rows)\n")

# Map probe -> gene_symbol via GPL22740. Query via GEOquery::getGEO("GPL22740") OR
# use Bioconductor's AgilentRGCEnvAccessions if available
cat("\n[", format(Sys.time()), "] Mapping probes to gene symbols\n", sep="")
# Try via GPL22740 download
gpl_url <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL22nnn/GPL22740/annot/GPL22740.annot.gz"
gpl_path <- "data/external/katsuda2019/GPL22740.annot.gz"
if (!file.exists(gpl_path)) {
  download.file(gpl_url, gpl_path, mode = "wb", quiet = TRUE)
}
if (file.exists(gpl_path) && file.size(gpl_path) > 100) {
  gpl_lines <- readLines(gzfile(gpl_path), n = 200)
  cat("  GPL header lines:\n"); cat(head(gpl_lines, 30), sep = "\n")
} else {
  cat("  GPL not available, falling back to no symbol mapping\n")
}
cat("[", format(Sys.time()), "] DONE Katsuda DE (probe-level)\n", sep="")
