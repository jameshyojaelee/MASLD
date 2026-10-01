# Shared helpers for ALT/patch-copy count correction (see memory
# reference-bulk-counts-alt-contig-undercount-2026-09-23).
suppressPackageStartupMessages(library(data.table))

PRIMARY_CHR <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

gene_table <- function(gtf) {
  g <- fread(cmd = paste("zcat", shQuote(gtf), "| awk -F'\\t' '$3 == \"gene\"'"),
             header = FALSE, sep = "\t", quote = "")
  g[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", V9)]
  g[, gene_name := sub('.*gene_name "([^"]+)".*', "\\1", V9)]
  g[, primary := V1 %in% PRIMARY_CHR]
  g[, .(gene_id, gene_name, primary)]
}

# Read a featureCounts table; columns renamed to run IDs taken from BAM file names.
read_featurecounts <- function(path) {
  f <- fread(path, skip = "Geneid")
  runs <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(names(f)[-(1:6)]))
  m <- as.matrix(f[, -(1:6)])
  dimnames(m) <- list(f$Geneid, runs)
  m
}

# Replace counts of primary genes that share a name with an ALT/patch copy by the
# multimapper-inclusive count of the primary gene plus its copies (rounded).
# Genes whose name maps to more than one primary gene are left unchanged.
# Returns list(counts = corrected matrix restricted to primary-chromosome genes,
#              replaced = gene_ids replaced).
fold_alt_copies <- function(orig, multi, genes) {
  runs <- intersect(colnames(orig), colnames(multi))
  stopifnot(length(runs) == ncol(orig))
  alt_names <- unique(genes[primary == FALSE, gene_name])
  one <- genes[primary == TRUE & gene_name %in% alt_names, .N, by = gene_name][N == 1, gene_name]
  fold <- genes[gene_name %in% one]
  ids <- intersect(fold$gene_id, rownames(multi))
  by_name <- rowsum(multi[ids, runs, drop = FALSE], fold[match(ids, gene_id), gene_name])
  tgt <- fold[primary == TRUE & gene_id %in% rownames(orig) & gene_name %in% rownames(by_name)]
  out <- orig[, runs, drop = FALSE]
  out[tgt$gene_id, ] <- round(by_name[tgt$gene_name, , drop = FALSE])
  keep <- rownames(out) %in% genes[primary == TRUE, gene_id]
  list(counts = out[keep, , drop = FALSE], replaced = tgt$gene_id)
}
