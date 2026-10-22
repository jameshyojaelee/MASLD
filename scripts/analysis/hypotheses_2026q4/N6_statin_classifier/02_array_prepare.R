#!/usr/bin/env Rscript
# N6: gene-level matrix from the deposited GSE130991 RMA values (HTA 2.0, GPL20265).
# Label-free: clusters -> first symbol of gene_assignment (category main), mean per symbol,
# keep the upper 60% of genes by mean expression across all arrays.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: 02_array_prepare.R series_matrix.txt.gz platform_table.txt samples.tsv out_dir")
mat_gz <- args[1]; plat_f <- args[2]; samp_f <- args[3]; out <- args[4]
if (dir.exists(out)) stop("Refusing to overwrite output dir")
dir.create(out, recursive = TRUE)
set.seed(20260923)

smp <- fread(samp_f, na.strings = "NA", check.names = FALSE)
expr <- fread(cmd = sprintf("zcat %s | sed -n '/^!series_matrix_table_begin/,/^!series_matrix_table_end/p' | sed '1d;$d'", shQuote(mat_gz)),
              sep = "\t", colClasses = list(character = 1L), quote = "\"")
setnames(expr, 1L, "ID")
stopifnot(identical(colnames(expr)[-1L], smp$gsm))
cat("clusters in matrix:", nrow(expr), " arrays:", ncol(expr) - 1L, "\n")

plat <- fread(cmd = sprintf("sed '1d;$d' %s", shQuote(plat_f)), sep = "\t", quote = "",
              select = c("ID", "category", "gene_assignment"), colClasses = list(character = 1:3))
plat <- plat[category %in% c("main", "main///normgene->exon")]
first <- sub(" /// .*$", "", plat$gene_assignment)
sym <- vapply(strsplit(first, " // ", fixed = TRUE), function(p) if (length(p) >= 2L) p[2L] else NA_character_, "")
sym[sym %in% c("---", "", "NULL")] <- NA_character_
plat[, symbol := sym]
map <- plat[!is.na(symbol), .(ID, symbol)]
j <- merge(expr[, "ID"], map, by = "ID")
cat("clusters with main category and a symbol:", nrow(j), "\n")

m <- as.matrix(expr[match(j$ID, expr$ID), -1L, with = FALSE])
storage.mode(m) <- "double"
S <- rowsum(m, group = j$symbol, reorder = TRUE)
cnt <- as.vector(table(factor(j$symbol, levels = rownames(S))))
G <- S / cnt
stopifnot(all(is.finite(G)))
gm <- rowMeans(G)
keep <- gm >= quantile(gm, 0.4)
Gk <- G[keep & apply(G, 1, sd) > 0, , drop = FALSE]
cat("symbols:", nrow(G), " kept (upper 60% by mean, nonzero variance):", nrow(Gk), "\n")
saveRDS(list(G = Gk, G_all = G, clusters_per_symbol = setNames(cnt, rownames(S))), file.path(out, "gene_matrix.rds"))
write_json(list(arrays = ncol(Gk), clusters_in_matrix = nrow(expr), clusters_mapped = nrow(j),
                symbols_total = nrow(G), symbols_kept = nrow(Gk), mean_cutoff_log2 = unname(quantile(gm, 0.4))),
           file.path(out, "prep_summary.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_prepare.txt"))
