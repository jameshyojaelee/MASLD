# H57 step 0 (prespec D3): control-only mean log2-CPM per gene from the refit DGEList.
# voom formula: log2((count + 0.5) / (lib.size * norm.factors + 1) * 1e6).
# Writes the all-unit mean too, so the Python step can check it against AveExpr (< 1e-8).
suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE)
out <- args[which(args == "--out") + 1]
stopifnot(length(out) == 1, dir.exists(out))

rds <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/merged_dge.rds"
expect <- "73f3c53fcb9c847e98b0e4f0918c52a04a6b0b21744e5b448cb45e5cbcfcad3b"
got <- strsplit(system2("sha256sum", rds, stdout = TRUE), " +")[[1]][1]
if (!identical(got, expect)) stop(sprintf("S1: sha256 mismatch for %s: %s", rds, got))

x <- readRDS(rds)
lib <- x$samples$lib.size * x$samples$norm.factors
lcpm <- log2(t((t(x$counts) + 0.5) / (lib + 1)) * 1e6)
is_ctl <- x$samples$group_binary == "Control"
cat("units:", ncol(x$counts), " control:", sum(is_ctl), " disease:", sum(!is_ctl), "\n")
print(table(x$samples$dataset, x$samples$group_binary))
stopifnot(sum(is_ctl) == 138, sum(!is_ctl) == 636)

res <- data.frame(gene = rownames(x$counts),
                  ave_all_units = rowMeans(lcpm),
                  ave_control = rowMeans(lcpm[, is_ctl, drop = FALSE]))
write.table(res, file.path(out, "control_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("wrote", nrow(res), "genes\n")
