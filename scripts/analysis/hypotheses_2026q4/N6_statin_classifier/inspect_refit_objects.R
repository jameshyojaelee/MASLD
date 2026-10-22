suppressPackageStartupMessages({library(edgeR); library(data.table)})
d <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
dge <- readRDS(file.path(d, "merged_dge.rds")); meta <- readRDS(file.path(d, "meta_matched.rds"))
cat("dge class:", class(dge), " dim:", dim(dge), "\n")
cat("samples cols:", colnames(dge$samples), "\n")
cat("genes cols:", if (is.null(dge$genes)) "NULL" else colnames(dge$genes), "\n")
cat("rownames head:", head(rownames(dge), 3), "\n")
print(head(dge$samples, 3)); print(summary(dge$samples$norm.factors))
cat("meta class:", class(meta), " dim:", dim(meta), "\n"); print(names(meta)); print(head(as.data.frame(meta), 3))
for (nm in intersect(names(meta), c("condition", "stage", "group_binary", "fibrosis", "nas", "histology"))) { cat("--", nm, "\n"); print(table(meta[[nm]], meta$dataset, useNA = "ifany")) }
print(table(meta$inferred_sex, meta$dataset, useNA = "ifany"))
dg <- fread(file.path(d, "deg_results.csv"), nrows = 3); print(dg)
cat("identical sample order dge vs meta:", identical(colnames(dge), meta$sample_id), "\n")
