# H79 audit: rank of the auxiliary cell-means design and GSE126848 sex. No expression read.
suppressPackageStartupMessages(library(data.table))
setwd('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
D <- 'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls'
a <- fread('RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv')
m <- as.data.table(readRDS(file.path(D, 'meta_matched.rds')))
x <- cbind(m, a[match(m$sample_id, a$run), .(PC1, PC2, PC3, PC4, possibly_mixed)])
x[, cell := interaction(dataset, group_binary, drop = TRUE)]
X <- model.matrix(~ 0 + cell + inferred_sex + PC1 + PC2 + PC3 + PC4, x)
cat('aux cell-means design rank', qr(X)$rank, 'of', ncol(X), '\n')
set.seed(20260923)
pm <- copy(x); idx <- unlist(lapply(split(seq_len(nrow(x)), x$cell), function(i) i[sample.int(length(i))]))
ord <- unlist(split(seq_len(nrow(x)), x$cell))
pm[ord, c('PC1','PC2','PC3','PC4') := x[idx, .(PC1, PC2, PC3, PC4)]]
cat('one within-cell permutation: max |cell-mean diff|', max(abs(as.matrix(pm[, lapply(.SD, mean), by = cell, .SDcols = PC1:PC4][, -1]) - as.matrix(x[, lapply(.SD, mean), by = cell, .SDcols = PC1:PC4][, -1]))), ' rank', qr(model.matrix(~ 0 + cell + inferred_sex + PC1 + PC2 + PC3 + PC4, pm))$rank, '\n')
cat('GSE126848 inferred_sex x group\n'); print(table(x[dataset == 'GSE126848', .(inferred_sex, group_binary)]))
cat('D2 n without possibly_mixed', sum(!x$possibly_mixed), '\n')
dge_cols <- colnames(readRDS(file.path(D, 'merged_dge.rds')))
cat('DGEList columns', length(dge_cols), 'same set as meta', setequal(dge_cols, m$sample_id), 'same order as meta', identical(dge_cols, m$sample_id), '\n')
