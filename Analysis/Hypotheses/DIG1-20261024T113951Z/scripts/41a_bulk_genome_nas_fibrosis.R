# DIG1 job 41a (exploratory; input for job 41b). Genome-wide per-gene association with NAS and fibrosis
# stage within bulk disease units, five-cohort fragment-count DGE (adopted 774-unit release).
# Partial Spearman: ranks of log CPM and of the score over pooled disease units, residualised on
# dataset + sex_final. NAS: GSE130970, GSE135251, GSE162694 (disease units with NAS). Fibrosis:
# disease units with a Kleiner stage (GSE130970, GSE135251, GSE162694, GSE213621). BH over genes.
# Genes: log CPM > 1 in at least half of all units (as job 36a). Also the pooled disease-vs-control
# logFC, padj from the adopted fit (deg_results.csv) joined by Ensembl ID.
# Usage: Rscript 41a_bulk_genome_nas_fibrosis.R <new_out_dir>
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(TRUE); stopifnot(length(args) == 1L, !dir.exists(args[1]))
out <- args[1]; dir.create(out, recursive = TRUE)
set.seed(20261024)
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
y <- readRDS(file.path(B, "merged_dge.rds")); m <- readRDS(file.path(B, "meta_matched.rds"))
m <- m[match(colnames(y), m$sample_id)]; stopifnot(identical(m$sample_id, colnames(y)))
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE, prior.count = 1)
lc <- lc[rowMeans(lc > 1) >= 0.5, ]
ens <- sub("\\..*$", "", rownames(lc))
pc_all <- function(sel, score) {
  z <- m[sel]; s <- as.numeric(z[[score]]); ok <- !is.na(s); z <- z[ok]; s <- s[ok]; X <- lc[, which(sel)[ok], drop = FALSE]
  D <- model.matrix(~ dataset + sex_final, data = z); q <- qr(D)
  rs <- qr.resid(q, rank(s)); RX <- qr.resid(q, apply(X, 1, rank))
  r <- as.numeric(crossprod(RX, rs)) / sqrt(colSums(RX^2) * sum(rs^2)); df <- nrow(D) - ncol(D) - 1
  p <- 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)
  data.table(ensembl = ens, n = nrow(D), rho = r, p = p, bh = p.adjust(p, "BH")) }
dis <- m$group_binary == "Disease"
nas <- pc_all(dis & m$dataset %in% c("GSE130970", "GSE135251", "GSE162694"), "nas_score")
fib <- pc_all(dis & m$dataset %in% c("GSE130970", "GSE135251", "GSE162694", "GSE213621"), "fibrosis_stage")
de <- fread(file.path(B, "deg_results.csv"))[, .(ensembl = sub("\\..*$", "", gene), symbol, logFC, t, padj)]
res <- merge(de, merge(nas, fib, by = "ensembl", suffixes = c("_nas", "_fib")), by = "ensembl", all = TRUE)
fwrite(res, file.path(out, "gene_disease_state.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("genes", nrow(res), "NAS n", nas$n[1], "fib n", fib$n[1], "\n")
print(res[symbol %in% c("HSPA4", "ACTG1", "KRT8", "PNPLA3", "TM6SF2", "HSD17B13", "GPAM", "MBOAT7")], digits = 3)
