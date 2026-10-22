# Recompute Fig 3D PCA with the stated recipe (read from scripts/figures/fig3d_endpoint_validation.R after the
# independent variants failed to match): TMM on the 392-unit subset, voom log-CPM, removeBatchEffect(dataset,
# inferred sex), top 2000 variable genes after adjustment, unscaled PCA. Roster is our own recomputation.
set.seed(20261017)
suppressPackageStartupMessages({library(edgeR); library(limma); library(matrixStats)})
base <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
O <- "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
dge_all <- readRDS(paste0(base, "merged_dge.rds"))
ros <- read.delim(paste0(O, "fig3d_recomputed_roster.tsv"))
dge <- calcNormFactors(dge_all[, ros$sample_id])
lib <- dge$samples$lib.size * dge$samples$norm.factors
E <- t(log2(t(dge$counts + 0.5) / (lib + 1) * 1e6))
sex_cov <- model.matrix(~ factor(ros$inferred_sex))[, -1, drop = FALSE]
adj <- removeBatchEffect(E, batch = factor(ros$dataset), covariates = sex_cov)
for (ntop in c(2000L, 1000L, 5000L)) {
  hvg <- order(rowVars(adj), decreasing = TRUE)[seq_len(ntop)]
  p <- prcomp(t(adj[hvg, ]), center = TRUE, scale. = FALSE)
  if (mean(p$x[ros$arm == "advanced", 1]) < mean(p$x[ros$arm == "strict_control", 1])) p$x[, 1] <- -p$x[, 1]
  pve <- 100 * p$sdev^2 / sum(p$sdev^2)
  cat("top", ntop, "PC1", round(pve[1], 2), "PC2", round(pve[2], 2), "\n")
  write.table(data.frame(sample_id = ros$sample_id, PC1 = p$x[, 1], PC2 = p$x[, 2]),
              paste0(O, "pca_variants/exactTMM392_top", ntop, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
# PC1 association with fibrosis among advanced units, exact stages only vs including GSE213621 bins
p <- read.delim(paste0(O, "pca_variants/exactTMM392_top2000.tsv")); d <- merge(p, ros)
a <- d[d$arm == "advanced", ]
cat("Spearman PC1~fibrosis, all advanced:", cor(a$PC1, a$fibrosis_stage, method = "spearman"), "n", nrow(a), "\n")
e <- a[a$dataset != "GSE213621", ]
cat("Spearman PC1~fibrosis, exact-stage advanced:", cor(e$PC1, e$fibrosis_stage, method = "spearman"), "n", nrow(e), "\n")
cat("mean PC1 by dataset x arm\n"); print(aggregate(PC1 ~ dataset + arm, d, mean))
cat("mean PC1 GSE213621 F3F4 bin vs exact F3 vs exact F4:", mean(a$PC1[a$dataset=="GSE213621"]), mean(e$PC1[e$fibrosis_stage==3]), mean(e$PC1[e$fibrosis_stage==4]), "\n")
writeLines(capture.output(sessionInfo()), paste0(O, "sessionInfo_12.txt"))
