# Fig 3D display-choice sensitivity. Same log-CPM (TMM on the 392 units, prior 0.5) and top-2000 unscaled PCA as the
# drawn panel, under three batch treatments: (a) drawn: removeBatchEffect(dataset, sex) with no protected design;
# (b) endpoint protected: design = ~endpoint passed to removeBatchEffect; (c) per-cohort centring on that cohort's
# strict controls (sex regressed). Reports PVE, AUC of PC1 advanced vs control, Spearman PC1 ~ fibrosis among advanced
# (all and exact-stage cohorts only), and per-cohort AUC.
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
adv <- ros$arm == "advanced"
auc <- function(x, y) { r <- rank(x); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
variants <- list()
variants$a_unprotected <- removeBatchEffect(E, batch = factor(ros$dataset), covariates = sex_cov)
variants$b_endpoint_protected <- removeBatchEffect(E, batch = factor(ros$dataset), covariates = sex_cov, design = model.matrix(~ adv))
s <- removeBatchEffect(E, covariates = sex_cov, design = model.matrix(~ factor(ros$dataset) + adv))
cc <- s; for (d in unique(ros$dataset)) { k <- ros$dataset == d; cc[, k] <- s[, k] - rowMeans(s[, k & !adv, drop = FALSE]) }
variants$c_control_centred <- cc
variants$d_no_batch_removal <- E
res <- list(); pcs <- list()
for (nm in names(variants)) {
  X <- variants[[nm]]; hvg <- order(rowVars(X), decreasing = TRUE)[1:2000]
  p <- prcomp(t(X[hvg, ]), center = TRUE, scale. = FALSE)
  if (mean(p$x[adv, 1]) < mean(p$x[!adv, 1])) p$x[, 1] <- -p$x[, 1]
  pve <- 100 * p$sdev^2 / sum(p$sdev^2)
  a <- ros[adv, ]; pc1a <- p$x[adv, 1]; ex <- a$dataset != "GSE213621"
  per <- sapply(sort(unique(ros$dataset)), function(d) { k <- ros$dataset == d; auc(p$x[k, 1], adv[k]) })
  r2 <- summary(lm(p$x[, 1] ~ factor(ros$dataset)))$r.squared
  res[[nm]] <- data.frame(variant = nm, PC1_pct = round(pve[1], 2), PC2_pct = round(pve[2], 2), AUC_PC1_adv_vs_ctrl = round(auc(p$x[, 1], adv), 3),
    spearman_fib_adv_all = round(cor(pc1a, a$fibrosis_stage, method = "spearman"), 3),
    spearman_fib_adv_exact = round(cor(pc1a[ex], a$fibrosis_stage[ex], method = "spearman"), 3),
    R2_PC1_on_cohort = round(r2, 3), t(round(per, 3)))
  pcs[[nm]] <- data.frame(sample_id = ros$sample_id, variant = nm, PC1 = p$x[, 1], PC2 = p$x[, 2])
}
R <- do.call(rbind, res); print(R, row.names = FALSE)
write.table(R, paste0(O, "fig3d_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, pcs), paste0(O, "fig3d_sensitivity_coords.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(O, "sessionInfo_18.txt"))
