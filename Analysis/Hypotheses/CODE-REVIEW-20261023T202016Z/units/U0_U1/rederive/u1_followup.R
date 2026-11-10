# Follow-up: endpoint fit with TMM (default calcNormFactors) on the 392 subset; PCA variants.
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table); library(matrixStats)})
args <- commandArgs(TRUE); OUT <- args[1]; U0 <- args[2]; U1 <- args[3]
dir.create(OUT, showWarnings = FALSE)
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
counts <- readRDS(file.path(R, "figures/candidates/resource-corrections-20261005-v1/verified_counts_r8_20261007/assembled_774_primary/counts_raw.rds"))
lab <- fread(file.path(U0, "my_774_labels.tsv"))[match(colnames(counts), sample_id)]
sx <- fread(file.path(U1, "sex_calls.tsv"))[match(colnames(counts), sample_id)]
mine <- fread(file.path(U1, "my_refit_774_deg.tsv"))
group <- factor(lab$group_binary, levels = c("Control", "Disease")); dataset <- factor(lab$dataset); sex <- factor(sx$expr_sex)
y <- DGEList(counts); keep <- filterByExpr(y, design = model.matrix(~0 + group)); y <- y[keep, , keep.lib.sizes = FALSE]
stopifnot(identical(rownames(y), mine$gene))
yR <- calcNormFactors(y, method = "RLE")
ep <- fread(file.path(U0, "my_endpoint_labels.tsv")); j <- match(ep$sample_id, colnames(y))
for (nf in c("TMM_on_subset", "RLE_on_subset")) {
  ys <- calcNormFactors(y[, j], method = sub("_on_subset", "", nf))
  eg <- factor(ep$endpoint, levels = c("StrictControl", "Advanced"))
  d2 <- model.matrix(~ factor(ep$dataset) + sex[j] + eg)
  t2 <- topTable(eBayes(lmFit(voomWithQualityWeights(ys, d2), d2)), coef = "egAdvanced", number = Inf, sort.by = "none")
  c2 <- t2$adj.P.Val < 0.05 & abs(t2$logFC) > 0.5
  cat(sprintf("[endpoint %s] rho=%.4f calls=%d recovered=%d of %d endpoint-only=%d opposite=%d\n", nf,
      cor(t2$logFC, mine$logFC, method = "spearman"), sum(c2), sum(c2 & mine$canon & sign(t2$logFC) == sign(mine$logFC)),
      sum(mine$canon), sum(c2 & !mine$canon), sum(c2 & mine$canon & sign(t2$logFC) != sign(mine$logFC))))
}
for (nm in c("RLE", "TMM")) {
  yy <- calcNormFactors(y, method = nm)
  for (wq in c("vqw", "voom")) {
    des <- model.matrix(~ dataset + sex + group)
    E <- if (wq == "vqw") voomWithQualityWeights(yy, des)$E else voom(yy, des)$E
    for (cov in c("batch2", "covariate")) {
      Ec <- if (cov == "batch2") removeBatchEffect(E, batch = dataset, batch2 = sex) else removeBatchEffect(E, batch = dataset, covariates = as.numeric(sex == "M"))
      for (sel in c("after", "before")) {
        vv <- if (sel == "after") rowVars(Ec) else rowVars(E)
        top <- order(vv, decreasing = TRUE)[1:2000]
        for (sc in c(FALSE, TRUE)) {
          p <- prcomp(t(Ec[top, ]), scale. = sc); pv <- 100 * p$sdev^2 / sum(p$sdev^2)
          cat(sprintf("[PCA %s %s %s hvg_%s scale=%s] PC1 %.2f PC2 %.2f\n", nm, wq, cov, sel, sc, pv[1], pv[2]))
        }
      }
    }
  }
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
