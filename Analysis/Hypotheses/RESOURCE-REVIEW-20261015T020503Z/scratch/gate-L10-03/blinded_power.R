t0 <- Sys.time()
suppressMessages({library(data.table); library(edgeR); library(limma)})
set.seed(20261015)
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
rel <- file.path(root, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls")
L <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L10-03"
dge_all <- readRDS(file.path(rel, "merged_dge.rds"))
meta <- as.data.table(readRDS(file.path(rel, "meta_matched.rds")))
st <- fread(file.path(root, "figures/candidates/fig3e-s3-774-20261011T151616Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv"),
            select = c("gene_id_versioned","gene_name","axis","contrast","logFC","FDR"))[axis == "fibrosis"]
g <- meta[dataset == "GSE213621" & group_binary == "Disease"]
g[, bin := factor(c("F0F1","F2","F3F4")[fibrosis_stage], levels = c("F0F1","F2","F3F4"))]
stopifnot(all(g$fibrosis_stage %in% 1:3)); print(table(g$bin))
setorder(g, sample_id)
dg <- dge_all[, match(g$sample_id, colnames(dge_all))]
rm(dge_all); invisible(gc())
hits <- list(binF3F4 = st[contrast == "F3_vs_F0" & FDR < 0.05], binF2 = st[contrast == "F2_vs_F0" & FDR < 0.05])
cat("x-axis hits (3E F3_vs_F0, F2_vs_F0):", sapply(hits, nrow), "\n")
fitp <- function() {
  lab <- sample(g$bin)
  des <- model.matrix(~ factor(inferred_sex) + lab, data.frame(inferred_sex = g$inferred_sex, lab = lab))
  colnames(des) <- c("Int","sexM","binF2","binF3F4")
  eBayes(lmFit(voom(dg, des), des))
}
res <- list(); k <- 0
f1 <- NULL
while (as.numeric(difftime(Sys.time(), t0, units = "secs")) < 200 && k < 40) {
  k <- k + 1; f <- fitp(); if (is.null(f1)) f1 <- f
  for (cf in names(hits)) {
    h <- hits[[cf]]; y <- f$coefficients[h$gene_id_versioned, cf]; x <- h$logFC
    res[[length(res)+1]] <- data.table(perm = k, contrast = cf, n_hits = length(x),
      slope = unname(coef(lm(y ~ x))[2]), conc = mean(sign(y) == sign(x)))
  }
}
res <- rbindlist(res); fwrite(res, file.path(L, "null_slopes_plainvoom.tsv"), sep = "\t")
z80 <- qnorm(0.8)
se <- data.table(gene = rownames(dg), F3F4 = f1$stdev.unscaled[, "binF3F4"] * sqrt(f1$s2.post),
                 F2 = f1$stdev.unscaled[, "binF2"] * sqrt(f1$s2.post), sigma = f1$sigma)
fwrite(se, file.path(L, "blinded_se_plainvoom.tsv.gz"), sep = "\t")
cat("perms:", k, " elapsed", format(Sys.time() - t0), "\n")
cat("Median sigma (all genes):", median(se$sigma), "\n")
for (cf in names(hits)) {
  h <- hits[[cf]]; s <- se[[sub("bin","",cf)]][match(h$gene_id_versioned, se$gene)]
  r <- res[contrast == cf]
  cat(sprintf("%s: n_hits %d; median SE %.4f; median per-gene MDE nominal0.05 %.3f, alpha0.005 %.3f, Bonferroni %.3f; median |disc logFC| %.3f\n",
    cf, nrow(h), median(s), median(2.8016*s), median((qnorm(1-0.0025)+z80)*s), median((qnorm(1-0.05/24115/2)+z80)*s), median(abs(h$logFC))))
  cat(sprintf("   null slope mean %.4f sd %.4f; null conc mean %.4f sd %.4f; MDE slope (one-sided 0.025, 80%%) %.4f; MDE conc %.4f\n",
    mean(r$slope), sd(r$slope), mean(r$conc), sd(r$conc), 2.8016*sd(r$slope), 2.8016*sd(r$conc)))
}
pc <- c("COL1A1","LUM","EPCAM","CYP2C19","GNMT")
gn <- unique(st[, .(gene_id_versioned, gene_name)])
pcid <- gn[gene_name %in% pc]
lc <- cpm(dg, log = TRUE, prior.count = 2)
out <- data.table(gene_name = pcid$gene_name, mean_logCPM = rowMeans(lc[pcid$gene_id_versioned, ]),
  frac_count_ge10 = rowMeans(dg$counts[pcid$gene_id_versioned, ] >= 10),
  SE_F3F4 = se$F3F4[match(pcid$gene_id_versioned, se$gene)], SE_F2 = se$F2[match(pcid$gene_id_versioned, se$gene)])
out[, `:=`(MDE_F3F4_nom = 2.8016*SE_F3F4, MDE_F3F4_bonf = (qnorm(1-0.05/24115/2)+z80)*SE_F3F4)]
print(out)
fwrite(out, file.path(L, "positive_control_blinded.tsv"), sep = "\t")
