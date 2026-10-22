# Recompute Fig 3D PCA on the 392-unit endpoint roster from the adopted 774-unit DGEList.
# Group-blind: removeBatchEffect(dataset, inferred sex) with no protected design. Several gene-set/scale variants
# are computed because the exact preprocessing is not stated in the legend; the variant matching PC1/PC2 % is reported.
set.seed(20261017)
suppressPackageStartupMessages({library(edgeR); library(limma)})
base <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
O <- "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
dge <- readRDS(paste0(base, "merged_dge.rds"))
cat("class", class(dge), "dim", dim(dge), "\n"); print(head(dge$samples, 3))
ros <- read.delim(paste0(O, "fig3d_recomputed_roster.tsv"))
meta <- readRDS(paste0(base, "meta_matched.rds"))
stopifnot(all(ros$sample_id %in% colnames(dge)))
variants <- list()
lc_all <- cpm(dge, log = TRUE, prior.count = 2)
# voom-style log-CPM
lib <- dge$samples$lib.size * dge$samples$norm.factors
E_voom <- t(log2(t(dge$counts + 0.5) / (lib + 1) * 1e6))
for (src in c("cpm2", "voomE")) {
  X0 <- if (src == "cpm2") lc_all else E_voom
  for (scope in c("subset_then_adjust", "adjust774_then_subset")) {
    if (scope == "adjust774_then_subset") {
      m <- meta[match(colnames(X0), meta$sample_id), ]
      Xa <- removeBatchEffect(X0, batch = m$dataset, batch2 = m$inferred_sex)
      Xa <- Xa[, ros$sample_id]
    } else {
      X <- X0[, ros$sample_id]
      Xa <- removeBatchEffect(X, batch = ros$dataset, batch2 = ros$inferred_sex)
    }
    v <- apply(Xa, 1, var)
    for (ntop in c(500, 1000, 2000, 5000, nrow(Xa))) {
      for (sc in c(FALSE, TRUE)) {
        g <- order(v, decreasing = TRUE)[seq_len(ntop)]
        p <- prcomp(t(Xa[g, ]), center = TRUE, scale. = sc)
        ve <- p$sdev^2 / sum(p$sdev^2)
        key <- sprintf("%s|%s|top%d|scale%s", src, scope, ntop, sc)
        variants[[key]] <- data.frame(variant = key, pc1_pct = 100 * ve[1], pc2_pct = 100 * ve[2])
        write.table(data.frame(sample_id = rownames(p$x), PC1 = p$x[, 1], PC2 = p$x[, 2], variant = key),
                    paste0(O, "pca_variants/", gsub("[|]", "_", key), ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      }
    }
  }
}
res <- do.call(rbind, variants); print(res, digits = 4)
write.table(res, paste0(O, "fig3d_pca_variant_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(O, "sessionInfo_09.txt"))
