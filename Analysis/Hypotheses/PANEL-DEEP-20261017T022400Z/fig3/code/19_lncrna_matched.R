# lncRNA inset: (1) gene-level logistic OR of a canonical call, lncRNA vs protein-coding, adjusted for abundance
# (AveExpr), variability (SD of log-CPM over 774 units), union-exon gene length (GENCODE v49) and cohort detection
# (number of 5 cohorts in which >=50% of units have a count); (2) 1:1 greedy nearest-neighbour matching of lncRNA calls
# to protein-coding calls on the same four z-scored covariates, 0.15-SD per-covariate caliper, random order (seeded);
# conditional OR of upregulation = discordant lncRNA-up / coding-up pairs.
set.seed(20261017)
suppressPackageStartupMessages({library(edgeR); library(splines); library(rtracklayer); library(GenomicRanges)})
base <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
O <- "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
dge <- readRDS(paste0(base, "merged_dge.rds")); meta <- readRDS(paste0(base, "meta_matched.rds"))
deg <- read.csv(paste0(base, "deg_results.csv"))
bt <- read.delim("figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv")[, c("gene_id_versioned", "gene_type")]
bt <- bt[!duplicated(bt$gene_id_versioned), ]
d <- merge(deg, bt, by.x = "gene", by.y = "gene_id_versioned"); d <- d[d$gene_type %in% c("protein_coding", "lncRNA"), ]
cnt <- dge$counts[d$gene, ]; lcpm <- cpm(dge[d$gene, ], log = TRUE, prior.count = 0.5)
d$variab <- apply(lcpm, 1, sd)
ds <- as.character(meta$dataset[match(colnames(cnt), meta$sample_id)])
if (any(is.na(ds))) ds <- as.character(dge$samples$dataset)
d$cohort_detect <- rowSums(sapply(unique(ds), function(k) rowMeans(cnt[, ds == k, drop = FALSE] > 0) >= 0.5))
gtf <- import("/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz", feature.type = "exon")
ex <- split(granges(gtf), gtf$gene_id); len <- sum(width(reduce(ex)))
d$len <- log10(len[d$gene]); cat("genes missing length:", sum(is.na(d$len)), "\n"); d <- d[!is.na(d$len), ]
d$call <- as.integer(d$padj < 0.05 & abs(d$logFC) > 0.5); d$lnc <- as.integer(d$gene_type == "lncRNA")
m <- glm(call ~ lnc + ns(AveExpr, 4) + ns(variab, 4) + ns(len, 4) + factor(cohort_detect), binomial, d)
ci <- exp(confint.default(m)["lnc", ]); cat("adjusted OR (spline abundance, variability, length; cohort-detection factor):", round(exp(coef(m)["lnc"]), 3), round(ci, 3), "n genes", nrow(d), "\n")
m2 <- glm(call ~ lnc + AveExpr + variab + len + cohort_detect, binomial, d)
ci <- exp(confint.default(m2)["lnc", ]); cat("adjusted OR (linear terms):", round(exp(coef(m2)["lnc"]), 3), round(ci, 3), "\n")
# matching among canonical calls
cv <- c("AveExpr", "variab", "len", "cohort_detect")
calls <- d[d$call == 1, ]; Z <- scale(calls[, cv]); calls$up <- calls$logFC > 0
L <- which(calls$lnc == 1); P <- which(calls$lnc == 0)
for (cal in c(0.15, 0.20)) {
  used <- rep(FALSE, nrow(calls)); pairs <- NULL
  for (i in sample(L)) {
    cand <- P[!used[P]]; dz <- abs(sweep(Z[cand, , drop = FALSE], 2, Z[i, ]))
    ok <- cand[apply(dz, 1, max) <= cal]; if (!length(ok)) next
    j <- ok[which.min(rowSums(sweep(Z[ok, , drop = FALSE], 2, Z[i, ])^2))]; used[j] <- TRUE; pairs <- rbind(pairs, c(i, j))
  }
  b <- sum(calls$up[pairs[, 1]] & !calls$up[pairs[, 2]]); c2 <- sum(!calls$up[pairs[, 1]] & calls$up[pairs[, 2]])
  smd <- (colMeans(calls[pairs[, 1], cv]) - colMeans(calls[pairs[, 2], cv])) / apply(calls[, cv], 2, sd)
  se <- sqrt(1 / b + 1 / c2)
  cat("caliper", cal, "pairs", nrow(pairs), "discordant lnc-up", b, "coding-up", c2, "OR", round(b / c2, 3), "CI", round(exp(log(b / c2) + c(-1.96, 1.96) * se), 3), "max |SMD|", round(max(abs(smd)), 3), "\n")
}
writeLines(capture.output(sessionInfo()), paste0(O, "sessionInfo_19.txt"))
