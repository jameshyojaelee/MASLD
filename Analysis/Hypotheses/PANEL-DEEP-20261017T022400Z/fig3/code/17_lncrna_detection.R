# lncRNA inset: expression- and detection-adjusted odds ratio of a canonical call (lncRNA vs protein-coding) from the
# 774-unit fit; detection = fraction of the 774 libraries with >=1 count. Natural splines (df 5) for both covariates.
set.seed(20261017)
suppressPackageStartupMessages({library(edgeR); library(splines)})
base <- "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
O <- "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
dge <- readRDS(paste0(base, "merged_dge.rds"))
deg <- read.csv(paste0(base, "deg_results.csv"))
bt <- read.delim("figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv")[, c("gene_id_versioned", "gene_type")]
bt <- bt[!duplicated(bt$gene_id_versioned), ]
d <- merge(deg, bt, by.x = "gene", by.y = "gene_id_versioned")
det <- rowMeans(dge$counts > 0); d$detect <- det[match(d$gene, rownames(dge$counts))]
d <- d[d$gene_type %in% c("protein_coding", "lncRNA"), ]
d$call <- as.integer(d$padj < 0.05 & abs(d$logFC) > 0.5); d$lnc <- as.integer(d$gene_type == "lncRNA")
show <- function(f, lab) { m <- glm(f, binomial, d); ci <- exp(confint.default(m)["lnc", ]); cat(lab, "OR", round(exp(coef(m)["lnc"]), 3), "CI", round(ci, 3), "\n") }
show(call ~ lnc, "unadjusted")
show(call ~ lnc + ns(AveExpr, 5), "AveExpr")
show(call ~ lnc + ns(AveExpr, 5) + ns(detect, 5), "AveExpr+detection")
show(call ~ lnc + AveExpr + detect, "linear AveExpr+detection")
cat("detection summary by type\n"); print(tapply(d$detect, d$gene_type, summary))
writeLines(capture.output(sessionInfo()), paste0(O, "sessionInfo_17.txt"))
