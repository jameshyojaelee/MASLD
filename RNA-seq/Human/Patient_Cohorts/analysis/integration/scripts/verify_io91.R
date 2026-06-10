#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
RDIR <- "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"
strip <- function(x) sub("\\.[0-9]+$","",x)

can <- fread(file.path(RDIR,"canonical_deg_results.csv")); can[, gid := strip(gene)]
io  <- fread(file.path(RDIR,"integration_only_C2/A2_integration_only_genes_C2.csv"))

cat("== IO file:", nrow(io), "integration-only genes (n_support==0) ==\n")
io91 <- io[abs(logFC) > 0.5]
cat("IO with |logFC|>0.5 :", nrow(io91), "\n\n")

m  <- match(io91$gid, can$gid); cj <- can[m]
cat("[1] present in canonical_deg_results.csv:", sum(!is.na(m)), "/", nrow(io91), "\n")
is_t1 <- cj$padj < 0.05 & abs(cj$logFC) > 0.5
cat("[2] qualify as canonical Tier-1 (padj<0.05 & |logFC|>0.5):", sum(is_t1, na.rm=TRUE), "/", nrow(io91), "\n")
cat(sprintf("[3] max |Δ logFC| = %.2e   max |Δ padj| = %.2e  (expect ~0)\n",
            max(abs(io91$logFC-cj$logFC),na.rm=TRUE), max(abs(io91$padj-cj$padj),na.rm=TRUE)))

t1 <- can[padj < 0.05 & abs(logFC) > 0.5]
cat("[4] canonical Tier-1 total:", nrow(t1),
    " | Tier-1 in full IO set:", sum(t1$gid %in% io$gid),
    " | identical set to io91:", setequal(t1$gid[t1$gid %in% io$gid], io91$gid), "\n\n")

ex <- c("GSTM1","HLA-DQB1","POSTN","ACACA","HNF4A","NR1H4","THRB","INSR")
exg <- can[symbol %in% ex, .(symbol, gid, logFC=round(logFC,3), padj=signif(padj,2))]
exg[, in_IO_all := gid %in% io$gid][, in_IO_91 := gid %in% io91$gid]
cat("== Manuscript Fig-1h example genes ==\n"); print(exg[order(symbol)])

io91s <- merge(io91, can[, .(gid, symbol)], by="gid", all.x=TRUE)[order(-abs(logFC))]
cat("\n== Top 25 of the 91 by |logFC| ==\n")
print(head(io91s[, .(symbol, gid, logFC=round(logFC,3), padj=signif(padj,2), n_concord)], 25))
cat("\nn_concord distribution among the 91 (sign-agreement count over 5 cohorts):\n")
print(table(io91s$n_concord))
fwrite(io91s[, .(gid, symbol, logFC, padj, n_concord)],
       file.path(RDIR,"integration_only_C2/io91_verified_against_canonical.csv"))
cat("\nwrote io91_verified_against_canonical.csv\n")
