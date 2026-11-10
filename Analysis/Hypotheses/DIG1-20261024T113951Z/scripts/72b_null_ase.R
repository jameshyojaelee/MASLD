# DIG1 job 72b. Reference-bias null for job 71b P1. Same het rule as 71b (depth >= 10, 0.1 < alt fraction < 0.9; reads
# >= 20). Per null SNV with >= 20 het samples: median and pooled alt fraction. PANX1 statistic: 3'UTR pooled alt fraction
# (0.454, 99 hets). Empirical p = (1 + #null pooled <= PANX1) / (1 + n null), one-sided (direction fixed by the eQTL).
# Same restricted to null SNVs of the same base-change classes as the PANX1 sites (G>C, C>G, T>G, A>C).
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); OUT <- a[1]
A <- fread(file.path(OUT, "allele_counts.tsv"), header = FALSE, col.names = c("s", "chr", "pos", "ref", "alt", "ad"))
T <- fread(file.path(OUT, "targets_annot.tsv"), header = FALSE, col.names = c("chr", "pos", "pos2", "gene", "eref", "ealt", "af"))
A <- merge(A, T[, .(chr, pos, gene, eref, ealt, af)], by = c("chr", "pos"))
A <- A[ref == eref]
cnt <- function(ad, alt, ealt) { v <- as.numeric(strsplit(ad, ",")[[1]]); al <- strsplit(alt, ",")[[1]]; k <- match(ealt, al); c(v[1], if (is.na(k)) 0 else v[k + 1]) }
cc <- t(mapply(cnt, A$ad, A$alt, A$ealt)); A[, `:=`(nref = cc[, 1], nalt = cc[, 2])][, dp := nref + nalt][, f := nalt / dp]
H <- A[dp >= 20 & f > 0.1 & f < 0.9]
N <- H[, .(n_het = .N, median_f = median(f), pooled_f = sum(nalt) / sum(dp), cls = paste0(eref[1], ">", ealt[1])), by = .(chr, pos, gene)][n_het >= 20]
fwrite(N, file.path(OUT, "null_site_summary.tsv"), sep = "\t")
U <- fread(a[2]); px <- sum(U$nalt) / sum(U$nref + U$nalt)
CL <- c("G>C", "C>G", "T>G", "A>C")
R <- data.table(set = c("all null SNVs", "same base-change classes"),
  n_sites = c(nrow(N), sum(N$cls %in% CL)),
  null_median_of_pooled = c(median(N$pooled_f), median(N[cls %in% CL]$pooled_f)),
  null_q05 = c(quantile(N$pooled_f, 0.05), quantile(N[cls %in% CL]$pooled_f, 0.05)),
  panx1_pooled = px,
  p_emp = c((1 + sum(N$pooled_f <= px)) / (1 + nrow(N)), (1 + sum(N[cls %in% CL]$pooled_f <= px)) / (1 + sum(N$cls %in% CL))))
fwrite(R, file.path(OUT, "panx1_vs_null.tsv"), sep = "\t"); print(R)
print(summary(N$pooled_f)); print(N[, .(.N, med = median(pooled_f)), by = cls][order(-N)])
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
