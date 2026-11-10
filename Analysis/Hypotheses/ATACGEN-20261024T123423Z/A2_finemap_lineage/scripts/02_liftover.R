#!/usr/bin/env Rscript
# A2 step 2: hg19 -> hg38 liftover of unique window positions; keep 1:1 same-chromosome maps.
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(rtracklayer)})
W <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage/work"
CHAIN <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/broadaway_eqtl/hg19ToHg38.over.chain"
b <- fread(file.path(W, "lift_in.bed"), header = FALSE, col.names = c("chrom", "start0", "end"))
gr <- GRanges(b$chrom, IRanges(b$end, width = 1))
ch <- import.chain(CHAIN)
res <- liftOver(gr, ch)
n <- as.integer(elementNROWS(res))
flat <- unlist(res, use.names = FALSE)
first <- head(cumsum(c(1L, n)), -1L)
one <- n == 1L
out <- data.table(chrom = b$chrom, pos19 = b$end, n_map = n,
                  chrom38 = NA_character_, pos38 = NA_integer_)
out[one, chrom38 := as.character(seqnames(flat))[first[one]]]
out[one, pos38 := start(flat)[first[one]]]
out[, status := fifelse(n == 0L, "unmapped", fifelse(n > 1L, "multi",
             fifelse(chrom38 != chrom, "other_chrom", "mapped")))]
print(out[, .N, by = status])
fwrite(out[, .N, by = status], file.path(W, "liftover_status_counts.tsv"), sep = "\t")
fwrite(out[status == "mapped", .(chrom, pos19, pos38)], file.path(W, "lift_out.tsv.gz"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(W, "..", "sessionInfo_liftover.txt"))
