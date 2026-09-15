#!/usr/bin/env Rscript
# F2/F3 step 1b: hg19 -> hg38 one-base liftover for the 1000G EUR common-SNV null panel.
# Same rule and same chain file as scripts/analysis/alphagenome_atlas/03_liftover_universeB.R:
# per position, the number of liftover mappings and the FIRST mapping's chromosome and position.
suppressPackageStartupMessages({ library(data.table); library(GenomicRanges); library(rtracklayer) })

args <- commandArgs(trailingOnly = TRUE)
in_path <- args[1]   # TSV with columns chrom_hg19 (with chr prefix), pos_hg19
out_path <- args[2]
chain_path <- args[3]
if (file.exists(out_path)) stop("refusing to overwrite: ", out_path)

x <- fread(in_path)
gr <- GRanges(seqnames = x$chrom_hg19, ranges = IRanges(start = x$pos_hg19, width = 1L))
lifted <- liftOver(gr, import.chain(chain_path))
n_map <- as.integer(elementNROWS(lifted))
flat <- unlist(lifted, use.names = FALSE)
first <- head(cumsum(c(1L, n_map)), -1L)
ok <- n_map > 0L
hg38_chrom <- rep(NA_character_, length(lifted))
hg38_pos <- rep(NA_integer_, length(lifted))
if (any(ok)) {
  hg38_chrom[ok] <- as.character(seqnames(flat))[first[ok]]
  hg38_pos[ok] <- start(flat)[first[ok]]
}
out <- data.table(chrom_hg19 = x$chrom_hg19, pos_hg19 = x$pos_hg19,
                  n_liftover_mappings = n_map, hg38_chrom = hg38_chrom, hg38_pos = hg38_pos)
fwrite(out, out_path, sep = "\t")
message(sprintf("liftover %s: %d positions, %d mapped once, %d unmapped, %d multimapped",
                basename(in_path), nrow(out), sum(n_map == 1L), sum(n_map == 0L), sum(n_map > 1L)))
