suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(VariantAnnotation)})
U <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15/S2K"
x <- fread(file.path(U, "S2K_unclassified_overlap.csv"))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr <- GRanges(paste0("chr", x$ch), IRanges(as.integer(x$pos), width = 1))
tx <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
intr <- unlist(intronsByTranscript(txdb))
x[, in_intron := overlapsAny(gr, intr, ignore.strand = TRUE)]
print(x[, .N, by = .(in_tx, in_intron, in_noncoding_tx_exon)])
i <- which(x$in_tx)[1]; print(x[i]); h <- tx[overlapsAny(tx, gr[i], ignore.strand = TRUE)]; print(h)
for (j in 1:3) { lv <- suppressWarnings(locateVariants(gr[j], txdb, IntronVariants())); cat(x$variant_key[j], "IntronVariants rows:", length(lv), "\n") }
lv <- suppressWarnings(locateVariants(gr, txdb, IntronVariants())); cat("IntronVariants rows on all 110:", length(lv), "\n")
gr2 <- gr; seqlevels(gr2) <- seqlevels(gr2)
cat("seqlevels:", paste(head(seqlevels(gr)), collapse=","), " isCircular:", paste(head(isCircular(gr)), collapse=","), "\n")
