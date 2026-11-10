# Why do 110 COLOC top variants get no locateVariants class? Check overlap with transcripts,
# exons of transcripts with and without CDS, and genes in UCSC hg19 knownGene. No randomness.
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(VariantAnnotation)})
U <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15/S2K"
cmp <- fread(file.path(U, "S2K_variant_class_compare.csv"))
x <- cmp[is.na(fine_class) | fine_class == ""]
cat("unclassified variants:", nrow(x), "; producer classes:", paste(names(table(x$ann_class)), table(x$ann_class)), "\n")
x[, c("ch", "pos") := tstrsplit(variant_key, ":")]
gr <- GRanges(paste0("chr", x$ch), IRanges(as.integer(x$pos), width = 1))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
cds_tx <- unique(as.integer(names(cdsBy(txdb, by = "tx"))))
ex <- exonsBy(txdb, by = "tx")
ex_nc <- unlist(ex[!(as.integer(names(ex)) %in% cds_tx)])
ex_cd <- unlist(ex[as.integer(names(ex)) %in% cds_tx])
g <- suppressMessages(genes(txdb, single.strand.genes.only = TRUE))
x[, in_tx := overlapsAny(gr, tx, ignore.strand = TRUE)]
x[, tx_has_gene_id := sapply(seq_along(gr), function(i) {h <- tx[overlapsAny(tx, gr[i], ignore.strand = TRUE)]; any(lengths(h$gene_id) > 0)})]
x[, in_noncoding_tx_exon := overlapsAny(gr, ex_nc, ignore.strand = TRUE)]
x[, in_coding_tx_exon := overlapsAny(gr, ex_cd, ignore.strand = TRUE)]
x[, in_gene := overlapsAny(gr, g, ignore.strand = TRUE)]
print(x[, .N, by = .(in_tx, tx_has_gene_id, in_noncoding_tx_exon, in_coding_tx_exon, in_gene)])
# re-run locateVariants on these variants with default IntergenicVariants to see what is returned
loc <- suppressWarnings(locateVariants(gr, txdb, AllVariants(promoter = PromoterVariants(upstream = 2000L, downstream = 200L))))
cat("locateVariants rows for these variants:", length(loc), "\n")
loc2 <- suppressWarnings(locateVariants(gr, txdb, IntergenicVariants()))
cat("IntergenicVariants rows:", length(loc2), "\n")
fwrite(x, file.path(U, "S2K_unclassified_overlap.csv"))
