# Second, locateVariants-free classification of the 1,263 S2K COLOC top variants.
# Direct overlaps with UCSC hg19 knownGene features, strand ignored (variants are unstranded):
# coding (CDS) > spliceSite (first/last 2 nt of an intron) > 5'UTR > 3'UTR >
# promoter (-2 kb/+200 bp of any transcript start) > intron > exon of a transcript without CDS
# ("ncRNA_exon", not a class in the legend) > intergenic.
# Also reruns locateVariants with ignore.strand = TRUE for comparison. No randomness.
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(VariantAnnotation)})
U <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15/S2K"
cmp <- fread(file.path(U, "S2K_variant_class_compare.csv"))
stopifnot(nrow(cmp) == 1263L)
cmp[, c("ch", "pos") := tstrsplit(variant_key, ":")]
gr <- GRanges(paste0("chr", cmp$ch), IRanges(as.integer(cmp$pos), width = 1))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

cds <- unlist(cdsBy(txdb, by = "tx"))
cds_tx <- unique(as.integer(names(cdsBy(txdb, by = "tx"))))
intr <- unlist(intronsByTranscript(txdb))
intr <- intr[width(intr) >= 4]
ss <- c(resize(intr, 2, fix = "start"), resize(intr, 2, fix = "end"))
u5 <- unlist(fiveUTRsByTranscript(txdb)); u3 <- unlist(threeUTRsByTranscript(txdb))
prom <- suppressWarnings(promoters(transcripts(txdb), upstream = 2000, downstream = 200))
ex <- exonsBy(txdb, by = "tx")
ex_nc <- unlist(ex[!(as.integer(names(ex)) %in% cds_tx)])

ov <- function(f) overlapsAny(gr, f, ignore.strand = TRUE)
flags <- data.table(coding = ov(cds), spliceSite = ov(ss), fiveUTR = ov(u5), threeUTR = ov(u3),
                    promoter = ov(prom), intron = ov(intr), ncRNA_exon = ov(ex_nc))
lev <- c(names(flags), "intergenic")
cmp[, direct_class := apply(as.matrix(flags), 1, function(r) if (any(r)) lev[which(r)[1]] else "intergenic")]

loc <- suppressWarnings(locateVariants(gr, txdb,
  AllVariants(promoter = PromoterVariants(upstream = 2000L, downstream = 200L)), ignore.strand = TRUE))
lv <- c("coding", "spliceSite", "fiveUTR", "threeUTR", "promoter", "intron", "intergenic")
ld <- data.table(q = loc$QUERYID, cls = as.character(loc$LOCATION))[cls %in% lv]
ld <- ld[, .(cls = lv[min(match(cls, lv))]), by = q]
cmp[, lv_ignore_strand := "no row"]
cmp[ld$q, lv_ignore_strand := ld$cls]

cat("Printed panel classes (producer annotation):\n"); print(cmp[, .N, by = ann_class][order(-N)])
cat("Direct-overlap classes:\n"); print(cmp[, .(N = .N, pct = round(100 * .N / nrow(cmp), 1)), by = direct_class][order(-N)])
cat("locateVariants ignore.strand=TRUE classes:\n"); print(cmp[, .N, by = lv_ignore_strand][order(-N)])
cat("Producer class x direct class:\n"); print(dcast(cmp[, .N, by = .(ann_class, direct_class)], ann_class ~ direct_class, fill = 0))
fwrite(cmp, file.path(U, "S2K_direct_overlap_classes.csv"))
