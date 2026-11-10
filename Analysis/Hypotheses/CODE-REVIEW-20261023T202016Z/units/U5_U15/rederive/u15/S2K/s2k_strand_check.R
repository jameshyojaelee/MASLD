# For producer-intronic and producer-intergenic-but-intronic variants, record the strands of the
# knownGene introns that contain them, and the transcript lengths. No randomness.
suppressPackageStartupMessages({library(data.table); library(GenomicRanges); library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)})
U <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15/S2K"
x <- fread(file.path(U, "S2K_direct_overlap_classes.csv"))
x <- x[direct_class == "intron"]
gr <- GRanges(paste0("chr", x$ch), IRanges(as.integer(x$pos), width = 1))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ib <- intronsByTranscript(txdb); intr <- unlist(ib)
h <- findOverlaps(gr, intr, ignore.strand = TRUE)
d <- data.table(i = queryHits(h), strand = as.character(strand(intr))[subjectHits(h)])
s <- d[, .(strands = paste(sort(unique(strand)), collapse = "")), by = i]
x[s$i, intron_strands := s$strands]
print(dcast(x[, .N, by = .(ann_class, intron_strands)], ann_class ~ intron_strands, fill = 0))
