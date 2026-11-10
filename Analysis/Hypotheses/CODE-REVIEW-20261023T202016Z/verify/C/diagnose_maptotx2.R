#!/usr/bin/env Rscript
# Verifier C step 3b: does the miss depend on which transcripts are in the subject list?
suppressPackageStartupMessages({ library(data.table); library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg19.knownGene) })
OUT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/verify/C/out"
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ibt <- intronsByTranscript(txdb)
miss <- fread(file.path(OUT, "S2K_no_locate_row_variants.tsv"))[fine_class == "intron"]
gm <- GRanges(paste0("chr", miss$chr), IRanges(miss$pos_hg19, width = 1))
ex <- GRanges("chr10", IRanges(101953705, width = 1))
hit_tx <- names(ibt)[overlapsAny(ibt, ex)]
cat("example: intron lists containing it:", paste(hit_tx, collapse = ","), "\n")
cat("full list                  :", length(mapToTranscripts(ex, ibt)), "\n")
cat("only its transcripts       :", length(mapToTranscripts(ex, ibt[hit_tx])), "\n")
on10 <- ibt[any(seqnames(ibt) == "chr10")]
cat("chr10 transcripts only     :", length(mapToTranscripts(ex, on10)), "\n")
nonempty <- ibt[lengths(ibt) > 0]
cat("non-empty intron lists only:", length(mapToTranscripts(ex, nonempty)), "\n")
# neighbourhood: transcripts overlapping +-1 Mb
win <- ibt[overlapsAny(ibt, resize(ex, 2e6, fix = "center"))]
cat("transcripts within 1 Mb    :", length(mapToTranscripts(ex, win)), " (n tx =", length(win), ")\n")
# all 99 intron no-row variants
cat("\nall 99: full list hits =", length(unique(mapToTranscripts(gm, ibt)$xHits)),
    "; per-variant own-transcript lists hits =",
    sum(vapply(seq_along(gm), function(i) length(mapToTranscripts(gm[i], ibt[overlapsAny(ibt, gm[i])])) > 0, TRUE)), "\n")
# overlap count of the unlisted introns at each missed variant
u <- unlist(ibt)
cat("findOverlaps(within) hits per missed variant (summary):\n"); print(summary(countOverlaps(gm, u, type = "within")))
# do hits disappear when there are many hits in total? compare with located intron variants
