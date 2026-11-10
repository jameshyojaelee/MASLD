#!/usr/bin/env Rscript
# Verifier C step 3: why IntronVariants (mapToTranscripts over intronsByTranscript) misses variants.
suppressPackageStartupMessages({ library(data.table); library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg19.knownGene) })
OUT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/verify/C/out"
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
ibt <- intronsByTranscript(txdb)
miss <- fread(file.path(OUT, "S2K_no_locate_row_variants.tsv"))
gm <- GRanges(paste0("chr", miss$chr), IRanges(miss$pos_hg19, width = 1))
cat("no-row variants by direct class:\n"); print(table(miss$fine_class))
ui <- unlist(ibt)
h <- findOverlaps(gm, ui)
d <- data.table(q = queryHits(h), tx = names(ui)[subjectHits(h)], strand = as.character(strand(ui))[subjectHits(h)])
d <- unique(d)
cat("\nintron-overlapping tx for no-row variants: strand table\n"); print(d[, table(strand)])
# per transcript: is the intron list ordered consistently, and does mapToTranscripts hit?
chk <- d[, {
  g <- ibt[[tx]]; s <- start(g)
  .(n_introns = length(g), ascending = !is.unsorted(s), descending = !is.unsorted(rev(s)),
    any_overlap_within_tx = any(countOverlaps(g, g) > 1L),
    m2t_hits = length(mapToTranscripts(gm[q], ibt[tx])))
}, by = .(q, tx, strand)]
print(chk[, .N, by = .(strand, ascending, descending, any_overlap_within_tx, m2t_hits)])
# control: located intron variants -> their transcripts
ex <- GRanges("chr10", IRanges(101953705, width = 1))
id <- names(ibt)[any(overlapsAny(ibt, ex))]
cat("\nexample tx ids:", id, "\n"); g <- ibt[[id[1]]]; print(g)
cat("mapToTranscripts as stored:", length(mapToTranscripts(ex, ibt[id[1]])), "\n")
gr_rev <- GRangesList(rev(g)); names(gr_rev) <- id[1]
cat("mapToTranscripts with introns reversed:", length(mapToTranscripts(ex, gr_rev)), "\n")
cat("findOverlaps type=within:", length(findOverlaps(ex, ibt[id[1]], type = "within")), "\n")
# all S2K-located intron hits: strand mix of their transcripts
fwrite(chk, file.path(OUT, "diagnose_maptotx.tsv"), sep = "\t")
