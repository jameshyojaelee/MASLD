#!/usr/bin/env Rscript
# Verifier C step 1: hg19 knownGene features for direct-overlap classification.
suppressPackageStartupMessages({
  library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg19.knownGene)
})
t0 <- Sys.time(); tick <- function(m) cat(sprintf("[%5.1fs] %s\n", as.numeric(Sys.time() - t0, units = "secs"), m))
OUT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/verify/C/out"
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx <- transcripts(txdb, columns = c("TXID", "TXNAME", "GENEID")); tick("transcripts")
tx$GENEID <- vapply(tx$GENEID, function(g) if (length(g)) g[1] else NA_character_, "")
cds_tx   <- cdsBy(txdb, "tx"); tick("cdsBy")
introns  <- unlist(intronsByTranscript(txdb)); tick("introns")
introns  <- introns[width(introns) >= 1]
ss_gr    <- c(resize(introns, pmin(2L, width(introns)), fix = "start", ignore.strand = TRUE),
              resize(introns, pmin(2L, width(introns)), fix = "end",   ignore.strand = TRUE))
five_gr  <- unlist(fiveUTRsByTranscript(txdb)); tick("5'UTR")
three_gr <- unlist(threeUTRsByTranscript(txdb)); tick("3'UTR")
prom_gr  <- suppressWarnings(trim(promoters(tx, upstream = 2000, downstream = 200)))
exon_gr  <- exons(txdb); tick("exons")
saveRDS(list(tx = tx, cds_tx = cds_tx, cds_gr = unlist(cds_tx), introns = introns, ss_gr = ss_gr,
             five_gr = five_gr, three_gr = three_gr, prom_gr = prom_gr, exon_gr = exon_gr),
        file.path(OUT, "features.rds"))
tick("saved")
