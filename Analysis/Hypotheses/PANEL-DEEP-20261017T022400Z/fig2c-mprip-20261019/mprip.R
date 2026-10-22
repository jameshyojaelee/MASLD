suppressPackageStartupMessages({library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db)})
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr <- GRanges("chr17", IRanges(17075180, width = 1))
tx <- transcripts(txdb, columns = c("TXID", "TXNAME", "GENEID"))
h <- tx[overlapsAny(tx, gr, ignore.strand = TRUE)]
print(data.frame(txname = h$TXNAME, geneid = sapply(h$GENEID, paste, collapse = ";"), start = start(h), end = end(h), strand = as.character(strand(h))))
intr <- intronsByTranscript(txdb, use.names = TRUE); iu <- unlist(intr)
near <- iu[overlapsAny(iu, resize(gr, 11, fix = "center"), ignore.strand = TRUE)]
print(data.frame(tx = names(near), start = start(near), end = end(near), width = width(near)))
cds <- cdsBy(txdb, by = "tx", use.names = TRUE); cu <- unlist(cds)
print(unique(names(cu)[overlapsAny(cu, gr, ignore.strand = TRUE)]))
z <- AnnotationDbi::select(org.Hs.eg.db, keys = "MPRIP", keytype = "SYMBOL", columns = "ENTREZID"); print(z)
