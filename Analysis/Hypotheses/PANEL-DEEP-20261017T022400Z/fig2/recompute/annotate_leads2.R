# Overlap-based hg19 consequence classes for COLOC lead variants (replaces locateVariants, which returned
# no rows for 81 of 934 lead variants, e.g. BRCA1 17:41225781 inside 13 transcripts).
# Priority: coding (CDS) or splice site (first/last 2 bp of an intron) > UTR > promoter (2 kb up, 200 bp down) > intron > intergenic.
suppressPackageStartupMessages({library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db)})
args <- commandArgs(trailingOnly = TRUE); inp <- read.delim(args[1], stringsAsFactors = FALSE); out <- args[2]
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
inp$variant <- sub("^([^:]+:[0-9]+).*", "\\1", inp$top_snp)
v <- unique(inp$variant)
gr <- GRanges(paste0("chr", sub(":.*", "", v)), IRanges(as.integer(sub(".*:", "", v)), width = 1)); names(gr) <- v
tx <- transcripts(txdb, columns = c("TXID", "GENEID")); tx2g <- setNames(as.character(unstrsplit(tx$GENEID, ";")), tx$TXID); tx2g[tx2g == ""] <- NA
feat_hits <- function(grl_or_gr, label) {
  empty <- data.frame(variant = character(0), loc = character(0), geneid = character(0), stringsAsFactors = FALSE)
  if (is(grl_or_gr, "GRangesList")) { f <- unlist(grl_or_gr); txid <- names(f) } else { f <- grl_or_gr; txid <- as.character(f$TXID) }
  h <- findOverlaps(gr, f, ignore.strand = TRUE); if (length(h) == 0) return(empty)
  data.frame(variant = names(gr)[queryHits(h)], loc = label, geneid = tx2g[txid[subjectHits(h)]], stringsAsFactors = FALSE)
}
cdsg <- cdsBy(txdb, by = "tx"); intr <- intronsByTranscript(txdb)
iu <- unlist(intr); iu <- iu[width(iu) >= 4]
spl <- c(resize(iu, 2, fix = "start"), resize(iu, 2, fix = "end"))
splr <- GRangesList(split(spl, names(spl)))
res <- rbind(feat_hits(cdsg, "coding"), feat_hits(splr, "spliceSite"),
             feat_hits(fiveUTRsByTranscript(txdb), "fiveUTR"), feat_hits(threeUTRsByTranscript(txdb), "threeUTR"),
             feat_hits(promoters(tx, upstream = 2000, downstream = 200), "promoter"), feat_hits(intr, "intron"))
rank <- c(coding = 1, spliceSite = 1, fiveUTR = 2, threeUTR = 2, promoter = 3, intron = 4)
lab <- c(`1` = "coding", `2` = "UTR", `3` = "promoter", `4` = "intron")
res$rank <- rank[res$loc]
best <- tapply(res$rank, res$variant, min)
inp$class <- ifelse(inp$variant %in% names(best), lab[as.character(best[inp$variant])], "intergenic")
cod <- res[res$rank == 1, ]
codg <- tapply(cod$geneid, cod$variant, function(x) paste(sort(unique(na.omit(x))), collapse = ";"))
inp$coding_geneids <- ifelse(inp$variant %in% names(codg), codg[inp$variant], "")
e2z <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(inp$ensembl), keytype = "ENSEMBL", columns = "ENTREZID"))
s2z <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(inp$symbol), keytype = "SYMBOL", columns = "ENTREZID"))
own <- mapply(function(ens, sym, ids) {
  z <- unique(c(e2z$ENTREZID[e2z$ENSEMBL == ens], s2z$ENTREZID[s2z$SYMBOL == sym])); z <- z[!is.na(z)]
  length(intersect(z, strsplit(ids, ";")[[1]])) > 0 }, inp$ensembl, inp$symbol, inp$coding_geneids)
inp$coding_owner <- ifelse(inp$class == "coding", ifelse(own, "this_gene", "other_gene"), "")
zs <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(na.omit(cod$geneid)), keytype = "ENTREZID", columns = "SYMBOL"))
inp$coding_gene_symbols <- sapply(inp$coding_geneids, function(ids) if (ids == "") "" else paste(zs$SYMBOL[match(strsplit(ids, ";")[[1]], zs$ENTREZID)], collapse = ";"))
inp$splice_only <- inp$variant %in% setdiff(res$variant[res$loc == "spliceSite"], res$variant[res$loc == "coding"])
write.table(inp, out, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), sub("\\.tsv$", "_sessionInfo.txt", out))
cat("classes:\n"); print(table(inp$class)); print(table(inp$coding_owner)); print(sum(inp$splice_only))
