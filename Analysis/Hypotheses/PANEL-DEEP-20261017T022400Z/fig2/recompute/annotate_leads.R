# Independent hg19 consequence annotation of COLOC lead variants (Fig 2B colour, Fig 2C partition).
# TxDb.Hsapiens.UCSC.hg19.knownGene + VariantAnnotation::locateVariants (AllVariants: promoter 2 kb up / 200 bp down).
# Input: TSV with columns ensembl, symbol, top_snp (chr:pos[:...]); output adds per-variant classes and whether
# a coding/splice hit belongs to the gene itself (Entrez via org.Hs.eg.db from Ensembl ID, then symbol).
suppressPackageStartupMessages({
  library(VariantAnnotation); library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db)
})
args <- commandArgs(trailingOnly = TRUE)
inp <- read.delim(args[1], stringsAsFactors = FALSE)
out <- args[2]
set.seed(20261017)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
v <- unique(sub("^([^:]+:[0-9]+).*", "\\1", inp$top_snp))
chr <- paste0("chr", sub(":.*", "", v)); pos <- as.integer(sub(".*:", "", v))
gr <- GRanges(chr, IRanges(pos, width = 1)); names(gr) <- v
gr <- keepSeqlevels(gr, intersect(seqlevels(gr), seqlevels(txdb)), pruning.mode = "coarse")
loc <- suppressWarnings(locateVariants(gr, txdb, AllVariants()))
df <- data.frame(variant = names(gr)[loc$QUERYID], location = as.character(loc$LOCATION),
                 geneid = as.character(loc$GENEID), stringsAsFactors = FALSE)
rank <- c(coding = 1, spliceSite = 1, fiveUTR = 2, threeUTR = 2, promoter = 3, intron = 4, intergenic = 5)
lab <- c(`1` = "coding", `2` = "UTR", `3` = "promoter", `4` = "intron", `5` = "intergenic")
df$rank <- rank[df$location]
best <- tapply(df$rank, df$variant, min)
codgenes <- tapply(df$geneid[df$rank == 1], df$variant[df$rank == 1], function(x) paste(sort(unique(na.omit(x))), collapse = ";"))
# map gene Ensembl -> Entrez, and symbol -> Entrez as fallback
e2z <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(inp$ensembl), keytype = "ENSEMBL", columns = "ENTREZID"))
s2z <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(inp$symbol), keytype = "SYMBOL", columns = "ENTREZID"))
inp$variant <- sub("^([^:]+:[0-9]+).*", "\\1", inp$top_snp)
inp$class <- lab[as.character(best[inp$variant])]
inp$class[is.na(inp$class)] <- "unannotated"
inp$coding_geneids <- ifelse(inp$variant %in% names(codgenes), codgenes[inp$variant], "")
own <- mapply(function(ens, sym, ids) {
  z <- unique(c(e2z$ENTREZID[e2z$ENSEMBL == ens], s2z$ENTREZID[s2z$SYMBOL == sym]))
  z <- z[!is.na(z)]
  ids <- strsplit(ids, ";")[[1]]
  length(intersect(z, ids)) > 0
}, inp$ensembl, inp$symbol, inp$coding_geneids)
inp$coding_owner <- ifelse(inp$class == "coding", ifelse(own, "this_gene", "other_gene"), "")
zsym <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(na.omit(df$geneid)), keytype = "ENTREZID", columns = "SYMBOL"))
inp$coding_gene_symbols <- sapply(inp$coding_geneids, function(ids) {
  if (ids == "") return(""); paste(zsym$SYMBOL[match(strsplit(ids, ";")[[1]], zsym$ENTREZID)], collapse = ";")
})
write.table(inp, out, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), sub("\\.tsv$", "_sessionInfo.txt", out))
cat("classes:\n"); print(table(inp$class)); print(table(inp$coding_owner))
