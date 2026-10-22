#!/usr/bin/env Rscript
# Blind check F2abce, step 3: hg19 consequence and coding-gene ownership of lead variants.
suppressPackageStartupMessages({
  library(data.table); library(VariantAnnotation)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db); library(AnnotationDbi)
})
set.seed(20260923)
args <- commandArgs(TRUE); CK <- args[1]; OUT <- file.path(CK, "out")
va <- fread(file.path(OUT, "variants_to_annotate.tsv"))
v <- unique(va[, .(top_snp)])
v[, c("chr", "pos") := tstrsplit(top_snp, ":", fixed = TRUE, keep = 1:2, type.convert = TRUE)]
gr <- GRanges(paste0("chr", v$chr), IRanges(v$pos, v$pos), vkey = v$top_snp)
loc <- locateVariants(gr, TxDb.Hsapiens.UCSC.hg19.knownGene,
                      AllVariants(promoter = PromoterVariants(upstream = 2000, downstream = 200)))
L <- data.table(vkey = mcols(gr)$vkey[loc$QUERYID], location = as.character(loc$LOCATION),
                geneid = as.character(loc$GENEID))
rank <- c(coding = 1L, spliceSite = 2L, fiveUTR = 3L, threeUTR = 4L, promoter = 5L, intron = 6L, intergenic = 7L)
L[, location := fifelse(location %chin% names(rank), location, "intergenic")]
cls <- L[, .(fine_class = names(rank)[min(rank[location])]), by = vkey]
cod <- L[location %chin% c("coding", "spliceSite") & !is.na(geneid), .(geneid = unique(geneid)), by = vkey]
sym <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(cod$geneid), keytype = "ENTREZID",
                                              columns = "SYMBOL"))
cod <- merge(cod, as.data.table(sym), by.x = "geneid", by.y = "ENTREZID", all.x = TRUE)
codagg <- cod[, .(coding_entrez = paste(sort(unique(geneid)), collapse = ";"),
                  coding_symbols = paste(sort(unique(SYMBOL)), collapse = ";")), by = vkey]
va <- merge(va, cls, by.x = "top_snp", by.y = "vkey", all.x = TRUE)
va[is.na(fine_class), fine_class := "intergenic"]
va <- merge(va, codagg, by.x = "top_snp", by.y = "vkey", all.x = TRUE)
va[is.na(coding_symbols), `:=`(coding_symbols = "", coding_entrez = "")]
# gene -> Entrez through official symbol and through alias (robust to symbol drift)
g2e <- rbind(
  suppressMessages(as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = unique(va$gene), keytype = "SYMBOL",
                                                     columns = "ENTREZID")))[, .(gene = SYMBOL, ENTREZID)],
  suppressMessages(as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys = intersect(unique(va$gene), keys(org.Hs.eg.db, "ALIAS")),
                                                     keytype = "ALIAS", columns = "ENTREZID")))[, .(gene = ALIAS, ENTREZID)])
g2e <- unique(g2e[!is.na(ENTREZID)])
va[, own_symbol := mapply(function(g, s) g %in% strsplit(s, ";")[[1]], gene, coding_symbols)]
va[, own_entrez := mapply(function(g, e) length(intersect(g2e[gene == g, ENTREZID], strsplit(e, ";")[[1]])) > 0,
                          gene, coding_entrez)]
va[, is_coding := fine_class %chin% c("coding", "spliceSite")]
va[, class3 := fifelse(!is_coding, "noncoding", fifelse(own_symbol | own_entrez, "coding_this_gene", "coding_other_gene"))]
fwrite(va, file.path(OUT, "lead_annotation.tsv"), sep = "\t")
for (p in unique(va$panel)) {
  cat("==", p, "\n"); print(va[panel == p, .N, by = class3])
  print(va[panel == p, .N, by = fine_class][order(-N)])
  cat("symbol-vs-entrez ownership disagreements:", va[panel == p & is_coding & own_symbol != own_entrez, .N], "\n")
}
cat("\n2C coding genes:\n"); print(va[panel == "2C" & is_coding, .(gene, top_snp, fine_class, coding_symbols, own_symbol, own_entrez)][order(own_symbol, top_snp)])
cat("2C coding distinct leads:", va[panel == "2C" & is_coding, uniqueN(top_snp)], "\n")
cat("\n2B susie coding points:\n"); print(va[panel == "2B_susie" & is_coding, .(gene, top_snp, fine_class, coding_symbols, class3)])
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo_annotate.txt"))
