# Shared helper (F2abce, review round 4).
# For a set of (gene symbol, lead variant hg19 chr:pos) pairs, return
#   txdb_class : strongest TxDb hg19 knownGene location class, same precedence
#                and promoter window as scripts/figures/fig2_coding_noncoding.R
#   coding_entrez / coding_symbols : genes whose CDS (coding) or splice site
#                contains the variant
#   own_gene_coding : TRUE when the COLOC gene itself is among those genes
suppressPackageStartupMessages({
  library(data.table); library(VariantAnnotation)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(GenomicRanges)
  library(org.Hs.eg.db); library(AnnotationDbi)
})

CLASS_RANK <- c(coding = 1L, spliceSite = 2L, fiveUTR = 3L,
                threeUTR = 4L, promoter = 5L, intron = 6L, intergenic = 7L)

annotate_leads <- function(gene, variant_key) {
  pairs <- unique(data.table(gene = gene, variant_key = variant_key))
  vk <- unique(pairs$variant_key)
  parts <- tstrsplit(vk, ":", fixed = TRUE)
  gr <- GRanges(seqnames = paste0("chr", parts[[1]]),
                ranges = IRanges(start = as.integer(parts[[2]]), width = 1),
                variant_key = vk)
  loc <- VariantAnnotation::locateVariants(
    gr, TxDb.Hsapiens.UCSC.hg19.knownGene,
    AllVariants(promoter = PromoterVariants(upstream = 2000, downstream = 200)))
  ld <- data.table(variant_key = mcols(gr)$variant_key[loc$QUERYID],
                   location = as.character(loc$LOCATION),
                   entrez = as.character(loc$GENEID))
  ld[, cls := fifelse(location %chin% names(CLASS_RANK), location, "intergenic")]
  ld[, rank := CLASS_RANK[cls]]
  cls <- ld[order(variant_key, rank)][, .(txdb_class = cls[1]), by = variant_key]
  cod <- ld[cls %chin% c("coding", "spliceSite") & !is.na(entrez) & entrez != "",
            .(coding_entrez = list(unique(entrez))), by = variant_key]
  out <- merge(pairs, cls, by = "variant_key", all.x = TRUE)
  out[is.na(txdb_class), txdb_class := "intergenic"]
  out <- merge(out, cod, by = "variant_key", all.x = TRUE)
  # COLOC gene symbol -> Entrez via official symbol, then alias.
  sym <- unique(out$gene)
  e1 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = sym, keytype = "SYMBOL",
                                                columns = "ENTREZID"))
  e2 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = sym, keytype = "ALIAS",
                                                columns = "ENTREZID"))
  names(e2)[1] <- "SYMBOL"
  emap <- unique(rbind(as.data.table(e1), as.data.table(e2))[!is.na(ENTREZID)])
  all_e <- unique(unlist(out$coding_entrez))
  all_e <- all_e[!is.na(all_e)]
  esym <- if (length(all_e)) as.data.table(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = all_e, keytype = "ENTREZID", columns = "SYMBOL"))) else
    data.table(ENTREZID = character(), SYMBOL = character())
  out[, own_gene_coding := mapply(function(g, ez) {
    if (is.null(ez) || !length(ez)) return(FALSE)
    any(ez %in% emap[SYMBOL == g, ENTREZID])
  }, gene, coding_entrez)]
  out[, coding_symbols := vapply(coding_entrez, function(ez) {
    if (is.null(ez) || !length(ez)) return("")
    paste(sort(unique(esym[ENTREZID %in% ez, SYMBOL])), collapse = ";")
  }, character(1))]
  out[, coding_entrez := vapply(coding_entrez, function(ez)
    if (is.null(ez)) "" else paste(ez, collapse = ";"), character(1))]
  out[]
}
