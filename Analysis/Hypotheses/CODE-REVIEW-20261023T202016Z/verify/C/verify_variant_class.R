#!/usr/bin/env Rscript
# Verifier C: Figure S2K and Figure 2C lead-variant classes.
# (1) reproduce the producer rule (locateVariants AllVariants, no row -> intergenic)
# (2) classify by direct overlap with hg19 knownGene features written here
# (3) diagnose why locateVariants returns no row for some variants.
# No randomness in this script.
suppressPackageStartupMessages({
  library(data.table)
  library(VariantAnnotation)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
})
set.seed(1)
R   <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT <- file.path(R, "Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/verify/C/out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ANN   <- file.path(R, "figures/candidates/fig2-coloc-r2-20260924T053244Z/variant_class/coloc_variant_annotation.csv")
COLOC <- file.path(R, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv")
TIER  <- file.path(R, "GWAS/finemapping/config/gwas_trait_tier.tsv")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

main <- fread(TIER)[placement == "main", study_name]
cat("main strata:", length(main), "\n")

# ---------------------------------------------------------------- features
# built by build_features.R from the same TxDb (transcripts, cdsBy, intronsByTranscript,
# five/threeUTRsByTranscript, promoters(-2000/+200), exons)
F <- readRDS(file.path(OUT, "features.rds"))
tx <- F$tx; cds_gr <- F$cds_gr; cds_tx <- F$cds_tx; introns <- F$introns; ss_gr <- F$ss_gr
five_gr <- F$five_gr; three_gr <- F$three_gr; prom_gr <- F$prom_gr; exon_gr <- F$exon_gr
MODE <- commandArgs(TRUE)[1]
lv_cache <- new.env(parent = emptyenv())   # shared locateVariants cache; does not change results
t0 <- Sys.time(); tick <- function(m) cat(sprintf("[%5.1fs] %s\n", as.numeric(Sys.time() - t0, units = "secs"), m))

direct_class <- function(gr) {
  hit <- function(f) overlapsAny(gr, f, ignore.strand = TRUE)
  h <- data.table(coding = hit(cds_gr), spliceSite = hit(ss_gr), fiveUTR = hit(five_gr),
                  threeUTR = hit(three_gr), promoter = hit(prom_gr), intron = hit(introns),
                  nc_exon = hit(exon_gr), in_tx = hit(tx))
  h[, fine_class := fcase(coding, "coding", spliceSite, "spliceSite", fiveUTR, "fiveUTR",
                          threeUTR, "threeUTR", promoter, "promoter", intron, "intron",
                          nc_exon, "nc_exon_only", default = "intergenic")]
  h
}

rank <- c(coding = 1, spliceSite = 2, fiveUTR = 3, threeUTR = 4, promoter = 5, intron = 6, intergenic = 7)
producer_class <- function(gr, keys) {
  loc <- locateVariants(gr, txdb,
           AllVariants(promoter = PromoterVariants(upstream = 2000, downstream = 200)), cache = lv_cache)
  d <- data.table(vk = keys[loc$QUERYID], cls = as.character(loc$LOCATION))
  d[!cls %in% names(rank), cls := "intergenic"]
  d[, r := rank[cls]]
  setorder(d, vk, r)
  d <- d[!duplicated(vk)]
  res <- data.table(vk = keys, any_row = keys %in% d$vk)
  res <- merge(res, d[, .(vk, producer = cls)], by = "vk", all.x = TRUE, sort = FALSE)
  res[is.na(producer), producer := "intergenic"]
  res[match(keys, vk)]
}

mk_gr <- function(chr, pos) GRanges(paste0("chr", chr), IRanges(pos, width = 1))
tab <- function(x, lev) { t <- table(factor(x, levels = lev)); data.table(class = lev, n = as.integer(t)) }

if (MODE == "s2k") {
# ================================================================= S2K
a <- fread(ANN)
a <- a[study %in% main & !is.na(gene_symbol) & gene_symbol != ""]
cat(sprintf("S2K input: %d gene-GWAS pairs, %d genes, %d unique lead variants\n",
            nrow(a), uniqueN(a$gene_symbol), uniqueN(a$variant_key)))

# independent rebuild of the lead set from the r2 release
cc <- fread(COLOC, select = c("gwas_name", "gene", "PP.H4.susie", "PP.H4.abf", "top_snp", "top_snp_PP"))
cc <- cc[gwas_name %in% main & !is.na(gene) & gene != "" & !is.na(top_snp) & top_snp != ""]
cc[, pp4 := suppressWarnings(pmax(as.numeric(PP.H4.susie), as.numeric(PP.H4.abf), na.rm = TRUE))]
cc <- cc[is.finite(pp4) & pp4 > 0.5]
setorder(cc, gene, gwas_name, -pp4, -top_snp_PP)
cc <- cc[!duplicated(cc, by = c("gene", "gwas_name"))]
cat(sprintf("Rebuilt from r2: %d pairs, %d genes, %d unique variants; pair-variant set identical to S2K input: %s\n",
            nrow(cc), uniqueN(cc$gene), uniqueN(cc$top_snp),
            identical(sort(paste(cc$gene, cc$gwas_name, cc$top_snp)),
                      sort(paste(a$gene_symbol, a$study, a$variant_key)))))

v <- unique(a[, .(variant_key, chr, pos_hg19, fine_class_file = fine_class)])
stopifnot(!anyDuplicated(v$variant_key))
gr <- mk_gr(v$chr, v$pos_hg19)
tick("start producer"); pc <- producer_class(gr, v$variant_key); tick("producer done")
dc <- direct_class(gr)
v[, producer := pc$producer]
v[, any_locate_row := pc$any_row]
v <- cbind(v, dc)
cat("Producer rule re-run equals fine_class in S2K input:", all(v$producer == v$fine_class_file), "\n")

lev <- c("coding", "spliceSite", "fiveUTR", "threeUTR", "promoter", "intron", "intergenic", "nc_exon_only")
s2k <- merge(tab(v$fine_class_file, lev)[, .(class, printed_rule_n = n)],
             tab(v$fine_class, lev)[, .(class, direct_n = n)], by = "class", sort = FALSE)
s2k[, `:=`(printed_rule_pct = round(100 * printed_rule_n / nrow(v), 1),
           direct_pct = round(100 * direct_n / nrow(v), 1))]
cat("\n--- S2K: unique lead variants (N =", nrow(v), ") ---\n"); print(s2k)
fwrite(s2k, file.path(OUT, "S2K_printed_vs_direct.tsv"), sep = "\t")
cat("\nS2K cross-tab producer x direct:\n"); print(v[, table(producer, fine_class)])
cat(sprintf("\nVariants with no locateVariants row: %d; of these inside a knownGene transcript: %d\n",
            v[any_locate_row == FALSE, .N], v[any_locate_row == FALSE & in_tx == TRUE, .N]))

# ================================================================= diagnostics
miss <- v[any_locate_row == FALSE]
fwrite(miss, file.path(OUT, "S2K_no_locate_row_variants.tsv"), sep = "\t")
gm <- mk_gr(miss$chr, miss$pos_hg19)
types <- list(coding = CodingVariants(), intron = IntronVariants(), fiveUTR = FiveUTRVariants(),
              threeUTR = ThreeUTRVariants(), spliceSite = SpliceSiteVariants(),
              promoter = PromoterVariants(upstream = 2000, downstream = 200),
              intergenic = IntergenicVariants())
cat("\n--- per-type locateVariants on the no-row variants ---\n")
for (nm in names(types)) {
  r <- tryCatch(locateVariants(gm, txdb, types[[nm]], cache = lv_cache), error = function(e) {cat(nm, "ERROR:", conditionMessage(e), "\n"); NULL})
  if (!is.null(r)) cat(sprintf("%-10s rows=%d variants=%d\n", nm, length(r), length(unique(r$QUERYID))))
}
# the same per-type calls on all S2K variants versus direct overlap of that feature
cat("\n--- per-type: locateVariants hits vs direct-overlap hits (all S2K variants) ---\n")
feat <- list(coding = cds_gr, intron = introns, fiveUTR = five_gr, threeUTR = three_gr,
             spliceSite = ss_gr, promoter = prom_gr)
for (nm in names(feat)) {
  r <- locateVariants(gr, txdb, types[[nm]], cache = lv_cache)
  lv_hit <- seq_along(gr) %in% r$QUERYID
  d_hit <- overlapsAny(gr, feat[[nm]], ignore.strand = TRUE)
  cat(sprintf("%-10s locateVariants=%d direct=%d both=%d direct_only=%d locate_only=%d\n",
              nm, sum(lv_hit), sum(d_hit), sum(lv_hit & d_hit), sum(d_hit & !lv_hit), sum(lv_hit & !d_hit)))
}
# IntergenicVariants on all variants: who gets an intergenic row
ig <- locateVariants(gr, txdb, IntergenicVariants(), cache = lv_cache)
cat(sprintf("\nIntergenicVariants on all S2K variants: %d variants get a row; of the no-row set: %d\n",
            length(unique(ig$QUERYID)), sum(which(!v$any_locate_row) %in% ig$QUERYID)))

# overlapping transcripts of the no-row variants: gene id present?
ov <- findOverlaps(gm, tx, ignore.strand = TRUE)
otx <- data.table(variant_key = miss$variant_key[queryHits(ov)],
                  TXNAME = tx$TXNAME[subjectHits(ov)], strand = as.character(strand(tx))[subjectHits(ov)],
                  GENEID = tx$GENEID[subjectHits(ov)])
cat("\nNo-row variants inside transcripts: overlapping transcripts with GENEID NA vs present:\n")
print(otx[, .(all_tx_geneid_na = all(is.na(GENEID))), by = variant_key][, .N, by = all_tx_geneid_na])
fwrite(otx, file.path(OUT, "S2K_no_row_overlapping_tx.tsv"), sep = "\t")
cat("\nGENEID NA share among all knownGene transcripts:", round(mean(is.na(tx$GENEID)), 3), "\n")
# located variants: share inside a gene-less transcript only
ov_all <- findOverlaps(gr, tx, ignore.strand = TRUE)
oa <- data.table(i = queryHits(ov_all), na = is.na(tx$GENEID[subjectHits(ov_all)]))[, .(all_na = all(na)), by = i]
v[, only_geneless_tx := FALSE]; v[oa$i[oa$all_na], only_geneless_tx := TRUE]
cat("\ncross-tab any_locate_row x (in_tx and only gene-less transcripts):\n")
print(v[in_tx == TRUE, table(any_locate_row, only_geneless_tx)])

# example 10:101953705
ex <- mk_gr(10, 101953705)
cat("\nExample 10:101953705 overlapping transcripts:\n")
print(as.data.table(tx[overlapsAny(tx, ex, ignore.strand = TRUE)])[, .(seqnames, start, end, strand, TXNAME, GENEID)])
for (nm in names(types)) {
  r <- locateVariants(ex, txdb, types[[nm]], cache = lv_cache)
  cat(sprintf("  %-10s rows=%d %s\n", nm, length(r), paste(unique(as.character(r$TXID)), collapse = ",")))
}
cat("\nVariantAnnotation internals used by IntronVariants:\n")
print(getMethod("locateVariants", c("GRanges", "TxDb", "IntronVariants")))
print(VariantAnnotation:::.makeResult)

}
if (MODE == "fig2c") {
# ================================================================= Fig 2C
f <- fread(COLOC, select = c("gwas_name", "gene", "PP.H4.susie", "PP.H4.abf", "top_snp"))
f <- f[gwas_name %in% main & !is.na(gene) & gene != ""]
f[, pp4 := suppressWarnings(pmax(as.numeric(PP.H4.susie), as.numeric(PP.H4.abf), na.rm = TRUE))]
f[!is.finite(pp4), pp4 := NA_real_]
setorder(f, gene, -pp4, gwas_name)
best <- f[, .SD[1], by = gene][is.finite(pp4) & pp4 > 0.5]
cat(sprintf("\nFig 2C: %d genes, %d distinct lead variants\n", nrow(best), uniqueN(best$top_snp)))
bv <- unique(best[, .(top_snp)])
bv[, c("chr", "pos") := tstrsplit(top_snp, ":", fixed = TRUE, type.convert = TRUE)]
g2 <- mk_gr(bv$chr, bv$pos)
bv[, producer := producer_class(g2, bv$top_snp)$producer]
bv <- cbind(bv, direct_class(g2)[, .(direct = fine_class)])
best <- merge(best, bv, by = "top_snp", sort = FALSE)
coarse <- function(x) fcase(x %in% c("coding", "spliceSite"), "coding",
                            x %in% c("fiveUTR", "threeUTR"), "UTR", default = x)
best[, `:=`(producer5 = coarse(producer), direct5 = coarse(direct))]
lev5 <- c("coding", "UTR", "promoter", "intron", "intergenic", "nc_exon_only")
f2c <- merge(tab(best$producer5, lev5)[, .(class, producer_rule_genes = n)],
             tab(best$direct5, lev5)[, .(class, direct_genes = n)], by = "class", sort = FALSE)
cat("--- Fig 2C genes per class ---\n"); print(f2c)
cat("Fig 2C cross-tab producer x direct (genes):\n"); print(best[, table(producer5, direct5)])
# coding owner: does the CDS hit belong to this gene?
cod <- best[direct5 == "coding"]
gc <- mk_gr(as.integer(sub(":.*", "", cod$top_snp)), as.integer(sub(".*:", "", cod$top_snp)))
h <- findOverlaps(gc, cds_tx, ignore.strand = TRUE)
egid <- tx$GENEID[match(as.integer(names(cds_tx))[subjectHits(h)], tx$TXID)]
sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, unique(na.omit(egid)), "SYMBOL", "ENTREZID"))
hd <- data.table(i = queryHits(h), sym = unname(sym[egid]))
own <- hd[, .(this_gene = any(sym == cod$gene[i[1]], na.rm = TRUE)), by = i]
cat(sprintf("coding, this gene = %d; coding, other gene = %d\n", sum(own$this_gene), sum(!own$this_gene)))
fwrite(best, file.path(OUT, "Fig2C_leads_classified.tsv"), sep = "\t")
fwrite(f2c, file.path(OUT, "Fig2C_producer_vs_direct.tsv"), sep = "\t")

}
writeLines(c(sprintf("%s  %s", tools::md5sum(c(ANN, COLOC, TIER)), c(ANN, COLOC, TIER))),
           file.path(OUT, "input_md5.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, paste0("sessionInfo_", MODE, ".txt")))
cat("\nDONE\n")
