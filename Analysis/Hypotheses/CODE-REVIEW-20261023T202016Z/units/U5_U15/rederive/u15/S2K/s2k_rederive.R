# Independent retabulation of Figure S2K from the adopted COLOC r2 master.
# Pairs: gene x stratum in the 35 primary strata with max(PP.H4 SuSiE, PP.H4 ABF) > 0.5.
# Variant: the master's top_snp for each pair; each unique variant counted once.
# Class: VariantAnnotation::locateVariants on UCSC hg19 knownGene, priority
# coding > spliceSite > fiveUTR > threeUTR > promoter (-2 kb/+200 bp) > intron > intergenic.
# No randomness.
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(VariantAnnotation)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
})
R <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT <- file.path(R, "Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U5_U15/rederive/u15/S2K")
MASTER <- file.path(R, "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv")
tier <- fread(file.path(R, "GWAS/finemapping/config/gwas_trait_tier.tsv"))
main <- tier[placement == "main", study_name]
stopifnot(length(main) == 35)

m <- fread(MASTER, select = c("gwas_name", "gene", "ensembl", "chr", "PP.H4.abf", "PP.H4.susie", "top_snp"))
cat("master rows", nrow(m), "\n")
m <- m[gwas_name %in% main]
cat("rows in 35 strata", nrow(m), "\n")
m[, best := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
p <- m[!is.na(best) & best > 0.5]
cat("rows with best PP.H4 > 0.5 (any gene label):", nrow(p), "\n")
cat("  of which empty/NA gene symbol:", p[is.na(gene) | gene == "", .N], "\n")
p <- p[!is.na(gene) & gene != ""]
cat("pairs (rows) with symbol:", nrow(p), "; unique gwas x gene:", uniqueN(p[, .(gwas_name, gene)]),
    "; genes:", uniqueN(p$gene), "; strata:", uniqueN(p$gwas_name), "\n")
cat("pairs with missing top_snp:", p[is.na(top_snp) | top_snp == "", .N], "\n")
cat("genes by method passing: SuSiE>0.5", uniqueN(p[PP.H4.susie > 0.5, gene]),
    "; ABF>0.5", uniqueN(p[PP.H4.abf > 0.5, gene]), "\n")

v <- unique(p[, .(top_snp)])
v[, c("vchr", "vpos") := tstrsplit(top_snp, ":", keep = 1:2)]
v[, vpos := as.integer(vpos)]
cat("unique top variants:", nrow(v), "; chromosomes:", paste(sort(unique(v$vchr)), collapse = ","), "\n")

gr <- GRanges(paste0("chr", v$vchr), IRanges(v$vpos, width = 1))
gr$QUERYID0 <- seq_along(gr)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
loc <- suppressWarnings(locateVariants(gr, txdb,
  AllVariants(promoter = PromoterVariants(upstream = 2000L, downstream = 200L))))
lev <- c("coding", "spliceSite", "fiveUTR", "threeUTR", "promoter", "intron", "intergenic")
ld <- data.table(q = loc$QUERYID, cls = as.character(loc$LOCATION))
ld <- ld[cls %in% lev]
ld[, rk := match(cls, lev)]
best <- ld[, .(cls = lev[min(rk)]), by = q]
v[, fine_class := NA_character_]
v[best$q, fine_class := best$cls]
cat("variants without any class:", v[is.na(fine_class), .N], "\n")
tab <- v[, .N, by = fine_class][, pct := round(100 * N / sum(N), 2)][order(-N)]
print(tab)

# per-class pair counts
p <- merge(p, v[, .(top_snp, fine_class)], by = "top_snp", all.x = TRUE)
print(p[, .(n_pairs = .N), by = fine_class][order(-n_pairs)])

# compare with the producer's per-variant annotation
ann <- fread(file.path(R, "figures/candidates/fig2-coloc-r2-20260924T053244Z/variant_class/coloc_variant_annotation.csv"))
ann <- unique(ann[study %in% main & !is.na(gene_symbol) & gene_symbol != "", .(variant_key, ann_class = fine_class)])
cmp <- merge(v[, .(variant_key = top_snp, fine_class)], ann, by = "variant_key", all = TRUE)
cat("variants: mine", nrow(v), "; producer annotation", nrow(ann), "; in both", cmp[!is.na(fine_class) & !is.na(ann_class), .N], "\n")
cat("class disagreements:", cmp[!is.na(fine_class) & !is.na(ann_class) & fine_class != ann_class, .N], "\n")
print(cmp[is.na(fine_class) | is.na(ann_class) | fine_class != ann_class][1:min(.N, 40)])
fwrite(tab, file.path(OUT, "S2K_rederived.csv"))
fwrite(cmp, file.path(OUT, "S2K_variant_class_compare.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
