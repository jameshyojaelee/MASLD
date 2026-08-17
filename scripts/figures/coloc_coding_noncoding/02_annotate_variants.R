#!/usr/bin/env Rscript
# 02_annotate_variants.R
# ----------------------------------------------------------------------------
# Classify each unique (chr, pos_hg19) into a fine class via
# VariantAnnotation::locateVariants() against TxDb.Hsapiens.UCSC.hg19.knownGene.
# Coarse rollup: coding (CDS or splice) vs non-coding.
#
# Reads:
#   RNA-seq/results/coloc_variant_classes/variants_long.csv
# Writes:
#   RNA-seq/results/coloc_variant_classes/variant_classification.csv
suppressPackageStartupMessages({
  library(data.table)
  library(VariantAnnotation)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(org.Hs.eg.db)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT  <- Sys.getenv(
  "FIG2_VARIANT_CLASS_DIR",
  file.path(BASE, "RNA-seq/results/coloc_variant_classes")
)

vlong <- fread(file.path(OUT, "variants_long.csv"))
cat(sprintf("Loaded variants_long.csv: %d rows | %d unique variant_keys\n",
            nrow(vlong), uniqueN(vlong$variant_key)))

# Build de-duplicated GRanges over (chr, pos_hg19)
uniq <- unique(vlong[, .(chr, pos_hg19, variant_key)])
uniq <- uniq[!is.na(chr) & !is.na(pos_hg19)]
input_variant_keys <- uniq$variant_key
# Keep the historical intergenic sentinel as a classifier fixture even when it
# is absent from the promoted release. It is removed before the output table is
# written, so the release universe remains exactly input-derived.
uniq <- unique(rbind(
  uniq,
  data.table(chr = 10L, pos_hg19 = 101193937L,
             variant_key = "10:101193937")
))
uniq[, seqnames := paste0("chr", chr)]
gr <- GRanges(seqnames = uniq$seqnames,
              ranges   = IRanges(start = uniq$pos_hg19, end = uniq$pos_hg19),
              variant_key = uniq$variant_key)
seqlevelsStyle(gr) <- "UCSC"
cat(sprintf("Built GRanges over %d unique variants\n", length(gr)))

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# locateVariants with all variant types incl. promoter (±2kb / 200 bp)
loc_all <- VariantAnnotation::locateVariants(
  gr, txdb,
  AllVariants(promoter = PromoterVariants(upstream = 2000, downstream = 200)))
cat(sprintf("locateVariants: %d annotation rows\n", length(loc_all)))

# Map LOCATION → fine class with precedence
# coding > spliceSite > fiveUTR > threeUTR > promoter > intron > intergenic
location_to_class <- function(loc) {
  s <- as.character(loc)
  fcase(
    s == "coding",     "coding",
    s == "spliceSite", "spliceSite",
    s == "fiveUTR",    "fiveUTR",
    s == "threeUTR",   "threeUTR",
    s == "promoter",   "promoter",
    s == "intron",     "intron",
    s == "intergenic", "intergenic",
    default            = "intergenic")
}
class_rank <- c("coding" = 1, "spliceSite" = 2, "fiveUTR" = 3,
                "threeUTR" = 4, "promoter" = 5, "intron" = 6,
                "intergenic" = 7)

ann <- data.table(
  variant_key = mcols(gr)$variant_key[loc_all$QUERYID],
  fine_class  = location_to_class(loc_all$LOCATION))
ann[, rank := class_rank[fine_class]]
setorder(ann, variant_key, rank)
class_one <- ann[!duplicated(variant_key), .(variant_key, fine_class)]
n_overlap <- ann[, .(n_overlap_classes = uniqueN(fine_class)), by = variant_key]
class_dt  <- merge(class_one, n_overlap, by = "variant_key", all.x = TRUE)

# Variants with no overlap at all → intergenic by default (most are >2kb from any feature)
unannot <- setdiff(uniq$variant_key, class_dt$variant_key)
if (length(unannot)) {
  class_dt <- rbind(class_dt, data.table(variant_key = unannot,
                                          fine_class = "intergenic",
                                          n_overlap_classes = 0L))
}

class_dt[, coarse_class := fifelse(fine_class %in% c("coding","spliceSite"),
                                    "coding", "non-coding")]

# Distance to nearest TSS
tss <- promoters(transcripts(txdb), upstream = 0, downstream = 1)
nearest <- distanceToNearest(gr, tss, ignore.strand = TRUE)
dist_dt <- data.table(
  variant_key = mcols(gr)$variant_key[queryHits(nearest)],
  distance_to_tss = mcols(nearest)$distance,
  tx_id = subjectHits(nearest))
# Get gene-symbol of nearest transcript
tx_to_gene <- AnnotationDbi::select(txdb, keys = as.character(seq_along(tss)),
                                    keytype = "TXID",
                                    columns = c("TXID","GENEID"))
setDT(tx_to_gene)
tx_to_gene[, TXID := as.integer(TXID)]
dist_dt <- merge(dist_dt, tx_to_gene, by.x = "tx_id", by.y = "TXID", all.x = TRUE)
sym_map <- tryCatch(
  setDT(AnnotationDbi::select(org.Hs.eg.db,
                              keys = unique(na.omit(dist_dt$GENEID)),
                              keytype = "ENTREZID",
                              columns = c("ENTREZID","SYMBOL"))),
  error = function(e) data.table(ENTREZID = character(), SYMBOL = character()))
dist_dt <- merge(dist_dt, sym_map, by.x = "GENEID", by.y = "ENTREZID", all.x = TRUE)
setnames(dist_dt, "SYMBOL", "nearest_gene_symbol")
dist_dt <- dist_dt[, .(variant_key, distance_to_tss, nearest_gene_symbol)]
dist_dt <- dist_dt[!duplicated(variant_key)]

class_dt <- merge(class_dt, dist_dt, by = "variant_key", all.x = TRUE)
class_dt[, distance_to_tss_kb := round(distance_to_tss / 1e3, 2)]

# predictCoding for variants with alleles (Definition B has them)
defB <- vlong[definition == "B_finemapped"]
if (nrow(defB)) {
  cs <- fread(file.path(BASE, "GWAS/finemapping/results/credible_sets.csv"),
              select = c("chromosome","position","allele1","allele2"))
  cs[, variant_key := paste0(chromosome, ":", position)]
  cs <- unique(cs[, .(variant_key, allele1, allele2)])
  defB_alleles <- merge(defB[, .(variant_key)], cs, by = "variant_key", all.x = TRUE)
  defB_alleles <- defB_alleles[!is.na(allele1) & !is.na(allele2)]
  defB_alleles <- unique(defB_alleles)
  if (nrow(defB_alleles)) {
    parts <- tstrsplit(defB_alleles$variant_key, ":", fixed = TRUE)
    gr_b <- GRanges(seqnames = paste0("chr", parts[[1]]),
                    ranges   = IRanges(start = as.integer(parts[[2]]),
                                       end   = as.integer(parts[[2]])),
                    variant_key = defB_alleles$variant_key)
    seqlevelsStyle(gr_b) <- "UCSC"
    # Upstream convention (00b_reformat_additional_gwas.R:48,67,86,105):
    # allele1 = effect_allele (= ALT), allele2 = other_allele (= REF).
    # predictCoding needs REF on the genome and varAllele = the substituted base.
    gr_b$REF <- DNAStringSet(defB_alleles$allele2)
    gr_b$varAllele <- DNAStringSet(defB_alleles$allele1)
    pc <- tryCatch(
      VariantAnnotation::predictCoding(gr_b, txdb,
                                        seqSource = NULL,
                                        varAllele = gr_b$varAllele),
      error = function(e) { cat("predictCoding skipped:", e$message, "\n"); NULL })
    if (!is.null(pc) && length(pc)) {
      pc_dt <- data.table(
        variant_key = mcols(gr_b)$variant_key[pc$QUERYID],
        coding_consequence = as.character(pc$CONSEQUENCE))
      pc_dt[, rank := fcase(
        coding_consequence == "nonsense",     1L,
        coding_consequence == "frameshift",   2L,
        coding_consequence == "nonsynonymous",3L,
        coding_consequence == "synonymous",   4L,
        default = 5L)]
      setorder(pc_dt, variant_key, rank)
      pc_dt <- pc_dt[!duplicated(variant_key), .(variant_key, coding_consequence)]
      class_dt <- merge(class_dt, pc_dt, by = "variant_key", all.x = TRUE)
    }
  }
}
if (!"coding_consequence" %in% names(class_dt)) class_dt[, coding_consequence := NA_character_]

# Verification — all spot-check coords MUST be present in the input;
# missing coords are a script-level failure (previous behaviour silently
# warned, which masked the fact that 3 of 5 historical coords were wrong).
verify <- function(key, expected_fine) {
  hit <- class_dt[variant_key == key]
  if (!nrow(hit)) {
    stop(sprintf("Verification target %s not in input variants — fix the spot-check coord", key))
  }
  ok <- hit$fine_class %in% expected_fine
  msg <- sprintf("  %s -> fine=%s coarse=%s (expected one of %s) %s\n",
                 key, hit$fine_class, hit$coarse_class,
                 paste(expected_fine, collapse="|"),
                 ifelse(ok, "OK", "MISMATCH"))
  cat(msg)
  if (!ok) stop("Verification mismatch for ", key)
}
cat("\n--- Verification spot-checks ---\n")
verify("22:44324727", c("coding"))      # PNPLA3 I148M (rs738409, hg19)
verify("19:19379549", c("coding"))      # TM6SF2 E167K (rs58542926, hg19)
verify("15:60883281", c("intron"))      # RORA region (hg19)
verify("10:102757934", c("promoter"))   # LZTS2 promoter sentinel
verify("10:101193937", c("intergenic")) # intergenic sentinel (3.4 kb from GOT1 TSS)

class_dt <- class_dt[variant_key %in% input_variant_keys]
fwrite(class_dt, file.path(OUT, "variant_classification.csv"))
cat("\n--- Class distribution ---\n")
print(class_dt[, .(n = .N), by = .(coarse_class, fine_class)][order(coarse_class, -n)])
cat(sprintf("\n[02_annotate_variants] DONE — wrote %d rows\n", nrow(class_dt)))
