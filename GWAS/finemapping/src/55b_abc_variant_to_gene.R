#!/usr/bin/env Rscript
# ============================================================================
# 55b_abc_variant_to_gene.R
# Tier-2 CROSS-VALIDATION layer: ABC (Activity-by-Contact) enhancer->gene
# variant-to-gene assignment from Nasser et al. 2021 (Nature).
#
# FRAMING: This is an EXTERNAL liver-regulatory cross-validation of the atlas's
# OWN MASLD-context scATAC (GSE244832; scripts 55-57), NOT primary evidence.
# HepG2 is a hepatocellular CARCINOMA line -- HepG2 ABC links are used only as
# an additional liver-lineage regulatory cell type for QC/overlap, NOT as MASLD
# evidence. The relevance gate is: a MASLD GWAS credible-set variant intersects
# an ABC enhancer that is predicted to regulate a target gene. ABC supplies a
# variant->gene assignment that is NOT nearest-gene (the headline value-add).
#
# Reuses the verified 55 / 58b frontend:
#   combined_finemapping.csv -> CS/PIP>0.1 filter -> dedup.
# ABC predictions (Nasser2021 AllPredictions.AvgHiC) are hg19 -- we therefore
# overlap the hg19 credible-set variants DIRECTLY against the hg19 ABC enhancers
# (matching ABC's genome build; no liftOver of the variant side needed). Target
# gene -> human_symbol is the ABC TargetGene column (already a gene symbol).
#
# Output: RNA-seq/results/multi_evidence/abc_atlas_columns.tsv (gene-level)
#         + GWAS/finemapping/results/gwas_atac/abc_variant_to_gene.csv (variant-level)
#         + an ABC enhancer subset BED (liver/HepG2) for reuse by 55c.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
  library(GenomicRanges)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR  <- file.path(BASE, "GWAS/finemapping")
OUT_DIR <- file.path(FM_DIR, "results/gwas_atac")
ABC_DIR <- file.path(BASE, "data/external/abc_liver")
ABC_RAW <- file.path(ABC_DIR, "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
ABC_URL <- "https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
ABC_LIVER_SUBSET <- file.path(ABC_DIR, "abc_liver_hepg2_enhancers.tsv.gz")  # reused by 55c
OUT_TSV <- file.path(BASE, "RNA-seq/results/multi_evidence/abc_atlas_columns.tsv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ABC_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("55b_abc_variant_to_gene.R  --  ABC enhancer->gene cross-validation\n")
cat("============================================================\n\n")

# ── 0. Download ABC predictions if absent ────────────────────────────────────
if (!file.exists(ABC_RAW)) {
  cat("Downloading Nasser2021 ABC predictions (~340 MB) ...\n")
  ok <- download.file(ABC_URL, ABC_RAW, method = "libcurl", quiet = FALSE, mode = "wb")
  if (ok != 0 || !file.exists(ABC_RAW)) stop("ABC download failed.")
}
cat("ABC predictions file: ", ABC_RAW, " (",
    round(file.info(ABC_RAW)$size / 1e6, 1), " MB)\n", sep = "")

# ── 1. 55 frontend: load -> filter CS/PIP>0.1 -> dedup (hg19) ────────────────
fm <- fread(file.path(FM_DIR, "results/combined_finemapping.csv"))
fm_cs <- fm %>% filter(
  (either_in_cs == TRUE) |
  (!is.na(recommended_pip) & recommended_pip > 0.1) |
  (is.na(recommended_pip) & !is.na(max_pip) & max_pip > 0.1))
variants <- as.data.table(fm_cs %>%
  group_by(chromosome, position, allele1, allele2) %>%
  summarise(max_pip = max(max_pip, na.rm = TRUE),
            max_rec_pip = if (any(!is.na(recommended_pip))) max(recommended_pip, na.rm = TRUE) else max(max_pip, na.rm = TRUE),
            n_studies = n_distinct(study),
            study_list = paste(sort(unique(study)), collapse = ";"),
            .groups = "drop"))
variants[, chr := paste0("chr", chromosome)]   # ABC uses chr-prefixed hg19
cat(sprintf("  credible-set variants (deduped, hg19): %d\n", nrow(variants)))

# ── 2. Stream ABC predictions; keep ONLY liver / HepG2 biosamples ────────────
# Nasser2021 CellType labels for liver lineage:
#   "liver-ENCODE" (primary liver tissue) and "HEPG2-Roadmap"/"HepG2" (carcinoma line).
# We discover the exact tokens from the file rather than hard-code (robust to
# label drift), matching the substrings "liver" or "hepg2" (case-insensitive).
cat("\n--- Streaming ABC, filtering to liver/HepG2 cell types ---\n")
# data.table::fread can read .gz directly; the file is ~340 MB compressed,
# ~2 GB uncompressed -> read full then filter (compute node has >=64 GB).
abc <- fread(ABC_RAW, showProgress = TRUE)
cat("  total ABC rows: ", nrow(abc), "\n", sep = "")
ct_all <- unique(abc$CellType)
liver_cts <- ct_all[grepl("liver|hepg2|hepatocyte", ct_all, ignore.case = TRUE)]
cat("  liver/HepG2 CellType tokens matched: ", paste(liver_cts, collapse = ", "), "\n", sep = "")
abc_liver <- abc[CellType %in% liver_cts]
rm(abc); gc(verbose = FALSE)
cat("  liver/HepG2 ABC enhancer-gene predictions: ", nrow(abc_liver),
    " (target genes: ", uniqueN(abc_liver$TargetGene), ")\n", sep = "")

# Save the liver/HepG2 enhancer subset for reuse by 55c (peak->gene via ABC).
fwrite(abc_liver[, .(chr, start, end, TargetGene, TargetGeneTSS, ABC.Score, CellType, isSelfPromoter)],
       ABC_LIVER_SUBSET, sep = "\t")
cat("  wrote liver ABC subset for 55c reuse: ", ABC_LIVER_SUBSET, "\n", sep = "")

# ── 3. Overlap credible-set variants (hg19) with ABC enhancers (hg19) ────────
# ABC enhancer coords are BED-style 0-based half-open -> 1-based GRanges.
cat("\n--- Overlapping credible-set variants with ABC enhancers ---\n")
gr_var <- GRanges(seqnames = variants$chr,
                  ranges = IRanges(start = variants$position, width = 1))
gr_enh <- GRanges(seqnames = abc_liver$chr,
                  ranges = IRanges(start = abc_liver$start + 1L, end = abc_liver$end))
hits <- findOverlaps(gr_var, gr_enh)
cat("  variant x ABC-enhancer overlaps: ", length(hits),
    "  (unique variants: ", uniqueN(queryHits(hits)), ")\n", sep = "")

if (length(hits) == 0) {
  cat("WARNING: no ABC overlaps. Writing empty output.\n")
  fwrite(data.table(human_symbol = character()), OUT_TSV, sep = "\t")
  quit(status = 0)
}

ovl <- data.table(
  variant_row = queryHits(hits),
  enh_row     = subjectHits(hits)
)
ovl[, `:=`(
  chromosome = variants$chromosome[variant_row],
  position   = variants$position[variant_row],
  allele1    = variants$allele1[variant_row],
  allele2    = variants$allele2[variant_row],
  max_rec_pip = variants$max_rec_pip[variant_row],
  study_list  = variants$study_list[variant_row],
  abc_target_gene = abc_liver$TargetGene[enh_row],
  abc_score       = abc_liver$ABC.Score[enh_row],
  abc_celltype    = abc_liver$CellType[enh_row],
  abc_is_self_promoter = abc_liver$isSelfPromoter[enh_row]
)]
ovl[, variant_id := paste(chromosome, position, allele1, allele2, sep = ":")]

# Variant-level export (one row per variant x ABC enhancer-gene link)
fwrite(ovl[, .(variant_id, chromosome, position, allele1, allele2,
               max_rec_pip, study_list,
               abc_target_gene, abc_score, abc_celltype, abc_is_self_promoter)],
       file.path(OUT_DIR, "abc_variant_to_gene.csv"))
cat("  wrote variant-level: ", file.path(OUT_DIR, "abc_variant_to_gene.csv"), "\n", sep = "")

# ── 4. Aggregate per target gene ─────────────────────────────────────────────
# Relevance gate: only count a gene if a MASLD credible-set variant lands in an
# ABC enhancer predicting that gene (which is exactly what `ovl` encodes).
agg <- ovl[, .(
  abc_max_score          = max(abc_score, na.rm = TRUE),
  abc_n_credset_variants = uniqueN(variant_id),
  abc_celltypes          = paste(sort(unique(abc_celltype)), collapse = ";"),
  # MASLD-disease relevance flag: any supporting variant from a disease endpoint
  abc_masld_gwas_driven  = any(grepl("NAFLD|NASH|cirrh|HCC|steato|hepatocell|ghodsian|nafl",
                                     study_list, ignore.case = TRUE))
), by = abc_target_gene]
setnames(agg, "abc_target_gene", "human_symbol")
agg <- agg[!is.na(human_symbol) & human_symbol != ""]
agg[, abc_v2g_hit := TRUE]
setcolorder(agg, c("human_symbol", "abc_v2g_hit", "abc_max_score",
                   "abc_n_credset_variants", "abc_celltypes", "abc_masld_gwas_driven"))
agg <- agg[order(-abc_max_score)]

fwrite(agg, OUT_TSV, sep = "\t")
cat(sprintf("\nWROTE %s : %d genes with an ABC credible-set variant->gene link\n",
            OUT_TSV, nrow(agg)))
cat(sprintf("  credible-set variants overlapping ABC liver enhancers: %d\n",
            uniqueN(ovl$variant_id)))
cat(sprintf("  ACTG1 present (ABC NAFLD positive control): %s\n",
            ifelse("ACTG1" %in% agg$human_symbol, "YES", "NO")))
if ("ACTG1" %in% agg$human_symbol) print(agg[human_symbol == "ACTG1"])
print(head(agg, 15))
cat("\n[", as.character(Sys.time()), "] 55b done.\n", sep = "")
