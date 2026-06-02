#!/usr/bin/env Rscript
# ============================================================================
# 55c_caqtl_credset.R
# Tier-2 CROSS-VALIDATION layer: bulk-liver caQTL (Currin & Mohlke 2025,
# Genome Research; 138 donors; doi:10.1101/gr.279741.124) credible-set overlap.
#
# FRAMING: external liver-regulatory cross-validation of the atlas's OWN
# MASLD-context scATAC (GSE244832; scripts 55-57), NOT primary evidence.
# Does a MASLD GWAS credible-set variant ALSO act as a SIGNIFICANT liver
# caQTL (i.e. alter chromatin accessibility of a liver peak)? Relevance gate =
# the variant is a credible-set variant AND a significant caQTL.
#
# Reuses the verified 55 / 58b frontend (CS/PIP>0.1 -> dedup -> liftOver hg19->hg38)
# and the Currin liftOver/orientation-match template from 56i.
# Currin `variant` is hg38 "chr:pos:REF:ALT"; we liftOver our hg19 credible-set
# variants to hg38 and match on (chr,pos38) with allele-set concordance.
#
# Peak->gene assignment cascade (NOT nearest-gene-only):
#   1. ABC liver/HepG2 enhancer predicting a gene (from 55b subset) at the
#      reconstructed caQTL peak center  ->  abc_target_gene
#   2. else nearest TSS (GENCODE v49)
# Currin files carry only a peak ID + distance_from_peakCenter, so the peak
# CENTER (hg38) is reconstructed as variant_pos - distance_from_peakCenter.
#
# Output: RNA-seq/results/multi_evidence/caqtl_atlas_columns.tsv (gene-level)
#         + GWAS/finemapping/results/gwas_atac/caqtl_variant_credset.csv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
  library(GenomicRanges); library(rtracklayer)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR  <- file.path(BASE, "GWAS/finemapping")
OUT_DIR <- file.path(FM_DIR, "results/gwas_atac")
CURRIN_DIR <- file.path(BASE, "data/external/currin_2025_caqtl")
ABC_LIVER_SUBSET <- file.path(BASE, "data/external/abc_liver/abc_liver_hepg2_enhancers.tsv.gz")
GTF <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
CAQTL_P <- 1e-5   # significant caQTL threshold (matches 56i "suggestive"); the
                  # nominal files carry no beta-permutation FDR, so 1e-5 is the
                  # honest significance call (same convention as 56i).
OUT_TSV <- file.path(BASE, "RNA-seq/results/multi_evidence/caqtl_atlas_columns.tsv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("55c_caqtl_credset.R  --  liver caQTL credible-set cross-validation\n")
cat("============================================================\n\n")

# ── 1. 55 frontend: load -> filter CS/PIP>0.1 -> dedup -> liftOver hg19->hg38 ─
fm <- fread(file.path(FM_DIR, "results/combined_finemapping.csv"))
fm_cs <- fm %>% filter(
  (either_in_cs == TRUE) |
  (!is.na(recommended_pip) & recommended_pip > 0.1) |
  (is.na(recommended_pip) & !is.na(max_pip) & max_pip > 0.1))
variants <- as.data.table(fm_cs %>%
  group_by(chromosome, position, allele1, allele2) %>%
  summarise(max_rec_pip = if (any(!is.na(recommended_pip))) max(recommended_pip, na.rm = TRUE) else max(max_pip, na.rm = TRUE),
            study_list = paste(sort(unique(study)), collapse = ";"),
            .groups = "drop"))
cat(sprintf("  credible-set variants (deduped, hg19): %d\n", nrow(variants)))

chain <- import.chain(file.path(BASE, "data/broadaway_eqtl/hg19ToHg38.over.chain"))
gr19 <- GRanges(paste0("chr", variants$chromosome),
                IRanges(variants$position, width = 1), strand = "*")
mcols(gr19) <- variants
lifted <- liftOver(gr19, chain)
gr38 <- unlist(lifted[lengths(lifted) == 1])
vh <- as.data.table(mcols(gr38))
vh[, `:=`(chr_hg38 = as.character(seqnames(gr38)), pos_hg38 = start(gr38))]
cat(sprintf("  lifted to hg38: %d\n", nrow(vh)))

# Build both-orientation hg38 variant keys to match Currin "chr:pos:REF:ALT".
vh[, var_id_AB := paste0(chr_hg38, ":", pos_hg38, ":", allele1, ":", allele2)]
vh[, var_id_BA := paste0(chr_hg38, ":", pos_hg38, ":", allele2, ":", allele1)]

# ── 2. Stream Currin per-chr, keep rows matching our credible-set variants ───
cat("\n--- Scanning Currin caQTL per chromosome ---\n")
currin_files <- list.files(CURRIN_DIR, pattern = "^caQTL_hg38_chr.*\\.txt\\.gz$",
                           full.names = TRUE)
stopifnot(length(currin_files) >= 22L)

chrs <- sort(unique(vh$chr_hg38)); chrs <- chrs[!is.na(chrs)]
hits_list <- list()
for (ch in chrs) {
  cf <- currin_files[grepl(paste0("caQTL_hg38_", ch, "_"), basename(currin_files), fixed = TRUE)]
  if (length(cf) != 1L) { cat("  [", ch, "] no Currin file; skip\n", sep = ""); next }
  vchr <- vh[chr_hg38 == ch]
  if (nrow(vchr) == 0L) next
  lookup <- unique(c(vchr$var_id_AB, vchr$var_id_BA)); lookup <- lookup[!is.na(lookup)]
  cdt <- fread(cf, sep = "\t", header = TRUE, showProgress = FALSE)
  cdt_hit <- cdt[variant %in% lookup]
  cat("  [", ch, "] Currin rows ", nrow(cdt), " | matched ", nrow(cdt_hit),
      " (variants ", uniqueN(cdt_hit$variant), "/", nrow(vchr), ")\n", sep = "")
  if (nrow(cdt_hit) > 0) hits_list[[ch]] <- cdt_hit
  rm(cdt, cdt_hit); gc(verbose = FALSE)
}
hits <- rbindlist(hits_list, fill = TRUE)
cat("\n  total Currin caQTL test rows for our variants: ", nrow(hits), "\n", sep = "")

if (nrow(hits) == 0) {
  cat("WARNING: no Currin matches. Writing empty output.\n")
  fwrite(data.table(human_symbol = character()), OUT_TSV, sep = "\t")
  quit(status = 0)
}

# Keep the most-significant caQTL test per (variant, peak): a variant tests
# against many peaks within 1 Mb; the peak it most strongly perturbs is the call.
hits_top <- hits[order(variant, pvalue), .SD[1], by = .(variant, peak)]
# Reconstruct hg38 peak CENTER from variant pos and distance_from_peakCenter.
hits_top[, c("v_chr", "v_pos", "v_ref", "v_alt") := tstrsplit(variant, ":", fixed = TRUE)]
hits_top[, v_pos := as.integer(v_pos)]
hits_top[, peak_center := v_pos - distance_from_peakCenter]

# Map each Currin hit back to our credible-set variant (track orientation).
ab <- vh[, .(var_id = var_id_AB, chromosome, position, allele1, allele2, max_rec_pip, study_list)]
ba <- vh[, .(var_id = var_id_BA, chromosome, position, allele1, allele2, max_rec_pip, study_list)]
v_key <- unique(rbind(ab, ba))
hits_top <- merge(hits_top, v_key, by.x = "variant", by.y = "var_id", all.x = TRUE)
cat("  unique credible-set variants with a Currin caQTL test: ",
    uniqueN(hits_top[, .(chromosome, position)]), "\n", sep = "")
cat("  ...with a SIGNIFICANT caQTL (p<", CAQTL_P, "): ",
    uniqueN(hits_top[pvalue < CAQTL_P, .(chromosome, position)]), "\n", sep = "")

# Keep only SIGNIFICANT caQTL peaks for gene assignment.
sig <- hits_top[pvalue < CAQTL_P]
if (nrow(sig) == 0) {
  cat("WARNING: no significant (p<", CAQTL_P, ") caQTLs. Writing empty output.\n", sep = "")
  fwrite(data.table(human_symbol = character()), OUT_TSV, sep = "\t")
  quit(status = 0)
}

# ── 3. Peak->gene: ABC liver enhancer at peak center, else nearest TSS ───────
cat("\n--- Assigning caQTL peaks to genes (ABC > nearest TSS) ---\n")
gr_peak <- GRanges(sig$v_chr, IRanges(sig$peak_center, width = 1))  # peak center, hg38

sig[, gene_assigned := NA_character_]
sig[, gene_source := NA_character_]

# 3a. ABC liver enhancer (hg19) -> need hg38; the 55b ABC subset is hg19, so we
# instead assign via ABC only when an ABC enhancer (lifted to hg38) covers the
# peak center. LiftOver the ABC enhancer subset hg19->hg38 once.
if (file.exists(ABC_LIVER_SUBSET)) {
  abc <- fread(ABC_LIVER_SUBSET)
  abc_gr19 <- GRanges(abc$chr, IRanges(abc$start + 1L, abc$end))
  mcols(abc_gr19)$TargetGene <- abc$TargetGene
  mcols(abc_gr19)$ABC.Score  <- abc$ABC.Score
  abc_lift <- liftOver(abc_gr19, chain)  # hg19 -> hg38
  keep <- lengths(abc_lift) == 1
  abc_gr38 <- unlist(abc_lift[keep])
  abc_meta <- as.data.table(mcols(abc_gr19)[keep, , drop = FALSE])
  ov <- findOverlaps(gr_peak, abc_gr38)
  if (length(ov) > 0) {
    od <- data.table(q = queryHits(ov),
                     gene = abc_meta$TargetGene[subjectHits(ov)],
                     sc = abc_meta$ABC.Score[subjectHits(ov)])
    best <- od[order(q, -sc)][, .SD[1], by = q]   # strongest ABC link per peak
    sig$gene_assigned[best$q] <- best$gene
    sig$gene_source[best$q]   <- "ABC"
  }
  cat("  peaks assigned via ABC liver enhancer: ", sum(sig$gene_source == "ABC", na.rm = TRUE), "\n", sep = "")
} else {
  cat("  (ABC liver subset not found -- run 55b first; falling back to nearest TSS only)\n")
}

# 3b. nearest TSS for the remainder (GENCODE v49, strand-aware 5' end).
genes <- rtracklayer::import(GTF, feature.type = "gene")
genes <- genes[seqnames(genes) %in% paste0("chr", c(1:22, "X", "Y"))]
tss <- resize(genes, width = 1L, fix = "start")
need <- which(is.na(sig$gene_assigned))
if (length(need) > 0) {
  ni <- nearest(gr_peak[need], tss, ignore.strand = TRUE)
  sig$gene_assigned[need] <- ifelse(!is.na(ni), genes$gene_name[ni], NA_character_)
  sig$gene_source[need]   <- "nearest_TSS"
}
cat("  peaks assigned via nearest TSS: ", sum(sig$gene_source == "nearest_TSS", na.rm = TRUE), "\n", sep = "")

# ── 4. Variant-level export ──────────────────────────────────────────────────
sig[, variant_id := paste(chromosome, position, allele1, allele2, sep = ":")]
fwrite(sig[, .(variant_id, chromosome, position, allele1, allele2, max_rec_pip, study_list,
               caqtl_peak = peak, caqtl_peak_center = peak_center,
               caqtl_pval = pvalue, caqtl_beta = beta, caqtl_maf = MAF,
               gene_assigned, gene_source)],
       file.path(OUT_DIR, "caqtl_variant_credset.csv"))
cat("\n  wrote variant-level: ", file.path(OUT_DIR, "caqtl_variant_credset.csv"), "\n", sep = "")

# ── 5. Aggregate per gene ────────────────────────────────────────────────────
agg <- sig[!is.na(gene_assigned) & gene_assigned != "", .(
  caqtl_credset_hit       = TRUE,
  caqtl_peak              = paste(sort(unique(peak)), collapse = ";"),
  caqtl_min_pval          = min(pvalue, na.rm = TRUE),
  caqtl_n_credset_variants = uniqueN(variant_id),
  caqtl_gene_source       = paste(sort(unique(gene_source)), collapse = ";"),
  caqtl_masld_gwas_driven = any(grepl("NAFLD|NASH|cirrh|HCC|steato|hepatocell|ghodsian|nafl",
                                      study_list, ignore.case = TRUE))
), by = .(human_symbol = gene_assigned)]
agg <- agg[order(caqtl_min_pval)]

fwrite(agg, OUT_TSV, sep = "\t")
cat(sprintf("\nWROTE %s : %d genes with a significant caQTL credible-set hit\n",
            OUT_TSV, nrow(agg)))
cat(sprintf("  credible-set variants that are significant caQTLs: %d\n",
            uniqueN(sig$variant_id)))
cat(sprintf("  ACTG1 present: %s\n", ifelse("ACTG1" %in% agg$human_symbol, "YES", "NO")))
print(head(agg, 15))
cat("\n[", as.character(Sys.time()), "] 55c done.\n", sep = "")
