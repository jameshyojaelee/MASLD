#!/usr/bin/env Rscript
# ============================================================================
# 55d_ccre_overlap.R
# Tier-2 CROSS-VALIDATION layer: ENCODE SCREEN cCRE (candidate cis-Regulatory
# Element) registry overlap (hg38). Does a MASLD GWAS credible-set variant fall
# in a registered cCRE, and of what class (PLS / pELS / dELS / CTCF-bound)?
#
# FRAMING: external regulatory-element cross-validation of the atlas's OWN
# MASLD-context scATAC (GSE244832; scripts 55-57), NOT primary evidence. The
# SCREEN registry is cell-type-agnostic (built across the whole ENCODE corpus);
# we OPTIONALLY add a HepG2 H3K27ac activity readout from an on-disk ENCODE
# bigWig for QC only (HepG2 = hepatocellular CARCINOMA line; NOT MASLD evidence).
# Relevance gate = credible-set variant intersects a cCRE.
#
# Reuses the verified 55 / 58b frontend (CS/PIP>0.1 -> dedup -> liftOver hg19->hg38).
# cCRE->gene assignment cascade (NOT nearest-gene-only): ABC liver enhancer
# (from 55b subset, lifted hg38) > nearest TSS (GENCODE v49).
#
# Output: RNA-seq/results/multi_evidence/ccre_atlas_columns.tsv (gene-level)
#         + GWAS/finemapping/results/gwas_atac/ccre_variant_overlap.csv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
  library(GenomicRanges); library(rtracklayer)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR  <- file.path(BASE, "GWAS/finemapping")
OUT_DIR <- file.path(FM_DIR, "results/gwas_atac")
CCRE    <- file.path(BASE, "data/external/screen_ccre_hepg2/GRCh38-cCREs.bed")
ABC_LIVER_SUBSET <- file.path(BASE, "data/external/abc_liver/abc_liver_hepg2_enhancers.tsv.gz")
GTF     <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
H3K27AC <- list.files(file.path(BASE, "data/external/encode_liver"),
                      pattern = "H3K27ac.*\\.bigWig$", full.names = TRUE)
OUT_TSV <- file.path(BASE, "RNA-seq/results/multi_evidence/ccre_atlas_columns.tsv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("55d_ccre_overlap.R  --  SCREEN cCRE credible-set cross-validation\n")
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

# ── 2. Load SCREEN cCRE registry (hg38 BED6) ─────────────────────────────────
ccre <- fread(CCRE, header = FALSE,
              col.names = c("chr", "start", "end", "acc_dnase", "acc_ccre", "ccre_class"))
cat("  SCREEN cCREs loaded: ", nrow(ccre), " (classes: ",
    paste(sort(unique(ccre$ccre_class)), collapse = ", "), ")\n", sep = "")

# ── 3. Overlap credible-set variants (hg38) with cCREs ───────────────────────
gr_var  <- GRanges(vh$chr_hg38, IRanges(vh$pos_hg38, width = 1))
gr_ccre <- GRanges(ccre$chr, IRanges(ccre$start + 1L, ccre$end))  # BED 0-based -> 1-based
ov <- findOverlaps(gr_var, gr_ccre)
cat("\n--- Variant x cCRE overlaps: ", length(ov),
    " (unique variants: ", uniqueN(queryHits(ov)), ") ---\n", sep = "")

if (length(ov) == 0) {
  cat("WARNING: no cCRE overlaps. Writing empty output.\n")
  fwrite(data.table(human_symbol = character()), OUT_TSV, sep = "\t")
  quit(status = 0)
}

hit <- data.table(
  variant_row = queryHits(ov),
  ccre_row    = subjectHits(ov))
hit[, `:=`(
  chromosome = vh$chromosome[variant_row],
  position   = vh$position[variant_row],
  allele1    = vh$allele1[variant_row],
  allele2    = vh$allele2[variant_row],
  chr_hg38   = vh$chr_hg38[variant_row],
  pos_hg38   = vh$pos_hg38[variant_row],
  max_rec_pip = vh$max_rec_pip[variant_row],
  study_list  = vh$study_list[variant_row],
  ccre_class  = ccre$ccre_class[ccre_row],
  ccre_acc    = ccre$acc_ccre[ccre_row],
  ccre_start  = ccre$start[ccre_row],
  ccre_end    = ccre$end[ccre_row])]
hit[, variant_id := paste(chromosome, position, allele1, allele2, sep = ":")]

# Collapse the SCREEN sub-classes into the canonical 4 + CA families.
# PLS=promoter-like, pELS=proximal-enhancer-like, dELS=distal-enhancer-like,
# CTCF=CTCF-bound (CA-CTCF), other CA-* / TF retained as "CA_other".
hit[, ccre_class_simple := fcase(
  grepl("^PLS", ccre_class), "PLS",
  grepl("^pELS", ccre_class), "pELS",
  grepl("^dELS", ccre_class), "dELS",
  grepl("CTCF", ccre_class), "CTCF",
  default = "CA_other")]

# ── 4. Optional HepG2 H3K27ac activity at the cCRE (QC only) ──────────────────
if (length(H3K27AC) >= 1) {
  cat("\n--- Scoring HepG2 H3K27ac activity at overlapped cCREs (QC only) ---\n")
  bw <- BigWigFile(H3K27AC[1])
  # Summarize mean H3K27ac signal over each overlapped cCRE element.
  uniq_ccre <- unique(hit[, .(chr_hg38, ccre_start, ccre_end, ccre_acc)])
  gr_q <- GRanges(uniq_ccre$chr_hg38, IRanges(uniq_ccre$ccre_start + 1L, uniq_ccre$ccre_end))
  sm <- tryCatch(
    rtracklayer::summary(bw, gr_q, type = "mean", size = 1L),
    error = function(e) { cat("  bigWig summary failed: ", conditionMessage(e), "\n", sep = ""); NULL })
  if (!is.null(sm)) {
    uniq_ccre[, h3k27ac_mean := vapply(sm, function(x) {
      v <- score(x); if (length(v) == 0 || all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE) }, numeric(1))]
    hit <- merge(hit, uniq_ccre[, .(chr_hg38, ccre_start, ccre_end, ccre_acc, h3k27ac_mean)],
                 by = c("chr_hg38", "ccre_start", "ccre_end", "ccre_acc"), all.x = TRUE)
    cat("  cCREs with HepG2 H3K27ac signal > 0: ",
        sum(hit$h3k27ac_mean > 0, na.rm = TRUE), " / ", nrow(hit), "\n", sep = "")
  } else {
    hit[, h3k27ac_mean := NA_real_]
  }
} else {
  cat("\n  (no HepG2 H3K27ac bigWig found -- skipping activity QC)\n")
  hit[, h3k27ac_mean := NA_real_]
}

# ── 5. cCRE->gene: ABC liver enhancer (hg38) > nearest TSS ───────────────────
cat("\n--- Assigning cCREs to genes (ABC > nearest TSS) ---\n")
gr_pos <- GRanges(hit$chr_hg38, IRanges(hit$pos_hg38, width = 1))
hit[, gene_assigned := NA_character_]
hit[, gene_source := NA_character_]

if (file.exists(ABC_LIVER_SUBSET)) {
  abc <- fread(ABC_LIVER_SUBSET)
  abc_gr19 <- GRanges(abc$chr, IRanges(abc$start + 1L, abc$end))
  mcols(abc_gr19)$TargetGene <- abc$TargetGene
  mcols(abc_gr19)$ABC.Score  <- abc$ABC.Score
  abc_lift <- liftOver(abc_gr19, chain)
  keep <- lengths(abc_lift) == 1
  abc_gr38 <- unlist(abc_lift[keep])
  abc_meta <- as.data.table(mcols(abc_gr19)[keep, , drop = FALSE])
  ab <- findOverlaps(gr_pos, abc_gr38)
  if (length(ab) > 0) {
    od <- data.table(q = queryHits(ab),
                     gene = abc_meta$TargetGene[subjectHits(ab)],
                     sc = abc_meta$ABC.Score[subjectHits(ab)])
    best <- od[order(q, -sc)][, .SD[1], by = q]
    hit$gene_assigned[best$q] <- best$gene
    hit$gene_source[best$q]   <- "ABC"
  }
  cat("  variants assigned via ABC liver enhancer: ", sum(hit$gene_source == "ABC", na.rm = TRUE), "\n", sep = "")
} else {
  cat("  (ABC liver subset not found -- run 55b first; nearest TSS only)\n")
}

genes <- rtracklayer::import(GTF, feature.type = "gene")
genes <- genes[seqnames(genes) %in% paste0("chr", c(1:22, "X", "Y"))]
tss <- resize(genes, width = 1L, fix = "start")
need <- which(is.na(hit$gene_assigned))
if (length(need) > 0) {
  ni <- nearest(gr_pos[need], tss, ignore.strand = TRUE)
  hit$gene_assigned[need] <- ifelse(!is.na(ni), genes$gene_name[ni], NA_character_)
  hit$gene_source[need]   <- "nearest_TSS"
}
cat("  variants assigned via nearest TSS: ", sum(hit$gene_source == "nearest_TSS", na.rm = TRUE), "\n", sep = "")

# ── 6. Variant-level export ──────────────────────────────────────────────────
fwrite(hit[, .(variant_id, chromosome, position, allele1, allele2, max_rec_pip, study_list,
               ccre_acc, ccre_class, ccre_class_simple, h3k27ac_mean,
               gene_assigned, gene_source)],
       file.path(OUT_DIR, "ccre_variant_overlap.csv"))
cat("\n  wrote variant-level: ", file.path(OUT_DIR, "ccre_variant_overlap.csv"), "\n", sep = "")

# ── 7. Aggregate per gene ────────────────────────────────────────────────────
# Class priority for the per-gene representative class: PLS > pELS > dELS > CTCF > CA_other
class_rank <- c(PLS = 1L, pELS = 2L, dELS = 3L, CTCF = 4L, CA_other = 5L)
hit[, class_pri := class_rank[ccre_class_simple]]
agg <- hit[!is.na(gene_assigned) & gene_assigned != ""][order(gene_assigned, class_pri), .(
  ccre_credset_hit       = TRUE,
  ccre_class             = ccre_class_simple[1],          # highest-priority class
  ccre_classes_all       = paste(sort(unique(ccre_class_simple)), collapse = ";"),
  ccre_n_credset_variants = uniqueN(variant_id),
  ccre_hepg2_active      = any(h3k27ac_mean > 0, na.rm = TRUE),   # QC flag only
  ccre_hepg2_h3k27ac_max = if (all(is.na(h3k27ac_mean))) NA_real_ else max(h3k27ac_mean, na.rm = TRUE),
  ccre_gene_source       = paste(sort(unique(gene_source)), collapse = ";"),
  ccre_masld_gwas_driven = any(grepl("NAFLD|NASH|cirrh|HCC|steato|hepatocell|ghodsian|nafl",
                                     study_list, ignore.case = TRUE))
), by = .(human_symbol = gene_assigned)]

fwrite(agg, OUT_TSV, sep = "\t")
cat(sprintf("\nWROTE %s : %d genes with a credible-set cCRE overlap\n", OUT_TSV, nrow(agg)))
cat(sprintf("  credible-set variants overlapping a cCRE: %d\n", uniqueN(hit$variant_id)))
cat("  cCRE class breakdown (per-gene representative):\n")
print(agg[, .N, by = ccre_class][order(-N)])
cat(sprintf("  ACTG1 present: %s\n", ifelse("ACTG1" %in% agg$human_symbol, "YES", "NO")))
print(head(agg, 15))
cat("\n[", as.character(Sys.time()), "] 55d done.\n", sep = "")
