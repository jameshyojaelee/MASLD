#!/usr/bin/env Rscript
# 55_gwas_atac_variant_overlap.R
# Intersects fine-mapped GWAS credible set variants with scATAC cell-type peaks
# Requires: 04_aggregate_results.R output (combined_finemapping.csv)
# Usage: Rscript 55_gwas_atac_variant_overlap.R

library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)
library(readr)

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
ATAC_DIR <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
OUT_DIR  <- file.path(FM_DIR, "results/gwas_atac")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("55_gwas_atac_variant_overlap.R\n")
cat("Intersecting fine-mapped variants with scATAC peaks\n")
cat("============================================================\n\n")

# ── 1. Load aggregated finemapping results ───────────────────────────────────
fm_file <- file.path(FM_DIR, "results/combined_finemapping.csv")
if (!file.exists(fm_file)) {
  stop("combined_finemapping.csv not found. Run 04_aggregate_results.R first.")
}
fm <- fread(fm_file)
cat("Loaded", nrow(fm), "variant-locus entries from combined finemapping\n")

# ── 2. Filter to credible set / high-PIP variants ───────────────────────────
# Keep variants in any credible set (either SuSiE or CARMA) OR PIP > 0.1
# Use recommended_pip (CARMA for EAS, max for EUR converged)
has_recommended <- "recommended_pip" %in% colnames(fm)

fm_cs <- fm %>%
  filter(
    (either_in_cs == TRUE) |
    (!is.na(recommended_pip) & recommended_pip > 0.1) |
    (is.na(recommended_pip) & !is.na(max_pip) & max_pip > 0.1)
  )
cat("Filtered to", nrow(fm_cs), "credible set / high-PIP variant-locus entries\n")

# ── 3. Deduplicate across studies ────────────────────────────────────────────
# Same variant can appear in multiple GWAS — keep max PIP and record studies
pip_col <- if (has_recommended) "recommended_pip" else "max_pip"
has_susie_clean <- "susie_pip_clean" %in% colnames(fm_cs)

variants <- fm_cs %>%
  group_by(chromosome, position, allele1, allele2) %>%
  summarise(
    max_pip      = max(max_pip, na.rm = TRUE),
    max_rec_pip  = if (has_recommended && any(!is.na(recommended_pip))) max(recommended_pip, na.rm = TRUE) else max(max_pip, na.rm = TRUE),
    n_studies    = n_distinct(study),
    study_list   = paste(sort(unique(study)), collapse = ";"),
    n_loci       = n_distinct(locus),
    best_method  = if (has_susie_clean && any(!is.na(susie_pip_clean) & susie_pip_clean > 0.1, na.rm = TRUE)) "susie" else "carma",
    in_susie_cs  = any(!is.na(susie_cs) & susie_cs > 0, na.rm = TRUE),
    in_carma_cs  = any(!is.na(carma_cs) & carma_cs > 0, na.rm = TRUE),
    .groups = "drop"
  )
cat("Deduplicated to", nrow(variants), "unique variants\n")

# ── 4. LiftOver hg19 → hg38 ─────────────────────────────────────────────────
cat("\n--- LiftOver hg19 → hg38 ---\n")

chain_file <- file.path(BASE_DIR, "data/broadaway_eqtl/hg19ToHg38.over.chain")
if (!file.exists(chain_file)) {
  # Try compressed version
  chain_gz <- paste0(chain_file, ".gz")
  if (file.exists(chain_gz)) {
    system(paste("gunzip -k", chain_gz))
  } else {
    stop("liftOver chain file not found at ", chain_file)
  }
}
chain <- import.chain(chain_file)

# Create GRanges in hg19
gr_hg19 <- GRanges(
  seqnames = paste0("chr", variants$chromosome),
  ranges   = IRanges(start = variants$position, width = 1),
  strand   = "*"
)
mcols(gr_hg19) <- variants

# LiftOver
lifted <- liftOver(gr_hg19, chain)
n_mapped <- lengths(lifted)

# Keep only 1:1 mappings
unique_idx <- which(n_mapped == 1)
n_failed <- sum(n_mapped != 1)
cat("LiftOver: ", length(unique_idx), "mapped, ", n_failed, "failed/multi-mapped\n")

gr_hg38 <- unlist(lifted[unique_idx])
variants_hg38 <- as.data.table(mcols(gr_hg38))
variants_hg38$chr_hg38 <- as.character(seqnames(gr_hg38))
variants_hg38$pos_hg38 <- start(gr_hg38)
variants_hg38$pos_hg19 <- variants_hg38$position  # preserve original

cat("Variants in hg38:", nrow(variants_hg38), "\n")

# Save unmapped variants for reference
if (n_failed > 0) {
  failed_idx <- which(n_mapped != 1)
  failed_variants <- variants[failed_idx, ]
  fwrite(failed_variants, file.path(OUT_DIR, "liftover_failed_variants.csv"))
  cat("Failed variants written to liftover_failed_variants.csv\n")
}

# ── 5. Load cell-type peaks ──────────────────────────────────────────────────
cat("\n--- Loading scATAC cell-type peaks ---\n")

peak_dir <- file.path(ATAC_DIR, "results/label_transfer/cell_type_peak_sets_v2")
peak_files <- list.files(peak_dir, pattern = "_peaks\\.bed$", full.names = TRUE)

if (length(peak_files) == 0) {
  stop("No peak BED files found in ", peak_dir)
}

cell_types <- gsub("_peaks\\.bed$", "", basename(peak_files))
cat("Found", length(peak_files), "cell-type peak sets:", paste(cell_types, collapse = ", "), "\n")

# Load all peak sets
peak_list <- list()
for (i in seq_along(peak_files)) {
  ct <- cell_types[i]
  # BED3 format: chr, start, end
  bed <- fread(peak_files[i], header = FALSE)

  if (ncol(bed) == 3) {
    colnames(bed) <- c("chr", "start", "end")
  } else if (ncol(bed) == 4) {
    # 4-col: ID, chr, start, end — skip first column
    colnames(bed) <- c("id", "chr", "start", "end")
    bed <- bed[, .(chr, start, end)]
  } else {
    cat("  WARNING: Unexpected column count in", basename(peak_files[i]), "— skipping\n")
    next
  }

  peak_list[[ct]] <- GRanges(
    seqnames = bed$chr,
    ranges   = IRanges(start = bed$start + 1L, end = bed$end)  # FIX (review): BED is 0-based half-open -> 1-based GRanges (matches 56e/56h)
  )
  cat("  ", ct, ":", length(peak_list[[ct]]), "peaks\n")
}

# ── 6. Intersect variants with peaks ─────────────────────────────────────────
cat("\n--- Intersecting variants with cell-type peaks ---\n")

# Create variant GRanges in hg38
gr_variants <- GRanges(
  seqnames = variants_hg38$chr_hg38,
  ranges   = IRanges(start = variants_hg38$pos_hg38, width = 1)
)

overlap_results <- list()

for (ct in names(peak_list)) {
  hits <- findOverlaps(gr_variants, peak_list[[ct]])
  if (length(hits) == 0) {
    cat("  ", ct, ": 0 overlaps\n")
    next
  }

  q_idx <- queryHits(hits)
  s_idx <- subjectHits(hits)

  ovl <- data.table(
    variant_row = q_idx,
    cell_type   = ct,
    peak_chr    = as.character(seqnames(peak_list[[ct]][s_idx])),
    peak_start  = start(peak_list[[ct]][s_idx]),
    peak_end    = end(peak_list[[ct]][s_idx])
  )
  ovl$peak_width <- ovl$peak_end - ovl$peak_start

  overlap_results[[ct]] <- ovl
  cat("  ", ct, ":", length(unique(q_idx)), "variants in", length(unique(s_idx)), "peaks\n")
}

if (length(overlap_results) == 0) {
  cat("\nWARNING: No variants overlap any peaks. Writing empty results.\n")
  fwrite(data.table(), file.path(OUT_DIR, "gwas_atac_variant_annotation.csv"))
  quit(save = "no", status = 0)
}

overlaps_all <- rbindlist(overlap_results)

# Merge variant metadata
overlaps_annotated <- overlaps_all %>%
  left_join(
    variants_hg38 %>% mutate(variant_row = 1:n()),
    by = "variant_row"
  )

cat("\nTotal variant-peak overlaps:", nrow(overlaps_annotated), "\n")

# Summary per variant: which cell types overlap
variant_ct_summary <- overlaps_annotated %>%
  group_by(chromosome, position, allele1, allele2) %>%
  summarise(
    n_cell_types_overlapping = n_distinct(cell_type),
    cell_types_overlapping   = paste(sort(unique(cell_type)), collapse = ","),
    overlaps_any_peak        = TRUE,
    .groups = "drop"
  )

# Add overlap summary to all variants (including non-overlapping)
variants_hg38 <- variants_hg38 %>%
  left_join(variant_ct_summary, by = c("chromosome", "position", "allele1", "allele2")) %>%
  mutate(
    overlaps_any_peak = replace_na(overlaps_any_peak, FALSE),
    n_cell_types_overlapping = replace_na(n_cell_types_overlapping, 0L),
    cell_types_overlapping = replace_na(cell_types_overlapping, "")
  )

n_in_peak <- sum(variants_hg38$overlaps_any_peak)
cat("Variants overlapping any peak:", n_in_peak, "/", nrow(variants_hg38),
    "(", round(100 * n_in_peak / nrow(variants_hg38), 1), "%)\n")

# ── 7. Annotate with hepatocyte DA ──────────────────────────────────────────
cat("\n--- Annotating with hepatocyte DA ---\n")

da_file <- file.path(ATAC_DIR, "results/snapatac2/scatac_da_corrected_hep.csv")
if (!file.exists(da_file)) {
  # Try alternative location
  da_file <- file.path(ATAC_DIR, "results/l8_annotated_corrected/scatac_da_corrected_hep.csv")
}

if (file.exists(da_file)) {
  da <- fread(da_file)
  cat("Loaded", nrow(da), "DA results\n")

  # Parse peak coordinates from feature name (e.g., "chr1:8026000-8026500")
  da[, c("da_chr", "da_coords") := tstrsplit(`feature name`, ":", fixed = TRUE)]
  da[, c("da_start", "da_end") := tstrsplit(da_coords, "-", fixed = TRUE, type.convert = TRUE)]

  # Create GRanges for DA peaks
  gr_da <- GRanges(
    seqnames = da$da_chr,
    ranges   = IRanges(start = da$da_start + 1L, end = da$da_end)  # FIX (review): BED 0-based half-open -> 1-based GRanges
  )
  mcols(gr_da)$hep_da_logFC <- da$`log2(fold_change)`
  mcols(gr_da)$hep_da_padj  <- da$`adjusted p-value`

  # Overlap DA peaks with variant-overlapping hepatocyte peaks
  hep_overlaps <- overlaps_annotated %>% filter(cell_type == "Hepatocytes")
  if (nrow(hep_overlaps) > 0) {
    gr_hep_peaks <- GRanges(
      seqnames = hep_overlaps$peak_chr,
      ranges   = IRanges(start = hep_overlaps$peak_start, end = hep_overlaps$peak_end)
    )

    hep_overlaps$hep_da_logFC <- NA_real_
    hep_overlaps$hep_da_padj  <- NA_real_
    da_hits <- findOverlaps(gr_hep_peaks, gr_da)
    if (length(da_hits) > 0) {
      hep_overlaps$hep_da_logFC[queryHits(da_hits)] <- mcols(gr_da)$hep_da_logFC[subjectHits(da_hits)]
      hep_overlaps$hep_da_padj[queryHits(da_hits)]  <- mcols(gr_da)$hep_da_padj[subjectHits(da_hits)]
      cat("  Matched", sum(!is.na(hep_overlaps$hep_da_logFC)), "hepatocyte peak-variant overlaps with DA\n")
    }
    # Merge DA info back
    overlaps_annotated <- overlaps_annotated %>%
      left_join(
        hep_overlaps %>%
          select(variant_row, cell_type, peak_start, peak_end, hep_da_logFC, hep_da_padj),
        by = c("variant_row", "cell_type", "peak_start", "peak_end")
      )
  }
} else {
  cat("WARNING: Hepatocyte DA file not found — skipping DA annotation\n")
  overlaps_annotated$hep_da_logFC <- NA_real_
  overlaps_annotated$hep_da_padj  <- NA_real_
}

# ── 8. Annotate with SCENIC+ enhancer-gene links ────────────────────────────
# Note: chromVAR TF motifs are per-TF, not per-peak — peak-specific motif
# disruption analysis is handled in Script 56 (motifbreakR + FIMO).
cat("\n--- Annotating with SCENIC+ enhancer-gene links ---\n")

enh_file <- file.path(ATAC_DIR, "scenic_plus/enhancer_gene_links.csv")
if (file.exists(enh_file)) {
  enh <- fread(enh_file)
  cat("Loaded", nrow(enh), "enhancer-gene links\n")

  # Create GRanges for enhancers
  gr_enh <- GRanges(
    seqnames = enh$enhancer_chr,
    ranges   = IRanges(start = enh$enhancer_start, end = enh$enhancer_end)
  )
  mcols(gr_enh)$scenic_target_gene   <- enh$target_gene
  mcols(gr_enh)$scenic_tf            <- enh$tf_name
  mcols(gr_enh)$scenic_rna_atac_corr <- enh$correlation_rna_atac
  mcols(gr_enh)$scenic_distance_tss  <- enh$distance_to_tss

  # Find which variant-peak overlaps fall in SCENIC+ enhancers
  gr_ovl_peaks <- GRanges(
    seqnames = overlaps_annotated$peak_chr,
    ranges   = IRanges(start = overlaps_annotated$peak_start, end = overlaps_annotated$peak_end)
  )

  enh_hits <- findOverlaps(gr_ovl_peaks, gr_enh)
  overlaps_annotated$scenic_target_gene   <- NA_character_
  overlaps_annotated$scenic_tf            <- NA_character_
  overlaps_annotated$scenic_rna_atac_corr <- NA_real_

  if (length(enh_hits) > 0) {
    # For multiple enhancer matches, take the one with highest RNA-ATAC correlation
    enh_dt <- data.table(
      ovl_row    = queryHits(enh_hits),
      target     = mcols(gr_enh)$scenic_target_gene[subjectHits(enh_hits)],
      tf         = mcols(gr_enh)$scenic_tf[subjectHits(enh_hits)],
      corr       = mcols(gr_enh)$scenic_rna_atac_corr[subjectHits(enh_hits)]
    )
    best_enh <- enh_dt[, .SD[which.max(corr)], by = ovl_row]

    overlaps_annotated$scenic_target_gene[best_enh$ovl_row]   <- best_enh$target
    overlaps_annotated$scenic_tf[best_enh$ovl_row]            <- best_enh$tf
    overlaps_annotated$scenic_rna_atac_corr[best_enh$ovl_row] <- best_enh$corr

    n_enh_linked <- sum(!is.na(overlaps_annotated$scenic_target_gene))
    cat("  Matched", n_enh_linked, "variant-peak overlaps to SCENIC+ enhancers\n")
  }
} else {
  cat("WARNING: Enhancer-gene links file not found — skipping\n")
  overlaps_annotated$scenic_target_gene   <- NA_character_
  overlaps_annotated$scenic_tf            <- NA_character_
  overlaps_annotated$scenic_rna_atac_corr <- NA_real_
}

# ── 10. Annotate with COLOC PP.H4 ───────────────────────────────────────────
cat("\n--- Annotating with COLOC PP.H4 ---\n")

coloc_file <- file.path(FM_DIR, "results/susie_coloc/gene_level_coloc.csv")
if (file.exists(coloc_file)) {
  coloc <- fread(coloc_file)
  # Remove rows with empty gene names (spurious matches)
  coloc <- coloc[gene != "" & !is.na(gene)]
  cat("Loaded", nrow(coloc), "gene-level COLOC results (after filtering empty genes)\n")
} else {
  cat("WARNING: gene_level_coloc.csv not found — skipping\n")
  coloc <- data.table()
}

# ── 11. Nearest gene annotation ─────────────────────────────────────────────
cat("\n--- Annotating with nearest gene ---\n")

gtf_file <- "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
if (file.exists(gtf_file)) {
  genes <- rtracklayer::import(gtf_file, feature.type = "gene")
  # Keep standard chromosomes
  genes <- genes[seqnames(genes) %in% paste0("chr", c(1:22, "X", "Y"))]
  cat("Loaded", length(genes), "genes from GENCODE v49\n")
  # FIX (review P1): distance_to_tss must be to the TSS, not the gene BODY (prior nearest()/distance()
  # used full gene ranges, so ~84% were 0 = inside a gene body). resize(fix="start") is strand-aware
  # (5' end = TSS). nearest_gene is now the nearest-TSS gene (standard regulatory annotation).
  tss <- resize(genes, width = 1L, fix = "start")

  # For each overlapping variant-peak pair, find nearest gene by TSS
  gr_var_pos <- GRanges(
    seqnames = overlaps_annotated$chr_hg38,
    ranges   = IRanges(start = overlaps_annotated$pos_hg38, width = 1)
  )

  nearest_idx <- nearest(gr_var_pos, tss, ignore.strand = TRUE)
  overlaps_annotated$nearest_gene    <- ifelse(!is.na(nearest_idx),
                                                genes$gene_name[nearest_idx], NA_character_)
  overlaps_annotated$distance_to_tss <- ifelse(!is.na(nearest_idx),
                                                distance(gr_var_pos, tss[nearest_idx]), NA_integer_)

  cat("  Annotated", sum(!is.na(overlaps_annotated$nearest_gene)), "variant-peak overlaps with nearest gene\n")
} else {
  cat("WARNING: GENCODE GTF not found — skipping nearest gene\n")
  overlaps_annotated$nearest_gene    <- NA_character_
  overlaps_annotated$distance_to_tss <- NA_integer_
}

# ── 12. Add COLOC for linked genes ──────────────────────────────────────────
if (nrow(coloc) > 0) {
  # Assign linked gene: SCENIC+ target > nearest gene
  overlaps_annotated <- overlaps_annotated %>%
    mutate(
      linked_gene = ifelse(!is.na(scenic_target_gene), scenic_target_gene, nearest_gene)
    )

  # Look up COLOC
  coloc_lookup <- coloc %>%
    select(gene, coloc_best_pp4, coloc_best_gwas, coloc_n_gwas_h4_05) %>%
    distinct(gene, .keep_all = TRUE)

  overlaps_annotated <- overlaps_annotated %>%
    left_join(coloc_lookup, by = c("linked_gene" = "gene"))

  n_coloc <- sum(!is.na(overlaps_annotated$coloc_best_pp4))
  cat("  Matched", n_coloc, "variant-peak overlaps with COLOC PP.H4\n")
} else {
  overlaps_annotated$linked_gene       <- overlaps_annotated$nearest_gene
  overlaps_annotated$coloc_best_pp4    <- NA_real_
  overlaps_annotated$coloc_best_gwas   <- NA_character_
  overlaps_annotated$coloc_n_gwas_h4_05 <- NA_integer_
}

# ── 13. Create variant ID and finalize ───────────────────────────────────────
overlaps_annotated <- overlaps_annotated %>%
  mutate(
    variant_id = paste(chromosome, position, allele1, allele2, sep = ":")
  ) %>%
  select(
    variant_id, chromosome, position, allele1, allele2,
    chr_hg38, pos_hg38, pos_hg19,
    max_pip, max_rec_pip, n_studies, study_list, best_method,
    in_susie_cs, in_carma_cs,
    cell_type, peak_chr, peak_start, peak_end, peak_width,
    hep_da_logFC, hep_da_padj,
    scenic_target_gene, scenic_tf, scenic_rna_atac_corr,
    linked_gene, nearest_gene, distance_to_tss,
    coloc_best_pp4, coloc_best_gwas, coloc_n_gwas_h4_05
  )

# ── 14. Write outputs ───────────────────────────────────────────────────────
cat("\n=== Writing outputs ===\n")

# Main annotation table: one row per variant × cell type overlap
fwrite(overlaps_annotated, file.path(OUT_DIR, "gwas_atac_variant_annotation.csv"))
cat("Variant annotation:", nrow(overlaps_annotated), "rows written\n")

# Variant-level summary (one row per unique variant)
variant_summary <- variants_hg38 %>%
  select(
    chromosome, position, allele1, allele2,
    chr_hg38, pos_hg38,
    max_pip, max_rec_pip, n_studies, study_list, best_method,
    overlaps_any_peak, n_cell_types_overlapping, cell_types_overlapping
  )
fwrite(variant_summary, file.path(OUT_DIR, "variant_overlap_summary.csv"))
cat("Variant summary:", nrow(variant_summary), "rows written\n")

# Quick stats
cat("\n=== Summary Statistics ===\n")
cat("Total unique variants (high-PIP / CS):", nrow(variants_hg38), "\n")
cat("  Successfully lifted to hg38:", nrow(variants_hg38), "\n")
cat("  Overlapping any peak:", n_in_peak, "(", round(100 * n_in_peak / nrow(variants_hg38), 1), "%)\n")
for (ct in names(peak_list)) {
  n_ct <- sum(overlaps_annotated$cell_type == ct)
  if (n_ct > 0) cat("    ", ct, ":", n_ct, "variant-peak overlaps\n")
}
cat("  With DA annotation:", sum(!is.na(overlaps_annotated$hep_da_logFC)), "\n")
cat("  With SCENIC+ enhancer link:", sum(!is.na(overlaps_annotated$scenic_target_gene)), "\n")
cat("  With COLOC PP.H4:", sum(!is.na(overlaps_annotated$coloc_best_pp4)), "\n")

cat("\n============================================================\n")
cat("Script 55 complete. Results in:", OUT_DIR, "\n")
cat("============================================================\n")
