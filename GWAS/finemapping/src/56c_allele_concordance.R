#!/usr/bin/env Rscript
# 56c_allele_concordance.R
# A2 step from the ATAC improvement plan.
#
# For every (variant, TF, peak, linked_gene) tuple in motif_disruption_scores.csv,
# evaluate 3-way allele-direction concordance between:
#   1. motif disruption    -- sign(alleleDiff) from motifbreakR (+ = ALT increases binding)
#   2. eQTL effect         -- Broadaway liver eQTL Beta (+ = ALT increases linked-gene expression)
#   3. peak DA in disease  -- hep_da_logFC from scATAC (+ = peak more accessible in disease)
#
# TF activator/repressor expectation: activator => motif-up should accompany expression-up.
# Repressor => motif-up should accompany expression-DOWN (flip eQTL expectation).
# Disease direction: motif-up + peak-up + (activator ? expr-up : expr-down) all align in disease.
#
# Inputs (read-only):
#   GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#   GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv
#   data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv  (hg19, primary lookup)
#   data/broadaway_eqtl/chr{1..22}_marginal_summary_results.tsv (hg19, fallback)
#   Analysis/ATAC/Human_Multiome/results/snapatac2/scatac_da_corrected_hep.csv  (hg38)
#   Analysis/ATAC/Human_Multiome/scenic_plus/disease_regulons.csv
#
# Outputs:
#   GWAS/finemapping/results/gwas_atac/allele_concordance.csv      (per-tuple)
#   GWAS/finemapping/results/gwas_atac/allele_concordance_summary.csv  (per-TF / per-gene)
#
# Conventions:
#   * Motif SNP_id is hg19 (e.g. "6:32191581:A:T"), motifbreakR seqnames/start are hg38.
#   * Broadaway Beta is in EA orientation -- we flip to ALT orientation explicitly.
#   * DA tiles are 500bp hg38 windows in scatac_da_corrected_hep.csv.
#   * linked_gene comes from gwas_atac_variant_annotation.csv (SCENIC+/HiC link); fall back to nearest_gene.
#   * "concordance" is in the *disease-progression* frame: a disease-increasing motif change
#     should match a disease-increasing accessibility and the expected expression direction
#     given activator vs repressor.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR     <- file.path(BASE_DIR, "GWAS/finemapping")
ATAC_DIR   <- file.path(BASE_DIR, "Analysis/ATAC/Human_Multiome")
EQTL_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
OUT_DIR    <- file.path(FM_DIR, "results/gwas_atac")

cat("============================================================\n")
cat("56c_allele_concordance.R  --  A2 3-way concordance\n")
cat("============================================================\n\n")

# ── 1. Load motif disruption + variant annotation ─────────────────────────────
md_file       <- file.path(OUT_DIR, "motif_disruption_scores.csv")
ann_file      <- file.path(OUT_DIR, "gwas_atac_variant_annotation.csv")
reg_file      <- file.path(ATAC_DIR, "scenic_plus/disease_regulons.csv")
da_hep_file   <- file.path(ATAC_DIR, "results/snapatac2/scatac_da_corrected_hep.csv")
da_all_file   <- file.path(ATAC_DIR, "results/snapatac2/scatac_da_results.csv")
leads_file    <- file.path(EQTL_DIR, "Liver_eQTL_Meta_Leads_ST3_20240530.tsv")

stopifnot(file.exists(md_file), file.exists(ann_file),
          file.exists(reg_file), file.exists(da_hep_file),
          file.exists(da_all_file), file.exists(leads_file))

md     <- fread(md_file)
ann    <- fread(ann_file)
regs   <- fread(reg_file)
da     <- fread(da_hep_file)    # Hep-only DA tiles (canonical)
da_all <- fread(da_all_file)    # Multi-cell-type DA tiles
cat("Loaded motif disruption rows: ", nrow(md), "\n", sep = "")
cat("Loaded variant annotations:   ", nrow(ann), "\n", sep = "")
cat("Loaded disease regulons:      ", nrow(regs), "\n", sep = "")
cat("Loaded Hep DA tiles:          ", nrow(da), "\n", sep = "")
cat("Loaded multi-CT DA tiles:     ", nrow(da_all), "\n", sep = "")

# ── 2. Build variant -> (linked_gene, hep_da_logFC) lookup ───────────────────
ann[, linked_gene_resolved := ifelse(linked_gene != "" & !is.na(linked_gene),
                                     linked_gene, nearest_gene)]
# Collapse to one row per variant (max PIP wins; SCENIC TF carried)
ann_uniq <- ann[order(-max_pip),
                .SD[1],
                by = .(variant_id)]
n_linked <- sum(ann_uniq$linked_gene_resolved != "" & !is.na(ann_uniq$linked_gene_resolved))
cat("Variants with linked_gene OR nearest_gene: ", n_linked, " / ", nrow(ann_uniq), "\n", sep = "")

# ── 3. Build chr/pos -> hg38 DA tile lookup via genomic overlap ──────────────
# DA tiles are 500bp hg38 windows (snapatac2). Match each MD variant's hg38 pos.
#
# We use TWO sources of DA values:
#  (a) scatac_da_corrected_hep.csv -- Hep-only canonical DA (1681 tiles), used
#      for the canonical hepatocyte-axis concordance test.
#  (b) scatac_da_results.csv      -- multi-cell-type DA (11686 tiles across
#      10 CTs). When (a) misses, we fall back to the DA value in the variant's
#      OWN cell type (carried via the variant annotation file).
#
# motif disruption rows are in hg38 coordinates (start col)
da[, c("chrom", "tile_start", "tile_end") := tstrsplit(`feature name`, "[:-]", fixed = FALSE)]
da[, tile_start := as.integer(tile_start)]
da[, tile_end   := as.integer(tile_end)]
setnames(da, "log2(fold_change)", "hep_da_logFC")
setnames(da, "adjusted p-value",  "hep_da_padj")

gr_md <- GRanges(seqnames = md$seqnames,
                 ranges   = IRanges(start = md$start, width = 1L))
gr_da <- GRanges(seqnames = da$chrom,
                 ranges   = IRanges(start = da$tile_start, end = da$tile_end))
ovl <- findOverlaps(gr_md, gr_da)
md[, hep_da_logFC := NA_real_]
md[, hep_da_padj  := NA_real_]
md[, hep_da_tile  := NA_character_]
if (length(ovl) > 0) {
  md$hep_da_logFC[queryHits(ovl)] <- da$hep_da_logFC[subjectHits(ovl)]
  md$hep_da_padj[queryHits(ovl)]  <- da$hep_da_padj[subjectHits(ovl)]
  md$hep_da_tile[queryHits(ovl)]  <- da$`feature name`[subjectHits(ovl)]
}
n_da_hits <- sum(!is.na(md$hep_da_logFC))
cat("Motif-disruption rows with Hep DA tile match: ", n_da_hits, " / ", nrow(md),
    "  (", round(100 * n_da_hits / nrow(md), 1), "%)\n", sep = "")

# Multi-cell-type DA -- fallback. Match by genomic overlap; keep ALL hits, then
# we'll pick the variant's own cell type later via the ann file's cell_type column.
da_all[, c("chrom", "tile_start", "tile_end") := tstrsplit(`feature name`, "[:-]", fixed = FALSE)]
da_all[, tile_start := as.integer(tile_start)]
da_all[, tile_end   := as.integer(tile_end)]
setnames(da_all, "log2(fold_change)", "ct_da_logFC")
setnames(da_all, "adjusted p-value",  "ct_da_padj")
setnames(da_all, "cell_type",         "ct_da_cell_type")

gr_da_all <- GRanges(seqnames = da_all$chrom,
                     ranges   = IRanges(start = da_all$tile_start, end = da_all$tile_end))
ovl_all <- findOverlaps(gr_md, gr_da_all)
# Build a per-MD-row -> list of (cell_type, logFC, padj) by aggregating overlaps.
if (length(ovl_all) > 0) {
  # Pivot: one row per (md_idx, cell_type), keep most significant per pair
  ovl_dt <- data.table(md_idx        = queryHits(ovl_all),
                       ct_cell_type  = da_all$ct_da_cell_type[subjectHits(ovl_all)],
                       ct_da_logFC   = da_all$ct_da_logFC[subjectHits(ovl_all)],
                       ct_da_padj    = da_all$ct_da_padj[subjectHits(ovl_all)],
                       ct_da_tile    = da_all$`feature name`[subjectHits(ovl_all)])
  ovl_dt <- ovl_dt[order(md_idx, ct_cell_type, ct_da_padj),
                   .SD[1], by = .(md_idx, ct_cell_type)]
  # Carry these as JSON-like collapsed columns (per md row, the union of CTs)
  ovl_summary <- ovl_dt[, .(ct_da_logFC_byct = paste(ct_cell_type, round(ct_da_logFC, 3), sep = ":", collapse = ";"),
                            ct_da_padj_byct  = paste(ct_cell_type, signif(ct_da_padj, 3), sep = ":", collapse = ";")),
                        by = md_idx]
  md[, ct_da_logFC_byct := NA_character_]
  md[, ct_da_padj_byct  := NA_character_]
  md$ct_da_logFC_byct[ovl_summary$md_idx] <- ovl_summary$ct_da_logFC_byct
  md$ct_da_padj_byct[ovl_summary$md_idx]  <- ovl_summary$ct_da_padj_byct
  cat("Motif-disruption rows with ANY multi-CT DA overlap: ",
      length(unique(ovl_dt$md_idx)), " / ", nrow(md),
      "  (", round(100 * length(unique(ovl_dt$md_idx)) / nrow(md), 1), "%)\n", sep = "")
} else {
  md[, ct_da_logFC_byct := NA_character_]
  md[, ct_da_padj_byct  := NA_character_]
  cat("No multi-CT DA overlaps found.\n")
}

# ── 4. Join SCENIC+ linked-gene & hep_da from variant annotation (preferred) ─
# Use the annotation file's hep_da_logFC where populated (variant-level call peaks);
# otherwise keep the tile-overlap value computed above.
md[ann_uniq, on = .(SNP_id = variant_id),
   `:=`(linked_gene_resolved = i.linked_gene_resolved,
        scenic_tf            = i.scenic_tf,
        scenic_target_gene   = i.scenic_target_gene,
        hep_da_logFC_ann     = i.hep_da_logFC,
        hep_da_padj_ann      = i.hep_da_padj,
        coloc_best_pp4       = i.coloc_best_pp4,
        variant_cell_type    = i.cell_type)]

# Prefer annotation hep_da (variant-level call peak overlap, snapATAC2-corrected)
# Fall back to tile-based overlap from corrected_hep.csv
md[, hep_da_logFC := ifelse(!is.na(hep_da_logFC_ann), hep_da_logFC_ann, hep_da_logFC)]
md[, hep_da_padj  := ifelse(!is.na(hep_da_padj_ann),  hep_da_padj_ann,  hep_da_padj)]

# ── 4b. Cell-type-matched DA + fallback DA
# Strategy:
#   (1) Try the variant's OWN cell-type DA in the multi-CT pool.
#   (2) If absent, use the most-significant CT-DA at that locus (any cell type).
#       Rationale: the DA file is sparse (only tested significant changes), and
#       a peak's accessibility change anywhere in liver immune/parenchymal cells
#       is informative for the cis-regulatory direction.
# Ann file uses "Hepatocytes" (plural), multi-CT DA uses "Hepatocyte" (singular).
ct_alias <- c("Hepatocytes" = "Hepatocyte",
              "Cholangiocytes" = "Cholangiocyte",
              "Macrophages" = "Macrophage",
              "Plasma_cells" = "Plasma_Cell",
              "Endothelial_cells" = "Endothelial",
              "T_cells" = "NK_T_Cell",
              "Resident_NK" = "NK_T_Cell",
              "Fibroblasts" = "Stellate_Cell",
              "B_cells" = "B_Cell")
md[, ct_da_match_key := ct_alias[variant_cell_type]]
md[is.na(ct_da_match_key), ct_da_match_key := variant_cell_type]

md[, matched_ct_da_logFC := NA_real_]
md[, matched_ct_da_padj  := NA_real_]
md[, matched_ct_da_ct    := NA_character_]
md[, anyct_da_logFC      := NA_real_]
md[, anyct_da_padj       := NA_real_]
md[, anyct_da_ct         := NA_character_]
if (exists("ovl_dt") && nrow(ovl_dt) > 0) {
  # (1) Own-CT match
  md_keys <- data.table(md_idx = seq_len(nrow(md)),
                        ct_da_match_key = md$ct_da_match_key)
  matched <- ovl_dt[md_keys, on = .(md_idx, ct_cell_type = ct_da_match_key),
                    nomatch = NULL]
  if (nrow(matched) > 0) {
    md$matched_ct_da_logFC[matched$md_idx] <- matched$ct_da_logFC
    md$matched_ct_da_padj[matched$md_idx]  <- matched$ct_da_padj
    md$matched_ct_da_ct[matched$md_idx]    <- matched$ct_cell_type
    cat("Own-cell-type DA hits: ", nrow(matched), "\n", sep = "")
  }
  # (2) Any-CT fallback: pick most significant CT per md_idx
  any_ct <- ovl_dt[order(md_idx, ct_da_padj), .SD[1], by = md_idx]
  md$anyct_da_logFC[any_ct$md_idx] <- any_ct$ct_da_logFC
  md$anyct_da_padj[any_ct$md_idx]  <- any_ct$ct_da_padj
  md$anyct_da_ct[any_ct$md_idx]    <- any_ct$ct_cell_type
  cat("Any-CT DA hits (most-sig per locus): ", nrow(any_ct), "\n", sep = "")
}

# Final DA value: prefer Hep canonical, then own-CT match, then any-CT fallback
md[, da_logFC_final := fifelse(!is.na(hep_da_logFC), hep_da_logFC,
                       fifelse(!is.na(matched_ct_da_logFC), matched_ct_da_logFC,
                       fifelse(!is.na(anyct_da_logFC), anyct_da_logFC, NA_real_)))]
md[, da_padj_final  := fifelse(!is.na(hep_da_padj), hep_da_padj,
                       fifelse(!is.na(matched_ct_da_padj), matched_ct_da_padj,
                       fifelse(!is.na(anyct_da_padj), anyct_da_padj, NA_real_)))]
md[, da_source      := fifelse(!is.na(hep_da_logFC), "Hep_canonical",
                       fifelse(!is.na(matched_ct_da_logFC),
                               paste0("CT_own_", matched_ct_da_ct),
                       fifelse(!is.na(anyct_da_logFC),
                               paste0("CT_any_", anyct_da_ct), NA_character_)))]
cat("Rows with final DA logFC (Hep+own-CT+any-CT): ",
    sum(!is.na(md$da_logFC_final)), " / ", nrow(md), "\n", sep = "")

# ── 5. Activator / repressor table for the 24 SCENIC+ disease regulons ──────
# Built from primary literature + regulon_activity_diff sign. Activators upregulate
# target gene when binding increases; repressors downregulate. Default = activator.
activator_table <- data.table(
  tf_name = c("HNF4A","NR1H4","THRB","RORA","CEBPA","CEBPB","FOXA1","FOXA2",
              "PPARA","PPARG","ESR1","AR","AHR","HNF1A","KLF15","RXRA",
              "MAX","SMAD4","MLXIPL","ARNT","RORC","NR1H3","XBP1","ETS1","ETS2",
              # Add common TFs the motifs hit
              "TBX3","MAFB","CTCF","TCF7L2","FOXO1","FOXO3","SREBF1","SREBF2",
              "STAT3","STAT1","NR2F1","NR2F2","RARA","REST"),
  is_activator = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                   TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                   TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                   # Context TFs
                   TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                   TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
  source = "literature"
)
# Merge in regulon-derived activator/repressor (sign of regulon_activity_diff
# tells us whether the TF tracks disease state -- positive diff => more active
# in MASLD => the TF is an activator in this context if its target gene is up
# in disease; we keep the literature call as primary).
md[, tf_upper := toupper(tf_name)]
md[activator_table, on = .(tf_upper = tf_name),
   is_activator := i.is_activator]
md[is.na(is_activator), is_activator := TRUE]   # default = activator

# ── 6. eQTL lookup -- leads first (light), then per-chr marginal (heavy) ────
# Broadaway is hg19; motif SNP_id is hg19 in the form "6:32191581:A:T".
# Convert to Broadaway variant id "6_32191581_A_T" (a1->NEA, a2->EA may differ).
md[, c("chr_h19", "pos_h19", "snp_a1", "snp_a2") :=
     tstrsplit(SNP_id, ":", fixed = TRUE)]
md[, pos_h19 := as.integer(pos_h19)]
# Broadaway variant key uses an unordered allele pair; we'll match in either order.
md[, var_key_a := paste(chr_h19, pos_h19, snp_a1, snp_a2, sep = "_")]
md[, var_key_b := paste(chr_h19, pos_h19, snp_a2, snp_a1, sep = "_")]

# FIX (review A11#1/B1#1): the eQTL Beta must be aligned to the GENOME ALT (motifbreakR's
# alleleDiff is signed REF->ALT), NOT to snp_a2 (the GWAS allele2 parsed from SNP_id, which
# is the genome ALT only when allele1 == genome REF). Script 56 now exports alt_genome.
if ("alt_genome" %in% names(md)) {
  md[, alt_align := toupper(alt_genome)]
} else {
  warning("56c: motif_disruption_scores.csv lacks alt_genome — re-run Script 56. ",
          "Falling back to GWAS allele2 (snp_a2) for eQTL beta alignment (may mis-sign).")
  md[, alt_align := toupper(snp_a2)]
}

# 6a. Leads (1 row per gene-signal):
leads <- fread(leads_file, na.strings = c("NA", ""))
leads <- leads[, .(Gene, Variant, EA, Marginal_Beta, Marginal_SE,
                   Marginal_Pvalue, FDR)]
leads[, chr := tstrsplit(Variant, "_", fixed = TRUE)[[1]]]

# Index: gene+variant key
md[, eqtl_beta_lead     := NA_real_]
md[, eqtl_se_lead       := NA_real_]
md[, eqtl_pval_lead     := NA_real_]
md[, eqtl_ea_lead       := NA_character_]
md[, eqtl_source        := NA_character_]

# Need linked_gene_resolved present and chr matched.
# NOTE: data.table's `.I` inside `dt[i, j]` with a filter returns positions
# within the FILTERED subset (1..n), NOT original-table indices. We must
# compute the original-table positions explicitly via which().
keep_idx <- which(!is.na(md$linked_gene_resolved) & md$linked_gene_resolved != "")
md_lookup <- data.table(idx                  = keep_idx,
                        chr_h19              = md$chr_h19[keep_idx],
                        pos_h19              = md$pos_h19[keep_idx],
                        snp_a1               = md$snp_a1[keep_idx],
                        snp_a2               = md$snp_a2[keep_idx],
                        alt_align            = md$alt_align[keep_idx],
                        linked_gene_resolved = md$linked_gene_resolved[keep_idx],
                        var_key_a            = md$var_key_a[keep_idx],
                        var_key_b            = md$var_key_b[keep_idx])
cat("Looking up eQTL via leads file: ", nrow(md_lookup), " rows...\n", sep = "")

# Build leads index: key = gene + variant
leads[, key_a := paste(Gene, Variant, sep = "|")]
setkey(leads, key_a)
md_lookup[, key_a := paste(linked_gene_resolved, var_key_a, sep = "|")]
md_lookup[, key_b := paste(linked_gene_resolved, var_key_b, sep = "|")]

# Join on key_a (snp_a1 = ref); also try key_b (snp_a2 = ref).
# Joins preserve all md_lookup columns (including idx and snp_a2).
hit_a <- leads[md_lookup, on = .(key_a = key_a), nomatch = NULL]
hit_b <- leads[md_lookup, on = .(key_a = key_b), nomatch = NULL]
n_lead_a <- nrow(hit_a)
n_lead_b <- nrow(hit_b)
cat("  leads matches (a1=ref orientation): ", n_lead_a, "\n", sep = "")
cat("  leads matches (a2=ref orientation): ", n_lead_b, "\n", sep = "")

apply_lead_hit <- function(hit) {
  if (nrow(hit) == 0) return(invisible(NULL))
  # FIX (review A11#1): align Beta to the GENOME ALT (alt_align), matching motifbreakR's
  # REF->ALT alleleDiff orientation. Keep sign if EA == genome ALT, else flip.
  beta_aligned <- ifelse(toupper(hit$EA) == hit$alt_align,
                         hit$Marginal_Beta, -hit$Marginal_Beta)
  md[hit$idx, `:=`(eqtl_beta_lead = beta_aligned,
                   eqtl_se_lead   = hit$Marginal_SE,
                   eqtl_pval_lead = hit$Marginal_Pvalue,
                   eqtl_ea_lead   = hit$EA,
                   eqtl_source    = "leads")]
}
apply_lead_hit(hit_a)
apply_lead_hit(hit_b)

n_with_lead <- sum(!is.na(md$eqtl_beta_lead))
cat("Variants with leads eQTL Beta: ", n_with_lead, " / ", nrow(md), "\n", sep = "")

# 6b. Per-chromosome marginal fallback for rows without a lead hit.
# Same .I caveat as above -- use which() for explicit original-row indices.
need_marg <- which(is.na(md$eqtl_beta_lead) &
                    !is.na(md$linked_gene_resolved) &
                    md$linked_gene_resolved != "")
need_marg_dt <- data.table(idx                  = need_marg,
                           chr_h19              = md$chr_h19[need_marg],
                           pos_h19              = md$pos_h19[need_marg],
                           snp_a1               = md$snp_a1[need_marg],
                           snp_a2               = md$snp_a2[need_marg],
                           alt_align            = md$alt_align[need_marg],
                           linked_gene_resolved = md$linked_gene_resolved[need_marg],
                           var_key_a            = md$var_key_a[need_marg],
                           var_key_b            = md$var_key_b[need_marg])
cat("Rows needing marginal fallback: ", nrow(need_marg_dt), "\n", sep = "")

# Only read each chr once -- group by chr
need_marg_dt[, chr_num := suppressWarnings(as.integer(chr_h19))]
need_marg_dt <- need_marg_dt[!is.na(chr_num) & chr_num >= 1 & chr_num <= 22]
cat("Rows needing marginal (autosomal only): ", nrow(need_marg_dt), "\n", sep = "")

if (nrow(need_marg_dt) > 0) {
  for (chr_i in sort(unique(need_marg_dt$chr_num))) {
    chr_file <- file.path(EQTL_DIR, sprintf("chr%d_marginal_summary_results.tsv", chr_i))
    if (!file.exists(chr_file)) next
    cat("  chr", chr_i, ": reading marginal file (", round(file.info(chr_file)$size / 1e9, 2),
        " GB)...\n", sep = "")
    # Only read needed columns
    marg <- fread(chr_file, select = c("Variant", "EA", "Beta", "SE", "PVAL",
                                       "GeneSymbol"),
                  showProgress = FALSE)
    # Build lookup key: gene|variant (chr_pos_a1_a2)
    marg[, key := paste(GeneSymbol, Variant, sep = "|")]
    sub <- need_marg_dt[chr_num == chr_i]
    sub[, key_a := paste(linked_gene_resolved, var_key_a, sep = "|")]
    sub[, key_b := paste(linked_gene_resolved, var_key_b, sep = "|")]

    hit_ma <- marg[sub, on = .(key = key_a), nomatch = NULL]
    hit_mb <- marg[sub, on = .(key = key_b), nomatch = NULL]

    if (nrow(hit_ma) > 0) {
      # FIX (review A11#1): align to genome ALT (alt_align), not GWAS allele2.
      beta_aligned <- ifelse(toupper(hit_ma$EA) == hit_ma$alt_align,
                             hit_ma$Beta, -hit_ma$Beta)
      md[hit_ma$idx, `:=`(eqtl_beta_lead = beta_aligned,
                          eqtl_se_lead   = hit_ma$SE,
                          eqtl_pval_lead = hit_ma$PVAL,
                          eqtl_ea_lead   = hit_ma$EA,
                          eqtl_source    = "marginal")]
    }
    if (nrow(hit_mb) > 0) {
      # FIX (review A11#1): align to genome ALT (alt_align), not GWAS allele2.
      beta_aligned <- ifelse(toupper(hit_mb$EA) == hit_mb$alt_align,
                             hit_mb$Beta, -hit_mb$Beta)
      md[hit_mb$idx, `:=`(eqtl_beta_lead = beta_aligned,
                          eqtl_se_lead   = hit_mb$SE,
                          eqtl_pval_lead = hit_mb$PVAL,
                          eqtl_ea_lead   = hit_mb$EA,
                          eqtl_source    = "marginal")]
    }
    cat("    chr", chr_i, " marginal a-orient hits: ", nrow(hit_ma),
        "  b-orient hits: ", nrow(hit_mb), "\n", sep = "")
    rm(marg, sub, hit_ma, hit_mb); gc(verbose = FALSE)
  }
}

n_with_eqtl <- sum(!is.na(md$eqtl_beta_lead))
cat("Total variants with eQTL Beta (leads + marginal): ", n_with_eqtl, " / ",
    nrow(md), "\n", sep = "")

# ── 7. Concordance assessment ────────────────────────────────────────────────
# Disease frame:
#   sign(motif_disruption_in_disease) = sign(alleleDiff)
#       (positive alleleDiff = ALT increases binding score)
#   sign(peak_in_disease)             = sign(hep_da_logFC)
#       (positive logFC = peak gains accessibility in MASLD vs healthy)
#   sign(expr_in_disease)             = sign(eqtl_beta_alt) * (activator ? +1 : -1)
#       (ALT-aligned Beta; if TF is repressor, flip the expected expression direction)
#
# 3-of-3 strict: all three signs equal in the disease frame.
# 2-of-3 lenient: at least two signs agree.

md[, sign_motif := sign(alleleDiff)]
# Use the final DA logFC (Hep canonical with CT-matched fallback).
md[, sign_da    := sign(da_logFC_final)]
md[, sign_eqtl  := sign(eqtl_beta_lead)]

# For repressor TFs, the expected co-direction in disease is FLIPPED for expression.
# We encode this by storing an "effective" expression sign that aligns with the
# motif-binding direction under the activator assumption.
md[, sign_eqtl_effective := ifelse(is_activator, sign_eqtl, -sign_eqtl)]

# Pairwise concordance booleans
md[, concord_motif_eqtl := !is.na(sign_motif) & !is.na(sign_eqtl_effective) &
                            sign_motif == sign_eqtl_effective]
md[, concord_motif_da   := !is.na(sign_motif) & !is.na(sign_da) &
                            sign_motif == sign_da]
md[, concord_da_eqtl    := !is.na(sign_da) & !is.na(sign_eqtl_effective) &
                            sign_da == sign_eqtl_effective]

# 3-of-3 strict (all three present AND all three concordant)
md[, n_signs_present := (!is.na(sign_motif)) + (!is.na(sign_da)) + (!is.na(sign_eqtl_effective))]
md[, concord_3of3 := n_signs_present == 3 &
                      sign_motif == sign_da & sign_motif == sign_eqtl_effective]
# 2-of-3 lenient
md[, n_pairs_concord := as.integer(concord_motif_eqtl) +
                         as.integer(concord_motif_da) +
                         as.integer(concord_da_eqtl)]
# When all three present, a 2-of-3 requires >=2 pair-concordances out of 3 pairs.
# When only two signs present, the single pair must agree.
md[, concord_2of3 := (n_signs_present == 3 & n_pairs_concord >= 2) |
                     (n_signs_present == 2 & n_pairs_concord >= 1)]

# Cascade label (descriptive)
md[, cascade := NA_character_]
md[concord_3of3 & sign_motif < 0, cascade := "LOF"]
md[concord_3of3 & sign_motif > 0, cascade := "GOF"]
md[!concord_3of3 & concord_2of3, cascade := "partial"]
md[n_signs_present == 3 & n_pairs_concord == 0, cascade := "incoherent"]

# ── 8. Write per-tuple output ────────────────────────────────────────────────
out_cols <- c("SNP_id", "seqnames", "start",
              "tf_name", "is_activator",
              "variant_cell_type",
              "linked_gene_resolved", "scenic_tf", "scenic_target_gene",
              "effect", "alleleDiff", "Refpvalue", "Altpvalue",
              "motif_in_disease_regulon", "regulon_activity_diff",
              "max_pip", "priority_score", "coloc_best_pp4",
              "hep_da_logFC", "hep_da_padj", "hep_da_tile",
              "matched_ct_da_logFC", "matched_ct_da_padj", "matched_ct_da_ct",
              "anyct_da_logFC", "anyct_da_padj", "anyct_da_ct",
              "ct_da_logFC_byct",
              "da_logFC_final", "da_padj_final", "da_source",
              "eqtl_beta_lead", "eqtl_se_lead", "eqtl_pval_lead",
              "eqtl_ea_lead", "eqtl_source",
              "sign_motif", "sign_da", "sign_eqtl", "sign_eqtl_effective",
              "concord_motif_eqtl", "concord_motif_da", "concord_da_eqtl",
              "n_signs_present", "n_pairs_concord",
              "concord_2of3", "concord_3of3", "cascade")
out_cols <- intersect(out_cols, names(md))
out <- md[, ..out_cols]
setorder(out, -concord_3of3, -concord_2of3, -priority_score)

out_file <- file.path(OUT_DIR, "allele_concordance.csv")
fwrite(out, out_file)
cat("\nWrote per-tuple concordance:", out_file, "\n  rows: ", nrow(out), "\n", sep = "")

# ── 9. Write per-TF / per-gene summary ───────────────────────────────────────
per_tf <- out[, .(n_total            = .N,
                  n_with_eqtl        = sum(!is.na(eqtl_beta_lead)),
                  n_with_da          = sum(!is.na(da_logFC_final)),
                  n_with_all_3       = sum(n_signs_present == 3),
                  n_concord_3of3     = sum(concord_3of3, na.rm = TRUE),
                  n_concord_2of3     = sum(concord_2of3, na.rm = TRUE),
                  n_LOF              = sum(cascade == "LOF",  na.rm = TRUE),
                  n_GOF              = sum(cascade == "GOF",  na.rm = TRUE),
                  is_activator       = unique(is_activator)[1],
                  in_disease_regulon = any(motif_in_disease_regulon, na.rm = TRUE)),
              by = tf_name][order(-n_concord_3of3, -n_concord_2of3)]
fwrite(per_tf, file.path(OUT_DIR, "allele_concordance_summary.csv"))

per_gene <- out[!is.na(linked_gene_resolved) & linked_gene_resolved != "",
                .(n_total       = .N,
                  n_tfs_distinct = length(unique(tf_name)),
                  n_concord_3of3 = sum(concord_3of3, na.rm = TRUE),
                  n_concord_2of3 = sum(concord_2of3, na.rm = TRUE),
                  tf_concord_3of3 = paste(unique(tf_name[concord_3of3]), collapse = ";")),
                by = linked_gene_resolved][order(-n_concord_3of3, -n_concord_2of3)]
fwrite(per_gene, file.path(OUT_DIR, "allele_concordance_per_gene.csv"))

# ── 10. Report ───────────────────────────────────────────────────────────────
cat("\n============================================================\n")
cat("3-of-3 concordant tuples: ", sum(out$concord_3of3, na.rm = TRUE), "\n", sep = "")
cat("2-of-3 concordant tuples: ", sum(out$concord_2of3, na.rm = TRUE), "\n", sep = "")
cat("Rows with all 3 signs present: ", sum(out$n_signs_present == 3), "\n", sep = "")
cat("Cascades:\n")
print(table(out$cascade, useNA = "ifany"))

cat("\nTop 10 TFs by 3-of-3 count:\n")
print(head(per_tf[n_concord_3of3 > 0,
                  .(tf_name, n_concord_3of3, n_concord_2of3, n_total,
                    in_disease_regulon, is_activator)], 10))

# Disease-regulon TFs specifically
disease_tfs <- toupper(regs$tf_name)
per_tf_disease <- per_tf[toupper(tf_name) %in% disease_tfs]
cat("\nDisease-regulon TF concordance:\n")
print(per_tf_disease[, .(tf_name, n_total, n_concord_3of3, n_concord_2of3, n_LOF, n_GOF)])

# HNF4A / THRB specifically (named in question)
for (tf in c("HNF4A", "THRB", "RORA", "NR1H4", "CEBPA")) {
  row <- per_tf[toupper(tf_name) == tf]
  if (nrow(row) > 0) {
    cat(sprintf("  %s: total=%d  3of3=%d  2of3=%d  LOF=%d  GOF=%d\n",
                tf, row$n_total, row$n_concord_3of3, row$n_concord_2of3,
                row$n_LOF, row$n_GOF))
  } else {
    cat(sprintf("  %s: no motif disruptions in input set\n", tf))
  }
}

cat("\nTop 20 genes by 3-of-3 count:\n")
print(head(per_gene[n_concord_3of3 > 0], 20))

cat("============================================================\n")
cat("Done. Outputs:\n")
cat("  ", file.path(OUT_DIR, "allele_concordance.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "allele_concordance_summary.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "allele_concordance_per_gene.csv"), "\n", sep = "")
cat("============================================================\n")
