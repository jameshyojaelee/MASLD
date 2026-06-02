#!/usr/bin/env Rscript
# 56i_currin_caqtl_concordance.R  --  B4 pivot (ATAC improvement plan)
#
# Cross-validate our scATAC motif-disrupted variants against an INDEPENDENT
# bulk-tissue liver caQTL panel from Currin & Mohlke 2025
# (138 donors; Genome Research; doi:10.1101/gr.279741.124; Zenodo 15025748).
#
# Inputs (read-only):
#   GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#   GWAS/finemapping/results/gwas_atac/allele_concordance.csv
#   data/external/currin_2025_caqtl/caQTL_hg38_chr{1..22}_*.txt.gz
#   data/external/currin_2025_caqtl/hg19ToHg38.over.chain.gz
#
# Outputs (write-only):
#   GWAS/finemapping/results/gwas_atac/currin_caqtl_concordance.csv
#   GWAS/finemapping/results/gwas_atac/currin_caqtl_concordance_summary.csv
#   GWAS/finemapping/results/gwas_atac/currin_caqtl_per_tf.csv
#
# Logic:
#   * Our SNP_id is hg19 (e.g. "1:155108287:G:A"); Currin variants are hg38
#     (e.g. "chr22:15539799:A:C"). LiftOver our variants to hg38, then match
#     on (chr, pos38, REF, ALT) -- allowing ref/alt swap.
#   * caQTL beta orientation: FastQTL beta is "ALT-allele effect on peak count"
#     after inverse-normal transform. We compare sign(beta_caQTL) to
#     sign(alleleDiff) (motif disruption: + = ALT increases binding).
#   * For each motif-disrupted variant, keep the most-significant caQTL hit
#     in cis (1Mb file).
#   * 4-way concordance: motif + eQTL + DA + caQTL. We re-use the existing
#     sign_motif / sign_da / sign_eqtl from allele_concordance.csv and add
#     sign_caqtl + 4-way pairwise concordance counts.
#
# Conventions:
#   * Currin caQTL p-value threshold: nominal p<1e-5 ("suggestive") and the
#     formal beta-FDR/beta-permutation threshold not in nominal file -- we
#     therefore report counts at p<1e-5 (suggestive), p<1e-3 (lenient).
#   * Directionality: "agree" means sign(beta_caQTL) == sign(alleleDiff).

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR    <- file.path(BASE_DIR, "GWAS/finemapping")
OUT_DIR   <- file.path(FM_DIR, "results/gwas_atac")
CURRIN_DIR <- file.path(BASE_DIR, "data/external/currin_2025_caqtl")
CHAIN_FILE <- file.path(CURRIN_DIR, "hg19ToHg38.over.chain.gz")

cat("============================================================\n")
cat("56i_currin_caqtl_concordance.R  --  B4 4-way concordance\n")
cat("============================================================\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load inputs
# ─────────────────────────────────────────────────────────────────────────────
md_file  <- file.path(OUT_DIR, "motif_disruption_scores.csv")
ac_file  <- file.path(OUT_DIR, "allele_concordance.csv")
stopifnot(file.exists(md_file), file.exists(ac_file), file.exists(CHAIN_FILE))

md <- fread(md_file)
ac <- fread(ac_file)
cat("Loaded motif_disruption_scores: ", nrow(md), " rows\n", sep = "")
cat("Loaded allele_concordance:      ", nrow(ac), " rows\n", sep = "")

# Use motif rows as the canonical set (1 row per SNP/TF/peak tuple).
# Unique variants for caQTL lookup:
md[, c("snp_chr_hg19", "snp_pos_hg19", "snp_ref", "snp_alt") :=
     tstrsplit(SNP_id, ":", fixed = TRUE)]
md[, snp_pos_hg19 := as.integer(snp_pos_hg19)]
md[, snp_chr := paste0("chr", snp_chr_hg19)]

variants_hg19 <- unique(md[, .(SNP_id, snp_chr, snp_pos_hg19, snp_ref, snp_alt)])
cat("Unique variants to lift over: ", nrow(variants_hg19), "\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# 2. LiftOver SNP_id hg19 -> hg38
# ─────────────────────────────────────────────────────────────────────────────
# rtracklayer::liftOver accepts a gzipped chain via gunzip-on-the-fly
chain_tmp <- tempfile(fileext = ".chain")
R.utils::gunzip(CHAIN_FILE, destname = chain_tmp, remove = FALSE, overwrite = TRUE)
chain <- import.chain(chain_tmp)
unlink(chain_tmp)

gr_hg19 <- GRanges(seqnames = variants_hg19$snp_chr,
                   ranges   = IRanges(start = variants_hg19$snp_pos_hg19, width = 1L),
                   SNP_id   = variants_hg19$SNP_id,
                   ref      = variants_hg19$snp_ref,
                   alt      = variants_hg19$snp_alt)

lifted <- liftOver(gr_hg19, chain)
n_lifted <- sum(elementNROWS(lifted) >= 1L)
cat("LiftOver hg19->hg38 success: ", n_lifted, " / ", nrow(variants_hg19),
    "  (", round(100 * n_lifted / nrow(variants_hg19), 1), "%)\n", sep = "")

# Take first mapping if multi-mapping
lifted_gr <- unlist(lifted[elementNROWS(lifted) >= 1L])
# Re-attach metadata: when a single hg19 position maps to multiple hg38 hits,
# elementNROWS > 1; we choose the first.
ll_first <- vapply(seq_along(lifted), function(i) {
  if (length(lifted[[i]]) == 0L) NA_integer_ else as.integer(start(lifted[[i]])[1])
}, integer(1))
ll_chr <- vapply(seq_along(lifted), function(i) {
  if (length(lifted[[i]]) == 0L) NA_character_ else as.character(seqnames(lifted[[i]])[1])
}, character(1))
variants_hg19[, snp_pos_hg38 := ll_first]
variants_hg19[, snp_chr_hg38 := ll_chr]

# Build hg38 variant_id matches in Currin format: "chr22:15539799:A:C"
# We will match in BOTH orientations (REF:ALT and ALT:REF) since Currin's
# REF/ALT assignment may differ from our liftOver source.
variants_hg19[, var_id_hg38_AB := paste0(snp_chr_hg38, ":", snp_pos_hg38, ":", snp_ref, ":", snp_alt)]
variants_hg19[, var_id_hg38_BA := paste0(snp_chr_hg38, ":", snp_pos_hg38, ":", snp_alt, ":", snp_ref)]
cat("Built hg38 variant_id keys (2 orientations per variant).\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Stream-load Currin caQTL files per chromosome and join
# ─────────────────────────────────────────────────────────────────────────────
currin_files <- list.files(CURRIN_DIR, pattern = "^caQTL_hg38_chr.*\\.txt\\.gz$",
                            full.names = TRUE)
stopifnot(length(currin_files) >= 22L)
cat("Currin caQTL files: ", length(currin_files), "\n", sep = "")

# For each chromosome, subset our variants and stream Currin for matching variant rows.
# Currin file size ~4M rows per chr; we filter using data.table on the variant col.
chrs <- sort(unique(variants_hg19$snp_chr_hg38))
chrs <- chrs[!is.na(chrs)]
cat("Chromosomes to scan: ", paste(chrs, collapse = ", "), "\n", sep = "")

hits_list <- vector("list", length(chrs))
for (i in seq_along(chrs)) {
  ch <- chrs[i]
  # Find corresponding Currin file
  pat <- paste0("caQTL_hg38_", ch, "_")
  cf <- currin_files[grepl(pat, basename(currin_files), fixed = TRUE)]
  if (length(cf) != 1L) {
    cat("  [", ch, "] no Currin file found; skipping.\n", sep = "")
    next
  }
  vchr <- variants_hg19[snp_chr_hg38 == ch]
  if (nrow(vchr) == 0L) next

  # Build a set of variant IDs to look up (both orientations)
  lookup_set <- unique(c(vchr$var_id_hg38_AB, vchr$var_id_hg38_BA))
  lookup_set <- lookup_set[!is.na(lookup_set)]

  # Stream Currin file using fread; filter on `variant %in% lookup_set`.
  # Currin schema: variant, peak, distance_from_peakCenter, pvalue, beta, varbeta, imputation_R2, MAF
  cdt <- fread(cf, sep = "\t", header = TRUE, showProgress = FALSE)
  if (nrow(cdt) == 0L) next
  cdt_hit <- cdt[variant %in% lookup_set]
  cat("  [", ch, "] Currin rows: ", nrow(cdt), "  matching our variants: ",
      nrow(cdt_hit), "  (unique variants matched: ",
      length(unique(cdt_hit$variant)), " / ", nrow(vchr), ")\n", sep = "")
  if (nrow(cdt_hit) == 0L) next

  hits_list[[i]] <- cdt_hit
  rm(cdt, cdt_hit); gc(verbose = FALSE)
}

hits <- rbindlist(hits_list, fill = TRUE)
cat("\nTotal Currin caQTL hits across all chromosomes: ", nrow(hits), "\n", sep = "")
cat("Unique variants with >=1 Currin caQTL test: ",
    length(unique(hits$variant)), "\n", sep = "")

if (nrow(hits) == 0L) {
  cat("WARNING: no Currin matches found. Writing empty output and exiting.\n")
  empty <- data.table()
  fwrite(empty, file.path(OUT_DIR, "currin_caqtl_concordance.csv"))
  fwrite(empty, file.path(OUT_DIR, "currin_caqtl_concordance_summary.csv"))
  fwrite(empty, file.path(OUT_DIR, "currin_caqtl_per_tf.csv"))
  quit(status = 0)
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. Map Currin hits back to our SNP_id (account for ref/alt orientation)
# ─────────────────────────────────────────────────────────────────────────────
# For each unique variant hit, keep the most-significant caQTL test (min p).
hits_top <- hits[order(variant, pvalue), .SD[1], by = variant]
cat("Most-significant caQTL per variant: ", nrow(hits_top), "\n", sep = "")

# Join to variants_hg19 via either orientation; track orientation flag.
v_AB <- variants_hg19[var_id_hg38_AB %in% hits_top$variant,
                       .(SNP_id, var_id_match = var_id_hg38_AB, orient = "AB")]
v_BA <- variants_hg19[var_id_hg38_BA %in% hits_top$variant &
                       !(SNP_id %in% v_AB$SNP_id),
                       .(SNP_id, var_id_match = var_id_hg38_BA, orient = "BA")]
v_lookup <- rbind(v_AB, v_BA)
cat("Variants matched (AB orientation): ", nrow(v_AB), "\n", sep = "")
cat("Variants matched (BA orientation): ", nrow(v_BA), "\n", sep = "")

setnames(hits_top, "variant", "var_id_match")
v_lookup <- merge(v_lookup, hits_top, by = "var_id_match", all.x = TRUE)
setnames(v_lookup,
         c("pvalue", "beta", "varbeta", "MAF", "peak", "distance_from_peakCenter"),
         c("currin_caqtl_pval", "currin_caqtl_beta", "currin_caqtl_varbeta",
           "currin_caqtl_maf", "currin_caqtl_peak", "currin_caqtl_dist"))
# Flip beta sign when matched in BA orientation (our ALT is Currin's REF)
v_lookup[orient == "BA", currin_caqtl_beta := -currin_caqtl_beta]
v_lookup[, currin_caqtl_se := sqrt(currin_caqtl_varbeta)]

cat("\nCurrin caQTL coverage:\n")
cat("  Variants with caQTL test:                       ", nrow(v_lookup), "\n", sep = "")
cat("  Variants with p < 1e-3:                         ",
    sum(v_lookup$currin_caqtl_pval < 1e-3, na.rm = TRUE), "\n", sep = "")
cat("  Variants with p < 1e-5 (suggestive):            ",
    sum(v_lookup$currin_caqtl_pval < 1e-5, na.rm = TRUE), "\n", sep = "")
cat("  Variants with p < 5e-8 (genome-wide):           ",
    sum(v_lookup$currin_caqtl_pval < 5e-8, na.rm = TRUE), "\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Join to allele_concordance and compute 4-way concordance
# ─────────────────────────────────────────────────────────────────────────────
# allele_concordance.csv is (SNP_id, tf_name, ...) -- can have multiple rows per SNP.
# We add Currin caQTL columns at the (SNP_id) level (per-variant assay).
ac_join <- merge(ac, v_lookup[, .(SNP_id, currin_caqtl_pval, currin_caqtl_beta,
                                   currin_caqtl_se, currin_caqtl_maf,
                                   currin_caqtl_peak, currin_caqtl_dist,
                                   currin_orient = orient)],
                  by = "SNP_id", all.x = TRUE)

# Compute sign_caqtl (+ = ALT increases peak accessibility, - = decreases).
# FIX (review F117): only SIGNIFICANT (p<1e-5) caQTLs get a sign; a p=0.48 noise caQTL must not
# count toward 4-way concordance. Unfiltered, genome-wide motif<->caQTL agreement was 46.3% (below
# the 50% chance null). NA sign propagates to concord_*_caqtl (computed below) and concord_4of4.
ac_join[, sign_caqtl := fifelse(!is.na(currin_caqtl_pval) & currin_caqtl_pval < 1e-5,
                                sign(currin_caqtl_beta), NA_real_)]

# caQTL "called" thresholds
ac_join[, currin_caqtl_called_p1e5 := !is.na(currin_caqtl_pval) & currin_caqtl_pval < 1e-5]
ac_join[, currin_caqtl_called_p1e3 := !is.na(currin_caqtl_pval) & currin_caqtl_pval < 1e-3]
ac_join[, currin_caqtl_tested      := !is.na(currin_caqtl_pval)]

# 4-way pairwise concordance (motif, eqtl, da, caqtl)
ac_join[, concord_motif_caqtl := sign_motif == sign_caqtl]
ac_join[, concord_da_caqtl    := sign_da    == sign_caqtl]
ac_join[, concord_eqtl_caqtl  := sign_eqtl_effective == sign_caqtl]

# 4-of-4 concordance: motif, da, eqtl (already concord_2of3==TRUE or
# concord_3of3==TRUE), AND caqtl agrees with motif (the primary regulatory axis).
ac_join[, concord_4of4 := concord_3of3 == TRUE &
                          !is.na(sign_caqtl) &
                          concord_motif_caqtl == TRUE]
ac_join[, concord_2of3_plus_caqtl := concord_2of3 == TRUE &
                                       !is.na(sign_caqtl) &
                                       concord_motif_caqtl == TRUE]

cat("\n4-way concordance counts:\n")
cat("  Total rows:                                     ", nrow(ac_join), "\n", sep = "")
cat("  Rows with Currin caQTL tested:                  ",
    sum(ac_join$currin_caqtl_tested), "\n", sep = "")
cat("  Rows where motif+caQTL signs agree:             ",
    sum(ac_join$concord_motif_caqtl, na.rm = TRUE), "\n", sep = "")
cat("  Rows where DA+caQTL signs agree:                ",
    sum(ac_join$concord_da_caqtl, na.rm = TRUE), "\n", sep = "")
cat("  Rows where eQTL+caQTL signs agree:              ",
    sum(ac_join$concord_eqtl_caqtl, na.rm = TRUE), "\n", sep = "")
cat("  3-of-3 (existing) + Currin caQTL motif-agree:   ",
    sum(ac_join$concord_4of4, na.rm = TRUE), "\n", sep = "")
cat("  2-of-3 (existing) + Currin caQTL motif-agree:   ",
    sum(ac_join$concord_2of3_plus_caqtl, na.rm = TRUE), "\n", sep = "")

# FIX (review F117): honest motif<->caQTL agreement among SIGNIFICANT caQTLs vs a binomial(0.5) null.
.sig <- ac_join[currin_caqtl_called_p1e5 == TRUE & !is.na(concord_motif_caqtl)]
if (nrow(.sig) > 0) {
  .n <- nrow(.sig); .k <- sum(.sig$concord_motif_caqtl)
  .bt <- binom.test(.k, .n, p = 0.5, alternative = "greater")
  cat(sprintf("  Motif<->caQTL agreement (SIGNIFICANT caQTLs p<1e-5): %d/%d = %.1f%%  binomial-vs-0.5 p=%.3g\n",
              .k, .n, 100 * .k / .n, .bt$p.value))
  cat(sprintf("  Honest 4-of-4 count (significant caQTLs only): %d\n", sum(.sig$concord_4of4, na.rm = TRUE)))
} else {
  cat("  No significant (p<1e-5) caQTLs with a motif sign — no honest 4-way agreement to report.\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. Subset headline columns + write
# ─────────────────────────────────────────────────────────────────────────────
keep_cols <- c("SNP_id", "tf_name", "linked_gene_resolved", "is_activator",
               "motif_in_disease_regulon", "variant_cell_type",
               "alleleDiff", "eqtl_beta_lead", "eqtl_pval_lead",
               "da_logFC_final", "da_padj_final", "coloc_best_pp4",
               "sign_motif", "sign_da", "sign_eqtl", "sign_eqtl_effective",
               "concord_motif_eqtl", "concord_motif_da", "concord_da_eqtl",
               "n_signs_present", "n_pairs_concord", "concord_2of3", "concord_3of3",
               "currin_caqtl_pval", "currin_caqtl_beta", "currin_caqtl_se",
               "currin_caqtl_maf", "currin_caqtl_peak", "currin_caqtl_dist",
               "currin_orient", "sign_caqtl",
               "currin_caqtl_tested", "currin_caqtl_called_p1e3",
               "currin_caqtl_called_p1e5",
               "concord_motif_caqtl", "concord_da_caqtl", "concord_eqtl_caqtl",
               "concord_4of4", "concord_2of3_plus_caqtl")
keep_cols <- intersect(keep_cols, names(ac_join))
out_main <- ac_join[, ..keep_cols]

out_path <- file.path(OUT_DIR, "currin_caqtl_concordance.csv")
fwrite(out_main, out_path)
cat("\nWrote: ", out_path, "  (", nrow(out_main), " rows x ",
    ncol(out_main), " cols)\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# 7. Summary stats (overall + per-2of3-subset)
# ─────────────────────────────────────────────────────────────────────────────
summarise_block <- function(dt, label) {
  data.table(
    subset                          = label,
    n_rows                          = nrow(dt),
    n_caqtl_tested                  = sum(dt$currin_caqtl_tested),
    pct_caqtl_tested                = round(100 * sum(dt$currin_caqtl_tested) / nrow(dt), 1),
    n_caqtl_p1e3                    = sum(dt$currin_caqtl_called_p1e3, na.rm = TRUE),
    n_caqtl_p1e5                    = sum(dt$currin_caqtl_called_p1e5, na.rm = TRUE),
    pct_motif_caqtl_agree           = round(100 * mean(dt$concord_motif_caqtl, na.rm = TRUE), 1),
    pct_da_caqtl_agree              = round(100 * mean(dt$concord_da_caqtl, na.rm = TRUE), 1),
    pct_eqtl_caqtl_agree            = round(100 * mean(dt$concord_eqtl_caqtl, na.rm = TRUE), 1)
  )
}

summ <- rbind(
  summarise_block(ac_join, "all_rows"),
  summarise_block(ac_join[concord_2of3 == TRUE], "concord_2of3"),
  summarise_block(ac_join[concord_3of3 == TRUE], "concord_3of3"),
  summarise_block(ac_join[motif_in_disease_regulon == TRUE], "disease_regulon"),
  summarise_block(ac_join[motif_in_disease_regulon == TRUE & concord_2of3 == TRUE],
                  "disease_regulon_concord_2of3")
)
summ_path <- file.path(OUT_DIR, "currin_caqtl_concordance_summary.csv")
fwrite(summ, summ_path)
cat("Wrote: ", summ_path, "\n", sep = "")
print(summ)

# ─────────────────────────────────────────────────────────────────────────────
# 8. Per-TF support table -- which disease regulon TFs survive 4-way?
# ─────────────────────────────────────────────────────────────────────────────
per_tf <- ac_join[, .(
  n_variants                  = .N,
  n_caqtl_tested              = sum(currin_caqtl_tested),
  n_caqtl_p1e3                = sum(currin_caqtl_called_p1e3, na.rm = TRUE),
  n_caqtl_p1e5                = sum(currin_caqtl_called_p1e5, na.rm = TRUE),
  n_motif_caqtl_agree         = sum(concord_motif_caqtl, na.rm = TRUE),
  pct_motif_caqtl_agree       = round(100 * mean(concord_motif_caqtl, na.rm = TRUE), 1),
  n_4of4                      = sum(concord_4of4, na.rm = TRUE),
  n_2of3_plus_caqtl           = sum(concord_2of3_plus_caqtl, na.rm = TRUE),
  any_disease_regulon         = any(motif_in_disease_regulon == TRUE, na.rm = TRUE)
), by = tf_name]
per_tf <- per_tf[order(-n_2of3_plus_caqtl, -pct_motif_caqtl_agree)]
per_tf_path <- file.path(OUT_DIR, "currin_caqtl_per_tf.csv")
fwrite(per_tf, per_tf_path)
cat("Wrote: ", per_tf_path, "  (", nrow(per_tf), " TFs)\n", sep = "")

# Top disease-regulon TFs surviving 4-way
cat("\nTop disease-regulon TFs by 4-way support (motif+eQTL+DA+caQTL):\n")
top_tf <- per_tf[any_disease_regulon == TRUE & n_motif_caqtl_agree > 0]
top_tf <- top_tf[order(-n_2of3_plus_caqtl, -pct_motif_caqtl_agree)][1:min(20, .N)]
print(top_tf)

cat("\n[", as.character(Sys.time()), "] Done.\n", sep = "")
