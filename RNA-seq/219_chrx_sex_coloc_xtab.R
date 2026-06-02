#!/usr/bin/env Rscript
# 219_chrx_sex_coloc_xtab.R — Cross-tab A9's chrX classification × chrX COLOC.
#
# B6 — chrX × MASLD COLOC pipeline (Team B).
# Inputs:
#   1) A9 chrX classification: RNA-seq/results/audit_sensitivity/sex_xci/chrx_stratified_classification.csv
#   2) B6 chrX ABF COLOC: GWAS/finemapping/results/chrx_coloc_abf/<gwas>/chrx_coloc_abf.csv
#
# Outputs:
#   RNA-seq/results/stratified_causal/chrx_sex_coloc.csv
#   RNA-seq/results/stratified_causal/chrx_sex_coloc_fisher.csv
#   RNA-seq/results/stratified_causal/chrx_x_inactive_male_biased_coloc.csv
#
# Tests:
#   * Fisher (chr_category × sex_class_revised × chrx_pp4>0.5) 3-way
#   * Focus on X_escape vs X_inactive enrichment for chrX-COLOC hits
#   * Sanity check on A9's X_inactive × Male_biased OR=2.17 — are these COLOC hits?

suppressPackageStartupMessages({
  library(data.table)
})

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

A9_FILE <- file.path(BASE_DIR,
  "RNA-seq/results/audit_sensitivity/sex_xci/chrx_stratified_classification.csv")
COLOC_ROOT <- file.path(BASE_DIR, "GWAS/finemapping/results/chrx_coloc_abf")
OUT_DIR <- file.path(BASE_DIR, "RNA-seq/results/stratified_causal")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PP4_THRESH <- 0.5  # canonical PP.H4 threshold from CLAUDE.md

cat("============================================================\n")
cat("B6 chrX × sex × COLOC cross-tabulation\n")
cat("============================================================\n")

# ---- Load A9 classification ------------------------------------------------
cat("\n[1/5] Loading A9 chrX classification...\n")
a9 <- fread(A9_FILE)
cat("  Rows:", nrow(a9), " cols:", ncol(a9), "\n")
cat("  chr_category distribution:\n")
print(table(a9$chr_category, useNA = "ifany"))

# Keep only what we need (the A9 file is wider than needed)
a9_use <- a9[, .(
  ensembl_base, gene_name, chromosome, chr_category,
  xci_status_oliva, xci_status_tukiainen,
  sex_class_canonical, sex_class_revised,
  logFC_M, padj_M, logFC_F, padj_F,
  interaction_logFC, interaction_padj,
  suspicious_flag
)]

# ---- Discover available COLOC outputs --------------------------------------
cat("\n[2/5] Discovering chrX COLOC outputs...\n")
coloc_dirs <- list.dirs(COLOC_ROOT, recursive = FALSE)
cat("  Found", length(coloc_dirs), "COLOC dirs:\n")
for (d in coloc_dirs) cat("   -", basename(d), "\n")

if (length(coloc_dirs) == 0) {
  stop("No chrX COLOC outputs found under ", COLOC_ROOT)
}

# Load each + tag with stratum
gwas_strata <- list()
for (d in coloc_dirs) {
  fp <- file.path(d, "chrx_coloc_abf.csv")
  if (file.exists(fp)) {
    tbl <- tryCatch(fread(fp), error = function(e) NULL)
    if (!is.null(tbl) && nrow(tbl) > 0) {
      strat <- basename(d)
      tbl[, stratum := strat]
      gwas_strata[[strat]] <- tbl
      cat("   loaded", strat, ":", nrow(tbl), "genes\n")
    }
  }
}

if (length(gwas_strata) == 0) {
  stop("All chrX COLOC outputs are empty.")
}

coloc_long <- rbindlist(gwas_strata, fill = TRUE)
cat("  Total COLOC (gene × stratum) rows:", nrow(coloc_long), "\n")

# ---- Reshape: per-gene best PP4 per stratum --------------------------------
cat("\n[3/5] Aggregating to per-gene best PP4 per stratum...\n")

# For each gene × stratum take MAX PP4 (in case multiple windows / runs)
coloc_long[, ensembl_base := sub("\\..*$", "", ensembl)]

agg <- coloc_long[, .(
  pp4_abf      = max(PP.H4.abf,    na.rm = TRUE),
  pp4_h3       = max(PP.H3.abf,    na.rm = TRUE),
  n_snps_used  = max(n_snps,       na.rm = TRUE),
  top_snp      = top_snp[which.max(PP.H4.abf)]
), by = .(ensembl_base, gene, stratum)]

# Replace -Inf (from max() on all-NA) with NA
agg[is.infinite(pp4_abf), pp4_abf := NA_real_]
agg[is.infinite(pp4_h3),  pp4_h3  := NA_real_]

# Wide: one column per stratum × metric (only pp4_abf needed for headline)
wide <- dcast(agg, ensembl_base + gene ~ stratum,
              value.var = c("pp4_abf"),
              fill = NA_real_)

# Rename: keep convention chrx_pp4_abf_<stratum>
strata_seen <- unique(agg$stratum)
for (s in strata_seen) {
  if (s %in% names(wide)) setnames(wide, s, paste0("chrx_pp4_abf_", tolower(s)))
}
cat("  Wide table dims:", nrow(wide), "x", ncol(wide), "\n")

# ---- Join with A9 ----------------------------------------------------------
cat("\n[4/5] Joining COLOC with A9 chrX classification...\n")
joined <- merge(a9_use, wide, by = "ensembl_base", all.x = TRUE,
                suffixes = c("", ".coloc"))
cat("  Joined rows:", nrow(joined), "\n")

# Compute "any" PP4 max across strata
pp4_cols <- grep("^chrx_pp4_abf_", names(joined), value = TRUE)
if (length(pp4_cols) > 0) {
  joined[, chrx_pp4_abf_any_max := do.call(pmax, c(.SD, list(na.rm = TRUE))),
         .SDcols = pp4_cols]
  joined[is.infinite(chrx_pp4_abf_any_max), chrx_pp4_abf_any_max := NA_real_]
}

# Hit flag (any stratum PP4 > threshold)
joined[, chrx_coloc_hit := !is.na(chrx_pp4_abf_any_max) &
                            chrx_pp4_abf_any_max > PP4_THRESH]

# Per-chrX-only summary
chrx_only <- joined[chromosome == "chrX"]
cat("  chrX rows:", nrow(chrx_only), "\n")
cat("  chrX genes tested by COLOC:",
    sum(!is.na(chrx_only$chrx_pp4_abf_any_max)), "\n")
cat("  chrX COLOC hits (PP4>", PP4_THRESH, "):",
    sum(chrx_only$chrx_coloc_hit, na.rm = TRUE), "\n")

# ---- Write joined output ---------------------------------------------------
out_main <- file.path(OUT_DIR, "chrx_sex_coloc.csv")
fwrite(joined, out_main)
cat("  Wrote:", out_main, "\n")

# ---- Fisher tests ---------------------------------------------------------
cat("\n[5/5] Fisher tests: chr_category × sex_class_revised × chrX-COLOC hit\n")

# Test (a): chr_category × COLOC hit (chrX strata only)
xchr_only <- joined[chromosome == "chrX" & !is.na(chrx_pp4_abf_any_max)]
cat("\nTest A: chr_category × chrx_coloc_hit (within chrX)\n")
tabA <- table(xchr_only$chr_category, xchr_only$chrx_coloc_hit)
print(tabA)
fisher_results <- list()
if (nrow(tabA) >= 2 && ncol(tabA) >= 2) {
  for (cat_level in rownames(tabA)) {
    if (cat_level == "X_inactive") next  # reference
    a <- tabA[cat_level, "TRUE"]
    b <- tabA[cat_level, "FALSE"]
    c <- tabA["X_inactive", "TRUE"]
    d <- tabA["X_inactive", "FALSE"]
    if (is.na(a) || is.na(b) || is.na(c) || is.na(d)) next
    fr <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    fisher_results[[length(fisher_results) + 1]] <- data.table(
      test       = "chr_category_vs_X_inactive",
      level      = cat_level,
      n_hit      = a, n_no_hit = b,
      n_ref_hit  = c, n_ref_no_hit = d,
      OR         = unname(fr$estimate),
      pval       = fr$p.value,
      ci_lower   = fr$conf.int[1],
      ci_upper   = fr$conf.int[2]
    )
  }
}

# Test (b): sex_class_revised × COLOC hit (chrX only)
cat("\nTest B: sex_class_revised × chrx_coloc_hit (within chrX)\n")
tabB <- table(xchr_only$sex_class_revised, xchr_only$chrx_coloc_hit)
print(tabB)
if (nrow(tabB) >= 2 && ncol(tabB) >= 2) {
  ref_lvl <- "Concordant"
  if (!ref_lvl %in% rownames(tabB)) ref_lvl <- rownames(tabB)[which.max(rowSums(tabB))]
  for (cls_level in rownames(tabB)) {
    if (cls_level == ref_lvl) next
    a <- tabB[cls_level, "TRUE"]
    b <- tabB[cls_level, "FALSE"]
    c <- tabB[ref_lvl,  "TRUE"]
    d <- tabB[ref_lvl,  "FALSE"]
    if (is.na(a) || is.na(b) || is.na(c) || is.na(d)) next
    if ((a + b) < 3) next  # skip very small classes
    fr <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    fisher_results[[length(fisher_results) + 1]] <- data.table(
      test       = paste0("sex_class_vs_", ref_lvl),
      level      = cls_level,
      n_hit      = a, n_no_hit = b,
      n_ref_hit  = c, n_ref_no_hit = d,
      OR         = unname(fr$estimate),
      pval       = fr$p.value,
      ci_lower   = fr$conf.int[1],
      ci_upper   = fr$conf.int[2]
    )
  }
}

# Test (c): X_escape vs X_inactive — within Female_biased only (the headline)
fb_only <- xchr_only[sex_class_revised == "Female_biased" |
                       grepl("escape_expected", sex_class_revised) |
                       grepl("Female_biased", sex_class_revised)]
cat("\nTest C: X_escape vs X_inactive among Female-biased calls\n")
if (nrow(fb_only) > 0) {
  tabC <- table(fb_only$chr_category, fb_only$chrx_coloc_hit)
  print(tabC)
  if ("X_escape" %in% rownames(tabC) && "X_inactive" %in% rownames(tabC) &&
      ncol(tabC) >= 2) {
    a <- tabC["X_escape",  "TRUE"]; b <- tabC["X_escape",  "FALSE"]
    c <- tabC["X_inactive","TRUE"]; d <- tabC["X_inactive","FALSE"]
    if (!is.na(a) && !is.na(b) && !is.na(c) && !is.na(d) && (a+b) > 0) {
      fr <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
      fisher_results[[length(fisher_results) + 1]] <- data.table(
        test = "Xescape_vs_Xinactive_within_FemaleBiased",
        level = "X_escape",
        n_hit = a, n_no_hit = b, n_ref_hit = c, n_ref_no_hit = d,
        OR = unname(fr$estimate), pval = fr$p.value,
        ci_lower = fr$conf.int[1], ci_upper = fr$conf.int[2]
      )
    }
  }
}

if (length(fisher_results) > 0) {
  fisher_dt <- rbindlist(fisher_results, fill = TRUE)
  fisher_dt[, pval_BH := p.adjust(pval, method = "BH")]
  fwrite(fisher_dt, file.path(OUT_DIR, "chrx_sex_coloc_fisher.csv"))
  cat("\n  Wrote Fisher:", file.path(OUT_DIR, "chrx_sex_coloc_fisher.csv"), "\n")
  print(fisher_dt)
} else {
  cat("\n  No Fisher tests produced (insufficient cell counts).\n")
}

# ---- Sanity check: X_inactive × Male_biased (A9 OR=2.17 enrichment) -------
cat("\n--- Sanity check: X_inactive × Male_biased COLOC status ---\n")
xim <- joined[chromosome == "chrX" &
                chr_category == "X_inactive" &
                sex_class_revised == "Male_biased"]
cat("  X_inactive × Male_biased gene count:", nrow(xim), "\n")
if (nrow(xim) > 0) {
  cat("  Tested by COLOC:", sum(!is.na(xim$chrx_pp4_abf_any_max)), "\n")
  cat("  COLOC hits (PP4>", PP4_THRESH, "):",
      sum(xim$chrx_coloc_hit, na.rm = TRUE), "\n")
  if (sum(xim$chrx_coloc_hit, na.rm = TRUE) > 0) {
    cat("  Hit gene names:\n")
    print(xim[chrx_coloc_hit == TRUE,
              .(gene_name, chr_category, sex_class_revised,
                chrx_pp4_abf_any_max)])
  }
}

out_sanity <- file.path(OUT_DIR, "chrx_x_inactive_male_biased_coloc.csv")
fwrite(xim, out_sanity)
cat("  Wrote:", out_sanity, "\n")

cat("\n============================================================\n")
cat("219: complete\n")
cat("============================================================\n")
