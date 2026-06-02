#!/usr/bin/env Rscript
# ============================================================================
# 343k_gated_fstage_audit.R
#
# Audit the gated F-stage axis lmm output (Andrews + Wang) to verify whether
# the gating from ungated 8082 -> gated 4887 is producing trustworthy LR-pair
# signals.
#
# Six audit dimensions:
#   1. Top-LR comparison (gated vs ungated already-replaced; we just inspect
#      gated, and re-fit Andrews-only as a within-gated robustness check).
#   2. Cell-type pair breakdown for the gated table.
#   3. Bulk concordance among gated F-stage top hits.
#   4. LOO replication rate among gated top hits.
#   5. Biological sanity of top-20 hits (textbook fibrosis ligands).
#   6. Andrews-only refit (drop Wang) to test cohort imbalance.
#
# Outputs:
#   gated_fstage_audit.tsv      (per-LR top-50 with andrews-only stats)
#   gated_fstage_audit_summary.txt
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(parallel)
})

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
options(mc.cores = N_CORES)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")

META_EXT <- file.path(OUT_DIR, "donor_metadata_extended.tsv")
MERGED_TSV <- file.path(OUT_DIR, "all_donor_lr_scores.tsv.gz")
FSTAGE_TSV <- file.path(OUT_DIR, "stage_lr_lmm_fstage.tsv")
COARSE_TSV <- file.path(OUT_DIR, "stage_lr_lmm_coarse.tsv")
BULK_TSV   <- file.path(OUT_DIR, "lr_bulk_concordance.tsv")
LOO_TSV    <- file.path(OUT_DIR, "loo_replication_rate_per_lr.tsv")

AUDIT_TSV  <- file.path(OUT_DIR, "gated_fstage_audit.tsv")
SUMMARY_TXT <- file.path(OUT_DIR, "gated_fstage_audit_summary.txt")

# Buffer for the human-readable summary
buf <- c()
addln <- function(...) {
  s <- paste0(...)
  buf <<- c(buf, s)
  cat(s, "\n", sep = "")
}

addln("====================================================================")
addln("Gated F-stage axis audit  (343k)")
addln(sprintf("Generated: %s", Sys.time()))
addln("====================================================================")

# ----- inputs ---------------------------------------------------------------
meta <- fread(META_EXT)
fstage <- fread(FSTAGE_TSV)
addln(sprintf("[input] gated stage_lr_lmm_fstage.tsv rows = %d", nrow(fstage)))
addln(sprintf("[input] donor metadata rows                = %d", nrow(meta)))

n_inf <- sum(!is.na(meta$F_stage_inferred))
n_doc <- sum(!is.na(meta$F_stage_documented))
addln(sprintf("[input] F_stage_inferred non-NA            = %d", n_inf))
addln(sprintf("[input] F_stage_documented non-NA          = %d", n_doc))
addln("[input] F_stage_source breakdown:")
src_tbl <- meta[, .N, by = F_stage_source][order(-N)]
for (i in seq_len(nrow(src_tbl))) {
  addln(sprintf("        %-30s %d", as.character(src_tbl$F_stage_source[i]), src_tbl$N[i]))
}
addln("[input] dataset breakdown (donors with F_stage_inferred):")
ds_tbl <- meta[!is.na(F_stage_inferred), .N, by = dataset][order(-N)]
for (i in seq_len(nrow(ds_tbl))) {
  addln(sprintf("        %-15s %d", ds_tbl$dataset[i], ds_tbl$N[i]))
}

# ----- score helper ---------------------------------------------------------
score_rank <- function(dt) {
  dt[, score := abs(Estimate) * (-log10(pmax(pval, 1e-300)))]
  setorder(dt, -score)
  dt[, rank := seq_len(.N)]
  dt
}
fstage <- score_rank(copy(fstage))

# ============================================================================
# DIMENSION 2 -- cell-type pair breakdown
# ============================================================================
addln("\n--- (2) cell-type pair breakdown of gated table ---")
ct_tbl <- fstage[, .N, by = ct_pair][order(-N)]
addln(sprintf("Total unique ct_pairs = %d", nrow(ct_tbl)))
addln("Top 20 ct_pairs by row count:")
for (i in seq_len(min(20, nrow(ct_tbl)))) {
  addln(sprintf("  %-50s %d", ct_tbl$ct_pair[i], ct_tbl$N[i]))
}

# 9 lineages mentioned by user
lineages_expected <- c("Hepatocytes","Cholangiocytes","Endothelial cells",
                       "Fibroblasts","Macrophages","T cells","B cells",
                       "Plasma cells","NK cells")
src_present <- unique(fstage$source)
tgt_present <- unique(fstage$target)
addln(sprintf("\nUnique sources represented = %d", length(src_present)))
addln(paste("  sources:", paste(sort(src_present), collapse = ", ")))
addln(sprintf("Unique targets represented = %d", length(tgt_present)))
addln(paste("  targets:", paste(sort(tgt_present), collapse = ", ")))

# ============================================================================
# DIMENSION 1 -- already gated; we just look at the gated top 50.
# ============================================================================
addln("\n--- (1) Gated top-50 (skipped vs documented-only as instructed) ---")
top50 <- fstage[1:50]
addln(sprintf("Gated top-50 |Estimate|*-log10(p) range = [%.3f, %.3f]",
              min(top50$score), max(top50$score)))
addln(sprintf("Gated top-50 unique ct_pairs = %d, unique LR pairs = %d",
              length(unique(top50$ct_pair)), length(unique(top50$lr_pair))))

# Reference: how many top-50 have nominal p<0.05 BH-within-ct or family-bonf<0.05?
addln(sprintf("  top-50 padj_within_ct < 0.05  : %d",
              sum(top50$padj_within_ct < 0.05, na.rm = TRUE)))
addln(sprintf("  top-50 family_bonferroni<0.05 : %d",
              sum(top50$family_bonferroni < 0.05, na.rm = TRUE)))

# ============================================================================
# DIMENSION 3 -- bulk concordance for gated top hits
# ============================================================================
addln("\n--- (3) bulk-concordance for gated F-stage top hits ---")
bulk <- fread(BULK_TSV)
# bulk has an axis col; look for F-stage axis there
addln(sprintf("Bulk table rows = %d, axis values:", nrow(bulk)))
ax_tbl <- bulk[, .N, by = axis][order(-N)]
for (i in seq_len(nrow(ax_tbl))) addln(sprintf("  %-25s %d", ax_tbl$axis[i], ax_tbl$N[i]))

fstage_bulk <- bulk[axis == "F_stage" | grepl("F_stage", axis)]
if (nrow(fstage_bulk) == 0) {
  # try term column instead
  fstage_bulk <- bulk[grepl("F_stage", term)]
}
addln(sprintf("F-stage rows in bulk table = %d", nrow(fstage_bulk)))

if (nrow(fstage_bulk) > 0) {
  key_top <- top50[, .(ct_pair, lr_pair)]
  fb_top <- merge(key_top, fstage_bulk, by = c("ct_pair","lr_pair"))
  addln(sprintf("Top-50 LR pairs matched into bulk table = %d", nrow(fb_top)))
  if (nrow(fb_top) > 0) {
    addln(sprintf("  ligand concordant   : %d / %d (%.1f%%)",
                  sum(fb_top$lig_concordant, na.rm=TRUE), nrow(fb_top),
                  100*mean(fb_top$lig_concordant, na.rm=TRUE)))
    addln(sprintf("  receptor concordant : %d / %d (%.1f%%)",
                  sum(fb_top$rec_concordant, na.rm=TRUE), nrow(fb_top),
                  100*mean(fb_top$rec_concordant, na.rm=TRUE)))
    addln(sprintf("  both concordant     : %d / %d (%.1f%%)",
                  sum(fb_top$both_concordant, na.rm=TRUE), nrow(fb_top),
                  100*mean(fb_top$both_concordant, na.rm=TRUE)))
  }
  # full F-stage axis concordance for reference
  addln(sprintf("Full F-stage axis bulk concordance (all %d rows):", nrow(fstage_bulk)))
  addln(sprintf("  ligand concordant   : %.1f%%", 100*mean(fstage_bulk$lig_concordant, na.rm=TRUE)))
  addln(sprintf("  receptor concordant : %.1f%%", 100*mean(fstage_bulk$rec_concordant, na.rm=TRUE)))
  addln(sprintf("  both concordant     : %.1f%%", 100*mean(fstage_bulk$both_concordant, na.rm=TRUE)))
}

# Coarse axis for comparison
coarse_bulk <- bulk[axis != "F_stage" & !grepl("F_stage", axis)]
if (nrow(coarse_bulk) > 0) {
  addln("\nCoarse-axis bulk concordance (for reference):")
  for (ax in unique(coarse_bulk$axis)) {
    sub <- coarse_bulk[axis == ax]
    addln(sprintf("  %-25s lig=%.1f%% rec=%.1f%% both=%.1f%% (n=%d)",
                  ax,
                  100*mean(sub$lig_concordant, na.rm=TRUE),
                  100*mean(sub$rec_concordant, na.rm=TRUE),
                  100*mean(sub$both_concordant, na.rm=TRUE),
                  nrow(sub)))
  }
}

# ============================================================================
# DIMENSION 4 -- LOO replication
# ============================================================================
addln("\n--- (4) LOO-dataset replication of gated top-50 ---")
loo <- fread(LOO_TSV)
addln(sprintf("LOO table rows = %d", nrow(loo)))
loo_top <- merge(top50[, .(ct_pair, lr_pair)], loo, by = c("ct_pair","lr_pair"))
addln(sprintf("Top-50 matched in LOO table = %d", nrow(loo_top)))
if (nrow(loo_top) > 0) {
  addln(sprintf("  mean replication rate            = %.3f", mean(loo_top$replication_rate, na.rm=TRUE)))
  addln(sprintf("  median replication rate          = %.3f", median(loo_top$replication_rate, na.rm=TRUE)))
  addln(sprintf("  >=0.5 replication                = %d (%.1f%%)",
                sum(loo_top$replication_rate >= 0.5, na.rm=TRUE),
                100*mean(loo_top$replication_rate >= 0.5, na.rm=TRUE)))
  addln(sprintf("  >=0.8 replication                = %d (%.1f%%)",
                sum(loo_top$replication_rate >= 0.8, na.rm=TRUE),
                100*mean(loo_top$replication_rate >= 0.8, na.rm=TRUE)))
}

# Compare to overall LOO distribution (across-all-LR baseline)
addln(sprintf("Overall LOO replication mean = %.3f, median = %.3f",
              mean(loo$replication_rate, na.rm=TRUE),
              median(loo$replication_rate, na.rm=TRUE)))

# ============================================================================
# DIMENSION 5 -- Biological sanity
# ============================================================================
addln("\n--- (5) biological sanity, top-20 gated F-stage hits ---")
top20 <- fstage[1:20]
# canonical fibrosis-related symbols
fibrosis_gene_set <- c(
  "COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2","COL6A1","COL6A2",
  "COL6A3","COL15A1","COL18A1","COL14A1",
  "TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2",
  "THBS1","THBS2","THBS4",
  "C3","C1QA","C1QB","C1QC","CFH","CFB",
  "PDGFA","PDGFB","PDGFC","PDGFD","PDGFRA","PDGFRB",
  "FN1","SPP1","DCN","BGN","LUM","POSTN","TIMP1","TIMP2","TIMP3","MMP2","MMP9","MMP14",
  "LOXL2","LOXL1","LOX","ITGAV","ITGB1","ITGB6","ITGB8",
  "CCR2","CCL2","CXCL10","CXCL12","CXCR4",
  "IGF1","IGFBP3","IGFBP7","HGF","EGF","EGFR",
  "APOA1","APOA2","APOB","APOC1","APOE","APOM","ALB","SERPINA1","SERPINE1",
  "CDH1","CDH11","CDH2","NRP1","NRP2","VEGFA","VEGFB",
  "PLAU","PLAUR","CD44","SPARC","SPARCL1",
  "NOTCH1","NOTCH2","NOTCH3","JAG1","DLL4","DLL1",
  "FGF1","FGF2","FGFR1","FGFR2"
)
top20[, lit_known := ligand_complex %in% fibrosis_gene_set |
                      receptor_complex %in% fibrosis_gene_set]
addln(sprintf("Top-20 with literature-known fibrosis member = %d / 20", sum(top20$lit_known)))
addln("Top 20 (rank, ct_pair, ligand__receptor, Estimate, p, lit_known):")
for (i in seq_len(nrow(top20))) {
  addln(sprintf("  %2d  %-50s  %-20s  est=%+.3f  p=%.2e  lit=%s",
                top20$rank[i], top20$ct_pair[i],
                paste0(top20$ligand_complex[i], "__", top20$receptor_complex[i]),
                top20$Estimate[i], top20$pval[i], top20$lit_known[i]))
}

# ============================================================================
# DIMENSION 6 -- Andrews-only refit (drop Wang)
# ============================================================================
addln("\n--- (6) Andrews-only refit (drop Wang) ---")
lr_long <- fread(MERGED_TSV)
addln(sprintf("Merged LR donor-rows = %d", nrow(lr_long)))

covar_cols <- c("sample","dataset","F_stage_documented","F_stage_inferred",
                "F_stage_source","age","sex_numeric")
covar_cols <- intersect(covar_cols, names(meta))
lr_long <- merge(lr_long, meta[, ..covar_cols], by = "sample", all.x = TRUE)

# Andrews-only = GSE202379 with F_stage_documented (the original documented set)
andrews_meta <- meta[dataset == "GSE202379" & !is.na(F_stage_documented)]
addln(sprintf("Andrews donors (GSE202379 with F_stage_documented) = %d",
              nrow(andrews_meta)))

lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long[, F_stage_numeric := as.numeric(F_stage_documented)]

lr_and <- lr_long[!is.na(F_stage_numeric) & dataset == "GSE202379"]
addln(sprintf("Andrews-only LR donor-rows after gating = %d", nrow(lr_and)))
addln(sprintf("Andrews-only unique donors                = %d",
              uniqueN(lr_and$sample)))

# Re-fit only the gated top-50 LR pairs on Andrews-only (single-dataset -> lm)
MIN_DONORS <- 20  # relaxed since Andrews-only is ~58 donors max
fit_andrews_one <- function(d) {
  if (nrow(d) < MIN_DONORS) return(NULL)
  rhs <- "F_stage_numeric + log10(pmax(n_source_cells,1)) + log10(pmax(n_target_cells,1))"
  fit <- try(stats::lm(stats::as.formula(paste("score ~", rhs)),
                       data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  co <- try(summary(fit)$coefficients, silent = TRUE)
  if (inherits(co, "try-error")) return(NULL)
  hits <- rownames(co)[startsWith(rownames(co), "F_stage_numeric")]
  if (!length(hits)) return(NULL)
  out <- as.data.table(co[hits, , drop = FALSE], keep.rownames = "term")
  setnames(out, c("term","Estimate_andrews","StdErr_andrews","tval_andrews","pval_andrews"))
  out[, n_donors_andrews := nrow(d)]
  out
}

addln("Re-fitting Andrews-only lm for each top-50 ct_pair x lr_pair ...")
# build key index for fast subset
setkey(lr_and, ct_pair, lr_pair)
top50[, key_ := paste(ct_pair, lr_pair, sep = "||")]
fit_list <- vector("list", nrow(top50))
for (i in seq_len(nrow(top50))) {
  d <- lr_and[ct_pair == top50$ct_pair[i] & lr_pair == top50$lr_pair[i]]
  r <- fit_andrews_one(d)
  if (!is.null(r)) {
    r[, ct_pair := top50$ct_pair[i]]
    r[, lr_pair := top50$lr_pair[i]]
    fit_list[[i]] <- r
  }
}
fit_and <- rbindlist(fit_list[!sapply(fit_list, is.null)], fill = TRUE)
addln(sprintf("Andrews-only fits succeeded for %d / 50 top hits", nrow(fit_and)))

# Merge full + Andrews-only estimates
top50_merged <- merge(top50, fit_and, by = c("ct_pair","lr_pair"), all.x = TRUE)
top50_merged[, sign_match := sign(Estimate) == sign(Estimate_andrews)]
addln(sprintf("Sign concordance gated-full vs Andrews-only = %d / %d (%.1f%%)",
              sum(top50_merged$sign_match, na.rm=TRUE),
              sum(!is.na(top50_merged$sign_match)),
              100*mean(top50_merged$sign_match, na.rm=TRUE)))
if (sum(!is.na(top50_merged$Estimate_andrews)) > 5) {
  cor_pe <- cor(top50_merged$Estimate, top50_merged$Estimate_andrews,
                use = "pairwise.complete.obs", method = "pearson")
  cor_sp <- cor(top50_merged$Estimate, top50_merged$Estimate_andrews,
                use = "pairwise.complete.obs", method = "spearman")
  addln(sprintf("Estimate correlation (gated vs Andrews-only): r=%.3f, rho=%.3f",
                cor_pe, cor_sp))
  addln(sprintf("  Andrews-only nominal p<0.05 among top-50 = %d (%.1f%%)",
                sum(top50_merged$pval_andrews < 0.05, na.rm=TRUE),
                100*mean(top50_merged$pval_andrews < 0.05, na.rm=TRUE)))
}

# Add lit_known flag to merged
top50_merged[, lit_known := ligand_complex %in% fibrosis_gene_set |
                            receptor_complex %in% fibrosis_gene_set]

# ----- final audit TSV ------------------------------------------------------
audit_cols <- c("rank","ct_pair","lr_pair","ligand_complex","receptor_complex",
                "Estimate","StdErr","pval","padj_within_ct","family_bonferroni",
                "n_donors","Estimate_andrews","StdErr_andrews","pval_andrews",
                "n_donors_andrews","sign_match","lit_known")
audit_cols <- intersect(audit_cols, names(top50_merged))
setorder(top50_merged, rank)
fwrite(top50_merged[, ..audit_cols], AUDIT_TSV, sep = "\t")
addln(sprintf("\n[output] %s (rows=%d)", AUDIT_TSV, nrow(top50_merged)))

# ============================================================================
# VERDICT
# ============================================================================
addln("\n====================================================================")
addln("VERDICT")
addln("====================================================================")
n_lit <- sum(top20$lit_known)
sign_concordance <- if (sum(!is.na(top50_merged$sign_match)) > 0) {
  mean(top50_merged$sign_match, na.rm=TRUE)
} else NA_real_
fstage_both <- if (nrow(fstage_bulk) > 0 && exists("fb_top") && nrow(fb_top) > 0) {
  mean(fb_top$both_concordant, na.rm=TRUE)
} else NA_real_
loo_mean <- if (exists("loo_top") && nrow(loo_top) > 0) {
  mean(loo_top$replication_rate, na.rm=TRUE)
} else NA_real_
addln(sprintf("  Top-20 lit-supported           : %d / 20",  n_lit))
addln(sprintf("  Sign concordance vs Andrews-only: %.1f%%",   100*sign_concordance))
addln(sprintf("  Top-50 bulk both-concordant    : %.1f%%",   100*fstage_both))
addln(sprintf("  Top-50 LOO mean replication    : %.3f",     loo_mean))

writeLines(buf, SUMMARY_TXT)
cat(sprintf("\n[output] %s\n", SUMMARY_TXT))
cat("[done] gated F-stage audit complete\n")
