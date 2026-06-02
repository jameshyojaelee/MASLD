#!/usr/bin/env Rscript
# 353_celltype_sex_coloc.R — Cell-type x sex_class x COLOC cross-tab + Fisher.
#
# For each of the 5 cell types, joins the per-cell-type sex_class (from Script 352)
# with the canonical SuSiE-COLOC table (per-GWAS PP4) and asks whether sex-biased
# genes are enriched for colocalization.
#
# Headline target: any (cell_type, gwas_family, sex_class) triple with OR > 2 at
# padj < 0.1 — i.e. localized sex x genetic-causal convergence that bulk doesn't
# show.
#
# Outputs:
#   RNA-seq/results/stratified_causal/sex_celltype_coloc.csv
#   RNA-seq/results/stratified_causal/sctwas_sex_celltype.csv
#   RNA-seq/results/stratified_causal/sex_celltype_coloc_summary.txt

suppressPackageStartupMessages({
  library(data.table)
})

proj_root <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
setwd(proj_root)

ddir_sc  <- "Analysis/SingleCell/results_gpu_v2/sex_celltype"
ddir_out <- "RNA-seq/results/stratified_causal"
dir.create(ddir_out, recursive = TRUE, showWarnings = FALSE)

PP4_THRESH <- 0.5

celltypes <- c("hepatocytes", "macrophages", "fibroblasts", "endothelial", "cholangiocytes")

# ---- Load COLOC: per-gene per-GWAS ----
coloc_long <- fread("GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")
cat("susie_coloc_all_gwas:", nrow(coloc_long), "rows;",
    length(unique(coloc_long$gene)), "genes;",
    length(unique(coloc_long$gwas_name)), "GWAS\n")

# Per-gene best PP4 (susie if present, else abf)
coloc_long[, best_pp4 := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
# Build GWAS-family map
gwas_to_family <- function(g) {
  g <- toupper(g)
  fcase(
    grepl("HCC|HEPATOCELLU", g), "hcc",
    grepl("CIRRHOSIS|CHIRHEP", g), "cirrhosis",
    grepl("ALT|AST|GGT", g), "enzyme",
    grepl("NAFLD|NASH|PDFF|STEATO", g), "steatosis",
    default = "other"
  )
}
coloc_long[, gwas_family := gwas_to_family(gwas_name)]
cat("GWAS family counts:\n"); print(coloc_long[, .N, by = gwas_family])

# Per-gene x family: max PP4
coloc_gf <- coloc_long[, .(best_pp4 = max(best_pp4, na.rm = TRUE),
                           n_gwas_tested = .N,
                           n_gwas_h4_05 = sum(best_pp4 > PP4_THRESH, na.rm = TRUE)),
                       by = .(gene, gwas_family)]
coloc_gf[is.infinite(best_pp4), best_pp4 := NA_real_]
# Also per-gene any (all families)
coloc_any <- coloc_long[, .(best_pp4 = max(best_pp4, na.rm = TRUE),
                            n_gwas_tested = .N,
                            n_gwas_h4_05 = sum(best_pp4 > PP4_THRESH, na.rm = TRUE)),
                        by = gene]
coloc_any[, gwas_family := "any"]
coloc_any[is.infinite(best_pp4), best_pp4 := NA_real_]
coloc_gf <- rbind(coloc_gf, coloc_any)

# ---- Fisher per (celltype, gwas_family, sex_class) ----
fisher_row <- function(n_class, n_coloc_in_class, n_other, n_coloc_in_other,
                       celltype, gwas_family, sex_class) {
  tbl <- matrix(c(n_coloc_in_class,
                  n_class - n_coloc_in_class,
                  n_coloc_in_other,
                  n_other - n_coloc_in_other),
                nrow = 2, byrow = TRUE)
  if (any(tbl < 0) || any(rowSums(tbl) == 0) || any(colSums(tbl) == 0)) {
    return(data.table(celltype = celltype, gwas_family = gwas_family,
                      sex_class = sex_class,
                      n_class = n_class, n_coloc_in_class = n_coloc_in_class,
                      n_other = n_other, n_coloc_in_other = n_coloc_in_other,
                      OR = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                      pvalue = NA_real_, direction = "ND"))
  }
  ft <- fisher.test(tbl)
  data.table(celltype = celltype, gwas_family = gwas_family,
             sex_class = sex_class,
             n_class = n_class, n_coloc_in_class = n_coloc_in_class,
             n_other = n_other, n_coloc_in_other = n_coloc_in_other,
             OR = unname(ft$estimate),
             ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
             pvalue = ft$p.value,
             direction = ifelse(ft$estimate > 1, "enriched", "depleted"))
}

results <- list()
for (ct in celltypes) {
  f <- file.path(ddir_sc, paste0(ct, "_sex_dream_results.csv"))
  if (!file.exists(f)) {
    cat("skip", ct, "(missing", f, ")\n")
    next
  }
  cls <- fread(f)
  if (nrow(cls) == 0) {
    cat("skip", ct, "— empty dream results\n")
    next
  }
  cat("\n=== ", ct, " — ", nrow(cls), " genes ===\n", sep = "")
  print(cls[, .N, by = sex_class])

  for (gf in unique(coloc_gf$gwas_family)) {
    cg <- coloc_gf[gwas_family == gf]
    merged <- merge(cls[, .(gene, sex_class)], cg[, .(gene, best_pp4)],
                    by = "gene", all.x = TRUE)
    merged[, is_coloc := !is.na(best_pp4) & best_pp4 > PP4_THRESH]
    for (sc in c("Female_biased", "Male_biased", "Divergent")) {
      n_class <- sum(merged$sex_class == sc, na.rm = TRUE)
      n_coloc_in_class <- sum(merged$sex_class == sc & merged$is_coloc, na.rm = TRUE)
      n_other <- sum(merged$sex_class == "Concordant", na.rm = TRUE)
      n_coloc_in_other <- sum(merged$sex_class == "Concordant" & merged$is_coloc, na.rm = TRUE)
      if (n_class == 0) next
      results[[length(results) + 1]] <- fisher_row(
        n_class, n_coloc_in_class, n_other, n_coloc_in_other,
        ct, gf, sc)
    }
  }
}

if (length(results) == 0) {
  cat("WARNING: no Fisher rows produced.\n")
  fwrite(data.table(), file.path(ddir_out, "sex_celltype_coloc.csv"))
  quit(status = 0)
}

res_dt <- rbindlist(results)
res_dt[, padj := p.adjust(pvalue, method = "BH")]
setorder(res_dt, pvalue)
fwrite(res_dt, file.path(ddir_out, "sex_celltype_coloc.csv"))
cat("\nWrote sex_celltype_coloc.csv:", nrow(res_dt), "rows\n")

# Headline hits
top <- res_dt[!is.na(OR) & OR > 2 & padj < 0.1]
cat("HEADLINE: (OR>2, padj<0.1):", nrow(top), "\n")
if (nrow(top) > 0) print(top)

# ---- Sc-TWAS x sex ----
sctwas_path <- "RNA-seq/results/causal_inference/sceqtl_twas/sceqtl_twas_all_results.csv"
if (file.exists(sctwas_path)) {
  sctwas <- fread(sctwas_path)
  cat("\nsctwas rows:", nrow(sctwas), "; celltypes:", paste(unique(sctwas$cell_type), collapse=","), "\n")
  # Significance: global_fdr < 0.1
  sctwas[, is_sig := !is.na(global_fdr) & global_fdr < 0.1]
  setnames(sctwas, "exposure", "gene", skip_absent = TRUE)
  # Map sctwas cell_type to our slug list
  ct_map <- c("hepatocyte" = "hepatocytes",
              "endothelial_cell" = "endothelial",
              "cholangiocyte" = "cholangiocytes",
              "stellate_cell" = "fibroblasts",
              "kupffer_cell" = "macrophages",
              "macrophage" = "macrophages")
  sctwas[, ct_slug := ct_map[cell_type]]
  out_sct <- list()
  for (ct in unique(sctwas$ct_slug)) {
    if (is.na(ct)) next
    f <- file.path(ddir_sc, paste0(ct, "_sex_dream_results.csv"))
    if (!file.exists(f)) next
    cls <- fread(f)
    if (nrow(cls) == 0) next
    sub <- sctwas[ct_slug == ct]
    merged <- merge(cls[, .(gene, sex_class)],
                    sub[, .(gene, is_sig, gwas)],
                    by = "gene")
    if (nrow(merged) == 0) next
    merged[, gwas_family := gwas_to_family(gwas)]
    for (gf in unique(merged$gwas_family)) {
      sub2 <- merged[gwas_family == gf]
      for (sc in c("Female_biased", "Male_biased", "Divergent")) {
        n_class <- sum(sub2$sex_class == sc, na.rm = TRUE)
        n_sig_in_class <- sum(sub2$sex_class == sc & sub2$is_sig, na.rm = TRUE)
        n_other <- sum(sub2$sex_class == "Concordant", na.rm = TRUE)
        n_sig_in_other <- sum(sub2$sex_class == "Concordant" & sub2$is_sig, na.rm = TRUE)
        if (n_class == 0) next
        tbl <- matrix(c(n_sig_in_class, max(0, n_class - n_sig_in_class),
                        n_sig_in_other, max(0, n_other - n_sig_in_other)),
                      nrow = 2, byrow = TRUE)
        if (any(rowSums(tbl) == 0) || any(colSums(tbl) == 0)) next
        ft <- fisher.test(tbl)
        out_sct[[length(out_sct) + 1]] <- data.table(
          celltype = ct, gwas_family = gf, sex_class = sc,
          n_class = n_class, n_sig_in_class = n_sig_in_class,
          n_other = n_other, n_sig_in_other = n_sig_in_other,
          OR = unname(ft$estimate),
          ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
          pvalue = ft$p.value,
          direction = ifelse(ft$estimate > 1, "enriched", "depleted"))
      }
    }
  }
  if (length(out_sct) > 0) {
    sct_dt <- rbindlist(out_sct)
    sct_dt[, padj := p.adjust(pvalue, method = "BH")]
    setorder(sct_dt, pvalue)
    fwrite(sct_dt, file.path(ddir_out, "sctwas_sex_celltype.csv"))
    cat("\nWrote sctwas_sex_celltype.csv:", nrow(sct_dt), "rows\n")
  } else {
    cat("\nNo sctwas Fisher rows produced.\n")
  }
} else {
  cat("\nNo sceqtl_twas_all_results.csv found at", sctwas_path, "\n")
}

# ---- Summary log ----
sink(file.path(ddir_out, "sex_celltype_coloc_summary.txt"))
cat("=== B4 Cell-type x Sex x COLOC summary ===\n")
cat("PP4 threshold:", PP4_THRESH, "\n")
cat("N Fisher rows:", nrow(res_dt), "\n\n")
cat("HEADLINE candidates (OR > 2, padj < 0.1):\n")
print(top)
cat("\nAll celltype x family enrichment medians:\n")
print(res_dt[, .(median_OR = median(OR, na.rm = TRUE),
                 min_padj = min(padj, na.rm = TRUE)),
             by = .(celltype, gwas_family)])
sink()

cat("\nDone.\n")
