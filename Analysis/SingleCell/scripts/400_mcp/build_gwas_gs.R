#!/usr/bin/env Rscript
# Build GWAS-anchored scDRS .gs file with side-by-side 1KG vs PolyFun panels.
#
# Output: 17 traits in GWAS_anchored_scdrs.gs
#   12 panel-dependent traits (6 trait classes × 2 LD panels: _1kg, _polyfun)
#    5 sceQTL traits (panel-independent — OneK1K LD, not affected by EUR LD swap)
#
# Both LD panels cover the same 17 EUR/Finnish/Ghouse GWAS so the comparison
# isolates the LD panel choice (1KG vs PolyFun-UKBB functional priors).
#
# Trait naming convention:
#   gwas_coloc_pooled_pp05_{1kg,polyfun}
#   gwas_coloc_pooled_pp08_{1kg,polyfun}
#   gwas_coloc_nafld_specific_{1kg,polyfun}
#   gwas_coloc_liver_enzymes_{1kg,polyfun}
#   gwas_coloc_pdff_{1kg,polyfun}
#   gwas_coloc_cirrhosis_hcc_{1kg,polyfun}
#   sceqtl_{hep,endo,chol,stell,merged_4ct}        (no panel suffix)

suppressPackageStartupMessages({
  library(data.table)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                           "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
setwd(PROJECT_ROOT)

OUT_DIR     <- file.path(PROJECT_ROOT,
                         "Analysis/SingleCell/results_gpu_v2/disease_signatures")
OUT_GS      <- file.path(OUT_DIR, "GWAS_anchored_scdrs.gs")
SUMMARY_CSV <- file.path(OUT_DIR, "GWAS_anchored_scdrs_summary.csv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LD_PANELS <- c("1kg", "polyfun")

PANEL_DIRS <- list(
  `1kg`     = "GWAS/finemapping/results/susie_coloc_1kg",
  polyfun   = "GWAS/finemapping/results/susie_coloc_polyfun"
)
PANEL_AGG_FILES <- list(
  `1kg`     = "GWAS/finemapping/results/susie_coloc_1kg/gene_level_coloc_1kg.csv",
  polyfun   = "GWAS/finemapping/results/susie_coloc_polyfun/gene_level_coloc_polyfun.csv"
)

# Phenotype-class → list of GWAS subdir names (must match per-panel dir contents)
PHENOTYPE_CLASSES <- list(
  nafld_specific = c(
    "2019_31311600_NAFLD_EUR",
    "2020_32298765_NAFLD_EUR",
    "2021_34841290_NAFLD_EUR",
    "2023_36280732_NAFLD_deCode_EUR",
    "2023_36280732_NAFLD_Intermountain_EUR",
    "2023_36280732_NAFLD_UKBB_EUR",
    "FinnGen_NAFLD",
    "FinnGen_NASH"
  ),
  liver_enzymes = c("UKBB_ALT", "UKBB_AST", "UKBB_GGT"),
  pdff = c(
    "2021_34128465_PDFF_EUR",
    "2021_34957434_PDFF_EUR",
    "2022_36402844_PDFF_EUR"
  ),
  cirrhosis_hcc = c("Ghouse_Cirrhosis", "Ghouse_HCC", "FinnGen_HCC")
)

SCEQTL_FILES <- list(
  sceqtl_hep   = "RNA-seq/results/causal_inference/sceqtl/coloc_results_hepatocyte.csv",
  sceqtl_endo  = "RNA-seq/results/causal_inference/sceqtl/coloc_results_endothelial_cell.csv",
  sceqtl_chol  = "RNA-seq/results/causal_inference/sceqtl/coloc_results_cholangiocyte.csv",
  sceqtl_stell = "RNA-seq/results/causal_inference/sceqtl/coloc_results_stellate_cell.csv"
)

SCEQTL_FALLBACK_TOPN  <- 25L
SCEQTL_MIN_FOR_STRICT <- 10L

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

format_geneset_row <- function(trait, genes, weights) {
  ok <- !is.na(genes) & nzchar(genes) & !is.na(weights) & is.finite(weights) & weights > 0
  genes <- genes[ok]; weights <- weights[ok]
  if (length(genes) == 0) {
    warning(sprintf("[%s] No valid genes", trait))
    return(NULL)
  }
  dt <- data.table(g = genes, w = weights)
  dt <- dt[, .(w = max(w)), by = g]
  setorder(dt, -w)
  paste0(trait, "\t", paste0(dt$g, ":", signif(dt$w, 6), collapse = ","))
}

write_gs_file <- function(rows, path) {
  rows <- Filter(Negate(is.null), rows)
  con <- file(path, "w")
  writeLines("TRAIT\tGENESET", con)
  for (r in rows) writeLines(r, con)
  close(con)
  cat(sprintf("Wrote %d traits to %s\n", length(rows), path))
}

SUMMARY_ROWS <- list()
record_summary <- function(trait, genes, weights, source_files) {
  ok <- !is.na(weights) & is.finite(weights) & weights > 0 &
        !is.na(genes) & nzchar(genes)
  weights <- weights[ok]; genes <- genes[ok]
  dt <- data.table(g = genes, w = weights)
  dt <- dt[, .(w = max(w)), by = g]
  SUMMARY_ROWS[[length(SUMMARY_ROWS) + 1L]] <<- data.table(
    trait        = trait,
    n_genes      = nrow(dt),
    weight_min   = if (nrow(dt)) min(dt$w) else NA_real_,
    weight_max   = if (nrow(dt)) max(dt$w) else NA_real_,
    weight_median= if (nrow(dt)) median(dt$w) else NA_real_,
    source_files = source_files
  )
  if (nrow(dt) < 10L) {
    warning(sprintf("[%s] < 10 genes (n=%d) — scDRS may be unstable",
                    trait, nrow(dt)))
  }
  cat(sprintf("  [%-44s]  n_genes=%5d  weight=[%.3f, %.3f]\n",
              trait, nrow(dt),
              if (nrow(dt)) min(dt$w) else NA_real_,
              if (nrow(dt)) max(dt$w) else NA_real_))
}

# Read per-GWAS susie_coloc_combined.csv → return data.table(gene, pp4)
load_per_gwas_combined <- function(panel, gwas_subdir) {
  fp <- file.path(PANEL_DIRS[[panel]], gwas_subdir, "susie_coloc_combined.csv")
  if (!file.exists(fp)) {
    warning(sprintf("[%s/%s] missing: %s", panel, gwas_subdir, fp))
    return(data.table(gene = character(0), pp4 = numeric(0)))
  }
  d <- fread(fp, select = c("gene", "PP.H4.susie", "PP.H4.abf"))
  # Prefer SuSiE PP4; fall back to ABF if SuSiE is NA
  d[, pp4 := fifelse(!is.na(PP.H4.susie), PP.H4.susie, PP.H4.abf)]
  d <- d[!is.na(gene) & nzchar(gene) & !is.na(pp4)]
  d[, .(pp4 = max(pp4)), by = gene]
}

# ----------------------------------------------------------------------------
# Build panel-dependent traits (12 total)
# ----------------------------------------------------------------------------

panel_rows <- list()

for (panel in LD_PANELS) {
  cat(sprintf("\n=== LD PANEL: %s ===\n", panel))

  # ------- Pooled (PP05 / PP08) from gene_level_coloc aggregated file -------
  agg_fp <- PANEL_AGG_FILES[[panel]]
  cat(sprintf("[%s] Pooled gene-level COLOC: %s\n", panel, agg_fp))
  coloc <- fread(agg_fp)
  coloc[, best_pp4 := fifelse(!is.na(coloc_best_susie_pp4),
                              coloc_best_susie_pp4,
                              coloc_best_pp4)]
  coloc <- coloc[!is.na(gene) & nzchar(gene)]
  coloc <- coloc[, .(best_pp4 = max(best_pp4, na.rm = TRUE)), by = gene]
  t05 <- coloc[best_pp4 > 0.5]
  t08 <- coloc[best_pp4 > 0.8]
  trait05 <- sprintf("gwas_coloc_pooled_pp05_%s", panel)
  trait08 <- sprintf("gwas_coloc_pooled_pp08_%s", panel)
  record_summary(trait05, t05$gene, t05$best_pp4, agg_fp)
  record_summary(trait08, t08$gene, t08$best_pp4, agg_fp)
  panel_rows[[length(panel_rows) + 1L]] <-
    format_geneset_row(trait05, t05$gene, t05$best_pp4)
  panel_rows[[length(panel_rows) + 1L]] <-
    format_geneset_row(trait08, t08$gene, t08$best_pp4)

  # ------- Phenotype-stratified (4 trait classes per panel) -----------------
  for (pheno in names(PHENOTYPE_CLASSES)) {
    gwas_list <- PHENOTYPE_CLASSES[[pheno]]
    cat(sprintf("[%s/%s] %d GWAS\n", panel, pheno, length(gwas_list)))
    sub <- rbindlist(lapply(gwas_list,
                            function(g) load_per_gwas_combined(panel, g)),
                     use.names = TRUE)
    sub <- sub[, .(pp4 = max(pp4)), by = gene]
    filt <- sub[pp4 > 0.5]
    src <- paste(file.path(PANEL_DIRS[[panel]], gwas_list,
                            "susie_coloc_combined.csv"), collapse = "; ")
    trait_name <- sprintf("gwas_coloc_%s_%s", pheno, panel)
    record_summary(trait_name, filt$gene, filt$pp4, src)
    panel_rows[[length(panel_rows) + 1L]] <-
      format_geneset_row(trait_name, filt$gene, filt$pp4)
  }
}

# ----------------------------------------------------------------------------
# Build sceQTL traits (5 total, panel-independent)
# ----------------------------------------------------------------------------

cat("\n=== sceQTL traits (panel-independent) ===\n")
sceqtl_rows <- list()
sceqtl_full_list <- list()

for (nm in names(SCEQTL_FILES)) {
  fp <- SCEQTL_FILES[[nm]]
  d <- fread(fp)
  stopifnot("PP.H4_masld" %in% names(d), "gene" %in% names(d))
  d <- d[!is.na(gene) & nzchar(gene) & !is.na(PP.H4_masld)]
  d <- d[, .(PP.H4_masld = max(PP.H4_masld)), by = gene]
  filt <- d[PP.H4_masld > 0.5]
  src <- fp
  if (nrow(filt) < SCEQTL_MIN_FOR_STRICT) {
    setorder(d, -PP.H4_masld)
    filt <- head(d, SCEQTL_FALLBACK_TOPN)
    src <- paste0(fp, " [top-", SCEQTL_FALLBACK_TOPN,
                  " fallback; strict PP>0.5 yielded n<10]")
    cat(sprintf("  [%s] strict filter yielded < 10; using top-%d\n",
                nm, SCEQTL_FALLBACK_TOPN))
  }
  record_summary(nm, filt$gene, filt$PP.H4_masld, src)
  sceqtl_rows[[length(sceqtl_rows) + 1L]] <-
    format_geneset_row(nm, filt$gene, filt$PP.H4_masld)
  sceqtl_full_list[[nm]] <- d
}

# ------- Trait: sceqtl_merged_4ct -----------------------------------------
merged_full <- rbindlist(sceqtl_full_list, use.names = TRUE)
merged_full <- merged_full[, .(PP.H4_masld = max(PP.H4_masld)), by = gene]
merged <- merged_full[PP.H4_masld > 0.5]
merged_src <- paste(unname(SCEQTL_FILES), collapse = "; ")
if (nrow(merged) < SCEQTL_MIN_FOR_STRICT) {
  setorder(merged_full, -PP.H4_masld)
  merged <- head(merged_full, SCEQTL_FALLBACK_TOPN)
  merged_src <- paste0(merged_src, " [top-", SCEQTL_FALLBACK_TOPN,
                       " merged-fallback]")
  cat(sprintf("  [merged] strict filter yielded < 10; using top-%d\n",
              SCEQTL_FALLBACK_TOPN))
}
record_summary("sceqtl_merged_4ct", merged$gene, merged$PP.H4_masld, merged_src)
sceqtl_rows[[length(sceqtl_rows) + 1L]] <-
  format_geneset_row("sceqtl_merged_4ct", merged$gene, merged$PP.H4_masld)

# ----------------------------------------------------------------------------
# Write final .gs (17 rows: 12 panel-dependent + 5 sceQTL)
# ----------------------------------------------------------------------------

all_rows <- c(panel_rows, sceqtl_rows)
write_gs_file(all_rows, OUT_GS)

summary_dt <- rbindlist(SUMMARY_ROWS)
fwrite(summary_dt, SUMMARY_CSV)
cat(sprintf("\nSummary written to %s\n", SUMMARY_CSV))

cat("\n=== Cross-panel sanity checks ===\n")
# Compare pooled-PP05 between 1KG and PolyFun
g_1kg <- summary_dt[trait == "gwas_coloc_pooled_pp05_1kg",     n_genes]
g_pf  <- summary_dt[trait == "gwas_coloc_pooled_pp05_polyfun", n_genes]
cat(sprintf("Pooled PP4>0.5: 1kg n=%d, polyfun n=%d, delta=%+d\n",
            g_1kg, g_pf, g_pf - g_1kg))

cat("\nDone.\n")
