#!/usr/bin/env Rscript
# 76b_govaere2026_crossmodal.R
# Cross-validate our existing pipeline signatures against Govaere 2026 Nat
# Genetics signatures (GeoMx spatial, snRNA-seq GPNMB+ macs, MetMac/LAM/IL32
# axis panels).
#
# Four checks:
#   1. Hepatocyte Progressor subtype signature vs Govaere GeoMx SH-vs-LS
#      (DE file from Spatial pipeline W3.A1, plus Sup Table 11 lookup).
#   2. Macrophage pseudotime late genes vs Govaere MetMac/LAM/GPNMB panels.
#   3. IL32 axis panel: per-gene atlas + Govaere signature column table.
#   4. Bulk dream DEGs vs Sup Table 11 (snRNA GPNMB+ MASH-vs-noMASH).
#
# Reads (relative to MASLD_PROJECT_ROOT):
#   Analysis/Spatial/results/govaere2026/geomx_de_sh_vs_ls.csv         (optional, W3.A1)
#   data/external/govaere2026_natgenetics/suptables/41588_2026_2600_MOESM4_ESM__Supplementary_Table_11.tsv
#   RNA-seq/results/govaere2026/govaere2026_signatures_wide.tsv
#   RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results_ashr.csv
#   Analysis/SingleCell/results_gpu_v2/pseudotime/pseudotime_corr_Macrophages.csv
#   Analysis/SingleCell/results_gpu_v2/pseudobulk_de/Hepatocytes_de.csv  (for sc_hepatocyte_logFC)
#
# Writes:
#   RNA-seq/results/govaere2026_crossmodal/progressor_signature_vs_geomx_sh_ls.csv
#   RNA-seq/results/govaere2026_crossmodal/progressor_signature_summary.txt
#   RNA-seq/results/govaere2026_crossmodal/macrophage_pseudotime_vs_metmac.csv
#   RNA-seq/results/govaere2026_crossmodal/macrophage_pseudotime_summary.txt
#   RNA-seq/results/govaere2026_crossmodal/il32_axis_atlas_lookup.tsv
#   RNA-seq/results/govaere2026_crossmodal/bulk_dream_vs_sup11.csv
#   RNA-seq/results/govaere2026_crossmodal/bulk_dream_vs_sup11_summary.txt
#   RNA-seq/results/govaere2026_crossmodal/README.md
#
# Env: rnaseq. Convention: data.table only (no dplyr).

suppressPackageStartupMessages({
  library(data.table)
})

PROJECT_ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/govaere2026_crossmodal")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[76b] project root: %s\n", PROJECT_ROOT))
cat(sprintf("[76b] output dir : %s\n", OUT_DIR))

# ---- Input paths ----------------------------------------------------------
GEOMX_DE_FILE <- file.path(PROJECT_ROOT,
  "Analysis/Spatial/results/govaere2026/geomx_de_sh_vs_ls.csv")
SUP11_FILE <- file.path(PROJECT_ROOT,
  "data/external/govaere2026_natgenetics/suptables",
  "41588_2026_2600_MOESM4_ESM__Supplementary_Table_11.tsv")
SIG_WIDE_FILE <- file.path(PROJECT_ROOT,
  "RNA-seq/results/govaere2026/govaere2026_signatures_wide.tsv")
ATLAS_FILE <- file.path(PROJECT_ROOT,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
DREAM_FILE <- file.path(PROJECT_ROOT,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "canonical_deg_results.csv")
MACRO_CORR_FILE <- file.path(PROJECT_ROOT,
  "Analysis/SingleCell/results_gpu_v2/pseudotime/pseudotime_corr_Macrophages.csv")
HEP_DE_FILE <- file.path(PROJECT_ROOT,
  "Analysis/SingleCell/results_gpu_v2/pseudobulk_de/Hepatocytes_de.csv")

# Progressor signature (canonical top markers per
# memory/govaere2026_integration.md + Hepatocyte Subclustering Pipeline 2026-04-10).
PROGRESSOR_GENES <- c("COL15A1", "CDH11", "LTBP2")
IL32_PANEL <- c("IL32", "GPNMB", "LPL", "AKR1B10", "FABP5", "HLA-DRA")

# Helper for safe Wilcoxon / Fisher
safe_wilcox <- function(x, y) {
  res <- tryCatch(
    suppressWarnings(wilcox.test(x, y, alternative = "two.sided")),
    error = function(e) NULL
  )
  if (is.null(res)) return(list(p = NA_real_, W = NA_real_))
  list(p = res$p.value, W = unname(res$statistic))
}

safe_fisher <- function(a, b, c, d) {
  m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  if (any(rowSums(m) == 0) || any(colSums(m) == 0))
    return(list(OR = NA_real_, p = NA_real_))
  res <- tryCatch(fisher.test(m), error = function(e) NULL)
  if (is.null(res)) return(list(OR = NA_real_, p = NA_real_))
  list(OR = unname(res$estimate), p = res$p.value)
}

# ===========================================================================
# Load shared resources
# ===========================================================================
cat("[76b] loading shared resources...\n")
sig_wide <- fread(SIG_WIDE_FILE)
cat(sprintf("  signatures_wide: %d genes x %d cols\n", nrow(sig_wide), ncol(sig_wide)))

sup11 <- fread(SUP11_FILE, skip = 1L)  # row 1 is title, row 2 is header
setnames(sup11, c("gene", "log2FC", "p_val", "p_val_adj"))
sup11 <- sup11[!is.na(log2FC)]
cat(sprintf("  sup11 GPNMB+ DE  : %d genes\n", nrow(sup11)))

# Atlas (large; we only need IL32 panel + dream cols, so don't shrink yet)
atlas <- fread(ATLAS_FILE)
stopifnot(all(c("bulk_padj", "bulk_logFC") %in% names(atlas)))
cat(sprintf("  atlas            : %d genes x %d cols\n",
            nrow(atlas), ncol(atlas)))

# Build sc_hepatocyte_logFC from pseudobulk Hepatocytes_de.csv (MASLD_vs_Healthy
# rows are symbol-keyed). Atlas itself does not carry this column.
hep_de_all <- fread(HEP_DE_FILE)
hep_de <- hep_de_all[contrast == "MASLD_vs_Healthy"
                     & !grepl("^ENSG", gene), .(gene, logFC, padj)]
setnames(hep_de, c("symbol", "sc_hepatocyte_logFC", "sc_hepatocyte_padj"))
cat(sprintf("  sc Hep DE rows   : %d (symbol-keyed, MASLD_vs_Healthy)\n",
            nrow(hep_de)))

# ===========================================================================
# CHECK 1: Hepatocyte Progressor signature vs Govaere GeoMx + Sup11
# ===========================================================================
cat("\n[76b] CHECK 1: Progressor signature vs GeoMx SH-vs-LS\n")

geomx_available <- file.exists(GEOMX_DE_FILE)
cat(sprintf("  GeoMx DE file present? %s\n", geomx_available))

prog_rows <- data.table(
  gene = PROGRESSOR_GENES,
  source = "progressor_top_marker"
)

# Pull GeoMx SH-vs-LS values: prefer fresh DE file; fall back to signatures_wide
# pre-canonicalized columns (these encode the paper's reported SH-vs-LS DE).
if (geomx_available) {
  geomx_de <- fread(GEOMX_DE_FILE)
  # Be tolerant of column naming until W3.A1 freezes its schema.
  sym_col <- intersect(c("gene", "symbol", "human_symbol"), names(geomx_de))[1]
  lfc_col <- intersect(c("logFC", "log2FC", "log2FoldChange", "lfc"), names(geomx_de))[1]
  pad_col <- intersect(c("padj", "p_val_adj", "adj.P.Val", "FDR"), names(geomx_de))[1]
  if (any(is.na(c(sym_col, lfc_col, pad_col)))) {
    cat("  WARNING: could not parse GeoMx DE schema; falling back to signatures_wide.\n")
    geomx_available <- FALSE
  } else {
    setnames(geomx_de, c(sym_col, lfc_col, pad_col),
             c("gene", "geomx_logFC", "geomx_padj"))
    geomx_de <- geomx_de[, .(gene, geomx_logFC, geomx_padj)]
  }
}
if (!geomx_available) {
  # TODO(W3.A1): rerun once Analysis/Spatial/results/govaere2026/geomx_de_sh_vs_ls.csv
  # is on disk. Until then, use paper-reported logFC stored in signatures_wide.
  cat("  TODO: GeoMx DE file not on disk yet. Using paper-reported SH-vs-LS\n")
  cat("        values from govaere2026_signatures_wide.tsv. Rerun 76b once\n")
  cat("        Spatial pipeline W3.A1 deposits geomx_de_sh_vs_ls.csv.\n")
  geomx_de <- sig_wide[, .(
    gene = human_symbol,
    geomx_logFC = as.numeric(signature_govaere2026_geomx_sh_vs_ls_logfc),
    geomx_padj  = as.numeric(signature_govaere2026_geomx_sh_vs_ls_padj)
  )][!is.na(geomx_logFC)]
}
cat(sprintf("  GeoMx SH-vs-LS rows used: %d\n", nrow(geomx_de)))

# Sup11 lookup for the Progressor genes (snRNA GPNMB+ MASH-vs-noMASH).
prog_sup11 <- sup11[gene %in% PROGRESSOR_GENES,
                    .(gene, sup11_log2FC = log2FC, sup11_padj = p_val_adj)]

prog_geomx <- geomx_de[gene %in% PROGRESSOR_GENES,
                       .(gene, geomx_logFC, geomx_padj)]

prog_out <- merge(prog_rows, prog_geomx, by = "gene", all.x = TRUE)
prog_out <- merge(prog_out, prog_sup11, by = "gene", all.x = TRUE)
prog_out[, geomx_sign_positive := geomx_logFC > 0]
fwrite(prog_out,
       file.path(OUT_DIR, "progressor_signature_vs_geomx_sh_ls.csv"))

# Wilcoxon: progressor genes' GeoMx logFC vs background (all other GeoMx genes).
fg_lfc <- geomx_de[gene %in% PROGRESSOR_GENES, geomx_logFC]
bg_lfc <- geomx_de[!gene %in% PROGRESSOR_GENES, geomx_logFC]
wx <- safe_wilcox(fg_lfc, bg_lfc)

mean_fg <- if (length(fg_lfc) > 0) mean(fg_lfc, na.rm = TRUE) else NA_real_
mean_bg <- if (length(bg_lfc) > 0) mean(bg_lfc, na.rm = TRUE) else NA_real_
n_fg_pos <- if (length(fg_lfc) > 0) sum(fg_lfc > 0, na.rm = TRUE) else 0L
frac_fg_pos <- if (length(fg_lfc) > 0) n_fg_pos / sum(!is.na(fg_lfc)) else NA_real_

prog_summary <- c(
  sprintf("Progressor signature (COL15A1, CDH11, LTBP2) vs Govaere GeoMx SH-vs-LS"),
  sprintf("(headline female-enriched hepatocyte subtype, 1.3%% -> 63.1%% expansion)"),
  "",
  sprintf("GeoMx DE source         : %s",
          if (geomx_available) "fresh DE (geomx_de_sh_vs_ls.csv)"
          else "TODO -- paper-reported (signatures_wide)"),
  sprintf("Genes tested            : %d / %d in GeoMx universe",
          nrow(prog_geomx), length(PROGRESSOR_GENES)),
  sprintf("Mean Progressor logFC   : %.3f", mean_fg),
  sprintf("Mean background logFC   : %.3f", mean_bg),
  sprintf("Fraction Progressor > 0 : %s (%d / %d)",
          if (is.na(frac_fg_pos)) "NA" else sprintf("%.1f%%", 100 * frac_fg_pos),
          n_fg_pos, sum(!is.na(fg_lfc))),
  sprintf("Wilcoxon p (fg vs bg)   : %s",
          if (is.na(wx$p)) "NA" else format.pval(wx$p, digits = 3, eps = 1e-300)),
  "",
  "Sup Table 11 (snRNA-seq GPNMB+ MASH-vs-noMASH) lookup for Progressor genes:",
  if (nrow(prog_sup11) == 0) "  (no matches in Sup11)"
  else paste0("  ",
              sprintf("%-10s log2FC=%+.3f padj=%.3g",
                      prog_sup11$gene, prog_sup11$sup11_log2FC,
                      prog_sup11$sup11_padj),
              collapse = "\n"),
  "",
  "Interpretation:",
  paste0("  Progressor signature is hepatocyte-stromal (COL15A1/CDH11/LTBP2 are",
         " stromal/EMT-like)."),
  paste0("  Paper's GeoMx SH lanes are hepatocyte-lobule, so these genes",
         " may be weak in GeoMx."),
  paste0("  Sup11 is GPNMB+ macrophage and reports macrophage genes only,",
         " so Progressor genes can be absent there."),
  paste0("  Report this for transparency; positive enrichment is plausible but",
         " not required.")
)
writeLines(prog_summary,
           file.path(OUT_DIR, "progressor_signature_summary.txt"))
cat(sprintf("  -> %d Progressor genes, mean lfc=%.3f, Wilcoxon p=%.3g\n",
            length(PROGRESSOR_GENES), mean_fg, wx$p))

# ===========================================================================
# CHECK 2: Macrophage pseudotime late genes vs Govaere MetMac/LAM/GPNMB
# ===========================================================================
cat("\n[76b] CHECK 2: macrophage pseudotime vs MetMac/LAM/GPNMB+\n")

macro_corr <- fread(MACRO_CORR_FILE)
setnames(macro_corr,
         c("gene", "spearman_rho", "pval", "padj"),
         c("symbol", "macro_rho", "macro_pval", "macro_padj"))
cat(sprintf("  macrophage pseudotime corr: %d genes\n", nrow(macro_corr)))

# Panel membership flags from sig_wide
panel_metmac <- sig_wide[!is.na(signature_govaere2026_metmac_member)
                         & signature_govaere2026_metmac_member != "",
                         unique(human_symbol)]
panel_lam <- sig_wide[!is.na(signature_govaere2026_lam_member)
                      & signature_govaere2026_lam_member != "",
                      unique(human_symbol)]
panel_gpnmb <- sig_wide[!is.na(signature_govaere2026_gpnmb_macrophage_member)
                        & signature_govaere2026_gpnmb_macrophage_member != "",
                        unique(human_symbol)]
gpnmb_logfc <- sig_wide[, .(symbol = human_symbol,
                            gpnmb_logfc = as.numeric(signature_govaere2026_gpnmb_logfc))][
                          !is.na(gpnmb_logfc)]
cat(sprintf("  panel sizes: MetMac=%d, LAM=%d, GPNMB+=%d, GPNMB logFC=%d\n",
            length(panel_metmac), length(panel_lam),
            length(panel_gpnmb), nrow(gpnmb_logfc)))

macro_corr[, is_metmac := symbol %in% panel_metmac]
macro_corr[, is_lam    := symbol %in% panel_lam]
macro_corr[, is_gpnmb  := symbol %in% panel_gpnmb]

# Top-quartile late-pseudotime macrophage genes = ρ in top quartile (most positive).
q3 <- quantile(macro_corr$macro_rho, probs = 0.75, na.rm = TRUE)
macro_corr[, in_late_q4 := macro_rho >= q3]

# Fisher: MetMac enrichment in top quartile
a <- macro_corr[in_late_q4 == TRUE & is_metmac == TRUE, .N]
b <- macro_corr[in_late_q4 == TRUE & is_metmac == FALSE, .N]
cc <- macro_corr[in_late_q4 == FALSE & is_metmac == TRUE, .N]
d <- macro_corr[in_late_q4 == FALSE & is_metmac == FALSE, .N]
fisher_metmac <- safe_fisher(a, b, cc, d)
# Same for LAM and GPNMB+ panels.
a2 <- macro_corr[in_late_q4 == TRUE & is_lam == TRUE, .N]
b2 <- macro_corr[in_late_q4 == TRUE & is_lam == FALSE, .N]
c2 <- macro_corr[in_late_q4 == FALSE & is_lam == TRUE, .N]
d2 <- macro_corr[in_late_q4 == FALSE & is_lam == FALSE, .N]
fisher_lam <- safe_fisher(a2, b2, c2, d2)
a3 <- macro_corr[in_late_q4 == TRUE & is_gpnmb == TRUE, .N]
b3 <- macro_corr[in_late_q4 == TRUE & is_gpnmb == FALSE, .N]
c3 <- macro_corr[in_late_q4 == FALSE & is_gpnmb == TRUE, .N]
d3 <- macro_corr[in_late_q4 == FALSE & is_gpnmb == FALSE, .N]
fisher_gpnmb <- safe_fisher(a3, b3, c3, d3)

# Spearman: per-gene macrophage pseudotime ρ vs paper's gpnmb logFC.
spearman_join <- merge(macro_corr[, .(symbol, macro_rho)], gpnmb_logfc,
                       by = "symbol")
spear_res <- if (nrow(spearman_join) >= 5) {
  s <- suppressWarnings(cor.test(spearman_join$macro_rho,
                                 spearman_join$gpnmb_logfc,
                                 method = "spearman"))
  list(rho = unname(s$estimate), p = s$p.value, n = nrow(spearman_join))
} else list(rho = NA_real_, p = NA_real_, n = nrow(spearman_join))

# Output table: per-panel-gene late-pseudotime status
macro_out <- macro_corr[is_metmac | is_lam | is_gpnmb,
                        .(symbol, macro_rho, macro_padj, in_late_q4,
                          is_metmac, is_lam, is_gpnmb)]
setorder(macro_out, -macro_rho)
fwrite(macro_out,
       file.path(OUT_DIR, "macrophage_pseudotime_vs_metmac.csv"))

macro_summary <- c(
  sprintf("Macrophage pseudotime vs Govaere MetMac/LAM/GPNMB+ panels"),
  sprintf("Source: pseudotime_corr_Macrophages.csv (Spearman rho vs fibrosis)"),
  sprintf("Top-quartile late threshold (rho >= %.3f, n = %d genes)",
          q3, macro_corr[in_late_q4 == TRUE, .N]),
  "",
  sprintf("Panel overlap with top-quartile late-pseudotime genes:"),
  sprintf("  MetMac : %d / %d in late-q4, Fisher OR=%.2f, p=%s",
          a, a + cc, fisher_metmac$OR,
          if (is.na(fisher_metmac$p)) "NA"
          else format.pval(fisher_metmac$p, digits = 3, eps = 1e-300)),
  sprintf("  LAM    : %d / %d in late-q4, Fisher OR=%.2f, p=%s",
          a2, a2 + c2, fisher_lam$OR,
          if (is.na(fisher_lam$p)) "NA"
          else format.pval(fisher_lam$p, digits = 3, eps = 1e-300)),
  sprintf("  GPNMB+ : %d / %d in late-q4, Fisher OR=%.2f, p=%s",
          a3, a3 + c3, fisher_gpnmb$OR,
          if (is.na(fisher_gpnmb$p)) "NA"
          else format.pval(fisher_gpnmb$p, digits = 3, eps = 1e-300)),
  "",
  sprintf("Per-gene Spearman (our macrophage rho vs paper gpnmb_logfc):"),
  sprintf("  rho = %.3f, p = %s, n = %d",
          spear_res$rho,
          if (is.na(spear_res$p)) "NA"
          else format.pval(spear_res$p, digits = 3, eps = 1e-300),
          spear_res$n),
  "",
  "Interpretation:",
  paste0("  Hypothesis (memory pseudotime_pipeline.md): macrophage pseudotime",
         " tracks fibrosis F3->F4 rho=0.39; late-pseudotime macrophages",
         " should enrich for the paper's MetMac/GPNMB+ panel."),
  paste0("  OR > 1 + nominal p < 0.05 = consistent with hypothesis."),
  paste0("  OR ~ 1 or < 1 means our pseudotime axis is not aligned with",
         " paper's metabolic-macrophage axis (possible if our trajectory",
         " captures Kupffer->scar-associated rather than metabolic).")
)
writeLines(macro_summary,
           file.path(OUT_DIR, "macrophage_pseudotime_summary.txt"))
cat(sprintf("  -> MetMac OR=%.2f p=%.3g; LAM OR=%.2f p=%.3g; GPNMB+ OR=%.2f p=%.3g; Spearman rho=%.3f n=%d\n",
            fisher_metmac$OR, fisher_metmac$p,
            fisher_lam$OR, fisher_lam$p,
            fisher_gpnmb$OR, fisher_gpnmb$p,
            spear_res$rho, spear_res$n))

# ===========================================================================
# CHECK 3: IL32 axis panel atlas lookup
# ===========================================================================
cat("\n[76b] CHECK 3: IL32 axis atlas + Govaere column lookup\n")

# Atlas needs to be subsetted on human_symbol; only requested cols.
atlas_sub <- atlas[human_symbol %in% IL32_PANEL,
                   .(symbol = human_symbol, bulk_logFC, bulk_padj)]

sig_sub <- sig_wide[human_symbol %in% IL32_PANEL,
                    .(symbol = human_symbol,
                      geomx_sh_vs_ls_logfc =
                        as.numeric(signature_govaere2026_geomx_sh_vs_ls_logfc),
                      geomx_sh_vs_ls_padj  =
                        as.numeric(signature_govaere2026_geomx_sh_vs_ls_padj),
                      gpnmb_logfc =
                        as.numeric(signature_govaere2026_gpnmb_logfc),
                      gpnmb_padj  =
                        as.numeric(signature_govaere2026_gpnmb_padj),
                      il32_axis_member =
                        signature_govaere2026_il32_axis_member)]

il32_panel_dt <- data.table(symbol = IL32_PANEL)
il32_out <- merge(il32_panel_dt, atlas_sub, by = "symbol", all.x = TRUE)
il32_out <- merge(il32_out, hep_de, by = "symbol", all.x = TRUE)
il32_out <- merge(il32_out, sig_sub, by = "symbol", all.x = TRUE)
setcolorder(il32_out, c("symbol", "bulk_logFC", "bulk_padj",
                        "sc_hepatocyte_logFC", "sc_hepatocyte_padj",
                        "geomx_sh_vs_ls_logfc", "geomx_sh_vs_ls_padj",
                        "gpnmb_logfc", "gpnmb_padj", "il32_axis_member"))

fwrite(il32_out, file.path(OUT_DIR, "il32_axis_atlas_lookup.tsv"), sep = "\t")

# Also drop a Markdown rendering for human review.
md_header <- "| symbol | bulk_logFC | bulk_padj | sc_hep_logFC | sc_hep_padj | geomx_SH_vs_LS_logFC | geomx_padj | gpnmb_logFC | gpnmb_padj | il32_axis |"
md_sep    <- "|---|---|---|---|---|---|---|---|---|---|"
fmt <- function(x, dig = 3) {
  if (is.na(x)) return("NA")
  if (is.numeric(x) && abs(x) < 1e-3 && x != 0)
    return(formatC(x, format = "e", digits = 2))
  formatC(x, digits = dig, format = "f")
}
md_rows <- vapply(seq_len(nrow(il32_out)), function(i) {
  r <- il32_out[i]
  sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
          r$symbol,
          fmt(r$bulk_logFC), fmt(r$bulk_padj),
          fmt(r$sc_hepatocyte_logFC), fmt(r$sc_hepatocyte_padj),
          fmt(r$geomx_sh_vs_ls_logfc), fmt(r$geomx_sh_vs_ls_padj),
          fmt(r$gpnmb_logfc), fmt(r$gpnmb_padj),
          ifelse(is.na(r$il32_axis_member) | r$il32_axis_member == "",
                 "-", as.character(r$il32_axis_member)))
}, character(1))
writeLines(c("# IL32-axis panel cross-modal lookup", "",
             md_header, md_sep, md_rows),
           file.path(OUT_DIR, "il32_axis_atlas_lookup.md"))
cat(sprintf("  -> wrote %d-row IL32-axis lookup (TSV + MD)\n", nrow(il32_out)))

# ===========================================================================
# CHECK 4: bulk dream DEGs vs Sup11 (GPNMB+ MASH-vs-noMASH)
# ===========================================================================
cat("\n[76b] CHECK 4: bulk dream vs Sup Table 11 concordance\n")

dream <- fread(DREAM_FILE)
dream <- dream[!is.na(symbol) & symbol != ""]
dream_min <- dream[, .(symbol,
                       bulk_logFC = logFC,
                       bulk_padj  = padj)]
merged <- merge(dream_min, sup11, by.x = "symbol", by.y = "gene")
cat(sprintf("  bulk genes: %d ; sup11 genes: %d ; intersect: %d\n",
            nrow(dream_min), nrow(sup11), nrow(merged)))

# Spearman on logFC.
if (nrow(merged) >= 10) {
  s4 <- suppressWarnings(cor.test(merged$bulk_logFC, merged$log2FC,
                                  method = "spearman"))
  spear4 <- list(rho = unname(s4$estimate), p = s4$p.value)
} else {
  spear4 <- list(rho = NA_real_, p = NA_real_)
}
# Sign concordance: same sign across both.
sign_match <- merged[!is.na(bulk_logFC) & !is.na(log2FC)
                     & bulk_logFC != 0 & log2FC != 0,
                     sign(bulk_logFC) == sign(log2FC)]
sign_pct <- 100 * mean(sign_match, na.rm = TRUE)

# Jaccard at padj < 0.05 (both directions).
bulk_sig <- dream_min[!is.na(bulk_padj) & bulk_padj < 0.05, symbol]
sup11_sig <- sup11[!is.na(p_val_adj) & p_val_adj < 0.05, gene]
jaccard <- length(intersect(bulk_sig, sup11_sig)) /
  max(1L, length(union(bulk_sig, sup11_sig)))

# Sign concordance restricted to genes sig in both.
both_sig <- intersect(bulk_sig, sup11_sig)
both_dt <- merged[symbol %in% both_sig]
both_sign_pct <- if (nrow(both_dt) > 0) {
  100 * mean(sign(both_dt$bulk_logFC) == sign(both_dt$log2FC), na.rm = TRUE)
} else NA_real_

fwrite(merged[, .(symbol, bulk_logFC, bulk_padj,
                  sup11_log2FC = log2FC, sup11_padj = p_val_adj)],
       file.path(OUT_DIR, "bulk_dream_vs_sup11.csv"))

bulk_summary <- c(
  "Bulk dream DEGs vs Govaere Sup Table 11 (snRNA-seq GPNMB+ MASH-vs-noMASH)",
  "",
  sprintf("Dream genes              : %d", nrow(dream_min)),
  sprintf("Sup11 genes              : %d", nrow(sup11)),
  sprintf("Intersect                : %d", nrow(merged)),
  sprintf("Spearman rho (logFC)     : %.3f (p = %s)",
          spear4$rho,
          if (is.na(spear4$p)) "NA"
          else format.pval(spear4$p, digits = 3, eps = 1e-300)),
  sprintf("Sign concordance (all)   : %.1f%%", sign_pct),
  sprintf("Sign concordance (both sig padj<0.05): %s",
          if (is.na(both_sign_pct)) "NA"
          else sprintf("%.1f%% (n=%d)", both_sign_pct, nrow(both_dt))),
  sprintf("Jaccard padj<0.05        : %.3f (intersect=%d, union=%d)",
          jaccard,
          length(intersect(bulk_sig, sup11_sig)),
          length(union(bulk_sig, sup11_sig))),
  "",
  "Interpretation:",
  paste0("  Bulk dream is a bulk-tissue Disease-vs-Control contrast spanning",
         " all liver cell types; Sup11 is GPNMB+ MACROPHAGE MASH-vs-noMASH."),
  paste0("  Disagreement is expected; this row documents the cross-resolution",
         " baseline. Macrophage-specific GPNMB+ signatures dilute when",
         " measured at bulk-tissue resolution."),
  paste0("  Treat low rho / low Jaccard as confirmation that bulk and",
         " GPNMB+ macrophage axes are non-redundant (motivates cell-type",
         " stratification, not a contradiction).")
)
writeLines(bulk_summary,
           file.path(OUT_DIR, "bulk_dream_vs_sup11_summary.txt"))
cat(sprintf("  -> Spearman rho=%.3f p=%.3g; sign %.1f%%; Jaccard=%.3f\n",
            spear4$rho, spear4$p, sign_pct, jaccard))

# ===========================================================================
# README
# ===========================================================================
readme <- c(
  "# Govaere 2026 cross-modal validation (76b)",
  "",
  sprintf("Generated %s by `RNA-seq/76b_govaere2026_crossmodal.R`.", Sys.time()),
  "",
  "Cross-validates four of our pipeline outputs against Govaere et al 2026",
  "(Nat Genetics) signatures.",
  "",
  "## Files",
  "",
  "| file | purpose |",
  "|---|---|",
  "| `progressor_signature_vs_geomx_sh_ls.csv` | Per-gene Progressor (COL15A1/CDH11/LTBP2) vs GeoMx SH-vs-LS + Sup11 |",
  "| `progressor_signature_summary.txt` | Wilcoxon p, mean logFC, sign concordance, interpretation |",
  "| `macrophage_pseudotime_vs_metmac.csv` | Per-panel-gene late-pseudotime overlap |",
  "| `macrophage_pseudotime_summary.txt` | Fisher OR + Spearman rho per panel |",
  "| `il32_axis_atlas_lookup.tsv` / `.md` | Atlas + Govaere values for IL32, GPNMB, LPL, AKR1B10, FABP5, HLA-DRA |",
  "| `bulk_dream_vs_sup11.csv` | Per-gene dream vs Sup11 merge |",
  "| `bulk_dream_vs_sup11_summary.txt` | Spearman rho, sign concordance, Jaccard |",
  "",
  "## Status",
  "",
  sprintf("- Check 1 (Progressor vs GeoMx): GeoMx DE %s",
          if (geomx_available) "loaded from fresh DE file"
          else "NOT YET ON DISK -- using paper-reported logFC from signatures_wide. TODO rerun once Spatial W3.A1 deposits `geomx_de_sh_vs_ls.csv`."),
  "- Check 2 (Macrophage pseudotime): completed against MetMac/LAM/GPNMB+ panels.",
  "- Check 3 (IL32 axis lookup): completed.",
  "- Check 4 (bulk dream vs Sup11): completed; documents cross-resolution baseline (disagreement expected).",
  "",
  "## Reproducibility",
  "",
  "```bash",
  "micromamba activate rnaseq",
  "Rscript RNA-seq/76b_govaere2026_crossmodal.R",
  "```"
)
writeLines(readme, file.path(OUT_DIR, "README.md"))

cat("\n[76b] DONE\n")
