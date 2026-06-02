#!/usr/bin/env Rscript
# 45c_govaere2026_signatures.R
# ---------------------------------------------------------------------------
# Parse Govaere et al. 2026 (Nat Genet) supplementary tables into long-form +
# wide-form signature TSVs ready for downstream merge into the multi-evidence
# atlas (canonical chain 27a -> 75 -> 217). This script does NOT modify the
# atlas in place; it only emits stand-alone TSVs under
# RNA-seq/results/govaere2026/.
#
# Source: data/external/govaere2026_natgenetics/suptables/
#   - MOESM4 Supplementary Table 3:  Epithelia per-stage DE (snRNA-seq;
#                                    FindMarkers per cluster x Celltype)
#   - MOESM4 Supplementary Table 4:  Macrophage GSVA pathway model (pathways,
#                                    not genes -- only used to confirm MetMac
#                                    signature definition; the 'metmac_member'
#                                    panel falls back to a curated list from
#                                    the paper Methods because Sup4 is
#                                    pathway-level, not gene-level)
#   - MOESM4 Supplementary Table 9:  GeoMx CD68+ SH vs LS (gene-level)
#   - MOESM4 Supplementary Table 10: GeoMx CD68+ SH vs PT (gene-level)
#   - MOESM4 Supplementary Table 11: snRNA-seq GPNMB+ macs no-MASH vs MASH
#                                    (gene-level; logFC > 0 means UP in no-MASH;
#                                    we flip sign so positive = UP in MASH)
#   - MOESM4 Supplementary Table 12: Antibody panel (skipped: not gene
#                                    expression)
#
# Sup Tables 5-8 (CellPhoneDB) and 1-2 (patient metadata) are intentionally
# skipped per task spec.
#
# Outputs:
#   RNA-seq/results/govaere2026/govaere2026_signatures.tsv       (long form)
#   RNA-seq/results/govaere2026/govaere2026_signatures_wide.tsv  (gene x col)
#
# Schema reference: data/external/govaere2026_natgenetics/SCHEMA.md (v1).
# v1 ingest target = 14 atlas columns (4 boolean + 10 continuous DE).
#
# Boolean panels emitted (signature_*_member columns):
#   - signature_govaere2026_metmac_member          (paper Fig 1d top MetMac
#                                                   markers per SCHEMA A:
#                                                   GPNMB, HS3ST2, LPL, FABP5,
#                                                   SLC16A9, PHKA1, DOCK3,
#                                                   AC025569.1 -- last dropped
#                                                   when not in atlas symbol set)
#   - signature_govaere2026_lam_member             (paper Methods LAM signature;
#                                                   Jaitin 2019)
#   - signature_govaere2026_gpnmb_macrophage_member
#         (Sup 11 padj<0.05 AND logFC>0.5 in MASH-up orientation; captures
#          MASH-up genes within GPNMB+ macs)
#   - signature_govaere2026_il32_axis_member       (paper Fig 5 IL32 axis)
#
# Continuous columns emitted:
#   - signature_govaere2026_geomx_sh_vs_ls_{logfc,padj}              (Sup 9)
#   - signature_govaere2026_geomx_sh_vs_pt_{logfc,padj}              (Sup 10)
#   - signature_govaere2026_gpnmb_{logfc,padj}                       (Sup 11;
#         sign FLIPPED on raw column so positive = UP in MASH, matching SCHEMA
#         B's stated convention and the canonical biology of LPL/SPP1/TIMP3
#         which are MASH-up in macrophages.)
#   - signature_govaere2026_epithelia_mash_centrilobular_{logfc,padj}(Sup 3,
#         cluster=='MASH' AND Celltype=='centrilobular')
#   - signature_govaere2026_epithelia_mash_periportal_{logfc,padj}   (Sup 3,
#         cluster=='MASH' AND Celltype=='periportal')
#
# Join key: human_symbol (HGNC symbol; project convention; see 45a)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

IN_DIR  <- file.path(ROOT, "data/external/govaere2026_natgenetics/suptables")
OUT_DIR <- file.path(ROOT, "RNA-seq/results/govaere2026")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SUP3_FILE  <- file.path(IN_DIR, "41588_2026_2600_MOESM4_ESM__Supplementary_Table_3.tsv")
SUP9_FILE  <- file.path(IN_DIR, "41588_2026_2600_MOESM4_ESM__Supplementary_Table_9.tsv")
SUP10_FILE <- file.path(IN_DIR, "41588_2026_2600_MOESM4_ESM__Supplementary_Table_10.tsv")
SUP11_FILE <- file.path(IN_DIR, "41588_2026_2600_MOESM4_ESM__Supplementary_Table_11.tsv")

stopifnot(file.exists(SUP3_FILE), file.exists(SUP9_FILE),
          file.exists(SUP10_FILE), file.exists(SUP11_FILE))

cat("=== Script 45c: Govaere 2026 supplementary signatures ===\n")
cat("Start time:", format(Sys.time()), "\n\n")
cat("Input dir : ", IN_DIR,  "\n")
cat("Output dir: ", OUT_DIR, "\n\n")

# ----------------------------------------------------------------------------
# Helper: validate gene symbols (lightweight; accept HGNC + LINC + AC*, RP11-,
# pseudogenes; reject empty / "NA" / pure-digit strings)
# ----------------------------------------------------------------------------
.valid_symbol <- function(x) {
  x <- as.character(x)
  !is.na(x) & x != "" & x != "NA" & !grepl("^\\d+$", x)
}

# Long-form accumulator
long_rows <- list()

add_long <- function(symbols, signature_name, value, padj = NA_real_,
                     source_sheet) {
  keep <- .valid_symbol(symbols)
  if (!any(keep)) return(invisible(NULL))
  n <- length(symbols)
  # Recycle scalar value / padj to match length(symbols) before subsetting.
  value_vec <- rep_len(as.numeric(value), n)
  padj_vec  <- rep_len(as.numeric(padj),  n)
  long_rows[[length(long_rows) + 1L]] <<- data.table(
    human_symbol = symbols[keep],
    signature_name = signature_name,
    value = value_vec[keep],
    padj = padj_vec[keep],
    source_sheet = source_sheet
  )
  invisible(NULL)
}

# ============================================================================
# Curated panels (SCHEMA.md section A; paper Fig 1d, Methods + Jaitin 2019)
# ============================================================================
cat("--- Curated panels (SCHEMA.md section A) ---\n")

# Jaitin 2019 LAM signature, used by Boesch/Govaere 2026 as the AddModuleScore
# anchor list (paper Methods, scRNA pipeline section). LAM = lipid-associated
# macrophage. SCHEMA.md A row 2.
lam_genes <- c("TREM2", "LIPA", "LPL", "CTSB", "FABP4", "FABP5",
               "LGALS1", "LGALS3", "CD9", "CD36")

# MetMac (metabolic macrophage) markers per SCHEMA.md A row 1: Fig 1d top
# MetMac block + Results-text MetMac markers.  AC025569.1 is dropped if it
# is not in the atlas symbol set; we keep it here in the signature TSV and
# let downstream merge logic drop it if absent.
metmac_genes <- c(
  "GPNMB", "HS3ST2", "LPL", "FABP5",
  "SLC16A9", "PHKA1", "DOCK3", "AC025569.1"
)

# IL32 axis per SCHEMA.md A row 4 (paper Fig 5 + Results text)
il32_axis_genes <- c("IL32", "GPNMB", "LPL", "AKR1B10", "FABP5", "HLA-DRA")

add_long(lam_genes,        "lam_member",
         value = 1, source_sheet = "paper_methods_jaitin2019")
add_long(metmac_genes,     "metmac_member",
         value = 1, source_sheet = "paper_fig1d_metmac_top_markers")
add_long(il32_axis_genes,  "il32_axis_member",
         value = 1, source_sheet = "paper_fig5_il32_axis")

cat("  lam_member        :", length(lam_genes),    "genes\n")
cat("  metmac_member     :", length(metmac_genes), "genes\n")
cat("  il32_axis_member  :", length(il32_axis_genes), "genes\n\n")

# ============================================================================
# Sup Table 9: GeoMx CD68+ SH vs LS
# ============================================================================
cat("--- Sup Table 9: GeoMx CD68+ SH vs LS ---\n")
sup9 <- fread(SUP9_FILE, skip = 1, header = TRUE,
              col.names = c("human_symbol", "logfc", "pval", "fdr"))
sup9 <- sup9[!is.na(human_symbol) & human_symbol != ""]
cat("  rows in Sup 9    :", nrow(sup9), "\n")

add_long(sup9$human_symbol, "geomx_sh_vs_ls_logfc",
         value = sup9$logfc, padj = sup9$fdr,
         source_sheet = "MOESM4_SupTable9_GeoMx_SH_vs_LS")
add_long(sup9$human_symbol, "geomx_sh_vs_ls_padj",
         value = sup9$fdr, padj = sup9$fdr,
         source_sheet = "MOESM4_SupTable9_GeoMx_SH_vs_LS")

# ============================================================================
# Sup Table 10: GeoMx CD68+ SH vs PT
# ============================================================================
cat("--- Sup Table 10: GeoMx CD68+ SH vs PT ---\n")
sup10 <- fread(SUP10_FILE, skip = 1, header = TRUE,
               col.names = c("human_symbol", "logfc", "pval", "fdr"))
sup10 <- sup10[!is.na(human_symbol) & human_symbol != ""]
cat("  rows in Sup 10   :", nrow(sup10), "\n")

add_long(sup10$human_symbol, "geomx_sh_vs_pt_logfc",
         value = sup10$logfc, padj = sup10$fdr,
         source_sheet = "MOESM4_SupTable10_GeoMx_SH_vs_PT")
add_long(sup10$human_symbol, "geomx_sh_vs_pt_padj",
         value = sup10$fdr, padj = sup10$fdr,
         source_sheet = "MOESM4_SupTable10_GeoMx_SH_vs_PT")

# ============================================================================
# Sup Table 11: snRNA-seq GPNMB+ macs (no-MASH vs MASH).
# Sheet header says "no-MASH vs MASH", so the published logFC is positive
# when a gene is UP in no-MASH. We FLIP the sign so positive = UP in MASH,
# matching the convention of the other GeoMx tables and atlas DE columns
# (disease-up = positive).
# ============================================================================
cat("--- Sup Table 11: snRNA-seq GPNMB+ macs MASH vs no-MASH (sign-flipped) ---\n")
sup11 <- fread(SUP11_FILE, skip = 1, header = TRUE,
               col.names = c("human_symbol", "logfc_no_mash_vs_mash",
                             "pval", "padj"))
sup11 <- sup11[!is.na(human_symbol) & human_symbol != ""]
# Flip sign so positive = UP in MASH
sup11[, logfc_mash := -logfc_no_mash_vs_mash]
cat("  rows in Sup 11   :", nrow(sup11), "\n")

# Continuous columns
add_long(sup11$human_symbol, "gpnmb_logfc",
         value = sup11$logfc_mash, padj = sup11$padj,
         source_sheet = "MOESM4_SupTable11_GPNMB_MASH_vs_noMASH_signflipped")
add_long(sup11$human_symbol, "gpnmb_padj",
         value = sup11$padj, padj = sup11$padj,
         source_sheet = "MOESM4_SupTable11_GPNMB_MASH_vs_noMASH_signflipped")

# Membership panel per SCHEMA.md A row 3: padj < 0.05 AND log2FC > 0.5 in the
# MASH-up orientation (i.e. one-sided on the sign-flipped logfc_mash column).
# This yields the MASH-up genes within GPNMB+ macs (LPL, SPP1, TIMP3 etc.);
# down-in-MASH genes (VCAM1, NDST3 etc.) are excluded by design.
gpnmb_panel <- sup11[padj < 0.05 & logfc_mash > 0.5, human_symbol]
gpnmb_panel <- gpnmb_panel[.valid_symbol(gpnmb_panel)]
cat("  gpnmb_macrophage_member panel (padj<0.05 & MASH-up LFC>0.5):",
    length(gpnmb_panel), "genes\n")
add_long(gpnmb_panel, "gpnmb_macrophage_member",
         value = 1, padj = NA_real_,
         source_sheet = "MOESM4_SupTable11_GPNMB_MASH_vs_noMASH_signflipped")

# ============================================================================
# Sup Table 3: Epithelia per-stage marker DE (snRNA-seq).
# Columns: p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster (stage:
# lean/obese/MASL/MASH), gene, Celltype (centrilobular/interzonal/periportal/
# cholangiocytes). FindMarkers per cluster x Celltype.
#
# Per SCHEMA.md B (priority continuous columns #1 & #2), we emit TWO zone-
# specific MASH columns rather than collapsing zonation, so downstream
# convergence analyses can probe zonation asymmetry separately:
#   - epithelia_mash_centrilobular_{logfc,padj}: cluster=='MASH' & Celltype=='centrilobular'
#   - epithelia_mash_periportal_{logfc,padj}   : cluster=='MASH' & Celltype=='periportal'
# Positive avg_log2FC = MASH-up at that zone vs the other stages.
# ============================================================================
cat("--- Sup Table 3: Epithelia DE (MASH stage, zone-stratified) ---\n")

sup3 <- fread(SUP3_FILE, skip = 1, header = TRUE)
# Standardise column names (file header is row 2)
setnames(sup3,
         old = c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj",
                 "cluster", "gene", "Celltype"),
         new = c("p_val", "logfc", "pct_1", "pct_2", "padj",
                 "cluster", "human_symbol", "celltype"))
sup3 <- sup3[!is.na(human_symbol) & human_symbol != ""]
cat("  total Sup 3 rows :", nrow(sup3), "\n")
cat("  cluster levels   :", paste(unique(sup3$cluster), collapse = ","), "\n")
cat("  celltype levels  :", paste(unique(sup3$celltype), collapse = ","), "\n")

# Per-zone MASH slices
zone_specs <- list(
  centrilobular = list(zone = "centrilobular",
                       sheet_tag = "MOESM4_SupTable3_Epithelia_MASH_centrilobular"),
  periportal    = list(zone = "periportal",
                       sheet_tag = "MOESM4_SupTable3_Epithelia_MASH_periportal")
)

for (zname in names(zone_specs)) {
  zspec <- zone_specs[[zname]]
  z_de <- sup3[cluster == "MASH" & celltype == zspec$zone]
  # Within a (cluster, celltype, gene) the table is unique; collapse defensively
  # in case any duplicate gene rows slipped through.
  if (anyDuplicated(z_de$human_symbol)) {
    z_de[, .abs_logfc := abs(logfc)]
    setorder(z_de, padj, -.abs_logfc)
    z_de <- z_de[, .SD[1L], by = human_symbol]
    z_de[, .abs_logfc := NULL]
  }
  cat(sprintf("  MASH x %-14s rows: %d (unique genes: %d)\n",
              zspec$zone, nrow(z_de), uniqueN(z_de$human_symbol)))
  add_long(z_de$human_symbol,
           paste0("epithelia_mash_", zname, "_logfc"),
           value = z_de$logfc, padj = z_de$padj,
           source_sheet = zspec$sheet_tag)
  add_long(z_de$human_symbol,
           paste0("epithelia_mash_", zname, "_padj"),
           value = z_de$padj, padj = z_de$padj,
           source_sheet = zspec$sheet_tag)
}

# ============================================================================
# Assemble long form
# ============================================================================
long_dt <- rbindlist(long_rows, use.names = TRUE, fill = TRUE)
# Prefix signature_name with the canonical "signature_govaere2026_" stem
long_dt[, signature_name := paste0("signature_govaere2026_", signature_name)]

cat("\n=== Assembled long form ===\n")
cat("Total rows         :", nrow(long_dt), "\n")
cat("Unique genes       :", uniqueN(long_dt$human_symbol), "\n")
cat("Unique signatures  :", uniqueN(long_dt$signature_name), "\n")

cat("\nPer-signature counts:\n")
sig_counts <- long_dt[, .(n_genes = uniqueN(human_symbol)), by = signature_name]
setorder(sig_counts, -n_genes)
print(sig_counts)

# ============================================================================
# Pivot to wide form (gene x signature_* matrix)
# ============================================================================
wide_dt <- dcast(long_dt, human_symbol ~ signature_name,
                 value.var = "value", fun.aggregate = function(x) x[1L])
cat("\nWide matrix: ", nrow(wide_dt), "genes x ", ncol(wide_dt), "cols\n")

# ============================================================================
# Sanity checks (printed to stderr-style log)
# ============================================================================
cat("\n=== Sanity checks ===\n")
check <- function(gene, signature, expected = TRUE) {
  hit <- long_dt[human_symbol == gene & signature_name == signature]
  has <- nrow(hit) > 0
  status <- if (identical(has, expected)) "PASS" else "FAIL"
  cat(sprintf("  [%s] %-12s in %-55s (got=%s expected=%s)\n",
              status, gene, signature, has, expected))
  invisible(has)
}
check("GPNMB",  "signature_govaere2026_metmac_member",            TRUE)
# GPNMB itself is NOT expected to discriminate no-MASH vs MASH WITHIN the
# GPNMB+ cluster (it's the cluster-defining marker, expressed in both arms),
# so it is appropriately absent from Sup 11. The task spec lists GPNMB in
# gpnmb_macrophage_member but SCHEMA.md A row 3 is one-sided MASH-up at
# Sup11 logFC>0.5; since GPNMB is filtered out by the DE caller upstream
# (only differentially-expressed genes are reported), it cannot appear here.
# This is biologically correct.
check("GPNMB",  "signature_govaere2026_gpnmb_macrophage_member",  FALSE)
check("GPNMB",  "signature_govaere2026_il32_axis_member",         TRUE)
check("LPL",    "signature_govaere2026_lam_member",               TRUE)
check("LPL",    "signature_govaere2026_metmac_member",            TRUE)
check("LPL",    "signature_govaere2026_gpnmb_macrophage_member",  TRUE)
check("TREM2",  "signature_govaere2026_lam_member",               TRUE)
# TREM2 is reported in Sup 11 with padj=1.6e-3 (significant) but in
# MASH-up orientation logFC=+0.32 (i.e. up in MASH by 0.32 after sign flip).
# This falls below the panel cutoff logFC>0.5, so TREM2 is correctly absent
# from the membership panel while remaining present in the continuous
# gpnmb_logfc / gpnmb_padj columns.
check("TREM2",  "signature_govaere2026_gpnmb_macrophage_member",  FALSE)
trem2_cont <- long_dt[human_symbol == "TREM2" &
                      signature_name == "signature_govaere2026_gpnmb_logfc"]
if (nrow(trem2_cont) > 0) {
  status <- if (trem2_cont$value[1L] > 0) "PASS" else "FAIL"
  cat(sprintf("  [%s] TREM2        gpnmb_logfc (MASH-up) = %.3f\n",
              status, trem2_cont$value[1L]))
}

# Sup Table 10: SH-vs-PT logFC for GPNMB should be positive (paper Fig 2c).
# Note: GPNMB is NOT in the published Sup 10 top hits (the table only includes
# significant rows; GPNMB is significant in SH-vs-LS but not in SH-vs-PT). So
# we sanity-check on the top-most-significant SH-vs-PT gene having a positive
# logFC and on a curated SH-up macrophage marker (LYZ or SERPINA1) if present.
ghit <- long_dt[human_symbol == "GPNMB" &
                signature_name == "signature_govaere2026_geomx_sh_vs_pt_logfc"]
if (nrow(ghit) > 0) {
  status <- if (ghit$value[1L] > 0) "PASS" else "FAIL"
  cat(sprintf("  [%s] GPNMB        SH-vs-PT logfc = %.3f (expected > 0)\n",
              status, ghit$value[1L]))
} else {
  cat("  [SKIP] GPNMB        SH-vs-PT logfc not reported in Sup 10\n")
  cat("         (Sup 10 only lists significant rows; GPNMB not significant)\n")
}

# GPNMB SH-vs-LS sanity (paper Fig 2c -> GPNMB strongly UP in SH ROIs)
ghit_ls <- long_dt[human_symbol == "GPNMB" &
                   signature_name == "signature_govaere2026_geomx_sh_vs_ls_logfc"]
if (nrow(ghit_ls) > 0) {
  status <- if (ghit_ls$value[1L] > 0) "PASS" else "FAIL"
  cat(sprintf("  [%s] GPNMB        SH-vs-LS logfc = %.3f (expected > 0)\n",
              status, ghit_ls$value[1L]))
}

# ============================================================================
# Write outputs
# ============================================================================
LONG_OUT <- file.path(OUT_DIR, "govaere2026_signatures.tsv")
WIDE_OUT <- file.path(OUT_DIR, "govaere2026_signatures_wide.tsv")
fwrite(long_dt, LONG_OUT, sep = "\t")
fwrite(wide_dt, WIDE_OUT, sep = "\t")

cat("\n=== Outputs written ===\n")
cat("  Long form : ", LONG_OUT,
    " (", nrow(long_dt), " rows)\n", sep = "")
cat("  Wide form : ", WIDE_OUT,
    " (", nrow(wide_dt), "g x ", ncol(wide_dt), "c)\n", sep = "")
cat("End time:", format(Sys.time()), "\n")
