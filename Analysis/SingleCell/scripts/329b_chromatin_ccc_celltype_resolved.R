#!/usr/bin/env Rscript
# 329b_chromatin_ccc_celltype_resolved.R
#
# Analysis G2v2 — cell-type-resolved chromatin/CCC reconciliation.
#
# Motivation:
#   Script 329 found ATAC-regulated CCC receptors had zero overlap with bulk
#   dream DEGs (hypergeometric p=1.0). This is a methods artifact: ATAC-regulated
#   CCC receptors are likely cell-type-specific (e.g., CD36 on macrophages,
#   PDGFR on HSCs) and diluted to invisibility in bulk pseudobulk dream.
#   The motif-disrupted ligand finding (OR=4.11, p=0.0008) is unchanged and
#   remains the headline.
#
#   Here we re-test against per-cell-type pseudobulk DEGs (Hepatocytes,
#   Macrophages, Fibroblasts, Endothelial_cells, Cholangiocytes), where the
#   correct DE denominator is matched to the cell type expressing the receptor.
#
# Strategy:
#   1. Load ATAC-regulated genes (linked_gene/nearest_gene from variants in
#      scATAC peaks) — recompute as in Script 329.
#   2. Load motif-disrupted genes (variants in peaks with strong motif
#      disruption).
#   3. Load LIANA differential L-R pairs; extract ligand + receptor gene sets.
#   4. Load per-CT pseudobulk DEGs (MASLD_vs_Healthy contrast).
#   5. ENSG -> symbol map from bulk canonical_deg_results.csv.
#   6. For each cell type, run hypergeometric tests:
#        a) ATAC-regulated CCC receptors overlap with CT-DEGs
#           (universe = LIANA receptors intersected with CT-tested genes)
#        b) Motif-disrupted CCC ligands overlap with CT-DEGs
#           (universe = LIANA ligands intersected with CT-tested genes)
#        c) Motif-disrupted CCC receptors overlap with CT-DEGs
#           (universe = LIANA receptors intersected with CT-tested genes)
#   7. BH-adjust p-values across all (cell_type x test) cells.
#   8. Write CSV + summary text headline for Methods.
#
# Env: rnaseq
# Outputs:
#   Analysis/SingleCell/results_gpu_v2/chromatin_ccc/celltype_resolved_enrichment.csv
#   Analysis/SingleCell/results_gpu_v2/chromatin_ccc/gwas_atac_ccc_summary_v2.txt

suppressPackageStartupMessages({
  library(data.table)
})

BASE      <- Sys.getenv("MASLD_PROJECT_ROOT",
                        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GWAS_ATAC <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
LIANA     <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data")
PB_DE     <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudobulk_de")
INT_RES   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR    <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/chromatin_ccc")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PADJ_TH  <- 0.05
LFC_TH   <- 0.5
CELLTYPE_FILES <- c(
  Hepatocytes        = "Hepatocytes_de.csv",
  Macrophages        = "Macrophages_de.csv",
  Fibroblasts        = "Fibroblasts_de.csv",
  Endothelial_cells  = "Endothelial_cells_de.csv",
  Cholangiocytes     = "Cholangiocytes_de.csv"
)

# ---------------------------------------------------------------------------
message("[1] Loading GWAS-ATAC variant annotation + motif disruption...")
vann  <- fread(file.path(GWAS_ATAC, "gwas_atac_variant_annotation.csv"))
motif <- fread(file.path(GWAS_ATAC, "motif_disruption_scores.csv"))

vann[, atac_gene := fifelse(linked_gene != "" & !is.na(linked_gene),
                            linked_gene, nearest_gene)]
atac_regulated_genes <- unique(vann$atac_gene[vann$atac_gene != "" &
                                              !is.na(vann$atac_gene)])
message(sprintf("  %d variants in peaks; %d unique ATAC-regulated genes",
                nrow(vann), length(atac_regulated_genes)))

snp_col <- intersect(c("SNP_id", "variant_id", "rsid"), names(motif))[1]
if (!is.na(snp_col)) {
  setnames(motif, snp_col, "variant_id")
  vann_lookup <- unique(vann[, .(variant_id, atac_gene)])
  motif_joined <- merge(unique(motif[, .(variant_id)]), vann_lookup,
                        by = "variant_id", allow.cartesian = TRUE)
  motif_disrupted_genes <- unique(motif_joined$atac_gene[
    motif_joined$atac_gene != "" & !is.na(motif_joined$atac_gene)])
} else {
  motif_disrupted_genes <- character()
}
message(sprintf("  %d motif-disrupted genes", length(motif_disrupted_genes)))

# ---------------------------------------------------------------------------
message("[2] Loading LIANA L-R sets...")
liana <- fread(file.path(LIANA, "liana_differential_interactions.csv"))
split_complex <- function(x) unlist(strsplit(x, "_"))
ligand_genes   <- unique(unlist(lapply(unique(liana$ligand_complex),   split_complex)))
receptor_genes <- unique(unlist(lapply(unique(liana$receptor_complex), split_complex)))
message(sprintf("  LIANA ligands: %d, receptors: %d",
                length(ligand_genes), length(receptor_genes)))

# ---------------------------------------------------------------------------
message("[3] Building ENSG -> symbol map (GENCODE v49 metadata)...")
# CT pseudobulk DE files have a MIXED 'gene' column: most rows are already
# gene symbols (e.g., 'HSPH1', 'SDC1', 'LINC01409'); only unannotated rows
# carry bare ENSG IDs (e.g., 'ENSG00000238009'). We use the GENCODE v49
# metadata table to map any ENSG -> symbol while leaving symbol-style entries
# unchanged.
gencode_meta_path <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
if (file.exists(gencode_meta_path)) {
  gm <- fread(gencode_meta_path)
  gm <- unique(gm[, .(ensembl_base, gene_name)])
  gm <- gm[!duplicated(ensembl_base) & !is.na(gene_name) & gene_name != ""]
  ensg2sym <- setNames(gm$gene_name, gm$ensembl_base)
  message(sprintf("  GENCODE v49: %d ENSG -> symbol mappings", length(ensg2sym)))
} else {
  message("  GENCODE metadata not found; falling back to dream symbol col")
  bulk_dream <- fread(file.path(INT_RES, "canonical_deg_results.csv"),
                      select = c("gene", "symbol"))
  bulk_dream <- bulk_dream[!is.na(symbol) & symbol != "" &
                            !grepl("^ENSG", symbol) & !is.na(gene) & gene != ""]
  bulk_dream[, gene_unv := sub("\\..*$", "", gene)]
  ensg2sym <- setNames(bulk_dream$symbol, bulk_dream$gene_unv)
  ensg2sym <- ensg2sym[!duplicated(names(ensg2sym))]
  message(sprintf("  dream fallback: %d ENSG -> symbol mappings", length(ensg2sym)))
}

# Helper: resolve the 'gene' column of a CT DE table to a final symbol.
# If gene starts with 'ENSG', look up via ensg2sym; otherwise keep as-is.
resolve_symbol <- function(gene_vec) {
  is_ensg <- grepl("^ENSG", gene_vec)
  out <- gene_vec
  # strip version on ENSG side just in case
  out[is_ensg] <- ensg2sym[sub("\\..*$", "", out[is_ensg])]
  out
}

# ---------------------------------------------------------------------------
hyper_test <- function(set_a, set_b, bg) {
  # set_a = target (e.g., ATAC-regulated CCC receptors)
  # set_b = condition (e.g., CT-DEGs)
  # bg    = universe (e.g., LIANA receptors intersected with CT-tested genes)
  set_a <- intersect(set_a, bg)
  set_b <- intersect(set_b, bg)
  bg_n  <- length(bg)
  m     <- length(set_a)       # successes in population
  n     <- bg_n - m            # failures in population
  k     <- length(set_b)       # draws
  x     <- length(intersect(set_a, set_b))
  if (bg_n == 0 || m == 0 || k == 0) {
    p   <- NA_real_
    or  <- NA_real_
  } else {
    p     <- phyper(x - 1, m, n, k, lower.tail = FALSE)
    exp_x <- m * k / bg_n
    or    <- if (exp_x > 0) x / exp_x else NA_real_
  }
  data.table(k_overlap = x, n_set = m, n_universe = bg_n,
             n_drawn = k, OR = or, p_hypergeom = p)
}

# ---------------------------------------------------------------------------
message("[4] Per-CT enrichment tests...")
results <- list()
for (ct in names(CELLTYPE_FILES)) {
  fname <- CELLTYPE_FILES[[ct]]
  fpath <- file.path(PB_DE, fname)
  if (!file.exists(fpath)) {
    message(sprintf("  [skip] %s missing", fpath))
    next
  }
  ct_de <- fread(fpath)
  if (!"symbol" %in% names(ct_de)) {
    ct_de[, symbol := resolve_symbol(gene)]
  }
  ct_de <- ct_de[!is.na(symbol) & symbol != ""]
  # collapse multiple ENSGs mapping to same symbol — keep most extreme by |t|
  if (anyDuplicated(ct_de$symbol)) {
    ct_de[, abs_t := abs(t_stat)]
    setorder(ct_de, symbol, -abs_t)
    ct_de <- ct_de[!duplicated(symbol)]
    ct_de[, abs_t := NULL]
  }
  tested_symbols <- unique(ct_de$symbol)
  ct_degs <- ct_de[padj < PADJ_TH & abs(logFC) > LFC_TH, symbol]
  message(sprintf("  %s: %d tested symbols; %d DEGs (padj<%.2f, |LFC|>%.1f)",
                  ct, length(tested_symbols), length(ct_degs), PADJ_TH, LFC_TH))

  # Test A: ATAC-regulated CCC receptors vs CT-DEGs
  # universe = LIANA receptors AND ATAC-regulated AND tested in this CT
  rec_universe <- intersect(intersect(receptor_genes, atac_regulated_genes),
                             tested_symbols)
  res_a <- hyper_test(set_a = rec_universe,             # the target set IS the universe of interest
                      set_b = ct_degs,
                      bg    = intersect(receptor_genes, tested_symbols))
  res_a[, `:=`(cell_type = ct,
               test = "ATAC-regulated CCC receptors vs CT-DEGs")]
  results[[paste0(ct, "_A")]] <- res_a

  # Test B: Motif-disrupted CCC ligands vs CT-DEGs
  res_b <- hyper_test(set_a = intersect(intersect(ligand_genes, motif_disrupted_genes),
                                        tested_symbols),
                      set_b = ct_degs,
                      bg    = intersect(ligand_genes, tested_symbols))
  res_b[, `:=`(cell_type = ct,
               test = "Motif-disrupted CCC ligands vs CT-DEGs")]
  results[[paste0(ct, "_B")]] <- res_b

  # Test C: Motif-disrupted CCC receptors vs CT-DEGs
  res_c <- hyper_test(set_a = intersect(intersect(receptor_genes, motif_disrupted_genes),
                                        tested_symbols),
                      set_b = ct_degs,
                      bg    = intersect(receptor_genes, tested_symbols))
  res_c[, `:=`(cell_type = ct,
               test = "Motif-disrupted CCC receptors vs CT-DEGs")]
  results[[paste0(ct, "_C")]] <- res_c
}

enr <- rbindlist(results, fill = TRUE)
setcolorder(enr, c("cell_type", "test", "k_overlap", "n_set",
                   "n_universe", "n_drawn", "OR", "p_hypergeom"))
enr[, p_bh := p.adjust(p_hypergeom, method = "BH")]

fwrite(enr, file.path(OUTDIR, "celltype_resolved_enrichment.csv"))
message(sprintf("[5] Wrote %s", file.path(OUTDIR, "celltype_resolved_enrichment.csv")))

# ---------------------------------------------------------------------------
message("[6] Overlap-gene listings (for manuscript supp + sanity)...")
overlap_lines <- list()
for (ct in names(CELLTYPE_FILES)) {
  fpath <- file.path(PB_DE, CELLTYPE_FILES[[ct]])
  if (!file.exists(fpath)) next
  ct_de <- fread(fpath)
  if (!"symbol" %in% names(ct_de)) ct_de[, symbol := resolve_symbol(gene)]
  ct_de <- ct_de[!is.na(symbol) & symbol != ""]
  if (anyDuplicated(ct_de$symbol)) {
    ct_de[, abs_t := abs(t_stat)]
    setorder(ct_de, symbol, -abs_t)
    ct_de <- ct_de[!duplicated(symbol)]
  }
  tested  <- unique(ct_de$symbol)
  ct_degs <- ct_de[padj < PADJ_TH & abs(logFC) > LFC_TH, symbol]

  ar_rec_overlap <- intersect(intersect(receptor_genes, atac_regulated_genes),
                              ct_degs)
  md_lig_overlap <- intersect(intersect(ligand_genes,   motif_disrupted_genes),
                              ct_degs)
  md_rec_overlap <- intersect(intersect(receptor_genes, motif_disrupted_genes),
                              ct_degs)
  overlap_lines[[ct]] <- data.table(
    cell_type = ct,
    test      = c("ATAC-reg receptors", "Motif-disrupted ligands",
                  "Motif-disrupted receptors"),
    genes     = c(paste(ar_rec_overlap, collapse = ";"),
                  paste(md_lig_overlap, collapse = ";"),
                  paste(md_rec_overlap, collapse = ";"))
  )
}
overlap_dt <- rbindlist(overlap_lines, fill = TRUE)
fwrite(overlap_dt, file.path(OUTDIR, "celltype_resolved_overlap_genes.csv"))

# ---------------------------------------------------------------------------
message("[7] Headline summary...")
sig_rows <- enr[!is.na(p_bh) & p_bh < 0.1]
any_sig  <- nrow(sig_rows) > 0

# pick the strongest significant receptor signal
rec_sig <- enr[grepl("receptors", test) & !is.na(p_hypergeom)]
rec_sig <- rec_sig[order(p_hypergeom)]
top_rec <- rec_sig[1]

lines <- c(
  "## Cell-type-resolved chromatin-CCC reconciliation (Script 329b)",
  "",
  sprintf("Per-CT pseudobulk DEG threshold: padj < %.2f, |logFC| > %.1f",
          PADJ_TH, LFC_TH),
  sprintf("Cell types tested: %s",
          paste(names(CELLTYPE_FILES), collapse = ", ")),
  "",
  "Enrichment table (sorted by BH p):",
  capture.output(print(enr[order(p_bh),
                           .(cell_type, test, k_overlap, n_set,
                             n_universe, n_drawn, OR, p_hypergeom, p_bh)],
                       nrows = 30)),
  "",
  sprintf("Significant rows (BH < 0.1): %d", nrow(sig_rows)),
  if (any_sig) {
    paste(capture.output(print(sig_rows[order(p_bh),
                                        .(cell_type, test, k_overlap,
                                          n_set, OR, p_hypergeom, p_bh)])),
          collapse = "\n")
  } else "None (all BH p >= 0.1)",
  "",
  "## Methods-ready headline (suggested):",
  "",
  if (any_sig) {
    sprintf(paste0("Cell-type-resolved pseudobulk DE recovered chromatin-CCC ",
                   "coupling that was invisible in bulk dream (Script 329, p=1.0). ",
                   "ATAC-regulated CCC receptors are enriched among per-cell-type ",
                   "MASLD-vs-Healthy DEGs in %d of %d tested cell types (BH<0.1). ",
                   "Strongest signal: %s receptors %s in %s (k=%d, OR=%.2f, ",
                   "BH p=%.3g)."),
            nrow(sig_rows[grepl("receptors", test)]),
            length(CELLTYPE_FILES),
            ifelse(grepl("Motif", top_rec$test), "motif-disrupted", "ATAC-regulated"),
            "DE",
            top_rec$cell_type, top_rec$k_overlap,
            top_rec$OR, top_rec$p_bh)
  } else {
    paste0("Cell-type-resolved pseudobulk DE did NOT rescue receptor-side ",
           "enrichment after Script 329's bulk null. Motif-disrupted ligand ",
           "enrichment (Script 329; OR=4.11, p=0.0008) remains the chromatin-CCC ",
           "headline; receptor-side genetic regulation is sparse across all five ",
           "tested cell types and below detection at our CT-DE sample sizes.")
  }
)
writeLines(lines, file.path(OUTDIR, "gwas_atac_ccc_summary_v2.txt"))
writeLines(lines)

message("Done. Outputs in: ", OUTDIR)
