#!/usr/bin/env Rscript
# 53_monogenic_convergence.R
# Annotate the multi-evidence atlas with the Okur et al. 2026 monogenic
# metabolic disease panel and compute common-variant × rare-variant convergence.
#
# Inputs:
#   - data/external/okur2026_monogenic/okur2026_panel_genes_partial.tsv
#       (30 of ~90 panel genes from PDF text + Table 3; full list pending Zenodo
#        embargo lift on 2026-05-01)
#   - data/external/okur2026_monogenic/okur2026_pl_variants_pdf.tsv
#       (31 P/LP variants in 31 participants, manually transcribed from PDF Table 3)
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
#
# Outputs:
#   - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv  (in-place; adds 3 cols)
#   - RNA-seq/results/monogenic_convergence/monogenic_x_coloc.csv  (per-gene overlay)
#   - RNA-seq/results/monogenic_convergence/monogenic_summary.txt   (overlap counts + ORs)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

PROJ <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

ATLAS_PATH    <- file.path(PROJ, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
PANEL_PATH    <- file.path(PROJ, "data/external/okur2026_monogenic/okur2026_panel_genes_partial.tsv")
VARIANTS_PATH <- file.path(PROJ, "data/external/okur2026_monogenic/okur2026_pl_variants_pdf.tsv")
OUT_DIR       <- file.path(PROJ, "RNA-seq/results/monogenic_convergence")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading inputs...\n")
panel    <- fread(PANEL_PATH, fill = TRUE)
variants <- fread(VARIANTS_PATH, fill = TRUE)
atlas    <- fread(ATLAS_PATH)

cat("  Panel genes (partial):    ", nrow(panel), "\n")
cat("  P/LP variants in cohort:  ", nrow(variants), " (in ", length(unique(variants$study_id)), " participants)\n", sep="")
cat("  Atlas dim:                ", nrow(atlas), "x", ncol(atlas), "\n")

# ----- Build per-gene annotation -----
genes_with_pl     <- unique(variants$gene)
genes_in_panel    <- unique(panel$gene)
genes_with_mafld  <- unique(variants[grepl("MAFLD|MASH", indication, ignore.case=TRUE), gene])

# Indication summary per gene (concatenate unique indications)
gene_ind <- variants[, .(monogenic_indication_okur2026 = paste(unique(indication), collapse=";")),
                     by = gene]

# Indication category from panel (HLD/HTG, T2DM, Lipodystrophy, Obesity, Syndromic_Obesity, Lipodystrophy+MAFLD, ...)
gene_cat <- panel[, .(monogenic_category_okur2026 = paste(unique(indication_category), collapse=";")),
                  by = gene]

# Merge per-gene info
mono <- merge(
  data.table(gene = unique(c(genes_in_panel, genes_with_pl))),
  gene_cat, by = "gene", all.x = TRUE)
mono <- merge(mono, gene_ind, by = "gene", all.x = TRUE)
mono[, monogenic_panel_okur2026  := gene %in% genes_in_panel]
mono[, monogenic_pl_okur2026     := gene %in% genes_with_pl]
mono[, monogenic_mafld_okur2026  := gene %in% genes_with_mafld]

cat("\n  Genes in panel (partial):       ", sum(mono$monogenic_panel_okur2026), "\n")
cat("  Genes with P/LP in cohort:      ", sum(mono$monogenic_pl_okur2026), "\n")
cat("  Genes with MAFLD/MASH P/LP:     ", sum(mono$monogenic_mafld_okur2026), "\n")

# ----- Annotate atlas (in-place) -----
# Atlas key: human_symbol
cat("\nAnnotating atlas...\n")
# Drop existing monogenic cols if rerunning
existing_mono_cols <- grep("^monogenic_.*_okur2026$", colnames(atlas), value=TRUE)
if (length(existing_mono_cols) > 0) {
  atlas[, (existing_mono_cols) := NULL]
  cat("  Dropped pre-existing monogenic cols:", length(existing_mono_cols), "\n")
}

atlas <- merge(atlas, mono, by.x = "human_symbol", by.y = "gene", all.x = TRUE)
# NA → FALSE for boolean cols
atlas[is.na(monogenic_panel_okur2026), monogenic_panel_okur2026 := FALSE]
atlas[is.na(monogenic_pl_okur2026),    monogenic_pl_okur2026    := FALSE]
atlas[is.na(monogenic_mafld_okur2026), monogenic_mafld_okur2026 := FALSE]

cat("  Atlas now:", nrow(atlas), "x", ncol(atlas), "\n")
cat("  Genes flagged monogenic_panel:    ", sum(atlas$monogenic_panel_okur2026), "\n")
cat("  Genes flagged monogenic_pl:       ", sum(atlas$monogenic_pl_okur2026), "\n")
cat("  Genes flagged monogenic_mafld:    ", sum(atlas$monogenic_mafld_okur2026), "\n")

# ----- Per-gene convergence overlay -----
cat("\nComputing common × rare convergence...\n")
top_susie <- atlas[!is.na(coloc_susie_best_pp4) & coloc_susie_best_pp4 > 0.5]
top_abf   <- atlas[!is.na(coloc_abf_best_pp4)   & coloc_abf_best_pp4   > 0.5]

cat("  SuSiE PP4 > 0.5 genes:           ", nrow(top_susie), "\n")
cat("  ABF PP4 > 0.5 genes:             ", nrow(top_abf), "\n")

overlap_panel_susie <- intersect(top_susie$human_symbol, genes_in_panel)
overlap_panel_abf   <- intersect(top_abf$human_symbol,   genes_in_panel)
overlap_pl_susie    <- intersect(top_susie$human_symbol, genes_with_pl)

cat("\n  Panel × SuSiE PP4>0.5 overlap:    ", length(overlap_panel_susie), "\n")
cat("    genes: ", paste(overlap_panel_susie, collapse=", "), "\n", sep="")
cat("  Panel × ABF PP4>0.5 overlap:      ", length(overlap_panel_abf), "\n")
cat("    genes: ", paste(overlap_panel_abf, collapse=", "), "\n", sep="")
cat("  P/LP-cohort × SuSiE PP4>0.5:      ", length(overlap_pl_susie), "\n")

# Fisher's exact test: panel-gene enrichment in SuSiE-priority set
# Only count genes present in atlas
panel_in_atlas <- intersect(genes_in_panel, atlas$human_symbol)
ct <- matrix(c(
  length(intersect(top_susie$human_symbol, panel_in_atlas)),       # panel & PP4>0.5
  length(panel_in_atlas) - length(intersect(top_susie$human_symbol, panel_in_atlas)),  # panel & PP4≤0.5
  nrow(top_susie) - length(intersect(top_susie$human_symbol, panel_in_atlas)),         # ¬panel & PP4>0.5
  nrow(atlas) - nrow(top_susie) - length(panel_in_atlas) +
    length(intersect(top_susie$human_symbol, panel_in_atlas))                          # ¬panel & PP4≤0.5
), nrow = 2, byrow = FALSE,
   dimnames = list(c("PP4>0.5", "PP4<=0.5"), c("Panel", "NotPanel")))

ft <- fisher.test(ct)
cat("\n  Enrichment 2x2 (atlas-restricted):\n")
print(ct)
cat("  Fisher OR =", round(ft$estimate, 3), " p =", signif(ft$p.value, 3),
    " 95%CI =", round(ft$conf.int[1], 2), "-", round(ft$conf.int[2], 2), "\n")

# ----- Per-gene overlay table for the figure -----
overlay <- atlas[
  human_symbol %in% panel_in_atlas | monogenic_pl_okur2026,
  .(human_symbol, ensembl_id, gene_biotype,
    coloc_susie_best_pp4, coloc_abf_best_pp4, coloc_susie_best_gwas,
    monogenic_panel_okur2026, monogenic_pl_okur2026, monogenic_mafld_okur2026,
    monogenic_category_okur2026, monogenic_indication_okur2026)
][order(-coloc_susie_best_pp4, -coloc_abf_best_pp4)]

# ----- Write outputs -----
overlay_path <- file.path(OUT_DIR, "monogenic_x_coloc.csv")
fwrite(overlay, overlay_path)
cat("\n  Wrote overlay (", nrow(overlay), " rows): ", overlay_path, "\n", sep="")

# Atlas back to disk (in-place; backup first to be safe)
atlas_bak <- file.path(PROJ, "RNA-seq/results/multi_evidence",
                      paste0("multi_evidence_atlas.bak_pre_monogenic_",
                             format(Sys.Date(), "%Y-%m-%d"), ".csv"))
if (!file.exists(atlas_bak)) {
  file.copy(ATLAS_PATH, atlas_bak)
  cat("  Backup written: ", atlas_bak, "\n", sep="")
}
fwrite(atlas, ATLAS_PATH)
cat("  Atlas updated: ", ATLAS_PATH, "  (", ncol(atlas), " cols)\n", sep="")

# Summary text
sumtxt <- file.path(OUT_DIR, "monogenic_summary.txt")
sink(sumtxt)
cat("Okur et al. 2026 monogenic convergence summary\n")
cat("Generated:", as.character(Sys.time()), "\n")
cat("Atlas version: ", ATLAS_PATH, "\n")
cat("Panel source:  ", PANEL_PATH, " (PARTIAL — Zenodo embargo until 2026-05-01)\n\n")
cat("INPUTS\n")
cat("  Atlas:                ", nrow(atlas), " genes\n", sep="")
cat("  Panel genes (partial):", length(genes_in_panel), "\n")
cat("  P/LP variant genes:   ", length(genes_with_pl), "\n")
cat("  MAFLD-relevant genes: ", length(genes_with_mafld), "\n\n")
cat("PRIORITY-SET OVERLAPS\n")
cat("  SuSiE PP4>0.5 atlas genes:", nrow(top_susie), "\n")
cat("  ABF PP4>0.5 atlas genes:  ", nrow(top_abf),   "\n\n")
cat("  panel ∩ SuSiE PP4>0.5: ", length(overlap_panel_susie), "  →  ",
    paste(overlap_panel_susie, collapse=", "), "\n", sep="")
cat("  panel ∩ ABF   PP4>0.5: ", length(overlap_panel_abf),   "  →  ",
    paste(overlap_panel_abf, collapse=", "), "\n\n", sep="")
cat("FISHER ENRICHMENT (panel-gene over-representation in SuSiE PP4>0.5)\n")
print(ct)
cat("\n  OR =", round(ft$estimate, 3), " p =", signif(ft$p.value, 3),
    " 95%CI =", round(ft$conf.int[1], 2), "-", round(ft$conf.int[2], 2), "\n")
sink()
cat("\n  Wrote summary: ", sumtxt, "\n", sep="")
cat("\nDone.\n")
