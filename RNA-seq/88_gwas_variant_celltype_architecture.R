#!/usr/bin/env Rscript
# 88_gwas_variant_celltype_architecture.R
#
# Analysis F2 (v1) — GWAS variant × cell-type architecture.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Strategy (v1, no individual-level genotypes available):
#   Use GWAS-ATAC variant annotation (per-variant cell-type peak coverage) to
#   build a cell-type architecture map of MASLD genetic risk:
#     1. For each MASLD GWAS top gene (PNPLA3, TM6SF2, MBOAT7, HSD17B13,
#        GCKR, TRIB1, etc.), enumerate which cell types have peaks containing
#        credible-set variants.
#     2. Test: are GWAS risk genes preferentially acting in specific cell types
#        (e.g., PNPLA3 variants in hepatocyte peaks only)?
#     3. Cross with A1 primary cell-type attribution for the GWAS risk genes.
#     4. Cross with LIANA CCC to identify if GWAS risk genes are ligands/
#        receptors (fueling the hypothesis that variants perturb signaling).
#
# v2 (with individual genotypes): variant-stratified cell-type composition,
# variant-stratified CCC topology.
#
# Env: rnaseq
# Outputs: RNA-seq/results/gwas_celltype/

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GWAS_ATAC <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
LIANA     <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data")
ATT       <- file.path(BASE, "RNA-seq/results/celltype_attribution")
OUTDIR    <- file.path(BASE, "RNA-seq/results/gwas_celltype")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Canonical MASLD GWAS risk genes (expanded list)
risk_genes <- c("PNPLA3","TM6SF2","MBOAT7","HSD17B13","GCKR","TRIB1",
                "APOC3","APOE","APOA5","LIPC","CETP","PPP1R3B","MLXIPL",
                "PNPLA2","GPAM","ENPP1","NCAN","LYPLAL1","MARC1","LPIN1",
                "SLC30A10","TRPC1","PARVB","TOR1B","IL32","ERLIN1")

message("[1] Loading GWAS-ATAC variant annotation...")
vann <- fread(file.path(GWAS_ATAC, "gwas_atac_variant_annotation.csv"))
message(sprintf("  variants: %d", nrow(vann)))

# Assign atac gene
vann[, atac_gene := fifelse(linked_gene != "" & !is.na(linked_gene),
                             linked_gene, nearest_gene)]
vann <- vann[atac_gene != "" & !is.na(atac_gene)]

message("[2] Per risk gene: cell-type coverage of credible-set variants...")
risk_ct <- vann[atac_gene %in% risk_genes,
                .(n_variants_in_peak = .N,
                  cell_types = paste(sort(unique(cell_type)), collapse = ","),
                  n_cell_types = uniqueN(cell_type),
                  max_pip = max(max_pip, na.rm = TRUE),
                  any_in_susie_cs = any(in_susie_cs, na.rm = TRUE),
                  any_in_carma_cs = any(in_carma_cs, na.rm = TRUE)),
                by = atac_gene][order(-n_variants_in_peak)]

fwrite(risk_ct, file.path(OUTDIR, "risk_gene_celltype_coverage.csv"))

message("[3] Cross with A1 attribution + LIANA ligand/receptor role...")
att <- fread(file.path(ATT, "celltype_primary_attribution.csv"))
liana <- fread(file.path(LIANA, "liana_differential_interactions.csv"))

risk_gene_roles <- data.table(atac_gene = risk_genes)
risk_gene_roles[, is_ligand   := atac_gene %in% liana$ligand_complex]
risk_gene_roles[, is_receptor := atac_gene %in% liana$receptor_complex]
risk_gene_roles <- merge(risk_gene_roles, att[, .(symbol, bulk_lfc, bulk_padj,
                                                   a1_primary = primary_celltype,
                                                   a1_class = attribution_class)],
                         by.x = "atac_gene", by.y = "symbol", all.x = TRUE)
risk_gene_roles <- merge(risk_gene_roles, risk_ct, by = "atac_gene", all.x = TRUE)

fwrite(risk_gene_roles, file.path(OUTDIR, "risk_gene_celltype_roles.csv"))

message("[4] Cell-type-specific GWAS enrichment (per cell type, how many risk genes)...")
all_ct_coverage <- vann[atac_gene %in% risk_genes,
                        .(atac_gene, cell_type, max_pip, in_susie_cs)]
ct_load <- all_ct_coverage[, .(n_risk_genes = uniqueN(atac_gene)),
                            by = cell_type][order(-n_risk_genes)]
fwrite(ct_load, file.path(OUTDIR, "celltype_risk_gene_load.csv"))

summary_lines <- c(
  sprintf("MASLD canonical risk genes: %d", length(risk_genes)),
  sprintf("  with credible-set variants in ATAC peaks: %d", nrow(risk_ct)),
  "",
  "=== Risk genes × cell-type coverage + A1/LIANA roles ===",
  capture.output(print(risk_gene_roles[, .(atac_gene, n_variants_in_peak,
                                            n_cell_types, max_pip,
                                            is_ligand, is_receptor, a1_primary,
                                            bulk_lfc, bulk_padj)],
                       nrows = 40)),
  "",
  "=== Cell-type risk gene load ===",
  capture.output(print(ct_load, nrows = 20)),
  "",
  "=== Top 15 variant-dense risk genes ===",
  capture.output(print(risk_ct[1:15], nrows = 15))
)
writeLines(summary_lines, file.path(OUTDIR, "gwas_celltype_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
