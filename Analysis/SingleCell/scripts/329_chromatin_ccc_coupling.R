#!/usr/bin/env Rscript
# 329_chromatin_ccc_coupling.R
#
# Analysis G2 — Chromatin-informed CCC coupling.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Strategy:
#   1. Load GWAS-ATAC variant annotation (3,724 variants in scATAC peaks,
#      from GWAS/finemapping/results/gwas_atac/).
#   2. Extract LIANA ligand + receptor gene sets (all unique).
#   3. Enrichment test: are ATAC-peak-linked genes more likely to be CCC
#      receptors or ligands than expected (vs background universe of all
#      expressed genes)?
#   4. For each LIANA MASLD-enriched LR pair, flag whether receptor is in
#      ATAC-regulated set AND variant is motif-disrupted — prioritizes
#      "genetically perturbed CCC axes" as therapeutic targets.
#   5. Compute: fraction of bulk DEGs that are (a) CCC receptors, (b) GWAS
#      ATAC-peak-linked, (c) both.
#
# Env: rnaseq
# Outputs: Analysis/SingleCell/results_gpu_v2/chromatin_ccc/

suppressPackageStartupMessages({
  library(data.table)
})

BASE    <- Sys.getenv("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GWAS_ATAC <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
LIANA   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data")
CCC_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc")
INT_RES <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/chromatin_ccc")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading GWAS-ATAC variant annotation...")
vann <- fread(file.path(GWAS_ATAC, "gwas_atac_variant_annotation.csv"))
motif <- fread(file.path(GWAS_ATAC, "motif_disruption_scores.csv"))
message(sprintf("  %d variants in peaks; %d motif disruption events",
                nrow(vann), nrow(motif)))

# ATAC-regulated genes: linked_gene or nearest_gene (prefer linked)
vann[, atac_gene := fifelse(linked_gene != "" & !is.na(linked_gene),
                             linked_gene, nearest_gene)]
atac_regulated_genes <- unique(vann$atac_gene[vann$atac_gene != "" & !is.na(vann$atac_gene)])
message(sprintf("  Unique ATAC-regulated genes: %d", length(atac_regulated_genes)))

# Join motif disruption to variant annotation via SNP_id / variant_id to get linked_gene
snp_col <- intersect(c("SNP_id","variant_id","rsid"), names(motif))[1]
if (!is.na(snp_col)) {
  setnames(motif, snp_col, "variant_id")
  vann_lookup <- unique(vann[, .(variant_id, atac_gene)])
  motif_joined <- merge(unique(motif[, .(variant_id)]),
                        vann_lookup, by = "variant_id", allow.cartesian = TRUE)
  motif_disrupted_genes <- unique(motif_joined$atac_gene[motif_joined$atac_gene != "" & !is.na(motif_joined$atac_gene)])
} else motif_disrupted_genes <- character()
message(sprintf("  Motif-disrupted genes: %d", length(motif_disrupted_genes)))

message("[2] Loading LIANA LR pairs + bulk DE...")
liana <- fread(file.path(LIANA, "liana_differential_interactions.csv"))
ligands   <- unique(liana$ligand_complex)
receptors <- unique(liana$receptor_complex)
# Complexes may have "_" separated subunits — split them
split_complex <- function(x) unlist(strsplit(x, "_"))
ligand_genes   <- unique(unlist(lapply(ligands, split_complex)))
receptor_genes <- unique(unlist(lapply(receptors, split_complex)))
message(sprintf("  LIANA unique ligand genes: %d, receptor genes: %d",
                length(ligand_genes), length(receptor_genes)))

bulk <- fread(file.path(INT_RES, "dream_results_ashr.csv"),
              select = c("gene","symbol","logFC","padj"))
bulk <- bulk[!is.na(symbol) & symbol != ""]
universe <- unique(bulk$symbol)
bulk_degs <- unique(bulk[padj < 0.05 & abs(logFC) > 0.5, symbol])
message(sprintf("  Bulk universe: %d symbols; DEGs: %d", length(universe), length(bulk_degs)))

message("[3] Enrichment tests (hypergeometric)...")
hyper_test <- function(set_a, set_b, bg) {
  set_a <- intersect(set_a, bg); set_b <- intersect(set_b, bg)
  k <- length(intersect(set_a, set_b))
  m <- length(set_a); n <- length(bg) - m; k_sample <- length(set_b)
  p <- phyper(k - 1, m, n, k_sample, lower.tail = FALSE)
  or <- (k / k_sample) / (m / length(bg))
  data.table(k = k, n_a = m, n_b = k_sample, n_bg = length(bg),
             odds_ratio = or, pvalue = p)
}

enr <- rbind(
  hyper_test(atac_regulated_genes, receptor_genes, universe)[, test := "ATAC-regulated enriched for LIANA receptors"],
  hyper_test(atac_regulated_genes, ligand_genes,   universe)[, test := "ATAC-regulated enriched for LIANA ligands"],
  hyper_test(motif_disrupted_genes, receptor_genes, universe)[, test := "Motif-disrupted enriched for LIANA receptors"],
  hyper_test(motif_disrupted_genes, ligand_genes,   universe)[, test := "Motif-disrupted enriched for LIANA ligands"],
  hyper_test(atac_regulated_genes, bulk_degs,       universe)[, test := "ATAC-regulated enriched for bulk DEGs"],
  hyper_test(intersect(atac_regulated_genes, receptor_genes),
             bulk_degs, universe)[, test := "ATAC-regulated CCC receptors enriched in bulk DEGs"],
  fill = TRUE
)
fwrite(enr, file.path(OUTDIR, "gwas_atac_ccc_enrichment.csv"))

message("[4] Prioritized receptor-disrupted CCC axes...")
# For each LIANA LR pair, flag receptor genes present in ATAC-regulated set
# (and motif-disrupted set for a higher tier)
rec_set <- data.table(receptor_gene = receptor_genes,
                      atac_regulated = receptor_genes %in% atac_regulated_genes,
                      motif_disrupted = receptor_genes %in% motif_disrupted_genes)
fwrite(rec_set, file.path(OUTDIR, "ligand_receptor_atac_annotation.csv"))

# Annotate LIANA differential pairs
liana_flagged <- copy(liana)
liana_flagged[, receptor_first_gene := sapply(strsplit(receptor_complex, "_"), `[`, 1)]
liana_flagged[, ligand_first_gene   := sapply(strsplit(ligand_complex,   "_"), `[`, 1)]
liana_flagged[, receptor_atac_regulated := receptor_first_gene %in% atac_regulated_genes]
liana_flagged[, ligand_atac_regulated   := ligand_first_gene   %in% atac_regulated_genes]
liana_flagged[, receptor_motif_disrupted := receptor_first_gene %in% motif_disrupted_genes]
liana_flagged[, ligand_motif_disrupted   := ligand_first_gene   %in% motif_disrupted_genes]

# Priority: MASLD-enriched + receptor ATAC-regulated + motif-disrupted
priority <- liana_flagged[score_diff > 0.1 &
                          (receptor_atac_regulated | ligand_atac_regulated) &
                          (receptor_motif_disrupted | ligand_motif_disrupted)]
priority <- priority[order(-score_diff)]
fwrite(priority[, .(source, target, ligand_complex, receptor_complex, score_diff,
                     receptor_atac_regulated, receptor_motif_disrupted,
                     ligand_atac_regulated,   ligand_motif_disrupted)],
       file.path(OUTDIR, "priority_variant_disrupted_CCC_axes.csv"))

message("[5] Summary...")
summary_lines <- c(
  sprintf("GWAS-ATAC variants in peaks: %d", nrow(vann)),
  sprintf("Unique ATAC-regulated genes: %d", length(atac_regulated_genes)),
  sprintf("Motif-disrupted genes: %d", length(motif_disrupted_genes)),
  sprintf("LIANA unique ligand genes: %d, receptor genes: %d",
          length(ligand_genes), length(receptor_genes)),
  "",
  "Enrichment tests:",
  capture.output(print(enr[, .(test, k, n_a, n_b, odds_ratio, pvalue)], nrows = 20)),
  "",
  sprintf("LIANA MASLD-enriched LR pairs with variant-disrupted receptor/ligand: %d",
          nrow(priority)),
  "Top 30 priority axes:",
  capture.output(print(priority[1:min(30, .N),
                                .(source, target, ligand_complex, receptor_complex,
                                  score_diff, receptor_motif_disrupted,
                                  ligand_motif_disrupted)], nrows = 30)))
writeLines(summary_lines, file.path(OUTDIR, "gwas_atac_ccc_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
