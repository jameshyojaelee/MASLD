#!/usr/bin/env Rscript

# Regenerate the canonical Figure 4 panel manifest from explicit generators and
# primary analysis products. This prevents a PDF from serving as its own source.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
FIG4 <- file.path(BASE, "figures/main/fig4_validation")

manifest <- data.table(
  callout = paste0("4", LETTERS[1:6]),
  filename = c(
    "panels/fig4a_overview_cascade.pdf",
    "panels/fig4b_protein_triage.pdf",
    "panels/fig4c_mrna_protein_composite.pdf",
    "panels/fig4d_snatac_accessibility.pdf",
    "panels/fig4e_spatial_program_maps.pdf",
    "panels/fig4f_multimodal_program_summary.pdf"
  ),
  source = c(
    "scripts/figures/fig4a_input_firewall.R",
    "scripts/figures/fig4b_protein_triage.R",
    "scripts/figures/composite_mrna_protein.R",
    "scripts/figures/fig4d_snatac_accessibility.R",
    "scripts/figures/fig4e_spatial_program_maps.py",
    "scripts/figures/fig4f_multimodal_program_summary.R"
  ),
  primary_inputs = c(
    "frozen_programs.tsv; fig4_input_contract.tsv",
    "protein_de_adjusted.tsv; multi_evidence_atlas.csv; prioritized_universe_FINAL.txt",
    "panel4c_adjusted_abundance.tsv; panel4c_mrna_protein.tsv; panel4c_histology_partial.tsv; panel4c_process_enrichment.tsv",
    "gene_level_gwas_atac.csv; multi_evidence_atlas.csv; prioritized_universe_FINAL.txt",
    "spatial_program_results.tsv; spatial_section_results.tsv; frozen_program_membership.tsv; two spatial_deconvolved h5ad files",
    "program_context_wide.tsv; module_protein_histology_burden_masld.tsv; spatial_program_results.tsv"
  ),
  status = c("accepted", "accepted", "accepted_cornerstone", "accepted_context", "accepted", "accepted_synthesis")
)

paths <- file.path(FIG4, manifest$filename)
if (any(!file.exists(paths))) {
  stop("Cannot build Figure 4 manifest; missing: ", paste(manifest$filename[!file.exists(paths)], collapse = ", "))
}
manifest[, `:=`(
  sha256 = vapply(paths, digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  bytes = as.numeric(file.info(paths)$size)
)]
setcolorder(manifest, c("callout", "filename", "source", "primary_inputs", "sha256", "bytes", "status"))
out <- file.path(FIG4, "CANONICAL_MAIN_PANELS.tsv")
fwrite(manifest, out, sep = "\t", quote = FALSE)
message("[fig4 manifest] wrote ", out, " with ", nrow(manifest), " canonical panels")
