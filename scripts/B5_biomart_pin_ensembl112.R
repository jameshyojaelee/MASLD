#!/usr/bin/env Rscript
# B5: Regenerate mouse↔human ortholog mapping pinned to Ensembl release 112
# (2024-10) and add gene symbols + biotype.
# Output: streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz (worktree copy)

suppressPackageStartupMessages({
  library(biomaRt)
  library(data.table)
})

WT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation"
out_dir <- file.path(WT, "streamlit_deg_explorer/data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

VER <- 112  # Ensembl release pinned 2024-10

cat("=== Connecting to Ensembl release", VER, "===\n")
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = VER)
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = VER)

cat("Querying mouse → human homologs...\n")
mh <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype",
                 "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_orthology_type",
                 "hsapiens_homolog_orthology_confidence",
                 "hsapiens_homolog_perc_id",
                 "hsapiens_homolog_perc_id_r1"),
  mart = mouse)
setDT(mh)
mh <- mh[hsapiens_homolog_ensembl_gene != ""]
cat("Mouse genes with human ortholog:", nrow(mh), "\n")

# Re-fetch human biotype
cat("Querying human gene metadata...\n")
hg <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = unique(mh$hsapiens_homolog_ensembl_gene),
  mart = human)
setDT(hg)
setnames(hg, c("human_ensembl_gene_id", "human_symbol", "human_biotype"))

setnames(mh,
  c("ensembl_gene_id", "external_gene_name", "gene_biotype",
    "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type", "hsapiens_homolog_orthology_confidence",
    "hsapiens_homolog_perc_id", "hsapiens_homolog_perc_id_r1"),
  c("mouse_ensembl_gene_id", "mouse_symbol", "mouse_biotype",
    "human_ensembl_gene_id", "human_symbol_homolog",
    "orthology_type", "orthology_confidence",
    "perc_id_mouse_in_human", "perc_id_human_in_mouse"))

ort <- merge(mh, hg, by = "human_ensembl_gene_id", all.x = TRUE)
ort[, ensembl_release := VER]
ort[, query_date := format(Sys.Date())]

# Order columns
col_order <- c("mouse_ensembl_gene_id", "human_ensembl_gene_id",
               "orthology_type", "orthology_confidence",
               "mouse_symbol", "human_symbol", "human_symbol_homolog",
               "mouse_biotype", "human_biotype",
               "perc_id_mouse_in_human", "perc_id_human_in_mouse",
               "ensembl_release", "query_date")
setcolorder(ort, intersect(col_order, names(ort)))

out_path <- file.path(out_dir, "mouse_human_orthologs.tsv.gz")
fwrite(ort, out_path, sep = "\t", compress = "gzip")
cat("Wrote", nrow(ort), "rows to", out_path, "\n")

# Summary table
cat("\nOrthology type counts:\n")
print(table(ort$orthology_type))
cat("\nConfidence:\n"); print(table(ort$orthology_confidence, useNA="always"))
