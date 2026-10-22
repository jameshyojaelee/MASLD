# H40: export the two covariate gene sets used to define the HNF4A-loss axis
# and the fibrosis (collagen) proxy score. Run once before any outcome is read.
#
# Usage: Rscript --vanilla export_covariate_gene_sets.R <out_dir>
#
# HNF4A axis: CollecTRI (Mueller-Dott 2023) HNF4A regulon via decoupleR /
#   OmnipathR, signed by `mor` (+1 activation, -1 repression).
# Fibrosis proxy: MSigDB REACTOME_COLLAGEN_FORMATION via msigdbr.
suppressPackageStartupMessages({
  library(decoupleR)
  library(msigdbr)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
stopifnot(dir.exists(out))

net <- as.data.table(decoupleR::get_collectri(organism = "human", split_complexes = FALSE))
hnf <- net[source == "HNF4A", .(source, target, mor)]
if (nrow(hnf) < 50L) stop("HNF4A regulon has fewer than 50 targets: ", nrow(hnf))
fwrite(hnf, file.path(out, "collectri_HNF4A_regulon.tsv"), sep = "\t")

ms <- as.data.table(msigdbr(species = "Homo sapiens"))
col <- unique(ms[gs_name == "REACTOME_COLLAGEN_FORMATION",
                 .(gs_name, gs_id, db_version, gene_symbol, ensembl_gene)])
if (nrow(col) < 30L) stop("REACTOME_COLLAGEN_FORMATION has fewer than 30 rows")
fwrite(col, file.path(out, "msigdb_reactome_collagen_formation.tsv"), sep = "\t")

cat(sprintf("HNF4A regulon: %d targets (%d activating, %d repressing)\n",
            nrow(hnf), sum(hnf$mor > 0), sum(hnf$mor < 0)))
cat(sprintf("Collagen formation: %d rows, %d unique Ensembl IDs\n",
            nrow(col), uniqueN(col$ensembl_gene)))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_gene_set_export.txt"))
