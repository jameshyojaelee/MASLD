#!/usr/bin/env Rscript
# Replace per-study loci counts in gwas_portfolio.csv with canonical numbers
# (PolyFun for UKBB EUR, 1KG for everything else). Also recompute per-cell
# unique COLOC gene counts so the mosaic doesn't double-count genes that hit
# multiple GWAS in the same (ancestry, trait) cell.
suppressPackageStartupMessages({ library(data.table) })

DAT <- "scripts/figures/sketches_fig3_intro/data"

ukbb_studies <- c("UKBB_ALT","UKBB_AST","UKBB_GGT",
                  "2023_36280732_NAFLD_UKBB_EUR",
                  "2021_34128465_PDFF_EUR",
                  "2021_34957434_PDFF_EUR",
                  "2022_36402844_PDFF_EUR")

# ---- Canonical per-study loci ----
canonical <- fread(file.path(DAT, "canonical_loci_per_study.csv"))

p <- fread(file.path(DAT, "gwas_portfolio.csv"))
p <- merge(p, canonical[, .(study = study, n_loci_canonical)],
           by.x = "study_name", by.y = "study", all.x = TRUE)
p[is.na(n_loci_canonical), n_loci_canonical := 0]
p[, susie_n_loci_converged := n_loci_canonical]
p[, n_loci_canonical := NULL]

# ---- Recompute UNIQUE COLOC genes per (ancestry, trait) cell ----
# Use the per-study canonical SuSiE-COLOC data (already aggregated).
read_per_study <- function(root) {
  out <- list()
  for (sd in list.dirs(root, recursive = FALSE, full.names = FALSE)) {
    files <- list.files(file.path(root, sd), pattern = "csv$", full.names = TRUE)
    if (!length(files)) next
    d <- rbindlist(lapply(files, fread), fill = TRUE, use.names = TRUE)
    if (!"method" %in% names(d)) next
    d <- d[method == "susie"]
    if (!nrow(d)) next
    keep <- intersect(c("gene","gwas_name","PP.H4.susie"), names(d))
    out[[sd]] <- d[, ..keep][, study := sd]
  }
  rbindlist(out, fill = TRUE)
}

reg <- fread("GWAS/finemapping/config/gwas_registry.tsv",
             select = c("study_name","ancestry"))
setnames(reg, "study_name", "study")

poly  <- read_per_study("GWAS/finemapping/results/susie_coloc_polyfun")
poly  <- poly[study %in% ukbb_studies]
non_ukbb_eur <- reg[ancestry == "EUR" & !study %in% ukbb_studies, study]
kg_eur <- read_per_study("GWAS/finemapping/results/susie_coloc_1kg")
kg_eur <- kg_eur[study %in% non_ukbb_eur]
kg_other <- read_per_study("GWAS/finemapping/results/susie_coloc")
kg_other <- kg_other[study %in% reg[ancestry %in% c("EAS","AFR","SAS"), study]]

all_susie <- rbindlist(list(poly, kg_eur, kg_other), fill = TRUE)
all_susie <- merge(all_susie, p[, .(study_name, ancestry, trait)],
                    by.x = "study", by.y = "study_name", all.x = TRUE)

# Per-(ancestry, trait) unique-gene COLOC counts
cell_genes <- all_susie[!is.na(PP.H4.susie) & PP.H4.susie >= 0.5,
                         .(coloc_n_genes_pp4_05_unique = uniqueN(gene)),
                         by = .(ancestry, trait)]
cell_genes[, ancestry := as.character(ancestry)]
cell_genes[, trait := as.character(trait)]

# Also recompute per-study unique counts (this is what the old portfolio had,
# but inconsistent with canonical setup so we re-derive)
study_genes <- all_susie[!is.na(PP.H4.susie) & PP.H4.susie >= 0.5,
                          .(coloc_n_genes_pp4_05_canonical = uniqueN(gene)),
                          by = study]
p <- merge(p, study_genes, by.x = "study_name", by.y = "study", all.x = TRUE)
p[is.na(coloc_n_genes_pp4_05_canonical), coloc_n_genes_pp4_05_canonical := 0]
p[, coloc_n_genes_pp4_05 := coloc_n_genes_pp4_05_canonical]
p[, coloc_n_genes_pp4_05_canonical := NULL]

fwrite(p, file.path(DAT, "gwas_portfolio.csv"))
fwrite(cell_genes, file.path(DAT, "cell_unique_genes.csv"))

cat("=== Portfolio refreshed with canonical loci/gene counts ===\n")
print(p[, .(study_name, ancestry, trait, susie_n_loci_converged, coloc_n_genes_pp4_05)])

cat("\n=== Per-cell UNIQUE COLOC genes (PP4>=0.5) ===\n")
print(cell_genes[order(ancestry, trait)])

cat("\nDone.\n")
