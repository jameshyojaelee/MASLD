#!/usr/bin/env Rscript
# SuSiE-COLOC genes per ancestry under canonical policy:
#   UKBB EUR  → PolyFun
#   All other → 1KG
#
# Approach: build a unified per-(gene, GWAS) SuSiE-COLOC result by:
#   - Reading per-study chr files for UKBB EUR studies from
#     susie_coloc_polyfun/{study}/susie_coloc_chr*.csv (method=="susie")
#   - Reading per-study chr files for all OTHER studies from
#     susie_coloc_1kg/{study}/susie_coloc_chr*.csv  (method=="susie")
# Then per (gene), keep best SuSiE PP.H4 across the union, attribute to the
# study whose run produced the best PP.H4 → ancestry rollup.
suppressPackageStartupMessages({ library(data.table) })

ukbb_studies <- c(
  "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
  "2023_36280732_NAFLD_UKBB_EUR",
  "2021_34128465_PDFF_EUR",
  "2021_34957434_PDFF_EUR",
  "2022_36402844_PDFF_EUR")

reg <- fread("GWAS/finemapping/config/gwas_registry.tsv",
             select = c("study_name","ancestry"))
setnames(reg, "study_name", "study")

read_per_study <- function(root) {
  out <- list()
  for (sd in list.dirs(root, recursive = FALSE, full.names = FALSE)) {
    files <- list.files(file.path(root, sd), pattern = "csv$", full.names = TRUE)
    if (!length(files)) next
    d <- rbindlist(lapply(files, fread), fill = TRUE,
                   use.names = TRUE)
    if (!"method" %in% names(d)) next
    d <- d[method == "susie"]
    if (!nrow(d)) next
    keep <- intersect(c("gene","ensembl","gwas_name","PP.H4.susie","PP.H4.abf"),
                      names(d))
    out[[sd]] <- d[, ..keep][, source := sd]
  }
  rbindlist(out, fill = TRUE)
}

cat("Reading PolyFun per-study (UKBB EUR studies only)...\n")
poly <- read_per_study("GWAS/finemapping/results/susie_coloc_polyfun")
poly <- poly[source %in% ukbb_studies]
cat("  rows:", nrow(poly), " studies:", uniqueN(poly$source), "\n")

# Non-UKBB EUR studies: 1KG-EUR (susie_coloc_1kg/)
cat("Reading 1KG per-study (non-UKBB EUR studies)...\n")
non_ukbb_eur <- reg[ancestry == "EUR" & !study %in% ukbb_studies, study]
kg_eur <- read_per_study("GWAS/finemapping/results/susie_coloc_1kg")
kg_eur <- kg_eur[source %in% non_ukbb_eur]
cat("  rows:", nrow(kg_eur), " studies:", uniqueN(kg_eur$source), "\n")

# EAS/AFR/SAS studies: ancestry-matched 1KG (live in susie_coloc/)
cat("Reading susie_coloc per-study (EAS/AFR/SAS — ancestry-matched 1KG)...\n")
non_eur_studies <- reg[ancestry %in% c("EAS","AFR","SAS"), study]
kg_other <- read_per_study("GWAS/finemapping/results/susie_coloc")
kg_other <- kg_other[source %in% non_eur_studies]
cat("  rows:", nrow(kg_other), " studies:", uniqueN(kg_other$source), "\n")

all_susie <- rbindlist(list(poly, kg_eur, kg_other), fill = TRUE)
setnames(all_susie, "source", "study")
all_susie <- merge(all_susie, reg, by = "study", all.x = TRUE)
all_susie[, pp4 := PP.H4.susie]

# ---------- Counts per ancestry × threshold ----------
thr <- c(0.5, 0.8, 0.9)
cnt <- list()
for (th in thr) {
  c1 <- all_susie[!is.na(pp4) & pp4 >= th,
                  .(n_genes_sig = uniqueN(gene)), by = ancestry]
  c1[, threshold := paste0("PP4>=", th)]
  cnt[[as.character(th)]] <- c1
}
ancestry_counts <- rbindlist(cnt)
ancestry_wide <- dcast(ancestry_counts, ancestry ~ threshold,
                        value.var = "n_genes_sig")
fwrite(ancestry_wide,
       "scripts/figures/sketches_fig3_intro/data/canonical_susie_coloc_per_ancestry.csv")

cat("\n=== Canonical SuSiE-COLOC genes per ancestry (PolyFun for UKBB EUR, 1KG otherwise) ===\n")
print(ancestry_wide)

# ---------- Per-study × threshold ----------
study_counts <- all_susie[!is.na(pp4),
  .(n_genes_pp4_05 = uniqueN(gene[pp4 >= 0.5]),
    n_genes_pp4_08 = uniqueN(gene[pp4 >= 0.8]),
    n_genes_pp4_09 = uniqueN(gene[pp4 >= 0.9])),
  by = .(study, ancestry)]
fwrite(study_counts,
       "scripts/figures/sketches_fig3_intro/data/canonical_susie_coloc_per_study.csv")

cat("\n=== Per-study (top 20 by PP4>=0.5 count) ===\n")
setorder(study_counts, -n_genes_pp4_05)
print(study_counts[1:20])

# ---------- Total unique genes ----------
total_genes <- all_susie[!is.na(pp4), .(
  total_pp4_05 = uniqueN(gene[pp4 >= 0.5]),
  total_pp4_08 = uniqueN(gene[pp4 >= 0.8]),
  total_pp4_09 = uniqueN(gene[pp4 >= 0.9]))]
cat("\n=== Total unique genes across all ancestries ===\n")
print(total_genes)

cat("\nDone.\n")
