#!/usr/bin/env Rscript
# Task list for the eQTL/GWAS sdY sensitivity of the ABF arm.
# Pairs: tier-1/2 studies, non-blank Ensembl, master PP.H4.abf in [0.1, 0.95]
# (the range from which a 0.5 crossing is plausible under a prior-SD change of
# <= 1.4x). One task per (study, chr); one gene filter per task.
suppressPackageStartupMessages(library(data.table))
A  <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive"
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping"
m <- fread(file.path(FM, "results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"),
           select = c("gwas_name", "ensembl", "chr", "PP.H4.abf"))
reg <- fread(file.path(FM, "config/gwas_registry.tsv"))
t12 <- reg$study_name[grepl("NAFLD|NASH|PDFF|ALT|AST|GGT", reg$study_name)]
stopifnot(length(t12) == 35L)
sdy <- fread(file.path(A, "sdy_by_study.tsv"))
EQTL_SDY <- sdy[arm == "eqtl", sdY_regression]        # coloc:::sdY.est method
x <- m[gwas_name %in% t12 & ensembl != "" & PP.H4.abf >= 0.1 & PP.H4.abf <= 0.95]
tasks <- unique(x[, .(gwas_name, chr)])[order(chr, gwas_name)]
tasks[, sdY_gwas := sdy$sdY_regression[match(gwas_name, sdy$study)]]
tasks[is.na(sdY_gwas), sdY_gwas := 1]                # binary studies: sdY not used
tasks[, sdY_eqtl := EQTL_SDY]
tasks[, filter := file.path(A, "abf_sdy/filters", sprintf("%s_chr%d.tsv", gwas_name, chr))]
for (i in seq_len(nrow(tasks)))
  fwrite(x[gwas_name == tasks$gwas_name[i] & chr == tasks$chr[i], .(gwas_name, chr, ensembl)],
         tasks$filter[i], sep = "\t")
fwrite(tasks, file.path(A, "abf_sdy/tasks.tsv"), sep = "\t")
cat("pairs", nrow(x), "| tasks", nrow(tasks), "| eQTL sdY", EQTL_SDY, "\n")
