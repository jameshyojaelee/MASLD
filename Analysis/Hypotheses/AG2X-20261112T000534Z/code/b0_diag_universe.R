#!/usr/bin/env Rscript
# Diagnostic only: which credible_sets.csv rows define the 6,579-variant substrate universe?
suppressPackageStartupMessages(library(data.table))
RP <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"; FM <- file.path(RP, "GWAS/finemapping/results")
akey <- function(c, p, a, b) paste(c, p, pmin(a, b), pmax(a, b), sep = ":")
sub <- fread(file.path(FM, "seqfunc/variant_substrate_hg38.tsv"))
sp <- tstrsplit(sub$variant_id_hg19, ":"); sub[, k := akey(sp[[1]], sp[[2]], sp[[3]], sp[[4]])]
cs <- fread(file.path(FM, "credible_sets.csv"))
cs[, k := akey(chromosome, position, allele1, allele2)]
cs[, in_sub := k %in% sub$k]
cat("rows", nrow(cs), "in_sub rows", sum(cs$in_sub), "unique in_sub", uniqueN(cs[in_sub == TRUE]$k), "of substrate", uniqueN(sub$k), "\n")
f <- function(x) print(cs[, .(rows = .N, uniq = uniqueN(k), uniq_in_sub = uniqueN(k[in_sub])), by = x][order(-rows)][1:min(.N, 60)])
cs[, susie_set := !is.na(susie_cs) & susie_cs > 0]
cs[, carma_set := !is.na(carma_cs) & carma_cs > 0]
f("susie_set"); f("carma_set"); f("in_joint_cs"); f("both_in_cs"); f("susie_reliable"); f("ancestry"); f("trait")
cs[, mp := cut(max_pip, c(-1, 0.001, 0.01, 0.1, 0.5, 0.9, 2))]; f("mp")
cs[, rp := !is.na(recommended_pip)]; f("rp")
f("study")
# per-variant: max over rows of max_pip vs in substrate
v <- cs[, .(mx = max(max_pip, na.rm = TRUE), anyrel = any(susie_reliable), anysusie = any(susie_set), anyjoint = any(in_joint_cs), in_sub = in_sub[1]), by = k]
print(v[, .N, by = .(in_sub, anysusie, anyjoint)][order(-N)])
print(v[, .(n = .N, sub = sum(in_sub)), by = .(b = cut(mx, c(-1, 0.001, 0.01, 0.05, 0.1, 0.5, 2)))][order(b)])
cons <- fread(file.path(FM, "credible_set_variant_consequences.csv")); cons[, k := akey(chr, pos, ref, alt)]
cat("consequence universe", nrow(cons), "in cs", sum(cons$k %in% cs$k), "\n")
print(summary(cons$max_pip))
