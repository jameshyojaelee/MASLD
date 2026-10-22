#!/usr/bin/env Rscript
# F2abce repair: re-derive the Fig 2B same-row x from the release files and
# measure how the x value depends on the rule used when the lead position has
# several fine-mapping rows (checker text problem on R02-04/R02-05).
#   y = max PP.H4.susie over the 35 main strata of COLOC r2 (first max in file order)
#   x = SuSiE PIP at y's GWAS and lead position, from combined_finemapping.csv and
#       combined_finemapping_coloc_targeted.csv; rule "max" (round-1 candidate)
#       versus rule "min" over all rows at that key.
suppressPackageStartupMessages(library(data.table))
set.seed(20260923)
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT  <- file.path(BASE, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce/out")
FM   <- file.path(BASE, "GWAS/finemapping/results")

MAIN <- fread(file.path(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv"))[placement == "main", study_name]
coloc <- fread(file.path(FM, "susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"),
               select = c("gwas_name", "gene", "PP.H4.susie", "top_snp"))
coloc <- coloc[gwas_name %in% MAIN & !is.na(gene) & gene != "" & !is.na(top_snp) & top_snp != ""]
coloc[, pp4 := suppressWarnings(as.numeric(PP.H4.susie))]
cat("main strata:", length(MAIN), "\n")
gy <- coloc[!is.na(pp4)][, .SD[which.max(pp4)], by = gene][, .(gene, y_gwas = gwas_name, y_top_snp = top_snp, pp4)]
gy[, c("chromosome", "position") := tstrsplit(y_top_snp, ":", fixed = TRUE)]
gy[, position := as.integer(position)]
cat("genes with finite SuSiE PP.H4:", nrow(gy), "\n")
keys <- unique(gy[, .(study = y_gwas, chromosome, position)])

rd <- function(f, src) {
  x <- fread(f, select = c("study", "chromosome", "position", "locus", "susie_pip"),
             colClasses = list(character = "chromosome"))
  x <- x[keys, on = .(study, chromosome, position), nomatch = 0L]
  x[, susie_pip := suppressWarnings(as.numeric(susie_pip))]
  x[, src := src]
  x
}
fm <- rbindlist(list(rd(file.path(FM, "combined_finemapping.csv"), "canonical"),
                     rd(file.path(FM, "combined_finemapping_coloc_targeted.csv"), "coloc_targeted")))
fwrite(fm, file.path(OUT, "fig2b_lead_finemapping_rows.tsv"), sep = "\t")

mx <- function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
mn <- function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
kx <- fm[, .(n_rows = .N, n_finite = sum(is.finite(susie_pip)),
             n_distinct_pip = uniqueN(round(susie_pip[is.finite(susie_pip)], 6)),
             n_loci = uniqueN(locus), sources = paste(sort(unique(src)), collapse = "+"),
             x_max = mx(susie_pip), x_min = mn(susie_pip)),
         by = .(study, chromosome, position)]
dt <- merge(gy, kx, by.x = c("y_gwas", "chromosome", "position"),
            by.y = c("study", "chromosome", "position"), all.x = TRUE)

roster <- c("HSD17B13", "TM6SF2", "MBOAT7", "GCKR", "RORA")
summ <- function(xcol) {
  x <- dt[[xcol]]
  pl <- is.finite(x) & (dt$pp4 > 0.05 | x >= 0.99 | dt$gene %in% roster)
  pl_nolab <- is.finite(x) & (dt$pp4 > 0.05 | x >= 0.99)
  p4 <- pmax(dt$pp4, 0)
  data.table(rule = xcol, finite_x = sum(is.finite(x)), plotted = sum(pl),
             plotted_without_label_rule = sum(pl_nolab),
             corner_pip_ge_0.9_pp4_lt_0.5 = sum(pl & x >= 0.9 & p4 < 0.5),
             topright_pip_ge_0.9_pp4_gt_0.5 = sum(pl & x >= 0.9 & p4 > 0.5),
             rora_x = round(dt[gene == "RORA"][[xcol]], 4))
}
res <- rbind(summ("x_max"), summ("x_min"))
dt[, plotted_max := is.finite(x_max) & (pp4 > 0.05 | x_max >= 0.99 | gene %in% roster)]
dt[, plotted_min := is.finite(x_min) & (pp4 > 0.05 | x_min >= 0.99 | gene %in% roster)]
extra <- data.table(
  genes_any_finemapping_row = dt[!is.na(n_rows), .N],
  genes_finite_susie_pip = dt[is.finite(x_max), .N],
  plotted_max_with_multi_rows = dt[plotted_max == TRUE & n_rows > 1, .N],
  plotted_max_with_differing_pips = dt[plotted_max == TRUE & n_distinct_pip > 1, .N],
  plotted_max_differing_pips_by_sources = dt[plotted_max == TRUE & n_distinct_pip > 1, paste(names(table(sources)), table(sources), sep = "=", collapse = ";")],
  plotted_max_differing_pips_multi_locus = dt[plotted_max == TRUE & n_distinct_pip > 1 & n_loci > 1, .N],
  median_abs_diff_differing = dt[plotted_max == TRUE & n_distinct_pip > 1, round(median(x_max - x_min), 4)],
  max_abs_diff_differing = dt[plotted_max == TRUE & n_distinct_pip > 1, round(max(x_max - x_min), 4)],
  changed_corner_or_gate = dt[plotted_max == TRUE & (plotted_min != plotted_max |
                               ((x_max >= 0.9) != (x_min >= 0.9))), paste(gene, collapse = ",")])
print(res); print(t(extra))
print(dt[gene %in% c(roster, "POLL"), .(gene, y_gwas, y_top_snp, pp4 = round(pp4, 4), n_rows, n_distinct_pip, n_loci, sources, x_max = round(x_max, 4), x_min = round(x_min, 4))])
fwrite(res, file.path(OUT, "fig2b_multirow_rule_summary.tsv"), sep = "\t")
fwrite(data.table(metric = names(extra), value = unlist(lapply(extra, as.character))),
       file.path(OUT, "fig2b_multirow_extra.tsv"), sep = "\t")
fwrite(dt[order(-plotted_max, gene)], file.path(OUT, "fig2b_same_row_x_by_rule.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_fig2b_multirow.txt"))
cat("done\n")
