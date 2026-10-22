# Leave-one-cohort-out for the Fig 4E four-cohort arm: REML + KH (published method), BH over 224.
suppressPackageStartupMessages({library(data.table); library(metafor)})
A <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_plus_GSE174478/"
ce <- fread(paste0(A, "cohort_effects.tsv"))[estimable == TRUE]
fitp <- function(y, s) { f <- tryCatch(rma(yi = y, sei = s, method = "REML", test = "knha"), error = function(e) NULL); if (is.null(f)) NA_real_ else f$pval }
for (mk in c("primary", "composition_pc")) for (drop in c("none", unique(ce$cohort))) {
  d <- ce[model_kind == mk & cohort != drop]
  r <- d[, .(p = fitp(.SD$beta, .SD$se_hc3), k = .N), by = .(feature_id, axis)]
  r[, q := p.adjust(p, "BH")]
  cat(sprintf("%-15s drop=%-10s k=%d fibrosis=%2d nas=%2d min_q=%.3g\n", mk, drop, max(r$k),
              r[axis == "fibrosis", sum(q < 0.05, na.rm = TRUE)], r[axis == "nas", sum(q < 0.05, na.rm = TRUE)], min(r$q, na.rm=TRUE)))
}
