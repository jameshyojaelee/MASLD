suppressPackageStartupMessages(library(data.table))
S <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-04-coverage/"
st <- readRDS(paste0(S, "state.rds")); base <- st$base; fm0 <- st$fm0
pb <- base[plotted == TRUE]
pb[, c("c","p") := tstrsplit(y_top_snp, ":", fixed = TRUE, type.convert = TRUE)]
for (m in c("susie_pip","carma_pip")) {
  f <- fm0[is.finite(get(m)), .(v = max(get(m))), by = .(study, chromosome, position)]
  x <- merge(pb[, .(gene, susie_gwas, c, p, causal_pip)], f, by.x = c("susie_gwas","c","p"), by.y = c("study","chromosome","position"), all.x = TRUE)
  cat(m, "at y GWAS+lead: of 695, x missing:", sum(is.na(x$v)), " x differs >1e-4 (incl missing):", sum(is.na(x$v) | abs(x$v - x$causal_pip) > 1e-4), "\n")
}
