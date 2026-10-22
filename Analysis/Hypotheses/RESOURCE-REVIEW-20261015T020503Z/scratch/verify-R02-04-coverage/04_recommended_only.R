suppressPackageStartupMessages(library(data.table))
S <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-04-coverage/"
st <- readRDS(paste0(S, "state.rds")); base <- st$base; fm0 <- st$fm0
roster <- c("HSD17B13","TM6SF2","MBOAT7","GCKR","RORA")
run <- function(expr, tag) {
  fm <- copy(fm0); fm[, cp := eval(expr)]; fm <- fm[is.finite(cp)]
  fm <- fm[order(study, chromosome, position, -cp, source_rank)]
  fu <- unique(fm, by = c("study","chromosome","position"))
  cm <- merge(st$coloc, fu[, .(study, chromosome, position, cp)], by.x = c("gwas_name","snp_chr","snp_pos"), by.y = c("study","chromosome","position"))
  gx <- cm[order(gene, -cp), .SD[1], by = gene][, .(gene, x = cp)]
  A <- merge(base[, .(gene, causal_pip, pp4, plotted)], gx, by = "gene", all.x = TRUE)
  A[, pl := !is.na(x) & ((pp4 > 0.05) | (x >= 0.99) | gene %in% roster)]
  Ap <- A[plotted == TRUE]
  cat(sprintf("%-28s changed among 695: %d | plotted: %d | corner(PIP>=0.9,PP.H4<0.5): %d | RORA x: %.4f\n", tag,
      sum(is.na(Ap$x) | abs(Ap$x - Ap$causal_pip) > 1e-4), sum(A$pl), sum(A$pl & A$x >= 0.9 & A$pp4 < 0.5), A[gene == "RORA", x]))
}
run(quote(recommended_pip), "recommended_pip only")
run(quote(fifelse(ancestry != "EUR" | !(susie_converged %in% TRUE), carma_pip, suppressWarnings(pmax(susie_pip, carma_pip, na.rm = TRUE)))), "pre-joint recommended")
run(quote(susie_pip), "SuSiE only")
run(quote(carma_pip), "CARMA only")
