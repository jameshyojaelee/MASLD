suppressPackageStartupMessages(library(data.table))
T4 <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction"
HO <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/gse193066_holdout/holdout-20260907T105813Z/results"
key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")
allz <- fread(file.path(T4, "program_interaction_edges.tsv.gz"),
              select = c("program_a","program_b","r_partial","z_partial_standardized","membership_entangled"))
coord <- fread(file.path(T4, "coordinated_edges.tsv"), select = c("program_a","program_b","r_partial","z_partial_standardized"))
allz[, ekey := key(program_a, program_b)]; coord[, ekey := key(program_a, program_b)]
allz[, is_coord := ekey %in% coord$ekey]
cat("pairs", nrow(allz), " coordinated", sum(allz$is_coord), " noncoord finite", allz[!is_coord & is.finite(r_partial), .N], "\n")
cat("coordinated r_partial quantiles:\n"); print(round(quantile(coord$r_partial, c(.05,.25,.5,.75,.95)),3))
cat("noncoord |r_partial| quantiles:\n"); print(round(quantile(abs(allz[!is_coord & is.finite(r_partial)]$r_partial), c(.05,.25,.5,.75,.95)),3))
cat("noncoord positive fraction:", round(allz[!is_coord & is.finite(r_partial), mean(r_partial > 0)],3), "\n")
pred <- function(r, zdisc, n_df, c) {
  pr <- pnorm(atanh(c * abs(r)) * sqrt(n_df - 3))           # replication sign = true sign
  pd <- pnorm(abs(zdisc))                                     # discovery sign = true sign
  pd * pr + (1 - pd) * (1 - pr)
}
nc <- allz[!is_coord & is.finite(r_partial)]
ho_df <- 105 - 8   # holdout design: intercept + sex + 4 fibrosis + 2 NAS = 8 columns
for (c in c(1, 0.8, 0.7, 0.6, 0.5)) {
  cat(sprintf("shrink c=%.1f  n=105: coord %.3f noncoord %.3f | n=258 (df 254): coord %.3f noncoord %.3f excess %.3f | n=129: coord %.3f noncoord %.3f\n",
    c, mean(pred(coord$r_partial, coord$z_partial_standardized, ho_df, c)),
       mean(pred(nc$r_partial, nc$z_partial_standardized, ho_df, c)),
       mean(pred(coord$r_partial, coord$z_partial_standardized, 254, c)),
       mean(pred(nc$r_partial, nc$z_partial_standardized, 254, c)),
       mean(pred(coord$r_partial, coord$z_partial_standardized, 254, c)) - mean(pred(nc$r_partial, nc$z_partial_standardized, 254, c)),
       mean(pred(coord$r_partial, coord$z_partial_standardized, 125, c)),
       mean(pred(nc$r_partial, nc$z_partial_standardized, 125, c))))
}
ho <- fread(file.path(HO, "A_edge_concordance.tsv.gz"))
m <- merge(ho, coord[, .(ekey, r_partial)], by = "ekey")
cat("holdout observed coord concordance", round(mean(sign(m$z_holdout) == sign(m$z_disc)),4), " n", nrow(m), "\n")
cat("holdout r vs discovery r_partial: slope (Deming-free OLS of tanh(z_holdout) on r_partial)",
    round(coef(lm(tanh(z_holdout) ~ r_partial, m))[2], 3), " mean ratio of means", round(mean(tanh(m$z_holdout))/mean(m$r_partial),3), "\n")
