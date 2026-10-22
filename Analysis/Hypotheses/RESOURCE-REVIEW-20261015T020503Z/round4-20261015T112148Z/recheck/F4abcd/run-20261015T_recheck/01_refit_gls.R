# Independent re-derivation of the H59 GLS numbers used in workstream F4abcd.
# Fits score ~ stage + dataset with varIdent(~1|dataset), REML, on the 60-donor
# universe (GSE189600 dropped); compares to H59 results_primary_family.tsv.
# No permutation is redone here (19,999 x 117 is a compute-node job); the
# permutation critical |t| is re-derived from the saved null t matrix.
suppressPackageStartupMessages({library(data.table); library(nlme)})
set.seed(20260923)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT <- file.path(ROOT, "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F4abcd/run-20261015T_recheck")
STAGES <- c("Healthy", "Steatosis", "Steatohepatitis")
sc <- fread(file.path(ROOT, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"))
meta <- fread(file.path(ROOT, "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"))
meta[, excl := tolower(as.character(exclude_stage_analysis)) %in% c("true", "t", "1")]
elig <- meta[!excl & disease_stage_coarse %in% STAGES, .(donor = sample, dataset, stg = disease_stage_coarse)]
elig[, stage := match(stg, STAGES) - 1]
d <- merge(sc[, .(cell_type, module, donor, score)], elig, by = "donor")
d <- d[is.finite(score)]
cat("donors in 62 roster:", uniqueN(d$donor), "\n")
d60 <- d[dataset != "GSE189600"]
cat("donors in 60 roster:", uniqueN(d60$donor), "\n")
# design table by lineage x dataset x stage (60 roster)
des <- dcast(d60[, .N, by = .(cell_type, dataset, stg, module)][, .(N = max(N)), by = .(cell_type, dataset, stg)],
             cell_type + dataset ~ stg, value.var = "N", fill = 0)
fwrite(des, file.path(OUT, "design_lineage_dataset_stage_60.tsv"), sep = "\t")
print(des)
# multi-stage datasets per lineage
ms <- d60[, .(nst = uniqueN(stage)), by = .(cell_type, module, dataset)]
fitd <- merge(d60, ms[nst > 1], by = c("cell_type", "module", "dataset"))
res <- fitd[, {
  df <- data.frame(y = score, stage = stage, dataset = factor(dataset))
  f <- gls(y ~ stage + dataset, data = df, weights = varIdent(form = ~1 | dataset), method = "REML",
           control = glsControl(maxIter = 200, msMaxIter = 200))
  b <- unname(coef(f)["stage"]); se <- sqrt(vcov(f)["stage", "stage"])
  ci <- tryCatch(intervals(f, which = "coef")$coef["stage", c(1, 3)], error = function(e) c(NA, NA))
  sl <- sapply(levels(df$dataset), function(ds) { s <- summary(lm(y ~ stage, df[df$dataset == ds, ]))$coef; s["stage", 1:2] })
  w <- 1 / sl[2, ]^2; pooled <- sum(w * sl[1, ]) / sum(w); Q <- sum(w * (sl[1, ] - pooled)^2)
  .(n_fit = .N, n_ds = nlevels(df$dataset), beta = b, se = se, lo = ci[1], hi = ci[2],
    slope1 = sl[1, 1], slope2 = sl[1, 2], Qp = pchisq(Q, ncol(sl) - 1, lower.tail = FALSE))
}, by = .(cell_type, module)]
h <- fread(file.path(ROOT, "Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/results_primary_family.tsv"))
m <- merge(res, h[, .(cell_type, module, h_n = n_fit, h_beta = gls_beta, h_se = gls_se, h_lo = gls_ci_low, h_hi = gls_ci_high,
                      h_s1 = slope_GSE202379, h_s2 = slope_GSE244832, h_Qp = cochran_Q_p, perm_p, perm_q, module_name)],
           by = c("cell_type", "module"))
cat("programs compared:", nrow(m), "\n")
cat("max|dn|", max(abs(m$n_fit - m$h_n)), " max|dbeta|", max(abs(m$beta - m$h_beta)), " max|dse|", max(abs(m$se - m$h_se)),
    " max|dlo|", max(abs(m$lo - m$h_lo)), " max|dhi|", max(abs(m$hi - m$h_hi)),
    " max|dslope|", max(abs(c(m$slope1 - m$h_s1, m$slope2 - m$h_s2))), " max|dQp|", max(abs(m$Qp - m$h_Qp)), "\n")
print(m[, .(n_fit = unique(n_fit), n_ds = unique(n_ds)), by = cell_type])
print(m[(cell_type == "hepatocytes" & module %in% c(20, 8)), .(module, module_name, n_fit, beta, lo, hi, slope1, slope2, Qp, perm_p, perm_q)])
fwrite(m, file.path(OUT, "gls_refit_vs_H59.tsv"), sep = "\t")
# Heatmap: dataset-centered donor score, median by stage, using the two multi-stage datasets
hm <- fitd[, cs := score - mean(score), by = .(cell_type, module, dataset)][, .(med = median(cs)), by = .(cell_type, module, stg)]
hm <- dcast(hm, cell_type + module ~ stg, value.var = "med")
fwrite(hm, file.path(OUT, "heatmap_medians_refit.tsv"), sep = "\t")
# permutation critical |t| from saved null GLS t
nt <- fread(file.path(ROOT, "Analysis/Hypotheses/H59-20261015T112432Z/out/null_gls_t.tsv.gz"))
cat("null t dims", dim(nt), "\n"); print(names(nt)[1:5])
fwrite(data.table(sessionInfo = capture.output(sessionInfo())), file.path(OUT, "sessionInfo_01.txt"))
