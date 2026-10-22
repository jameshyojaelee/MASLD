# Fig 4C and 4E independent refits from upstream analysis outputs (light models only).
# 4C: H59 GLS (score ~ stage + dataset, varIdent by dataset, REML) on GSE202379 + GSE244832
#     informative donors, from frozen donor scores and 2026-10-14 corrected stages.
# 4E: three-cohort REML random-effects meta with Knapp-Hartung from per-cohort HC3 estimates;
#     BH over 224 (112 programs x 2 axes). Also REML z and KH floored at REML variants.
suppressPackageStartupMessages({library(data.table); library(nlme); library(metafor)})
set.seed(20261017)
root <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out <- file.path(root, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute")

## ---------------- 4C
S <- fread(file.path(root, "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/donor_program_scores_primary.tsv"))
M <- fread(file.path(root, "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"))
S <- merge(S, M[, .(donor = sample, stage_corr = disease_stage_coarse, excl_corr = exclude_stage_analysis)], by = "donor")
S <- S[dataset %in% c("GSE202379", "GSE244832") & stage_corr %in% c("Healthy", "Steatosis", "Steatohepatitis") & excl_corr == FALSE]
S[, stage := match(stage_corr, c("Healthy", "Steatosis", "Steatohepatitis")) - 1]
fit_one <- function(d) {
  d <- as.data.frame(d)
  fit <- tryCatch(gls(score ~ stage + dataset, data = d, weights = varIdent(form = ~ 1 | dataset), method = "REML",
                      control = glsControl(maxIter = 200, msMaxIter = 200)), error = function(e) e)
  if (inherits(fit, "error")) return(data.table(cell_type = d$cell_type[1], module = d$module[1], n = nrow(d), beta = NA_real_, se = NA_real_,
                                                ci_t_lo = NA_real_, ci_t_hi = NA_real_, ci_z_lo = NA_real_, ci_z_hi = NA_real_, df = NA_real_, err = conditionMessage(fit)))
  tt <- summary(fit)$tTable; b <- tt["stage", 1]; se <- tt["stage", 2]; df <- fit$dims$N - fit$dims$p
  data.table(cell_type = d$cell_type[1], module = d$module[1], n = nrow(d), beta = b, se = se,
             ci_t_lo = b - qt(.975, df) * se, ci_t_hi = b + qt(.975, df) * se,
             ci_z_lo = b - qnorm(.975) * se, ci_z_hi = b + qnorm(.975) * se, df = df, err = NA_character_)
}
res4c <- rbindlist(lapply(split(S, by = c("cell_type", "module")), fit_one))
cat("4C fits failed:", sum(is.na(res4c$beta)), "\n"); if (any(is.na(res4c$beta))) print(unique(res4c$err))
H <- fread(file.path(root, "Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/results_primary_family.tsv"))
cmp <- merge(res4c, H[, .(cell_type, module, n_fit, gls_beta, gls_se, gls_ci_low, gls_ci_high)], by = c("cell_type", "module"))
cat("4C refit: programs", nrow(cmp), "; max|dbeta|", max(abs(cmp$beta - cmp$gls_beta), na.rm = TRUE), "; max|dse|", max(abs(cmp$se - cmp$gls_se), na.rm = TRUE),
    "; max|dn|", max(abs(cmp$n - cmp$n_fit)), "\n")
cat("4C CI: max|d| t-CI", max(abs(cmp$ci_t_lo - cmp$gls_ci_low), na.rm = TRUE), "; z-CI", max(abs(cmp$ci_z_lo - cmp$gls_ci_low), na.rm = TRUE), "\n")
fwrite(cmp, file.path(out, "r4c_gls_refit_vs_H59.tsv"), sep = "\t")

## ---------------- 4E
D <- file.path(root, "figures/candidates/fig4e-program-axis-774-20261012T010830Z/analysis/new774_three_cohorts")
ce <- fread(file.path(D, "cohort_effects.tsv"))[model_kind == "primary" & estimable == TRUE]
cat("4E cohort rows", nrow(ce), "; programs", uniqueN(ce$feature_id), "; cohorts", paste(unique(ce$cohort), collapse = ","), "\n")
ce_n <- unique(ce[, .(cohort, n)]); print(ce_n)
meta_one <- function(d) {
  f_kh <- rma(yi = d$beta, sei = d$se_hc3, method = "REML", test = "knha", control = list(stepadj = 0.5, maxiter = 1000))
  f_z <- rma(yi = d$beta, sei = d$se_hc3, method = "REML", control = list(stepadj = 0.5, maxiter = 1000))
  # KH with the variance multiplier floored at 1 (not allowed below the REML z variance)
  s2 <- max(1, f_kh$se^2 / f_z$se^2)
  se_fl <- f_z$se * sqrt(s2); t_fl <- f_z$b[1] / se_fl
  data.table(feature_id = d$feature_id[1], axis = d$axis[1], k = nrow(d), beta = f_kh$b[1], se_kh = f_kh$se, p_kh = f_kh$pval,
             se_z = f_z$se, p_z = f_z$pval, p_khfloor = 2 * pt(-abs(t_fl), df = nrow(d) - 1), tau2 = f_kh$tau2, QE = f_kh$QE,
             n_same = max(sum(d$beta > 0), sum(d$beta < 0)), n_units = sum(d$n))
}
meta <- rbindlist(lapply(split(ce, by = c("feature_id", "axis")), meta_one))
meta[, q_kh := p.adjust(p_kh, "BH")]
meta[, q_z := p.adjust(p_z, "BH")]
meta[, q_khfloor := p.adjust(p_khfloor, "BH")]
meta[, mlq := -log10(q_kh)]
cat("4E: rows", nrow(meta), "; calls q<0.05 KH by axis:\n"); print(meta[, .(n = .N, kh = sum(q_kh < .05), z = sum(q_z < .05), khfloor = sum(q_khfloor < .05),
                                                                   min_q = min(q_kh), max_mlq = max(mlq), max_abs_beta = max(abs(beta))), by = axis])
ref <- fread(file.path(D, "meta_effects.tsv"))[model_kind == "primary"]
cm <- merge(meta, ref[, .(feature_id, axis, beta_meta, se_meta, q_value)], by = c("feature_id", "axis"))
cat("4E vs analysis meta_effects: max|dbeta|", max(abs(cm$beta - cm$beta_meta)), "; max|dse|", max(abs(cm$se_kh - cm$se_meta)),
    "; max|dq|", max(abs(cm$q_kh - cm$q_value)), "\n")
fwrite(meta, file.path(out, "r4e_meta_recomputed.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_r4ce.txt"))

## ---------------- 4E per-cohort refit from participant program scores (independent of cohort_effects.tsv)
H <- fread(file.path(root, "RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z/programs/hotspot/hotspot_participant_scores.tsv.gz"))
man <- fread(file.path(root, "figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv"), na.strings = "NA")
man <- man[dataset %in% c("GSE130970", "GSE135251", "GSE162694") & !is.na(fibrosis_stage) & !is.na(nas_score)]
cat("4E roster from manifest:\n"); print(man[, .N, by = dataset])
hc3_fit <- function(y, X) {
  XtXi <- solve(crossprod(X)); b <- XtXi %*% crossprod(X, y); e <- as.vector(y - X %*% b)
  h <- rowSums((X %*% XtXi) * X); meat <- crossprod(X * (e / (1 - h)))
  V <- XtXi %*% meat %*% XtXi; list(b = as.vector(b), se = sqrt(diag(V)))
}
pc <- list()
for (sc in c("outcome_z", "score_raw")) {
  for (ds in unique(man$dataset)) {
    m <- man[dataset == ds]
    m[, `:=`(fib_c = fibrosis_stage - mean(fibrosis_stage), nas_c = nas_score - mean(nas_score), male = as.numeric(inferred_sex == "M"))]
    X <- cbind(1, m$fib_c, m$nas_c, m$male)
    hh <- H[dataset == ds & sample_id %in% m$sample_id]
    W <- dcast(hh, sample_id ~ program_uid, value.var = sc)
    W <- W[match(m$sample_id, W$sample_id)]
    stopifnot(nrow(W) == nrow(m), all(W$sample_id == m$sample_id))
    for (pu in setdiff(names(W), "sample_id")) {
      y <- W[[pu]]; if (any(!is.finite(y)) || sd(y) == 0) next
      f <- hc3_fit(y, X)
      pc[[length(pc) + 1]] <- data.table(score = sc, cohort = ds, feature_id = pu, n = nrow(m),
                                         b_fib = f$b[2], se_fib = f$se[2], b_nas = f$b[3], se_nas = f$se[3])
    }
  }
}
pc <- rbindlist(pc)
ce2 <- fread(file.path(D, "cohort_effects.tsv"))[model_kind == "primary"]
for (sc in c("outcome_z", "score_raw")) {
  x <- merge(pc[score == sc], ce2[axis == "fibrosis", .(cohort, feature_id, beta, se_hc3)], by = c("cohort", "feature_id"))
  cat("4E per-cohort refit (", sc, ") fibrosis: rows", nrow(x), "; max|dbeta|", max(abs(x$b_fib - x$beta)), "; max|dse|", max(abs(x$se_fib - x$se_hc3)), "\n")
}
fwrite(pc, file.path(out, "r4e_cohort_refit_hc3.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_r4ce.txt"))
