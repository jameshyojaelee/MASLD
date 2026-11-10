# U8a: independent observed fits for H59 (Fig 4C) and descriptive numbers.
# Usage: Rscript --vanilla 01_observed.R <outdir>
args <- commandArgs(trailingOnly = TRUE)
here <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U8a/rederive"
source(file.path(here, "u8a_lib.R"))
out <- if (length(args)) args[1] else file.path(here, "out/observed")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(nlme))

u62 <- load_donor_scores()
u60 <- u62[dataset != "GSE189600"]
fit <- u60[dataset %in% FIT_DS]
say("60-universe donors %d; fitted donors %d", uniqueN(u60$donor), uniqueN(fit$donor))
print(dcast(unique(u62[, .(donor, dataset, stage_label)]), dataset ~ stage_label, fun.aggregate = length))

# donors per lineage in the fitted set
lin_n <- fit[, .(n = uniqueN(donor)), by = .(cell_type, program)][, .(n_min = min(n), n_max = max(n)), by = cell_type]
print(lin_n)

reg <- fread(REGISTRY_FILE)
reg[, program := paste0(cell_type, "_", module)]
reg <- reg[, .(program, module_name, stability_median)]
say("registry rows %d; stability_median >= 0.5: %d", nrow(reg), sum(reg$stability_median >= 0.5))

progs <- sort(unique(fit$program))
stopifnot(length(progs) == 117L)

res <- rbindlist(lapply(progs, function(p) {
  d <- fit[program == p]
  g <- gls_fit(d$score, d$stage, d$dataset)
  nl <- gls(score ~ stage + dataset, data = d, weights = varIdent(form = ~ 1 | dataset),
            method = "REML", control = glsControl(maxIter = 200, msMaxIter = 200))
  tt <- summary(nl)$tTable
  # per-dataset slopes
  sl <- lapply(FIT_DS, function(l) {
    e <- d[dataset == l]
    m <- summary(lm(score ~ stage, data = e))$coefficients
    c(b = m["stage", 1], se = m["stage", 2], p = m["stage", 4], n = nrow(e))
  })
  w <- 1 / c(sl[[1]]["se"], sl[[2]]["se"])^2
  bb <- c(sl[[1]]["b"], sl[[2]]["b"])
  biv <- sum(w * bb) / sum(w)
  Q <- sum(w * (bb - biv)^2)
  # cell-mean-centred SD per dataset (D2)
  sdcm <- sapply(FIT_DS, function(l) {
    e <- d[dataset == l]
    r <- e$score - ave(e$score, e$stage)
    sqrt(sum(r^2) / (nrow(e) - uniqueN(e$stage)))
  })
  sxx <- sapply(FIT_DS, function(l) { s <- d[dataset == l]$stage; sum((s - mean(s))^2) })
  sxy <- sapply(FIT_DS, function(l) { e <- d[dataset == l]; sum((e$stage - mean(e$stage)) * (e$score - mean(e$score))) })
  sig <- c(g$sigma1, g$sigma2)
  data.table(program = p, cell_type = d$cell_type[1], n_fit = g$n,
    beta = g$beta, se = g$se, t = g$t, sigma_202379 = g$sigma1, sigma_244832 = g$sigma2,
    nlme_beta = tt["stage", 1], nlme_t = tt["stage", 3],
    param_p = 2 * pt(-abs(g$t), g$n - 3),
    ci_lo = g$beta - qt(0.975, g$n - 3) * g$se, ci_hi = g$beta + qt(0.975, g$n - 3) * g$se,
    slope_202379 = sl[[1]]["b"], se_202379 = sl[[1]]["se"], p_202379 = sl[[1]]["p"], n_202379 = sl[[1]]["n"],
    slope_244832 = sl[[2]]["b"], se_244832 = sl[[2]]["se"], p_244832 = sl[[2]]["p"], n_244832 = sl[[2]]["n"],
    Q = Q, Q_p = pchisq(Q, 1, lower.tail = FALSE),
    sd_ratio_cm = sdcm[2] / sdcm[1], sd_ratio_gls = g$sigma2 / g$sigma1,
    info_gls_d2_202379 = (sxx[1] / sdcm[1]^2) / sum(sxx / sdcm^2),
    info_gls_fit_202379 = (sxx[1] / sig[1]^2) / sum(sxx / sig^2),
    ols_num_share_244832 = sxy[2] / sum(sxy))
}))

# OLS families on 62 and 60 donors (all datasets, dataset intercepts)
ols_family <- function(dat, tag) {
  rbindlist(lapply(progs, function(p) {
    d <- dat[program == p]
    o <- ols_fit(d$score, d$stage, d$dataset)
    data.table(program = p, n = nrow(d), beta = o$beta, t = o$t, df = o$df,
      p = 2 * pt(-abs(o$t), o$df), p_hc3 = 2 * pt(-abs(o$beta / o$se_hc3), o$df), maxlev = o$maxlev)
  }))[, `:=`(q = p.adjust(p, "BH", n = 117), q_hc3 = p.adjust(p_hc3, "BH", n = 117))][, set := tag][]
}
o62 <- ols_family(u62, "ols62")
o60 <- ols_family(u60, "ols60")

res <- merge(res, reg, by = "program", all.x = TRUE)
stopifnot(nrow(res) == 117L, !anyNA(res$stability_median))
fwrite(res, file.path(out, "observed_gls_117.tsv"), sep = "\t")
fwrite(rbind(o62, o60), file.path(out, "ols_62_60.tsv"), sep = "\t")

say("max |own t - nlme t| = %.3g; max |own beta - nlme beta| = %.3g",
    max(abs(res$t - res$nlme_t)), max(abs(res$beta - res$nlme_beta)))
say("opposite-sign dataset slopes: %d", sum(sign(res$slope_202379) != sign(res$slope_244832)))
say("Cochran Q BH smallest q: %.4f; n q<0.05 %d", min(p.adjust(res$Q_p, "BH")), sum(p.adjust(res$Q_p, "BH") < 0.05))
say("SD ratio (cell-mean) median %.3f range %.3f-%.3f", median(res$sd_ratio_cm), min(res$sd_ratio_cm), max(res$sd_ratio_cm))
say("SD ratio (GLS fit)   median %.3f range %.3f-%.3f", median(res$sd_ratio_gls), min(res$sd_ratio_gls), max(res$sd_ratio_gls))
named <- c("hepatocytes_20", "hepatocytes_8")
print(res[program %in% named, .(program, module_name, n_fit, beta, se, t, ci_lo, ci_hi, slope_202379, slope_244832,
  p_244832, Q_p, sd_ratio_cm, sd_ratio_gls, info_gls_d2_202379, info_gls_fit_202379, ols_num_share_244832)])
print(o62[program %in% named]); print(o60[program %in% named])
say("ols62 calls q<0.05: %d (HC3 %d); ols60 calls %d", sum(o62$q < 0.05), sum(o62$q_hc3 < 0.05, na.rm = TRUE), sum(o60$q < 0.05))
say("ols62 HC3 NA rows: %d", sum(is.na(o62$p_hc3)))

# S4A rows: 62-donor OLS HC3 for the two named programs; per-dataset and omission fits
s4a <- rbindlist(lapply(named, function(p) {
  d <- u62[program == p]
  rows <- list(pooled = d, GSE202379 = d[dataset == "GSE202379"], GSE244832 = d[dataset == "GSE244832"])
  for (l in sort(unique(d$dataset))) rows[[paste0("omit_", l)]] <- d[dataset != l]
  rbindlist(lapply(names(rows), function(r) {
    e <- rows[[r]]; o <- ols_fit(e$score, e$stage, e$dataset)
    data.table(program = p, row = r, n = nrow(e), beta = o$beta,
      hc3_lo = o$beta - qt(0.975, o$df) * o$se_hc3, hc3_hi = o$beta + qt(0.975, o$df) * o$se_hc3)
  }))
}))
fwrite(s4a, file.path(out, "s4a_rows.tsv"), sep = "\t")
print(s4a)

# Compare with the H59 producer outputs
pr <- fread(file.path(H59, "out/outcome/results_primary_family.tsv"))
cmp <- merge(res[, .(program, beta, se, t, n_fit)], pr[, .(program, p_beta = gls_beta, p_se = gls_se, p_t = gls_t, p_n = n_fit)], by = "program")
say("vs producer: max|dbeta| %.3g max|dse| %.3g max|dt| %.3g n mismatches %d", max(abs(cmp$beta - cmp$p_beta)),
    max(abs(cmp$se - cmp$p_se)), max(abs(cmp$t - cmp$p_t)), sum(cmp$n_fit != cmp$p_n))
cmp2 <- merge(o62[, .(program, q, q_hc3)], pr[, .(program, ols62_q, ols62_q_hc3)], by = "program")
say("ols62 vs producer: max|dq| %.3g max|dq_hc3| %.3g", max(abs(cmp2$q - cmp2$ols62_q)), max(abs(cmp2$q_hc3 - cmp2$ols62_q_hc3), na.rm = TRUE))
cmp3 <- merge(o60[, .(program, q)], pr[, .(program, ols60_q)], by = "program")
say("ols60 vs producer: max|dq| %.3g", max(abs(cmp3$q - cmp3$ols60_q)))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
system2("sha256sum", c(SCORE_FILE, META_FILE, REGISTRY_FILE, file.path(H59, "out/outcome/results_primary_family.tsv")),
        stdout = file.path(out, "inputs.sha256"))
