#!/usr/bin/env Rscript
# H35 step 2: Alu A-to-I editing index versus fibrosis stage (see PRESPEC_H35.md).
# Usage: 02_analysis.R PERSON_TABLE.tsv OUT_DIR N_PERM
# Everything below follows PRESPEC_H35.md; no model choice depends on the results.
suppressPackageStartupMessages({
  library(data.table)
  library(aod)
  library(metafor)
  library(parallel)
  library(numDeriv)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "lib_bb.R"))
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 3L) stop("Usage: 02_analysis.R person_table out_dir n_perm")
person_file <- argv[[1]]; out <- argv[[2]]; n_perm <- as.integer(argv[[3]])
if (dir.exists(out)) stop("Refusing to overwrite output dir: ", out)
dir.create(out, recursive = TRUE)
SEED <- 20260923L
set.seed(SEED)
n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
COMP <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz"
PRIMARY <- c("GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE240729")
ISG <- c("ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "MX1", "MX2",
         "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "XAF1", "USP18", "HERC5", "CMPK2")
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

# ------------------------------------------------------------------ data
d <- fread(person_file, na.strings = c("NA", ""))
for (cn in c("stage", "alu_ref", "alu_alt", "alu_sites_cov10", "alu_ts_ref", "alu_ts_alt", "other_base_alu",
             "age", "rpg_libsize")) d[, (cn) := suppressWarnings(as.numeric(get(cn)))]
d[, sex := fifelse(inferred_sex %in% c("M", "F"), inferred_sex,
                   fifelse(expr_sex %in% c("M", "F"), expr_sex, "U"))]
d[, depth := log10(alu_ref + alu_alt)]
d[, n_reads := alu_ref + alu_alt]
d[, disease := fifelse(group_binary == "Disease", 1, fifelse(group_binary == "Control", 0, NA_real_))]
# ISG score: mean of cohort-wise z-scored log2CPM; zero-variance genes dropped per cohort
isg_cols <- paste0("lcpm_", ISG)
dropped_isg <- list()
d[, isg := NA_real_]
for (co in unique(d$cohort)) {
  ix <- which(d$cohort == co)
  m <- as.matrix(d[ix, ..isg_cols])
  sdv <- apply(m, 2, sd)
  keep <- is.finite(sdv) & sdv > 0
  dropped_isg[[co]] <- ISG[!keep]
  z <- scale(m[, keep, drop = FALSE])
  set(d, ix, "isg", rowMeans(z))
}
# exclusions (PRESPEC: fewer than 1,000 Alu sites at depth >= 10)
excl <- d[alu_sites_cov10 < 1000]
d <- d[alu_sites_cov10 >= 1000]
fwrite(data.table(excluded_low_coverage = nrow(excl)), file.path(out, "exclusions.tsv"), sep = "\t")
comp <- fread(cmd = paste("zcat", shQuote(COMP)), select = c("sample_id", "Hepatocytes"))
d <- merge(d, comp, by = "sample_id", all.x = TRUE, sort = FALSE)
d[, hep_logit := qlogis(pmin(pmax(Hepatocytes, 0.001), 0.999))]

# sex levels with < 5 persons in a cohort are merged into U (within analysis data)
prep_cohort <- function(dd) {
  dd <- copy(dd)
  tb <- table(dd$sex)
  small <- names(tb)[tb < 5]
  dd[sex %in% small, sex := "U"]
  dd[, sex := factor(sex)]
  dd
}

# ------------------------------------------------------------------ model helpers
dl_meta <- function(est, se) {
  ok <- is.finite(est) & is.finite(se) & se > 0
  est <- est[ok]; se <- se[ok]; k <- length(est)
  if (k < 2L) return(list(est = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_, tau2 = NA_real_, I2 = NA_real_,
                          Q = NA_real_, Qp = NA_real_, k = k))
  w <- 1 / se^2
  mu_fe <- sum(w * est) / sum(w)
  Q <- sum(w * (est - mu_fe)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (se^2 + tau2)
  mu <- sum(ws * est) / sum(ws)
  sem <- sqrt(1 / sum(ws))
  z <- mu / sem
  list(est = mu, se = sem, z = z, p = 2 * pnorm(-abs(z)), tau2 = tau2,
       I2 = max(0, (Q - (k - 1)) / Q) * 100, Q = Q, Qp = pchisq(Q, k - 1, lower.tail = FALSE), k = k)
}
terms_primary <- function(dd, extra = character()) {
  tm <- c("stage")
  if (nlevels(droplevels(dd$sex)) > 1L) tm <- c(tm, "sex")
  c(tm, "lcpm_ADAR", "lcpm_ADARB1", "isg", "depth", extra)
}
# one fit; method in betabin (own ML, lib_bb.R) | aod (aod::betabin as prespecified) | quasibin | ols
# returns coefficient table and metadata. betabin and aod fall back to quasibin when they do not converge.
fit_one <- function(dd, terms, y_succ = "alu_alt", y_fail = "alu_ref", method = "betabin", maxit = 2000L, center = TRUE) {
  dd <- dd[complete.cases(dd[, c(y_succ, y_fail, setdiff(terms, "sex")), with = FALSE])]
  if ("sex" %in% terms) dd <- dd[!is.na(sex)]
  dd <- droplevels(dd)
  if ("sex" %in% terms && nlevels(dd$sex) < 2L) terms <- setdiff(terms, "sex")
  # Centering numeric covariates re-parameterises the intercept only: slopes, SEs and the likelihood are unchanged.
  if (center) for (tt in setdiff(terms, "sex")) set(dd, j = tt, value = dd[[tt]] - mean(dd[[tt]]))
  f <- as.formula(paste0("cbind(", y_succ, ", ", y_fail, ") ~ ", paste(terms, collapse = " + ")))
  used <- method; note <- ""
  tab <- NULL; logL <- NA_real_; bb_code <- NA_integer_; phi <- NA_real_
  if (method == "betabin") {
    X <- model.matrix(reformulate(terms), data = as.data.frame(dd))
    yv <- dd[[y_succ]]; nv <- yv + dd[[y_fail]]
    bf <- tryCatch(bb_fit(X, yv, nv), error = function(e) list(converged = FALSE, why = conditionMessage(e)))
    if (isTRUE(bf$converged)) {
      k <- ncol(X)
      tab <- cbind(Estimate = bf$par[seq_len(k)], `Std. Error` = bf$se[seq_len(k)]); rownames(tab) <- colnames(X)
      logL <- bf$logL; bb_code <- 0L; phi <- bf$rho
    } else {
      note <- paste("betabin ML not converged:", bf$why); used <- "quasibin"; bb_code <- 1L
    }
  }
  if (method == "aod") {
    w <- character()
    fit <- withCallingHandlers(
      tryCatch(aod::betabin(f, ~ 1, data = as.data.frame(dd), link = "logit", control = list(maxit = maxit)), error = function(e) e),
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    if (inherits(fit, "error")) {
      note <- paste("aod betabin error:", conditionMessage(fit)); used <- "quasibin"
    } else {
      tab <- as.matrix(summary(fit)@Coef)
      logL <- fit@logL; bb_code <- as.integer(fit@code)
      phi <- tryCatch(as.numeric(fit@random.param[1]), error = function(e) NA_real_)
      if (fit@code != 0 || any(!is.finite(tab[, "Std. Error"]))) {
        note <- paste("aod betabin code", fit@code, paste(w, collapse = "; ")); used <- "quasibin"; tab <- NULL
      } else {
        note <- paste(w, collapse = "; "); used <- "aod_betabin"
      }
    }
  }
  if (is.null(tab) && (method == "quasibin" || used == "quasibin")) {
    g <- glm(f, family = quasibinomial(), data = as.data.frame(dd))
    tab <- coef(summary(g))[, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"), drop = FALSE]
    used <- "quasibin"
  }
  if (method == "ols") {
    dd[, y := qlogis((get(y_succ) + 0.5) / (get(y_succ) + get(y_fail) + 1))]
    g <- lm(as.formula(paste("y ~", paste(terms, collapse = " + "))), data = as.data.frame(dd))
    tab <- coef(summary(g))[, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"), drop = FALSE]
    used <- "ols_logit"
  }
  list(tab = tab, n = nrow(dd), used = used, note = gsub("\\s+", " ", note), terms = terms,
       logL = logL, code = bb_code, phi = phi)
}
row_of <- function(res, term, cohort, label) {
  if (is.null(res$tab) || !term %in% rownames(res$tab)) {
    return(data.table(cohort = cohort, analysis = label, term = term, n = res$n, estimate = NA_real_, se = NA_real_,
                      z = NA_real_, p = NA_real_, method = res$used, note = res$note,
                      logL = res$logL, bb_code = res$code, phi = res$phi))
  }
  e <- res$tab[term, "Estimate"]; s <- res$tab[term, "Std. Error"]
  data.table(cohort = cohort, analysis = label, term = term, n = res$n, estimate = e, se = s, z = e / s,
             p = 2 * pnorm(-abs(e / s)), method = res$used, note = res$note,
             logL = res$logL, bb_code = res$code, phi = res$phi)
}

# ------------------------------------------------------------------ primary cohorts
dp <- lapply(PRIMARY, function(co) prep_cohort(d[cohort == co & stage_resolution == "exact" & !is.na(stage)]))
names(dp) <- PRIMARY
cohort_n <- rbindlist(lapply(PRIMARY, function(co) {
  x <- dp[[co]]
  data.table(cohort = co, n_persons = nrow(x), stages = paste(sort(unique(x$stage)), collapse = ","),
             median_sites_cov10 = median(x$alu_sites_cov10), median_alu_reads = median(x$n_reads),
             sex_levels = paste(names(table(x$sex)), table(x$sex), sep = ":", collapse = " "),
             isg_genes_dropped = paste(dropped_isg[[co]], collapse = ","))
}))
fwrite(cohort_n, file.path(out, "primary_cohort_summary.tsv"), sep = "\t")
log_line("primary persons:", sum(cohort_n$n_persons))

run_family <- function(label, y_succ = "alu_alt", y_fail = "alu_ref", method = "betabin", terms_fun = terms_primary,
                       cohorts = PRIMARY, data_list = dp, target = "stage", extra = character(), transform = identity,
                       maxit = 2000L, center = TRUE) {
  rows <- list()
  for (co in cohorts) {
    dd <- transform(data_list[[co]])
    if (is.null(dd) || nrow(dd) < 20L) next
    tm <- terms_fun(dd, extra)
    res <- fit_one(dd, tm, y_succ, y_fail, method, maxit, center)
    rows[[co]] <- rbindlist(lapply(unique(c(target, "lcpm_ADAR")), function(tt) row_of(res, tt, co, label)))
  }
  rbindlist(rows)
}
meta_of <- function(tab, term_name, label) {
  x <- tab[tab$term == term_name & !is.na(tab$estimate)]
  m <- dl_meta(x$estimate, x$se)
  data.table(analysis = label, term = term_name, k = m$k, estimate = m$est, se = m$se, z = m$z, p = m$p,
             OR_per_unit = exp(m$est), OR_lo = exp(m$est - 1.96 * m$se), OR_hi = exp(m$est + 1.96 * m$se),
             tau2 = m$tau2, I2 = m$I2, Q = m$Q, Q_p = m$Qp)
}

log_line("primary fits")
t0 <- Sys.time()
prim <- run_family("primary_betabin")
log_line("primary done in", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s")
fwrite(prim, file.path(out, "primary_per_cohort.tsv"), sep = "\t")
meta_stage <- meta_of(prim, "stage", "primary_betabin")
meta_adar <- meta_of(prim, "lcpm_ADAR", "primary_betabin")
# cross-check of the hand-written DL against metafor
chk <- prim[term == "stage" & !is.na(estimate)]
rm_chk <- rma(yi = chk$estimate, sei = chk$se, method = "DL")
stopifnot(abs(as.numeric(rm_chk$beta) - meta_stage$estimate) < 1e-8, abs(rm_chk$tau2 - meta_stage$tau2) < 1e-8)

# positive control: ADAR coefficient, one-sided
pos_p <- pnorm(meta_adar$z, lower.tail = FALSE)
pos_pass <- is.finite(meta_adar$estimate) && meta_adar$estimate > 0 && pos_p < 0.05

# negative control: non-A-to-I mismatch rate at the same sites
neg_tab <- {
  rows <- list()
  for (co in PRIMARY) {
    dd <- copy(dp[[co]])
    dd[, other_n := other_base_alu][, ref_n := alu_ref + alu_alt]
    res <- fit_one(dd, terms_primary(dd), "other_n", "ref_n", "betabin")
    rows[[co]] <- row_of(res, "stage", co, "negative_control_other_base")
  }
  rbindlist(rows)
}
fwrite(neg_tab, file.path(out, "negative_control_per_cohort.tsv"), sep = "\t")
meta_neg <- meta_of(neg_tab, "stage", "negative_control_other_base")
neg_pass <- is.finite(meta_neg$p) && meta_neg$p >= 0.05

# ------------------------------------------------------------------ permutation (stage within cohort)
set.seed(SEED)
perm_list <- lapply(PRIMARY, function(co) {
  n <- nrow(dp[[co]])
  sapply(seq_len(n_perm), function(b) sample.int(n))
})
names(perm_list) <- PRIMARY
perm_one <- function(b) {
  est <- se <- setNames(rep(NA_real_, length(PRIMARY)), PRIMARY)
  n_fallback <- 0L
  for (co in PRIMARY) {
    dd <- copy(dp[[co]])
    dd[, stage := stage[perm_list[[co]][, b]]]
    r <- tryCatch(fit_one(dd, terms_primary(dd), "alu_alt", "alu_ref", "betabin"), error = function(e) NULL)
    if (!is.null(r) && !is.null(r$tab) && "stage" %in% rownames(r$tab)) {
      est[co] <- r$tab["stage", "Estimate"]; se[co] <- r$tab["stage", "Std. Error"]
      if (r$used != "betabin") n_fallback <- n_fallback + 1L
    }
  }
  m <- dl_meta(est, se)
  c(z = m$z, k = m$k, n_fallback = n_fallback, setNames(est / se, paste0("z_", PRIMARY)))
}
perm_res <- NULL
if (n_perm > 0L) {
  log_line("permutations:", n_perm, "on", n_cores, "cores")
  t1 <- Sys.time()
  perm_out <- mclapply(seq_len(n_perm), perm_one, mc.cores = n_cores, mc.set.seed = FALSE)
  bad <- vapply(perm_out, function(x) inherits(x, "try-error") || length(x) != 3L + length(PRIMARY), logical(1))
  if (any(bad)) log_line("permutations that failed:", sum(bad))
  perm_out[bad] <- list(c(z = NA_real_, k = NA_real_, n_fallback = NA_real_, setNames(rep(NA_real_, length(PRIMARY)), paste0("z_", PRIMARY))))
  perm_res <- as.data.table(do.call(rbind, perm_out))
  log_line("permutations done in", round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1), "s")
  fwrite(perm_res, file.path(out, "permutation_meta_z.tsv"), sep = "\t")
}
perm_p <- if (!is.null(perm_res)) (1 + sum(abs(perm_res$z) >= abs(meta_stage$z), na.rm = TRUE)) / (1 + sum(!is.na(perm_res$z))) else NA_real_

# ------------------------------------------------------------------ Family B: per-cohort slopes, BH over 5
famB <- prim[term == "stage"]
famB[, p_BH := p.adjust(p, method = "BH")]
fwrite(famB, file.path(out, "familyB_per_cohort_stage_BH.tsv"), sep = "\t")

# ------------------------------------------------------------------ Family C: secondary cohorts and contrasts
mk <- function(co, filt) prep_cohort(d[cohort == co][filt(d[cohort == co])])
c_rows <- list()
d_prj <- mk("PRJNA512027", function(x) x$stage_resolution == "exact" & x$stage %in% c(3, 4))
res <- fit_one(d_prj, terms_primary(d_prj), method = "betabin")
c_rows$prj <- row_of(res, "stage", "PRJNA512027", "familyC_F4_vs_F3")
d_213 <- mk("GSE213621", function(x) !is.na(x$stage))
res <- fit_one(d_213, terms_primary(d_213), method = "betabin")
c_rows$g213 <- row_of(res, "stage", "GSE213621", "familyC_group_code_slope")
dvc <- rbindlist(lapply(c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621", "PRJNA512027"), function(co) {
  x <- prep_cohort(d[cohort == co & !is.na(disease)])
  if (nrow(x) < 20L || length(unique(x$disease)) < 2L || min(table(x$disease)) < 5L) return(NULL)
  tm <- c("disease", if (nlevels(droplevels(x$sex)) > 1L) "sex", "lcpm_ADAR", "lcpm_ADARB1", "isg", "depth")
  r <- fit_one(x, tm, method = "betabin")
  row_of(r, "disease", co, "disease_vs_control")
}))
meta_dvc <- meta_of(dvc, "disease", "familyC_disease_vs_control_DL")
famC <- rbindlist(list(c_rows$prj, c_rows$g213,
                       data.table(cohort = "meta6", analysis = "familyC_disease_vs_control_DL", term = "disease",
                                  n = sum(dvc$n), estimate = meta_dvc$estimate, se = meta_dvc$se,
                                  z = meta_dvc$z, p = meta_dvc$p, method = "DL", note = "")), fill = TRUE)
famC[, p_BH := p.adjust(p, method = "BH")]
fwrite(famC, file.path(out, "familyC_secondary_BH.tsv"), sep = "\t")
fwrite(dvc, file.path(out, "disease_vs_control_per_cohort.tsv"), sep = "\t")
fwrite(meta_dvc, file.path(out, "disease_vs_control_meta.tsv"), sep = "\t")

# ------------------------------------------------------------------ sensitivities S1-S8 (descriptive)
sens <- list(primary = meta_stage)
s1 <- run_family("S1_quasibin", method = "quasibin"); sens$S1 <- meta_of(s1, "stage", "S1_quasibin")
s2 <- run_family("S2_ols_logit", method = "ols"); sens$S2 <- meta_of(s2, "stage", "S2_ols_logit")
dp_ts <- lapply(dp, function(x) { y <- copy(x); y[, ts_n := alu_ts_ref + alu_ts_alt]; y[ts_n > 0 & !is.na(ts_n)] })
s3 <- run_family("S3_strand_aware", y_succ = "alu_ts_alt", y_fail = "alu_ts_ref", data_list = dp_ts,
                 transform = function(x) { x <- copy(x); x[, depth := log10(alu_ts_ref + alu_ts_alt)]; x })
sens$S3 <- meta_of(s3, "stage", "S3_strand_aware")
s4 <- run_family("S4_unadjusted", terms_fun = function(dd, extra) c("stage", if (nlevels(droplevels(dd$sex)) > 1L) "sex"))
sens$S4 <- meta_of(s4, "stage", "S4_unadjusted_stage_sex")
age_ok <- sapply(dp, function(x) mean(!is.na(x$age)) >= 0.9)
s5 <- run_family("S5_plus_age", cohorts = PRIMARY[age_ok], terms_fun = function(dd, extra) c(terms_primary(dd), "age"))
sens$S5 <- if (nrow(s5)) meta_of(s5, "stage", paste0("S5_plus_age_k", sum(age_ok))) else NULL
chkp <- prim[term == "stage" & !is.na(estimate)]
rm_hk <- rma(yi = chkp$estimate, sei = chkp$se, method = "REML", test = "knha")
sens$S6 <- data.table(analysis = "S6_REML_HartungKnapp", term = "stage", k = rm_hk$k, estimate = as.numeric(rm_hk$beta),
                      se = rm_hk$se, z = as.numeric(rm_hk$zval), p = rm_hk$pval, OR_per_unit = exp(as.numeric(rm_hk$beta)),
                      OR_lo = exp(rm_hk$ci.lb), OR_hi = exp(rm_hk$ci.ub), tau2 = rm_hk$tau2, I2 = rm_hk$I2, Q = rm_hk$QE, Q_p = rm_hk$QEp)
hep_co <- intersect(PRIMARY, c("GSE130970", "GSE135251", "GSE162694"))
s7 <- run_family("S7_plus_hepatocyte_fraction", cohorts = hep_co, terms_fun = function(dd, extra) c(terms_primary(dd), "hep_logit"))
sens$S7 <- meta_of(s7, "stage", "S7_plus_hepatocyte_fraction_3cohorts")
sens$S7_ref <- meta_of(prim[cohort %in% hep_co], "stage", "S7_reference_primary_same_3cohorts")
# S9, S10: the prespecified fitter aod::betabin (default optimizer settings), with and without centered covariates
s9 <- run_family("S9_aod_betabin_centered", method = "aod", center = TRUE); sens$S9 <- meta_of(s9, "stage", "S9_aod_betabin_centered")
s10 <- run_family("S10_aod_betabin_prespecified_defaults", method = "aod", center = FALSE)
sens$S10 <- meta_of(s10, "stage", "S10_aod_betabin_prespecified_defaults")
mk_audit <- function(x, tag) x[term %in% c("stage", "lcpm_ADAR"),
  setNames(.(cohort, term, estimate, se, method, logL, bb_code), c("cohort", "term", paste0(c("est_", "se_", "method_", "logL_", "code_"), tag)))]
audit <- Reduce(function(a, b) merge(a, b, by = c("cohort", "term")),
                list(mk_audit(prim, "ml"), mk_audit(s9, "aod_centered"), mk_audit(s10, "aod_uncentered"),
                     s1[term %in% c("stage", "lcpm_ADAR"), .(cohort, term, est_quasibin = estimate)]))
audit[, logL_gain_ml_vs_aod_uncentered := logL_ml - logL_aod_uncentered]
fwrite(audit, file.path(out, "convergence_audit.tsv"), sep = "\t")
loo <- lapply(PRIMARY, function(co) meta_of(prim[cohort != co], "stage", paste0("S8_leave_out_", co)))
sens$S8 <- rbindlist(loo)
fwrite(rbindlist(c(list(s1, s2, s3, s4, s5, s7, s9, s10), list(prim)), fill = TRUE), file.path(out, "sensitivity_per_cohort.tsv"), sep = "\t")
sens_tab <- rbindlist(sens, fill = TRUE)
fwrite(sens_tab, file.path(out, "sensitivity_meta.tsv"), sep = "\t")

# ------------------------------------------------------------------ power note
pw <- prim[term == "stage", .(cohort, n, se, min_detectable_OR_per_step_80pct = exp(2.8 * se))]
pw <- rbind(pw, data.table(cohort = "meta_DL", n = sum(prim[term == "stage"]$n), se = meta_stage$se,
                           min_detectable_OR_per_step_80pct = exp(2.8 * meta_stage$se)))
fwrite(pw, file.path(out, "power_note.tsv"), sep = "\t")

# ------------------------------------------------------------------ descriptive: pooled index by stage (after all tests)
desc <- rbindlist(lapply(PRIMARY, function(co) {
  x <- dp[[co]]
  x[, .(n = .N, median_index = median(alu_alt / (alu_alt + alu_ref)), median_ADAR_log2cpm = median(lcpm_ADAR),
        median_ISG = median(isg)), by = .(stage)][, cohort := co][]
}))
fwrite(desc, file.path(out, "descriptive_by_stage.tsv"), sep = "\t")

# ------------------------------------------------------------------ decision
same_sign <- sum(sign(prim[term == "stage"]$estimate) == sign(meta_stage$estimate), na.rm = TRUE)
dec <- list(
  n_persons_primary = sum(cohort_n$n_persons),
  meta_stage = as.list(meta_stage),
  perm_p = perm_p, n_perm = n_perm,
  cohorts_same_sign_as_meta = same_sign, k_cohorts = sum(!is.na(prim[term == "stage"]$estimate)),
  positive_control = list(adar_est = meta_adar$estimate, one_sided_p = pos_p, pass = pos_pass),
  negative_control = list(est = meta_neg$estimate, p = meta_neg$p, pass = neg_pass),
  rule = list(a_meta_p_lt_0.05 = meta_stage$p < 0.05, b_same_sign_ge_4 = same_sign >= 4L,
              c_perm_p_lt_0.05 = isTRUE(perm_p < 0.05), d_negative_control_pass = neg_pass),
  supported = isTRUE(meta_stage$p < 0.05 && same_sign >= 4L && isTRUE(perm_p < 0.05) && neg_pass && pos_pass),
  uninformative_if_positive_control_fails = !pos_pass,
  verdict = if (!pos_pass) "uninformative (positive control failed)" else
    if (isTRUE(meta_stage$p < 0.05 && same_sign >= 4L && isTRUE(perm_p < 0.05) && neg_pass)) "supported" else "not supported",
  seed = SEED
)
writeLines(jsonlite::toJSON(dec, auto_unbox = TRUE, pretty = TRUE, na = "string", digits = 8),
           file.path(out, "decision.json"))
fwrite(rbindlist(list(meta_stage, meta_adar, meta_neg), fill = TRUE), file.path(out, "meta_primary_controls.tsv"), sep = "\t")

# ------------------------------------------------------------------ post hoc diagnostics (NOT prespecified; cannot change the verdict)
sp <- function(a, b) suppressWarnings(cor(a, b, method = "spearman", use = "complete.obs"))
ph_cor <- rbindlist(lapply(PRIMARY, function(co) {
  x <- dp[[co]]
  idx <- x$alu_alt / (x$alu_alt + x$alu_ref)
  r2 <- summary(lm(lcpm_ADAR ~ lcpm_ADARB1 + isg + stage + depth, data = x))$r.squared
  data.table(cohort = co, n = nrow(x),
             rho_index_ADAR = sp(idx, x$lcpm_ADAR), rho_index_ADARB1 = sp(idx, x$lcpm_ADARB1), rho_index_ISG = sp(idx, x$isg),
             rho_index_stage = sp(idx, x$stage), rho_index_depth = sp(idx, x$depth),
             rho_ADAR_ISG = sp(x$lcpm_ADAR, x$isg), rho_ADAR_ADARB1 = sp(x$lcpm_ADAR, x$lcpm_ADARB1),
             rho_ADAR_stage = sp(x$lcpm_ADAR, x$stage), rho_ISG_stage = sp(x$isg, x$stage),
             VIF_ADAR_in_primary_model = 1 / (1 - r2))
}))
fwrite(ph_cor, file.path(out, "posthoc_correlations.tsv"), sep = "\t")
sex_term <- function(dd) if (nlevels(droplevels(dd$sex)) > 1L) "sex"
p2 <- run_family("P2_ADAR_only_sex_depth", terms_fun = function(dd, extra) c(sex_term(dd), "lcpm_ADAR", "depth"), target = "lcpm_ADAR")
p3 <- run_family("P3_ADAR_ISG_sex_depth", terms_fun = function(dd, extra) c(sex_term(dd), "lcpm_ADAR", "isg", "depth"), target = "lcpm_ADAR")
ph_meta <- rbindlist(list(meta_of(p2, "lcpm_ADAR", "P2_ADAR_only_sex_depth"), meta_of(p3, "lcpm_ADAR", "P3_ADAR_ISG_sex_depth"),
                          meta_of(prim, "lcpm_ADAR", "primary_model_ADAR")))
ph_meta[, p_one_sided_positive := pnorm(z, lower.tail = FALSE)]
fwrite(ph_meta, file.path(out, "posthoc_ADAR_meta.tsv"), sep = "\t")
fwrite(rbindlist(list(p2, p3)), file.path(out, "posthoc_ADAR_per_cohort.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
log_line("done")
print(meta_stage); print(prim[term == "stage"]); print(meta_adar); print(meta_neg); print(audit); print(ph_cor); print(ph_meta)
cat("perm_p", perm_p, " same_sign", same_sign, "\n")
