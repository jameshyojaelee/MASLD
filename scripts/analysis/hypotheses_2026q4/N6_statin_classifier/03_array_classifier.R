#!/usr/bin/env Rscript
# N6: statin-exposure classifier on GSE130991 (Affymetrix HTA 2.0, deposited RMA log2).
# Native methods for this modality: limma on log2 intensities for gene selection,
# ridge logistic regression on a small gene set; a within-sample rank score is the
# platform-agnostic form. Unit: patient (one array each). Gene selection happens
# inside every training fold. See PRESPEC_N6.md.
suppressPackageStartupMessages({library(data.table); library(limma); library(glmnet); library(jsonlite); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: 03_array_classifier.R gene_matrix.rds samples.tsv out_dir")
out <- args[3]
if (dir.exists(out)) stop("Refusing to overwrite output dir")
dir.create(out, recursive = TRUE)
SEED <- 20260923L
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
NREP <- 20L; NPERM <- 200L; NBOOT <- 2000L
set.seed(SEED)

gm <- readRDS(args[1]); G <- gm$G
smp <- fread(args[2], na.strings = "NA", check.names = FALSE)
stopifnot(identical(colnames(G), smp$gsm))
y <- as.integer(smp$statin_treatment == "Statin")
stopifnot(!anyNA(y), sum(y) >= 30L)                       # gate G1
cov <- data.frame(age = as.numeric(smp$age), sex = factor(smp$gender), bmi = as.numeric(smp$bmi))
cov$diab <- as.integer(!is.na(smp[["diabetes_treatment"]]) & smp[["diabetes_treatment"]] == "Yes")
cov$hba1c <- as.numeric(smp[["hba1c_(%)"]])
cov$hba1c[is.na(cov$hba1c)] <- median(cov$hba1c, na.rm = TRUE)
stopifnot(!anyNA(cov$age), !anyNA(cov$bmi), !anyNA(cov$sex))
RK <- apply(G, 2, rank) / nrow(G)                          # within-array percentile rank, label-free
form_primary <- ~ y + age + sex + bmi
form_diab <- ~ y + age + sex + bmi + diab + hba1c

auc_fn <- function(score, yy) {
  r <- rank(score); n1 <- sum(yy == 1L); n0 <- sum(yy == 0L)
  (sum(r[yy == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
make_folds <- function(yy, k, seed) {
  set.seed(seed); f <- integer(length(yy))
  for (cl in 0:1) { ix <- which(yy == cl); f[ix] <- sample(rep_len(seq_len(k), length(ix))) }
  f
}
select_sig <- function(tr, yv, K, covform) {
  d <- cov[tr, , drop = FALSE]; d$y <- yv[tr]
  fit <- eBayes(lmFit(G[, tr, drop = FALSE], model.matrix(covform, d)))
  tt <- fit$t[, "y"]
  idx <- order(abs(tt), decreasing = TRUE)[seq_len(K)]
  list(idx = idx, sgn = sign(tt[idx]))
}
score_test <- function(tr, te, yv, K, covform) {
  s <- select_sig(tr, yv, K, covform)
  cvf <- cv.glmnet(t(G[s$idx, tr, drop = FALSE]), yv[tr], family = "binomial", alpha = 0,
                   nfolds = 10, type.measure = "deviance")
  ridge <- as.numeric(predict(cvf, newx = t(G[s$idx, te, drop = FALSE]), s = "lambda.min", type = "link"))
  up <- s$idx[s$sgn > 0]; dn <- s$idx[s$sgn < 0]
  mm <- function(ix) if (length(ix)) colMeans(RK[ix, te, drop = FALSE]) else rep(0, length(te))
  list(ridge = ridge, rank = mm(up) - mm(dn))
}
run_cv <- function(yv, K, covform, fold_seed, inner_seed) {
  f <- make_folds(yv, 5L, fold_seed); set.seed(inner_seed)
  o_ridge <- o_rank <- numeric(length(yv)); fr <- fk <- numeric(5L)
  for (k in 1:5) {
    te <- which(f == k); tr <- which(f != k)
    s <- score_test(tr, te, yv, K, covform)
    o_ridge[te] <- s$ridge; o_rank[te] <- s$rank
    fr[k] <- auc_fn(s$ridge, yv[te]); fk[k] <- auc_fn(s$rank, yv[te])
  }
  list(ridge = o_ridge, rank = o_rank, auc_ridge = auc_fn(o_ridge, yv), auc_rank = auc_fn(o_rank, yv),
       fold_ridge = mean(fr), fold_rank = mean(fk))
}
pmc <- function(X, FUN) {
  r <- mclapply(X, FUN, mc.cores = NCORES, mc.set.seed = FALSE)
  bad <- vapply(r, function(z) inherits(z, "try-error") || is.null(z), TRUE)
  if (any(bad)) stop("worker failure: ", paste(unique(unlist(lapply(r[bad], as.character))), collapse = "; "))
  r
}

# ---- repeated CV: K = 20, 50, 100 (primary K = 50) and the diabetes-adjusted sensitivity
tasks <- rbind(data.frame(rep = rep(seq_len(NREP), 3L), K = rep(c(20L, 50L, 100L), each = NREP), form = "primary"),
               data.frame(rep = seq_len(NREP), K = 50L, form = "diab"))
res <- pmc(seq_len(nrow(tasks)), function(i) {
  tk <- tasks[i, ]
  run_cv(y, tk$K, if (tk$form == "primary") form_primary else form_diab,
         fold_seed = SEED + tk$rep, inner_seed = SEED + 1000L * tk$K + tk$rep)
})
auc_rep <- rbindlist(lapply(seq_len(nrow(tasks)), function(i)
  data.table(tasks[i, ], auc_ridge = res[[i]]$auc_ridge, auc_rank = res[[i]]$auc_rank,
             fold_mean_ridge = res[[i]]$fold_ridge, fold_mean_rank = res[[i]]$fold_rank)))
fwrite(auc_rep, file.path(out, "cv_auc_by_repeat.tsv"), sep = "\t")

avg_scores <- function(fm, KK, what) {
  ix <- which(tasks$form == fm & tasks$K == KK)
  Reduce(`+`, lapply(ix, function(i) res[[i]][[what]])) / length(ix)
}
boot_ci_headline <- function(idx_list, what, B, seed) {
  set.seed(seed); n <- length(y)
  b <- replicate(B, { s <- sample.int(n, n, replace = TRUE)
    if (length(unique(y[s])) < 2L) NA_real_ else mean(vapply(idx_list, function(i) auc_fn(res[[i]][[what]][s], y[s]), 0)) })
  quantile(b, c(0.025, 0.975), na.rm = TRUE)
}
boot_ci <- function(score, yy, B, seed) {
  set.seed(seed); n <- length(yy)
  b <- replicate(B, { s <- sample.int(n, n, replace = TRUE); if (length(unique(yy[s])) < 2L) NA_real_ else auc_fn(score[s], yy[s]) })
  quantile(b, c(0.025, 0.975), na.rm = TRUE)
}
summ <- list()
for (fm in c("primary", "diab")) for (KK in unique(tasks$K[tasks$form == fm])) for (mth in c("ridge", "rank")) {
  a <- auc_rep[form == fm & K == KK][[paste0("auc_", mth)]]
  sc <- avg_scores(fm, KK, mth)
  ci <- boot_ci(sc, y, NBOOT, SEED + 7L)
  ixs <- which(tasks$form == fm & tasks$K == KK)
  chd <- boot_ci_headline(ixs, if (mth == "ridge") "ridge" else "rank", NBOOT, SEED + 8L)
  fw <- auc_rep[form == fm & K == KK][[paste0("fold_mean_", mth)]]
  summ[[length(summ) + 1L]] <- data.table(adjust = fm, K = KK, method = mth, n = length(y), n_statin = sum(y),
    mean_auc_over_repeats = mean(a), headline_boot_lo = chd[[1]], headline_boot_hi = chd[[2]], sd_over_repeats = sd(a),
    mean_foldwise_auc = mean(fw), auc_of_averaged_oof_score = auc_fn(sc, y), boot_lo = ci[[1]], boot_hi = ci[[2]])
}
summ <- rbindlist(summ)
fwrite(summ, file.path(out, "cv_summary.tsv"), sep = "\t")
oof <- data.table(gsm = smp$gsm, statin = y, ridge_K50 = avg_scores("primary", 50L, "ridge"),
                  rank_K50 = avg_scores("primary", 50L, "rank"))
fwrite(oof, file.path(out, "oof_scores_K50.tsv.gz"), sep = "\t")

# ---- permuted-label nulls (whole pipeline, one 5-fold CV each, K = 50)
age_t <- cut(cov$age, quantile(cov$age, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = FALSE)
strata <- interaction(cov$sex, cov$diab, age_t, drop = TRUE)
perm_fun <- function(stratified) function(i) {
  set.seed(SEED + (if (stratified) 90000L else 50000L) + i)
  yp <- y
  if (stratified) { for (s in levels(strata)) { ix <- which(strata == s); yp[ix] <- y[ix][sample.int(length(ix))] } }
  else yp <- sample(y)
  r <- run_cv(yp, 50L, form_primary, fold_seed = SEED + 777L + i, inner_seed = SEED + 888L + i)
  data.table(perm = i, auc_ridge = r$auc_ridge, auc_rank = r$auc_rank)
}
null_simple <- rbindlist(pmc(seq_len(NPERM), perm_fun(FALSE)))
null_strat <- rbindlist(pmc(seq_len(NPERM), perm_fun(TRUE)))
fwrite(null_simple, file.path(out, "perm_null_unrestricted.tsv"), sep = "\t")
fwrite(null_strat, file.path(out, "perm_null_stratified.tsv"), sep = "\t")
obs <- summ[adjust == "primary" & K == 50L]
emp_p <- function(nullv, o) (1 + sum(nullv >= o)) / (length(nullv) + 1)
perm_summary <- data.table(
  method = c("ridge", "rank"),
  observed_mean_auc = c(obs[method == "ridge", mean_auc_over_repeats], obs[method == "rank", mean_auc_over_repeats]),
  null_unrestricted_mean = c(mean(null_simple$auc_ridge), mean(null_simple$auc_rank)),
  null_unrestricted_max = c(max(null_simple$auc_ridge), max(null_simple$auc_rank)),
  p_unrestricted = c(emp_p(null_simple$auc_ridge, obs[method == "ridge", mean_auc_over_repeats]),
                     emp_p(null_simple$auc_rank, obs[method == "rank", mean_auc_over_repeats])),
  null_stratified_mean = c(mean(null_strat$auc_ridge), mean(null_strat$auc_rank)),
  null_stratified_max = c(max(null_strat$auc_ridge), max(null_strat$auc_rank)),
  p_stratified = c(emp_p(null_strat$auc_ridge, obs[method == "ridge", mean_auc_over_repeats]),
                   emp_p(null_strat$auc_rank, obs[method == "rank", mean_auc_over_repeats])))
fwrite(perm_summary, file.path(out, "perm_summary.tsv"), sep = "\t")

# ---- confound checks
num <- function(nm) as.numeric(smp[[nm]])
vars <- list(age = num("age"), bmi = num("bmi"), height_cm = num("height_(cm)"), waist_cm = num("waist_circonference_(cm)"),
  sbp = num("systolic_blood_pressure_(mmhg)"), dbp = num("diastolic_blood_pressure_(mmhg)"),
  glucose_mmol_l = num("fasting_blood_glucose_(mmol/l)"), insulin_mui_l = num("fasting_insulin_(mui/l)"),
  hba1c_pct = num("hba1c_(%)"), chol_total = num("cholest_total_(mmol/l)"), ldl = num("ldl_cholesterol_(mmol/l)"),
  hdl = num("hdl_cholesterol_(mmol/l)"), tg = num("triglycerides_(mmol/l)"), c_peptide = num("fasting_c-peptide_(ngml)"),
  array_index = suppressWarnings(as.numeric(smp$array_index)))
bal <- rbindlist(lapply(names(vars), function(nm) {
  v <- vars[[nm]]; x1 <- v[y == 1L & !is.na(v)]; x0 <- v[y == 0L & !is.na(v)]
  data.table(variable = nm, type = "continuous", n_statin = length(x1), n_nostatin = length(x0),
             statin = sprintf("%.2f (%.2f)", mean(x1), sd(x1)), nostatin = sprintf("%.2f (%.2f)", mean(x0), sd(x0)),
             smd = (mean(x1) - mean(x0)) / sqrt((var(x1) + var(x0)) / 2), p_nominal = t.test(x1, x0)$p.value)
}))
mk_bin <- function(nm, v) {
  tb <- table(statin = y, v = v)
  data.table(variable = nm, type = "binary", n_statin = sum(y == 1L), n_nostatin = sum(y == 0L),
             statin = sprintf("%d (%.1f%%)", tb["1", "1"], 100 * tb["1", "1"] / sum(tb["1", ])),
             nostatin = sprintf("%d (%.1f%%)", tb["0", "1"], 100 * tb["0", "1"] / sum(tb["0", ])),
             smd = NA_real_, p_nominal = fisher.test(tb)$p.value)
}
bal <- rbind(bal, mk_bin("male", as.integer(cov$sex == "Male")), mk_bin("diabetes_treatment", cov$diab))
bal[, q_BH_15 := p.adjust(p_nominal, "BH")]
fwrite(bal, file.path(out, "balance_table.tsv"), sep = "\t")

folds_rep <- lapply(seq_len(NREP), function(r) make_folds(y, 5L, SEED + r))
clin_cv <- function(cols) {
  mean(vapply(seq_len(NREP), function(r) {
    f <- folds_rep[[r]]; p <- numeric(length(y))
    for (k in 1:5) {
      tr <- which(f != k); te <- which(f == k)
      d <- cbind(cov[, cols, drop = FALSE], y = y)
      fit <- glm(y ~ ., data = d[tr, ], family = binomial())
      p[te] <- predict(fit, newdata = d[te, ], type = "link")
    }
    auc_fn(p, y)
  }, 0))
}
ai <- suppressWarnings(as.numeric(smp$array_index))
ok_ai <- !is.na(ai)
tp <- factor(smp$title_pattern)
tp_tab <- table(statin = y, pattern = tp)
tp_p <- chisq.test(tp_tab)$p.value
tp_cv <- mean(vapply(seq_len(NREP), function(r) {
  f <- folds_rep[[r]]; p <- numeric(length(y))
  for (k in 1:5) { tr <- which(f != k); te <- which(f == k)
    pr <- tapply(y[tr], tp[tr], mean); p[te] <- pr[as.character(tp[te])] }
  auc_fn(p, y) }, 0))
fwrite(as.data.table(as.data.frame.matrix(tp_tab), keep.rownames = "statin"), file.path(out, "title_pattern_by_statin.tsv"), sep = "\t")
baselines <- data.table(
  baseline = c("clinical_age_sex_bmi", "clinical_age_sex_bmi_diab_hba1c", "deposit_index_only_abs(auc-0.5)+0.5_nonmissing", "title_pattern_only_cv"),
  cv_auc = c(clin_cv(c("age", "sex", "bmi")), clin_cv(c("age", "sex", "bmi", "diab", "hba1c")),
             abs(auc_fn(ai[ok_ai], y[ok_ai]) - 0.5) + 0.5, tp_cv))
baselines[, extra := NA_character_]
baselines[3, extra := sprintf("n=%d; spearman(rank score, index)=%.3f", sum(ok_ai), cor(oof$rank_K50[ok_ai], ai[ok_ai], method = "spearman"))]
baselines[4, extra := sprintf("chi-square p=%.3g across patterns", tp_p)]
fwrite(baselines, file.path(out, "baselines.tsv"), sep = "\t")

sc_z <- as.numeric(scale(oof$rank_K50))
d_all <- cbind(cov, y = y, sc = sc_z)
m0 <- glm(y ~ age + sex + bmi + diab + hba1c, data = d_all, family = binomial())
m1 <- glm(y ~ sc + age + sex + bmi + diab + hba1c, data = d_all, family = binomial())
co <- summary(m1)$coefficients["sc", ]
comb <- data.table(or_per_sd_rank_score = unname(exp(co["Estimate"])), or_lo = unname(exp(co["Estimate"] - 1.96 * co["Std. Error"])),
                   or_hi = unname(exp(co["Estimate"] + 1.96 * co["Std. Error"])),
                   lrt_p = pchisq(m0$deviance - m1$deviance, 1L, lower.tail = FALSE), n = length(y))
fwrite(comb, file.path(out, "combined_logistic_rank_score.tsv"), sep = "\t")

strat_def <- list(sex_male = cov$sex == "Male", sex_female = cov$sex == "Female",
                  diabetes_treated = cov$diab == 1L, diabetes_not_recorded = cov$diab == 0L,
                  matched_subset = !is.na(smp[["matched_cohort_statin_treatment"]]),
                  age_tertile1 = age_t == 1L, age_tertile2 = age_t == 2L, age_tertile3 = age_t == 3L)
sa <- rbindlist(lapply(names(strat_def), function(nm) {
  ix <- which(strat_def[[nm]])
  rbindlist(lapply(c("ridge_K50", "rank_K50"), function(m) {
    ci <- boot_ci(oof[[m]][ix], y[ix], 1000L, SEED + 11L)
    data.table(stratum = nm, method = m, n = length(ix), n_statin = sum(y[ix]),
               auc = auc_fn(oof[[m]][ix], y[ix]), boot_lo = ci[[1]], boot_hi = ci[[2]])
  }))
}))
fwrite(sa, file.path(out, "strata_auc.tsv"), sep = "\t")

# ---- full-data limma (descriptive), final signature, positive control
d_full <- cov; d_full$y <- y
fit_full <- eBayes(lmFit(G, model.matrix(form_primary, d_full)))
tt_full <- topTable(fit_full, coef = "y", number = Inf, sort.by = "none")
tt_full$symbol <- rownames(tt_full)
fwrite(as.data.table(tt_full), file.path(out, "limma_full_statin_vs_none.tsv.gz"), sep = "\t")
n_bh <- sum(tt_full$adj.P.Val < 0.05)
top <- order(abs(fit_full$t[, "y"]), decreasing = TRUE)[1:50]
sig <- data.table(symbol = rownames(G)[top], logFC = fit_full$coefficients[top, "y"], t = fit_full$t[top, "y"],
                  P = fit_full$p.value[top, "y"], sign = sign(fit_full$t[top, "y"]))
sig[, adj_P_BH_all_genes := tt_full$adj.P.Val[match(symbol, tt_full$symbol)]]
fwrite(sig, file.path(out, "signature_K50_all_arrays.tsv"), sep = "\t")

pc_genes <- c("HMGCR", "HMGCS1", "LDLR", "INSIG1", "SQLE", "IDI1", "FDFT1", "LSS", "CYP51A1", "MVK", "PMVK", "MVD",
              "FDPS", "DHCR7", "DHCR24", "SC5D", "NSDHL", "HSD17B7", "TM7SF2", "PCSK9")
present <- intersect(pc_genes, rownames(G)); absent <- setdiff(pc_genes, rownames(G))
tv <- fit_full$t[, "y"]
set.seed(SEED + 31L)
nullm <- replicate(10000L, mean(tv[sample.int(length(tv), length(present))]))
pc <- data.table(n_present = length(present), absent = paste(absent, collapse = ","), mean_t = mean(tv[present]),
                 n_positive_t = sum(tv[present] > 0), p_random_sets = (1 + sum(nullm >= mean(tv[present]))) / 10001)
fwrite(pc, file.path(out, "positive_control.tsv"), sep = "\t")
fwrite(data.table(symbol = present, t = tv[present], logFC = fit_full$coefficients[present, "y"]),
       file.path(out, "positive_control_genes.tsv"), sep = "\t")

# ---- gate G2
r50 <- obs[method == "rank"]; d50 <- obs[method == "ridge"]
gate <- list(G1_statin_n = sum(y), G1_pass = sum(y) >= 30L,
             G2_rank_K50_mean_auc = r50$mean_auc_over_repeats, G2_rank_pass = r50$mean_auc_over_repeats >= 0.70,
             ridge_K50_mean_auc = d50$mean_auc_over_repeats, ridge_pass = d50$mean_auc_over_repeats >= 0.70,
             rank_boot_ci = c(r50$boot_lo, r50$boot_hi), ridge_boot_ci = c(d50$boot_lo, d50$boot_hi),
             n_genes_BH05_full_data = n_bh, seed = SEED, repeats = NREP, permutations = NPERM, bootstraps = NBOOT)
write_json(gate, file.path(out, "gate.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(summ); print(perm_summary); print(pc); print(baselines); print(comb); print(sa)
cat(toJSON(gate, auto_unbox = TRUE, pretty = TRUE, digits = 6), "\n")
