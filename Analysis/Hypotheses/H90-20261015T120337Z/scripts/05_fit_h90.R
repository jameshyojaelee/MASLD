# H90 job 6302 (b): primary test, fixed-sequence secondary, G5, PC4, NC1-NC3, E1-E10,
# and RESULT_H90.md with the verdict computed here.
# Unit: donor (one library per donor). Model: score ~ stage + sex + log10(n_hep), OLS,
# HC3 SE (closed form), Freedman-Lane permutation p (19,999 draws, seed 20260923).
suppressPackageStartupMessages(library(edgeR))
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN <- file.path(ROOT, "Analysis/Hypotheses/H90-20261015T120337Z")
IN <- Sys.getenv("H90_IN", file.path(RUN, "outputs"))
OUT <- Sys.getenv("H90_OUT", file.path(RUN, "outputs"))
TEST <- Sys.getenv("H90_TEST", "0") == "1"
RESULT_MD <- if (TEST) file.path(OUT, "RESULT_H90.md") else file.path(RUN, "RESULT_H90.md")
SEED <- 20260923
NPERM <- 19999
setwd(ROOT)

if (!TEST) {
  chk <- readLines(file.path(RUN, "CHECKED_BEFORE_OUTCOME.txt"))
  if (!any(grepl("PROCEED: yes", chk, fixed = TRUE))) stop("STOP: CHECKED_BEFORE_OUTCOME.txt not PROCEED: yes")
  ps <- system2("sha256sum", c(file.path(RUN, "PRESPEC_H90.md"),
                 "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H90.md"), stdout = TRUE)
  if (!all(substr(ps, 1, 64) == "ce61a0eb5c6820a5b1f0bb3a0b7cd6d3a63358d1cb8b82dc5b12815e8d62da27")) stop("STOP: prespec hash differs")
  st <- system2("sha256sum", c("-c", "--quiet", file.path(RUN, "inputs.sha256")), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(st, "status")) && attr(st, "status") != 0) stop("STOP: input hash check failed: ", paste(st, collapse = "; "))
  cat("prespec and 25 input hashes OK\n")
}
for (f in c("results_primary.tsv")) if (file.exists(file.path(OUT, f))) stop("refusing to overwrite ", f)

rd <- function(f, ...) read.delim(file.path(IN, f), check.names = FALSE, stringsAsFactors = FALSE, ...)
wr <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- fitting helpers ------------------------------------------------------------
# Closed-form OLS with HC3 for coefficient j; Y may be a matrix (one column per outcome)
hc3_fit <- function(Y, X, j) {
  Y <- as.matrix(Y)
  XtXi <- solve(crossprod(X))
  Cm <- XtXi %*% t(X)                     # p x n
  B <- Cm %*% Y
  E <- Y - X %*% B
  h <- rowSums((X %*% XtXi) * X)
  W <- E^2 / (1 - h)^2
  se_hc3 <- sqrt(as.vector((Cm[j, ]^2) %*% W))
  df <- nrow(X) - ncol(X)
  sigma <- sqrt(colSums(E^2) / df)
  se_ols <- sigma * sqrt(XtXi[j, j])
  list(beta = B[j, ], se = se_hc3, t = B[j, ] / se_hc3, df = df, sigma = sigma, se_ols = se_ols, h = h, E = E)
}
fl_null <- function(y, X, j, P) {
  X0 <- X[, -j, drop = FALSE]
  f0 <- X0 %*% solve(crossprod(X0), crossprod(X0, y))
  r0 <- as.vector(y - f0)
  Ys <- matrix(r0[t(P)], nrow = length(y)) + as.vector(f0)   # column k = f0 + r0[P[k, ]]
  hc3_fit(Ys, X, j)$t
}
perm_p <- function(t, tnull) (1 + sum(abs(tnull) >= abs(t))) / (length(tnull) + 1)
full_test <- function(y, X, j, P, label) {
  f <- hc3_fit(y, X, j)
  tn <- fl_null(y, X, j, P)
  qt975 <- qt(0.975, f$df)
  list(row = data.frame(test = label, n = length(y), beta = f$beta, hc3_se = f$se, hc3_t = f$t,
                        ci_low = f$beta - qt975 * f$se, ci_high = f$beta + qt975 * f$se,
                        resid_sd = f$sigma, beta_sd_units = f$beta / f$sigma,
                        ci_low_sd = (f$beta - qt975 * f$se) / f$sigma, ci_high_sd = (f$beta + qt975 * f$se) / f$sigma,
                        perm_p = perm_p(f$t, tn), ols_p_printed_not_used = 2 * pt(-abs(f$beta / f$se_ols), f$df),
                        hc3_p_printed_not_used = 2 * pt(-abs(f$t), f$df), df = f$df, n_perm = length(tn)),
       tnull = tn, fit = f)
}
draw_perm <- function(n) { set.seed(SEED); t(replicate(NPERM, sample.int(n))) }

# ---- inputs ---------------------------------------------------------------------
sc <- rd("donor_scores.tsv")
cv <- rd("donor_covariates_hep.tsv")
cv <- cv[order(cv$gsm), ]
cv$sex_f <- as.numeric(cv$sex == "female")
cv$log10_n_hep <- log10(cv$n_hep)
n <- nrow(cv)
getscore <- function(ct, mod, gsm) { s <- sc[sc$cell_type == ct & sc$module == mod, ]; s$score[match(gsm, s$gsm)] }
y <- getscore("hepatocytes", 20, cv$gsm)
if (any(is.na(y))) stop("primary score missing for some included donors")
X <- cbind(1, stage = cv$stage, sex_f = cv$sex_f, log10_n_hep = cv$log10_n_hep)
J <- 2
P <- draw_perm(n)
gz <- gzfile(file.path(OUT, "permutation_matrix.tsv.gz"), "w"); write.table(P, gz, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); close(gz)

# ---- G5 permutation integrity and HC3 check against sandwich --------------------
prim <- full_test(y, X, J, P, "primary: Ductular injury (hepatocytes 20)")
X0 <- X[, -J, drop = FALSE]; f0 <- X0 %*% solve(crossprod(X0), crossprod(X0, y)); r0 <- as.vector(y - f0)
t_id <- hc3_fit(as.vector(f0) + r0[seq_len(n)], X, J)$t
g5_perm <- abs(t_id - prim$fit$t) < 1e-10
tmpd <- file.path(OUT, "g5_tmp"); dir.create(tmpd, showWarnings = FALSE)
write.table(data.frame(y = y, X[, -1]), file.path(tmpd, "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c('d <- read.delim(commandArgs(TRUE)[1]); library(sandwich)',
             'f <- lm(y ~ stage + sex_f + log10_n_hep, d)',
             'v <- vcovHC(f, type = "HC3")',
             'cat(sprintf("%.15e\\t%s\\n", sqrt(v["stage","stage"]), as.character(packageVersion("sandwich"))))'),
           file.path(tmpd, "sandwich_check.R"))
sw <- system2("/gpfs/commons/home/jameslee/micromamba/envs/celltype_bio/bin/Rscript",
              c("--vanilla", file.path(tmpd, "sandwich_check.R"), file.path(tmpd, "design.tsv")), stdout = TRUE)
sw <- strsplit(tail(sw, 1), "\t")[[1]]
se_sw <- as.numeric(sw[1])
g5_sw <- abs(se_sw - prim$fit$se) < 1e-10
g5 <- data.frame(check = c("identity permutation reproduces observed HC3 t", "closed-form HC3 SE equals sandwich::vcovHC HC3"),
                 value = c(sprintf("%.3e", abs(t_id - prim$fit$t)), sprintf("%.3e (sandwich %s)", abs(se_sw - prim$fit$se), sw[2])),
                 tolerance = "1e-10", pass = c(g5_perm, g5_sw))
wr(g5, "G5_hc3_check.tsv")
if (!all(g5$pass)) {
  writeLines("untestable: G5 permutation integrity or HC3 check failed", file.path(OUT, "STOP_G5.txt"))
  stop("STOP G5")
}
gz <- gzfile(file.path(OUT, "null_t_primary.tsv.gz"), "w"); writeLines(format(prim$tnull, digits = 15), gz); close(gz)

# ---- PC4 detection on the real residual structure -------------------------------
d4 <- rd("D4_bridge.tsv")
d4p <- d4[d4$score == "projection" & d4$dataset == "n-weighted", ]
ns <- as.vector(table(factor(cv$stage, levels = 0:2)))
xx <- rep(0:2, ns); mux <- rep(c(0, d4p$m1, d4p$m2), ns)
implied <- cov(xx, mux) / var(xx)
q95 <- quantile(abs(prim$tnull), 0.95, names = FALSE)
sdr <- prim$fit$sigma
pc4 <- do.call(rbind, lapply(c(gate_MDE = 0.44, fixed_smaller = 0.28, discovery_implied = implied), function(delta) {
  det <- vapply(1:1000, function(k) {
    sp <- cv$stage[P[k, ]]
    Xk <- X; Xk[, "stage"] <- sp
    f <- hc3_fit(y + delta * sp * sdr, Xk, J)
    f$beta > 0 && abs(f$t) > q95
  }, logical(1))
  data.frame(delta_sd_per_step = delta, detection = mean(det))
}))
pc4$label <- c("gate MDE per step", "fixed smaller effect (about 70% of discovery-implied slope)", "discovery-implied slope (D4 projection means at realized n)")
pc4$null_abs_t_q95 <- q95
wr(pc4, "PC4_detection.tsv")
det_implied <- pc4$detection[3]

# ---- NC1, NC2, NC3 ----------------------------------------------------------------
y2 <- getscore("hepatocytes", 2, cv$gsm)
nc1 <- if (all(!is.na(y2))) full_test(y2, X, J, P, "NC1: Sinusoidal endothelial (hepatocytes 2)")$row else data.frame(test = "NC1", note = "not testable")
nc1$call <- if (!is.null(nc1$perm_p)) nc1$perm_p < 0.05 else NA
wr(nc1, "NC1_null_program.tsv")
nc2_frac <- mean(2 * pt(-abs(prim$tnull), prim$fit$df) < 0.05)
nc2 <- data.frame(fraction_null_draws_hc3_p_below_0.05 = nc2_frac, predicted = "0.04 to 0.08",
                  hc3_p_label = ifelse(nc2_frac > 0.08, "anticonservative", "not anticonservative"))
wr(nc2, "NC2_calibration.tsv")
nc3s <- as.matrix(read.delim(file.path(IN, "NC3_set_scores.tsv.gz"), check.names = FALSE))
nc3s <- nc3s[, cv$gsm]
t_sets <- hc3_fit(t(nc3s), X, J)$t
nc3_pct <- 100 * mean(t_sets < prim$fit$t)
nc3 <- data.frame(n_sets = length(t_sets), observed_t = prim$fit$t, percentile_of_observed = nc3_pct,
                  set_t_median = median(t_sets), set_t_q95 = quantile(t_sets, 0.95, names = FALSE))
wr(nc3, "NC3_matched_sets.tsv")

# ---- primary result and fixed-sequence secondary ----------------------------------
pr <- prim$row
verdict <- if (pr$beta > 0 && pr$perm_p < 0.05) "pass" else if (pr$beta < 0 && pr$perm_p < 0.05) "opposite direction" else "fail"
cov3 <- rd("D3_coverage.tsv")
g3s <- with(cov3[cov3$cell_type == "hepatocytes" & cov3$module == 8, ], n_scored >= 8 & retained_weight >= 0.80)
y8 <- getscore("hepatocytes", 8, cv$gsm)
sec <- if (g3s && all(!is.na(y8))) full_test(y8, X, J, P, "secondary: Stromal ECM (hepatocytes 8)") else NULL
sec_status <- if (is.null(sec)) "untestable (G3 coverage)" else if (verdict != "pass") "not tested (primary did not pass); descriptive only (E9)" else
  if (sec$row$beta > 0 && sec$row$perm_p < 0.05) "pass" else "fail"
res <- pr; res$verdict <- verdict
if (!is.null(sec)) { s <- sec$row; s$verdict <- sec_status; res <- rbind(res, s) }
wr(res, "results_primary.tsv")

# ---- descriptive arms E1-E10 -------------------------------------------------------
E <- list()
E$E1 <- full_test(y, cbind(1, stage = cv$stage), 2, P, "E1 no covariates: score ~ stage")$row
cv$log10_med_umi <- log10(cv$median_umi_hep)
E$E2 <- full_test(y, cbind(1, stage = cv$stage, sex_f = cv$sex_f, log10_med_umi = cv$log10_med_umi), 2, P, "E2 covariate swap: log10 median UMI")$row
masl <- as.numeric(cv$group == "MASL"); mash <- as.numeric(cv$group %in% c("eMASH", "aMASH")); masld <- masl + mash
cvs <- cbind(sex_f = cv$sex_f, log10_n_hep = cv$log10_n_hep)
E$E3 <- rbind(
  full_test(y, cbind(1, masld = masld, mash = mash, cvs), 3, P, "E3 MASH vs MASL")$row,
  full_test(y, cbind(1, masl = masl, mash = mash, cvs), 2, P, "E3 MASL vs ctrl")$row,
  full_test(y, cbind(1, masl = masl, mash = mash, cvs), 3, P, "E3 MASH vs ctrl")$row,
  full_test(y, cbind(1, masld = masld, cvs), 2, P, "E3 MASLD vs ctrl")$row)
fl <- which(cv$series == "FL")
P_fl <- draw_perm(length(fl))   # A7: own matrix, seed set immediately before, sorted GSM order
E$E4 <- full_test(y[fl], X[fl, ], J, P_fl, "E4 FL series only, ordinal slope")$row
stage4 <- c(ctrl = 0, MASL = 1, eMASH = 2, aMASH = 3)[cv$group]
E$E5 <- full_test(y, cbind(1, stage4 = stage4, cvs), 2, P, "E5 four-level coding ctrl0 MASL1 eMASH2 aMASH3")$row
gm <- do.call(rbind, lapply(c("ctrl", "MASL", "eMASH", "aMASH"), function(g) { v <- y[cv$group == g]; m <- mean(v); se <- sd(v) / sqrt(length(v))
  data.frame(group = g, n = length(v), mean = m, ci_low = m - qt(.975, length(v) - 1) * se, ci_high = m + qt(.975, length(v) - 1) * se) }))
wr(gm, "E5_group_means.tsv")
e6 <- rd("E6_pseudobulk_counts.tsv"); e6 <- e6[match(cv$gsm, e6$gsm), ]
sw6 <- rd("scored_genes_weights.tsv"); sw6 <- sw6[sw6$cell_type == "hepatocytes" & sw6$module == 20, ]
cnt <- t(as.matrix(e6[, sw6$gene]))
lc <- cpm(cnt, lib.size = e6$lib_size, log = TRUE, prior.count = 2)
zc <- t(scale(t(lc)))
y6 <- as.vector(t(zc) %*% sw6$weight)
E$E6 <- full_test(y6, X, J, P, "E6 pseudobulk log-CPM score")$row
lodo <- t(vapply(seq_len(n), function(i) { f <- hc3_fit(y[-i], X[-i, ], J); c(f$beta, f$t, 2 * pt(-abs(f$t), f$df)) }, numeric(3)))
E$E7 <- data.frame(test = "E7 leave-one-donor-out", n = n - 1, beta_min = min(lodo[, 1]), beta_max = max(lodo[, 1]),
                   t_min = min(lodo[, 2]), t_max = max(lodo[, 2]), n_deletions_hc3_p_ge_0.05 = sum(lodo[, 3] >= 0.05),
                   dropped_donor_at_min_t = cv$gsm[which.min(lodo[, 2])])
wr(data.frame(gsm = cv$gsm, beta = lodo[, 1], hc3_t = lodo[, 2], hc3_p = lodo[, 3]), "E7_lodo.tsv")
e10 <- rd("E10_scores.tsv"); k10 <- match(e10$gsm, cv$gsm)
X10 <- cbind(1, stage = cv$stage[k10], sex_f = cv$sex_f[k10], log10_n_hep = log10(e10$n_hep_kept))
P10 <- if (length(k10) == n) P else draw_perm(length(k10))
E$E10 <- full_test(e10$score, X10, 2, P10, "E10 hepatocyte nuclei with any KRT7/KRT19 count removed")$row
E$E9 <- if (!is.null(sec)) { s <- sec$row; s$test <- paste("E9", s$test, if (verdict != "pass") "(descriptive only)" else "(also the secondary)"); s } else NULL
rb <- function(L) { cols <- unique(unlist(lapply(L, names))); do.call(rbind, lapply(L, function(d) { for (c in setdiff(cols, names(d))) d[[c]] <- NA; d[cols] })) }
wr(rb(E[c("E1", "E2", "E3", "E4", "E5", "E6", "E9", "E10")]), "E_descriptive_arms.tsv")
wr(E$E7, "E7_summary.tsv")

# ---- E8: 117 programs in their own lineage ----------------------------------------
ln <- rd("lineage_donor_nuclei.tsv")
dq <- rd("donor_qc.tsv")
pe <- read.delim("figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/primary_effects.tsv")
reg <- read.delim("figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv")
progs <- unique(reg[, c("cell_type", "module", "module_name")])
pcache <- list()
e8 <- do.call(rbind, lapply(seq_len(nrow(progs)), function(i) {
  ct <- progs$cell_type[i]; md <- progs$module[i]
  s <- sc[sc$cell_type == ct & sc$module == md, ]
  base <- data.frame(cell_type = ct, module = md, module_name = progs$module_name[i], testable = FALSE, n = nrow(s),
                     beta = NA, hc3_t = NA, perm_p = 1, fig4c_beta = pe$beta[pe$cell_type == ct & pe$module == md])
  if (nrow(s) < 10) return(base)
  s <- s[order(s$gsm), ]
  l <- ln[ln$cell_type == ct, ]; nn <- l$n_nuclei[match(s$gsm, l$gsm)]
  d <- dq[match(s$gsm, dq$gsm), ]
  if (length(unique(d$stage)) < 2) return(base)
  Xe <- cbind(1, stage = d$stage, sex_f = as.numeric(d$sex == "female"), log10_n = log10(nn))
  key <- paste(s$gsm, collapse = ",")
  Pe <- if (identical(s$gsm, cv$gsm)) P else { if (is.null(pcache[[key]])) pcache[[key]] <<- draw_perm(nrow(s)); pcache[[key]] }
  r <- full_test(s$score, Xe, 2, Pe, paste(ct, md))$row
  base$testable <- TRUE; base$beta <- r$beta; base$hc3_t <- r$hc3_t; base$perm_p <- r$perm_p; base
}))
e8$bh_q_n117 <- p.adjust(e8$perm_p, "BH", n = 117)
e8$sign_agrees_fig4c <- ifelse(e8$testable, sign(e8$beta) == sign(e8$fig4c_beta), NA)
wr(e8, "h90_117_descriptive.tsv")
e8_agree <- sum(e8$sign_agrees_fig4c, na.rm = TRUE); e8_test <- sum(e8$testable)

# ---- RESULT_H90.md ------------------------------------------------------------------
gate4 <- rd("gate_G4_bridge.tsv")
rho <- as.numeric(gate4$value[grepl("G4c", gate4$check)])
g13 <- rd("gate_G1_G3.tsv", header = TRUE)
fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")
fp <- function(p) if (p < 1e-3) formatC(p, digits = 2, format = "e") else fmt(p, 4)
nstr <- sprintf("%d control, %d MASL, %d MASH donors", ns[1], ns[2], ns[3])
e3m <- E$E3[E$E3$test == "E3 MASH vs MASL", ]
main_sentence <- switch(verdict,
  "pass" = sprintf("The Ductular injury score was higher in later-stage donors (cross-sectional) in held-out GSE289173 (%s; beta %s SD per stage step, 95%% CI %s to %s; Freedman-Lane permutation p %s).",
                   nstr, fmt(pr$beta_sd_units, 2), fmt(pr$ci_low_sd, 2), fmt(pr$ci_high_sd, 2), fp(pr$perm_p)),
  "fail" = sprintf("In held-out GSE289173 the Ductular injury score did not show a stage-associated increase (%s; beta %s SD per stage step, 95%% CI %s to %s; Freedman-Lane permutation p %s).",
                   nstr, fmt(pr$beta_sd_units, 2), fmt(pr$ci_low_sd, 2), fmt(pr$ci_high_sd, 2), fp(pr$perm_p)),
  "opposite direction" = sprintf("In held-out GSE289173 the Ductular injury score was lower in later-stage donors (cross-sectional), the opposite direction to discovery (%s; beta %s SD per stage step, 95%% CI %s to %s; Freedman-Lane permutation p %s).",
                   nstr, fmt(pr$beta_sd_units, 2), fmt(pr$ci_low_sd, 2), fmt(pr$ci_high_sd, 2), fp(pr$perm_p)))
pc4_sentence <- if (verdict == "pass") NULL else if (det_implied >= 0.80)
  sprintf("At the discovery-implied slope (%s SD per step) the cohort had adequate detection (%s on the real residual structure) for the discovery effect and did not find it.", fmt(implied, 2), fmt(det_implied, 2)) else
  sprintf("Did not replicate; this cohort could not exclude the discovery effect (detection %s at the discovery-implied slope of %s SD per step).", fmt(det_implied, 2), fmt(implied, 2))
caveats <- c(
  sprintf("Control source: 20 of 23 GEO controls come from a separately named series (ctrl_xx), so the stage slope mixes stage with control source. MASH versus MASL (all FL series, E3): beta %s (%s SD), 95%% CI %s to %s, permutation p %s.",
          fmt(e3m$beta), fmt(e3m$beta_sd_units, 2), fmt(e3m$ci_low), fmt(e3m$ci_high), fp(e3m$perm_p)),
  sprintf("Score: the tested score is the frozen-weight projection of Ductular injury, not the Fig 4C Hotspot module score; within-dataset Spearman correlation of the two in discovery (G4c) = %s.%s",
          fmt(rho, 2), if (rho < 0.5) " The projection agrees weakly with the Fig 4C score." else ""),
  sprintf("NC1 (Sinusoidal endothelial, hepatocyte module 2): beta %s, permutation p %s (%s).", fmt(nc1$beta), fp(nc1$perm_p), ifelse(isTRUE(nc1$call), "called; a global shift that moves an unrelated program", "no call")),
  sprintf("NC3 (1,000 matched random gene sets): observed t at the %s percentile.", fmt(nc3_pct, 1)))
spec_flag <- verdict == "pass" && (isTRUE(nc1$call) || nc3_pct < 95)
md <- c("# RESULT_H90", "",
  sprintf("Run: `%s`. Prespec sha256 ce61a0eb5c6820a5b1f0bb3a0b7cd6d3a63358d1cb8b82dc5b12815e8d62da27. Verdict computed in `scripts/05_fit_h90.R`.", basename(RUN)), "",
  sprintf("**Primary verdict (Ductular injury, hepatocyte module 20): %s**", toupper(verdict)), "",
  main_sentence, if (!is.null(pc4_sentence)) c("", pc4_sentence), "",
  if (spec_flag) c("The rise is not specific to the program (NC1 call or NC3 percentile below 95).", "") else NULL,
  "## Qualifiers", "", paste("-", caveats), "",
  "## Primary and secondary", "",
  sprintf("- Primary model `score ~ stage + sex + log10(hepatocyte nuclei)`, n = %d donors (%s). beta %s (score units), HC3 SE %s, HC3 t %s, residual df %d, permutation p %s (19,999 Freedman-Lane draws, seed 20260923). Family: one test, alpha 0.05, no correction.",
          pr$n, nstr, fmt(pr$beta, 4), fmt(pr$hc3_se, 4), fmt(pr$hc3_t, 3), pr$df, fp(pr$perm_p)),
  sprintf("- Secondary (Stromal ECM, hepatocyte module 8, fixed sequence): %s.%s", sec_status,
          if (!is.null(sec)) sprintf(" beta %s SD per step, permutation p %s.", fmt(sec$row$beta_sd_units, 2), fp(sec$row$perm_p)) else ""),
  sprintf("- PC4 detection: %s.", paste(sprintf("%s at %s SD/step", fmt(pc4$detection, 3), fmt(pc4$delta_sd_per_step, 3)), collapse = "; ")),
  sprintf("- NC2: %s of null draws had parametric HC3 p < 0.05 (%s).", fmt(nc2_frac, 4), nc2$hc3_p_label), "",
  "## Descriptive (outside the family; no p here is called significant)", "",
  sprintf("- %s: beta %s (%s SD), permutation p %s", E$E1$test, fmt(E$E1$beta), fmt(E$E1$beta_sd_units, 2), fp(E$E1$perm_p)),
  sprintf("- %s: beta %s (%s SD), permutation p %s", E$E3$test, fmt(E$E3$beta), fmt(E$E3$beta_sd_units, 2), sapply(E$E3$perm_p, fp)),
  sprintf("- %s: n %d, beta %s (%s SD), permutation p %s", E$E4$test, E$E4$n, fmt(E$E4$beta), fmt(E$E4$beta_sd_units, 2), fp(E$E4$perm_p)),
  sprintf("- %s: beta %s, permutation p %s", E$E6$test, fmt(E$E6$beta), fp(E$E6$perm_p)),
  sprintf("- E7 leave-one-donor-out: beta %s to %s; %d of %d deletions with HC3 p >= 0.05.", fmt(E$E7$beta_min), fmt(E$E7$beta_max), E$E7$n_deletions_hc3_p_ge_0.05, n),
  sprintf("- %s: beta %s, permutation p %s", E$E10$test, fmt(E$E10$beta), fp(E$E10$perm_p)),
  sprintf("- E8: %d of 117 programs testable in their lineage; %d of %d testable programs share the sign of the Fig 4C 62-donor beta; %d with BH q < 0.05 over 117 (description of this cohort, not a replication call).",
          e8_test, e8_agree, e8_test, sum(e8$bh_q_n117 < 0.05)),
  sprintf("- Group means of the primary score: %s.", paste(sprintf("%s %s (n %d)", gm$group, fmt(gm$mean), gm$n), collapse = "; ")), "",
  "## Reading under the H59 table", "",
  "Apply the PRESPEC_H90 'Dependence on H59' table with the H59 verdict when it arrives; the H90 column is the primary verdict above. A pass replicates one prespecified program's direction, not the 117-program selection.", "",
  "## Sealed status", "",
  "GSE289173 is registered `withheld_sealed` in `Analysis/MASLD_Model_Benchmark/config/datasets/gse289173.toml`. Two events ended that status: on 2026-10-15 the feasibility gate read disease group and sex for all 48 samples from the GEO series matrix (beyond the 2 samples in the benchmark's 2026-08-21 incident), and on 2026-10-15 this run downloaded and opened all 48 expression matrices. This run did not edit the registry; the owner should update it.", "",
  "## Not verified", "",
  "- Whether aMASH denotes advanced fibrosis and how MASL and MASH were assigned (paper not read) [Unverified].",
  "- Hepatocyte calls are CellTypist on unintegrated nuclei of one study, not the harmonized discovery atlas (D2).")
if (TEST) md <- c("# TEST RUN ON SYNTHETIC SCORES - NOT A RESULT", md)
writeLines(md, RESULT_MD)
writeLines(c(capture.output(sessionInfo()), paste("sandwich (celltype_bio env, G5 check):", sw[2])), file.path(OUT, "sessionInfo.txt"))
cat("verdict:", verdict, "\n")
