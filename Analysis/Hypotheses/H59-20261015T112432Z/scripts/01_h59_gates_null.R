# H59 stage A: input hashes, gates G0-G6, diagnostics D1-D8, controls PC1-PC4,
# NC1, NC2, and permutation null matrices. Writes NO observed stage statistic.
RUN <- Sys.getenv("H59_RUN_DIR", "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H59-20261015T112432Z")
source(file.path(RUN, "scripts/h59_lib.R"))
NULLDATA <- Sys.getenv("H59_NULLDATA", "0") == "1"
if (!dir.exists(OUT)) dir.create(OUT, recursive = FALSE) else if (length(list.files(OUT))) stop("OUT not empty: ", OUT)
msg("stage A; B =", B, "cores =", NCORES, "OUT =", OUT, "NULLDATA =", NULLDATA)

gate <- data.table(gate = character(), pass = logical(), detail = character())
add_gate <- function(g, pass, detail) { gate <<- rbind(gate, data.table(gate = g, pass = pass, detail = detail)); msg("GATE", g, pass, detail) }
stop_rule <- function(text) {
  wtsv(gate, "gates.tsv")
  writeLines(text, file.path(OUT, "STOP_RULE.txt"))
  msg("STOP:", text); quit(save = "no", status = 0)
}

hashes <- check_hashes()
wtsv(hashes, "inputs.sha256.tsv")
dat <- load_data()
if (NULLDATA) {   # debug only: shuffle stage within dataset per donor (outcome-free code test)
  set.seed(1); e <- dat$elig
  e[, stage := sample(stage), by = dataset]
  m <- e[, .(donor, stage)]
  for (nm in c("p", "eq")) { x <- dat[[nm]]; x[, stage := NULL]; dat[[nm]] <- merge(x, m, by = "donor") }
}
programs <- dat$programs
p62 <- dat$p
p60 <- p62[dataset != DROP_DATASET]
pfit <- p60[dataset %in% FIT_DATASETS]

# ------------------------------------------------------------------ G0 / PC1 / PC2
exprs <- parse(inp("producer_512"))
fail <- function(...) stop(..., call. = FALSE)
direction_of <- function(x) fifelse(is.na(x), NA_character_, fifelse(x > 0, "positive", fifelse(x < 0, "negative", "zero")))
EXPECTED_PROGRAMS <- 117L
for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                     as.character(e[[2]]) %in% c("fit_fixed_ols", "fit_program_family")) eval(e)
prog_keys <- programs[, .(cell_type, module)]
prod_in <- function(x) { y <- copy(x); y[, stage_ordinal := stage]; y }
a62 <- fit_program_family(prod_in(p62), prog_keys, "stage_ordinal", TRUE, TRUE)
a60 <- fit_program_family(prod_in(p60), prog_keys, "stage_ordinal", TRUE, TRUE)
pe <- fread(inp("primary_effects"))
c1 <- merge(a62, pe, by = c("cell_type", "module"), suffixes = c("", "_ref"))
d1 <- c(beta = max(abs(c1$beta - c1$beta_ref)), q = max(abs(c1$qvalue - c1$qvalue_ref)),
        hc3q = max(abs(c1$hc3_qvalue - c1$hc3_qvalue_ref)), n = max(abs(c1$n_donors - c1$n_donors_ref)))
duct62 <- a62[cell_type == "hepatocytes" & module == 20]; strom62 <- a62[cell_type == "hepatocytes" & module == 8]
g0a <- nrow(c1) == 117 && d1["beta"] <= 1e-8 && d1["q"] <= 1e-8 && d1["hc3q"] <= 1e-8 && d1["n"] == 0 &&
  abs(duct62$beta - 0.402086) < 5e-7 && abs(duct62$qvalue - 0.037633) < 5e-7 && abs(strom62$qvalue - 0.085458) < 5e-7
v60 <- fread(inp("verifier_60"))
c2 <- merge(a60, v60, by = c("cell_type", "module"), suffixes = c("", "_ver"))
amb <- fread(inp("ambient_effects"))[analysis_universe == "complete_case_common_universe"]
c3 <- merge(a60, amb[, .(cell_type, module, raw_qvalue, n_model_donors)], by = c("cell_type", "module"))
d2 <- c(q_ver = max(abs(c2$qvalue - c2$q)), n_ver = max(abs(c2$n_donors - c2$n)),
        q_amb = max(abs(c3$qvalue - c3$raw_qvalue)), n_amb = max(abs(c3$n_donors - c3$n_model_donors)))
duct60 <- a60[cell_type == "hepatocytes" & module == 20]
g0b <- nrow(c2) == 117 && nrow(c3) == 117 && d2["q_ver"] <= 1e-8 && d2["q_amb"] <= 1e-8 &&
  d2["n_ver"] == 0 && d2["n_amb"] == 0 && abs(duct60$qvalue - 0.098523) < 5e-7
pc12 <- data.table(check = c(names(d1), names(d2), "ductular_beta62", "ductular_q62", "stromal_q62",
                             "ductular_q60", "calls_q05_62", "calls_q05_60"),
                   value = c(d1, d2, duct62$beta, duct62$qvalue, strom62$qvalue, duct60$qvalue,
                             sum(a62$qvalue < .05, na.rm = TRUE), sum(a60$qvalue < .05, na.rm = TRUE)))
wtsv(pc12, "PC1_PC2_G0_reproduction.tsv")
add_gate("G0a", g0a, sprintf("max|dbeta|=%.3g max|dq|=%.3g max|dHC3q|=%.3g max|dn|=%g; Ductular beta %.6f q %.6f; Stromal q %.6f",
                             d1["beta"], d1["q"], d1["hc3q"], d1["n"], duct62$beta, duct62$qvalue, strom62$qvalue))
add_gate("G0b", g0b, sprintf("verifier max|dq|=%.3g, ambient raw max|dq|=%.3g, max|dn| %g/%g; Ductular q60 %.6f; calls60 %d",
                             d2["q_ver"], d2["q_amb"], d2["n_ver"], d2["n_amb"], duct60$qvalue, sum(a60$qvalue < .05)))
if (!(g0a && g0b) && !NULLDATA) stop_rule("untestable: pipeline does not reproduce Fig 4C")

# ------------------------------------------------------------------ D1, G1
D1 <- p62[, .N, by = .(cell_type, module, dataset, stage)][order(cell_type, module, dataset, stage)]
wtsv(D1, "D1_donors_by_program_dataset_stage_62.tsv")
excl <- unique(rbind(
  dat$meta[!donor %in% dat$elig$donor, .(donor, dataset, reason = paste0("outside 62-donor universe: stage ", disease_stage_coarse, ", exclude_stage_analysis ", exclude_stage_analysis))],
  p62[dataset == DROP_DATASET, .(donor, dataset, reason = "GSE189600 removed: one donor per stage, zero residual df, no ambient correction")],
  p60[!dataset %in% FIT_DATASETS, .(donor, dataset, reason = "single-stage dataset: no information on stage coefficient under dataset intercept and variance")]))
wtsv(excl, "D1_donors_excluded.tsv")
cnt <- pfit[, .N, by = .(cell_type, module, dataset, stage)]
g1 <- programs[, .(cell_type, module)][, {
  ok <- sapply(FIT_DATASETS, function(ds) {
    x <- cnt[cell_type == .BY$cell_type & module == .BY$module & dataset == ds]
    sum(x$N) >= 8 && sum(x$N >= 2) >= 2 })
  .(g1_pass = all(ok), n_202379 = cnt[cell_type == .BY$cell_type & module == .BY$module & dataset == "GSE202379", sum(N)],
    n_244832 = cnt[cell_type == .BY$cell_type & module == .BY$module & dataset == "GSE244832", sum(N)])
}, by = .(cell_type, module)]
wtsv(g1, "G1_min_n.tsv")
named_g1 <- merge(NAMED, g1, by = c("cell_type", "module"))$g1_pass
add_gate("G1", all(named_g1), sprintf("%d of 117 programs pass; named programs pass: %s", sum(g1$g1_pass), paste(named_g1, collapse = ",")))
testable <- g1[g1_pass == TRUE, .(cell_type, module)]

# ------------------------------------------------- D2, D3, D4, D5/G2, D6, D8
zc <- qnorm(1 - 0.05 / 117 / 2) + qnorm(0.80); z05 <- qnorm(0.975) + qnorm(0.80)
diag <- rbindlist(lapply(seq_len(nrow(programs)), function(i) {
  d <- prog_data(pfit, programs$cell_type[i], programs$module[i])
  sg <- sigma_cells(d); sx <- sxx_by(d)
  r <- cell_centred(d)
  info_gls <- (sx / sg^2) / sum(sx / sg^2); info_ols <- sx / sum(sx)
  se_exp <- 1 / sqrt(sum(sx / sg^2))
  arr <- sapply(FIT_DATASETS, function(ds) { s <- d$stage[d$dataset == ds]; exp(lfactorial(length(s)) - sum(lfactorial(table(s)))) })
  data.table(cell_type = programs$cell_type[i], module = programs$module[i], n_fit = nrow(d),
             sd_GSE202379 = sg[1], sd_GSE244832 = sg[2], sd_ratio_244832_202379 = sg[2] / sg[1],
             bf_p_between_datasets = brown_forsythe(r, d$dataset),
             sxx_GSE202379 = sx[1], sxx_GSE244832 = sx[2],
             info_gls_GSE202379 = info_gls[1], info_gls_GSE244832 = info_gls[2],
             info_ols_GSE202379 = info_ols[1], info_ols_GSE244832 = info_ols[2],
             se_expected = se_exp, mde_lone_call = zc * se_exp, mde_alpha05 = z05 * se_exp,
             arrangements_GSE202379 = round(arr[1]), arrangements_GSE244832 = round(arr[2]))
}))
wtsv(diag, "D2_D3_D4_D5_D6_diagnostics.tsv")
d2sum <- diag[, .(stat = c("min", "median", "max"), sd_ratio = c(min(sd_ratio_244832_202379), median(sd_ratio_244832_202379), max(sd_ratio_244832_202379)))]
wtsv(d2sum, "D2_sd_ratio_summary.tsv")
msg("D2 ratio min/median/max", paste(signif(d2sum$sd_ratio, 4), collapse = " / "),
    "; D3 count BF p<0.05:", sum(diag$bf_p_between_datasets < 0.05, na.rm = TRUE))
g2 <- merge(NAMED, diag[, .(cell_type, module, mde_lone_call, mde_alpha05)], by = c("cell_type", "module"))
g2[, g2_pass := mde_lone_call <= printed_beta]
wtsv(g2, "G2_mde_named.tsv")
add_gate("G2", all(g2$g2_pass), paste(sprintf("%s MDE %.4f vs printed %.6f", g2$label, g2$mde_lone_call, g2$printed_beta), collapse = "; "))
D8 <- rbindlist(lapply(seq_len(nrow(NAMED)), function(i) {
  d <- prog_data(pfit, NAMED$cell_type[i], NAMED$module[i])
  rbindlist(lapply(FIT_DATASETS, function(ds) { x <- d[dataset == ds]
    data.table(label = NAMED$label[i], dataset = ds, bf_p_across_stages = brown_forsythe(x$score, x$stage),
               sd_stage0 = sd(x$score[x$stage == 0]), sd_stage1 = sd(x$score[x$stage == 1]), sd_stage2 = sd(x$score[x$stage == 2])) }))
}))
wtsv(D8, "D8_named_within_dataset_stage_variance.tsv")

# ------------------------------------------------------------------ keys
K <- key_matrix(B)
msg("key matrix", nrow(K), "x", ncol(K), "first key", sprintf("%.10f", K[1, 1]))
donors62 <- dat$donors62
Kid <- matrix(seq_len(62) / 63, nrow = 1)   # identity: keys increase in donor-ID order

# ------------------------------------------------------- null GLS (primary)
null_one <- function(i, scores, sel_rows = NULL) {
  ct <- programs$cell_type[i]; mod <- programs$module[i]
  d <- prog_data(scores, ct, mod, FIT_DATASETS)
  if (nrow(d) == 0) return(NULL)
  obs <- fit_gls(d$score, d$stage, d$dataset)
  idp <- perm_stage(d, Kid, donors62)[1, ]
  idf <- fit_gls(d$score, idp, d$dataset)
  P <- perm_stage(d, K, donors62)
  res <- t(vapply(seq_len(nrow(P)), function(k) fit_gls(d$score, P[k, ], d$dataset), numeric(6)))
  list(obs_converged = is.finite(obs["t"]), identity_diff = abs(obs["t"] - idf["t"]),
       identity_diff_t3 = abs(obs["t3"] - idf["t3"]), identity_stage_equal = identical(as.numeric(idp), as.numeric(d$stage)),
       null = res, n = nrow(d))
}
msg("primary GLS null start")
t0 <- Sys.time()
NL <- mclapply(seq_len(nrow(programs)), null_one, scores = pfit, mc.cores = NCORES, mc.preschedule = FALSE)
msg("primary GLS null done in", format(Sys.time() - t0))
bad <- which(vapply(NL, function(x) inherits(x, "try-error") || is.null(x), TRUE))
if (length(bad)) stop("worker errors in programs: ", paste(bad, collapse = ","))
pid <- paste(programs$cell_type, programs$module, sep = "_")
Tn <- sapply(NL, function(x) x$null[, "t"]); T3n <- sapply(NL, function(x) x$null[, "t3"])
Bn <- sapply(NL, function(x) x$null[, "beta"])
if (!is.matrix(Tn)) { Tn <- matrix(Tn, nrow = 1); T3n <- matrix(T3n, nrow = 1); Bn <- matrix(Bn, nrow = 1) }
colnames(Tn) <- colnames(T3n) <- colnames(Bn) <- pid
nfit <- sapply(NL, `[[`, "n")
fwrite(as.data.table(Tn), file.path(OUT, "null_gls_t.tsv.gz"), sep = "\t")
fwrite(as.data.table(T3n), file.path(OUT, "null_wls_hc3_t3.tsv.gz"), sep = "\t")
fwrite(as.data.table(Bn), file.path(OUT, "null_gls_beta.tsv.gz"), sep = "\t")

# G4, G5, D6
g4 <- data.table(program = pid, identity_diff_t = sapply(NL, `[[`, "identity_diff"),
                 identity_diff_t3 = sapply(NL, `[[`, "identity_diff_t3"),
                 identity_stage_equal = sapply(NL, `[[`, "identity_stage_equal"),
                 obs_converged = sapply(NL, `[[`, "obs_converged"),
                 n_failed_draws = colSums(!is.finite(Tn)), n_failed_t3 = colSums(!is.finite(T3n)))
g4[, frac_failed := n_failed_draws / nrow(Tn)]
wtsv(g4, "G4_G5_D6_identity_convergence.tsv")
g4pass <- all(g4$identity_stage_equal) && all(is.finite(g4$identity_diff_t)) && max(g4$identity_diff_t) <= 1e-10
add_gate("G4", g4pass, sprintf("max identity |dt| = %.3g; stage vectors identical %d/117", max(g4$identity_diff_t), sum(g4$identity_stage_equal)))
g5pass <- sum(!g4$obs_converged) <= 5 && max(g4$frac_failed) <= 0.01
add_gate("G5", g5pass, sprintf("observed non-converged %d; max failed-draw fraction %.4g; programs with 0-1%% failures %d",
                               sum(!g4$obs_converged), max(g4$frac_failed), sum(g4$n_failed_draws > 0)))
if (!g4pass && !NULLDATA) stop_rule("untestable: G4 identity permutation does not reproduce the observed t")
if (!g5pass && !NULLDATA) stop_rule("untestable: G5 convergence failure")

# ------------------------------------------------------- null OLS 60 and 62
ols_null <- function(scores) {
  sapply(seq_len(nrow(programs)), function(i) {
    d <- prog_data(scores, programs$cell_type[i], programs$module[i])
    ols_t_mat(d$score, perm_stage(d, K, donors62), d$dataset)$t
  })
}
O60 <- ols_null(p60); O62 <- ols_null(p62)
if (!is.matrix(O60)) { O60 <- matrix(O60, nrow = 1); O62 <- matrix(O62, nrow = 1) }
colnames(O60) <- colnames(O62) <- pid
rdf60 <- sapply(seq_len(nrow(programs)), function(i) nrow(prog_data(p60, programs$cell_type[i], programs$module[i])) - uniqueN(prog_data(p60, programs$cell_type[i], programs$module[i])$dataset) - 1)
rdf62 <- sapply(seq_len(nrow(programs)), function(i) nrow(prog_data(p62, programs$cell_type[i], programs$module[i])) - uniqueN(prog_data(p62, programs$cell_type[i], programs$module[i])$dataset) - 1)
fwrite(as.data.table(O60), file.path(OUT, "null_ols60_t.tsv.gz"), sep = "\t")
fwrite(as.data.table(O62), file.path(OUT, "null_ols62_t.tsv.gz"), sep = "\t")
wtsv(data.table(program = pid, n_fit_gls = nfit, rdf_ols60 = rdf60, rdf_ols62 = rdf62), "null_design_df.tsv")

# ------------------------------------------------------------------ NC1 / G6
nc1_draws <- min(2000L, nrow(Tn))
A <- abs(Tn)
pnull <- matrix(1, nc1_draws, ncol(A))
for (j in seq_len(ncol(A))) {
  v <- A[, j]; ok <- is.finite(v); s <- sort(v[ok]); ns <- length(s)
  for (k in seq_len(nc1_draws)) if (ok[k]) {
    ge <- ns - findInterval(v[k], s, left.open = TRUE) - 1   # others with |t*| >= |t*_k|
    pnull[k, j] <- (1 + ge) / ns                              # ns - 1 others, + 1
  }
}
for (j in seq_len(ncol(A))) if (!g1$g1_pass[j]) pnull[, j] <- 1
calls <- apply(pnull, 1, function(p) sum(bh117(p) < 0.05))
nc1 <- data.table(n_draws = nc1_draws, frac_draws_any_call = mean(calls > 0), max_calls = max(calls))
wtsv(nc1, "NC1_null_family.tsv")
add_gate("G6", nc1$frac_draws_any_call <= 0.10, sprintf("NC1: %.4f of %d null draws have >=1 BH call (prediction <= 0.06)", nc1$frac_draws_any_call, nc1_draws))
if (nc1$frac_draws_any_call > 0.10 && !NULLDATA) stop_rule("untestable: G6 permutation null gives too many calls")

# ------------------------------------------------------------------ NC2
nc2 <- data.table(program = pid,
                  gls_rate = colMeans(2 * pt(abs(Tn), rep(nfit - 3, each = nrow(Tn)), lower.tail = FALSE) < 0.05, na.rm = TRUE),
                  ols60_rate = colMeans(2 * pt(abs(O60), rep(rdf60, each = nrow(O60)), lower.tail = FALSE) < 0.05, na.rm = TRUE),
                  ols62_rate = colMeans(2 * pt(abs(O62), rep(rdf62, each = nrow(O62)), lower.tail = FALSE) < 0.05, na.rm = TRUE),
                  sd_ratio = diag$sd_ratio_244832_202379)
wtsv(nc2, "NC2_parametric_calibration.tsv")
nc2_anticons <- sum(nc2$gls_rate > 0.075) > 10
nc2s <- data.table(gls_median = median(nc2$gls_rate), n_gls_above_0.075 = sum(nc2$gls_rate > 0.075),
                   gls_parametric_label_anticonservative = nc2_anticons,
                   ols60_median_sdratio_gt1.5 = median(nc2[sd_ratio > 1.5, ols60_rate]),
                   ols60_median_sdratio_le1.5 = median(nc2[sd_ratio <= 1.5, ols60_rate]),
                   ductular_ols60 = nc2[program == "hepatocytes_20", ols60_rate],
                   ductular_ols62 = nc2[program == "hepatocytes_20", ols62_rate],
                   ductular_gls = nc2[program == "hepatocytes_20", gls_rate])
wtsv(nc2s, "NC2_summary.tsv")

# ------------------------------------------------------------------ PC3
dsn <- c(rep("GSE202379", 37), rep("GSE244832", 18))
stg <- c(rep(0, 4), rep(1, 7), rep(2, 26), rep(0, 5), rep(1, 4), rep(2, 9))
sdv <- ifelse(dsn == "GSE202379", 0.33, 0.95)
nsim <- if (B < B_FULL) 50L else 2000L
set.seed(SEED)
sims <- lapply(c(0, 0.30), function(b) lapply(seq_len(nsim), function(r) b * stg + rnorm(55) * sdv))
pc3 <- rbindlist(lapply(1:2, function(s) {
  tb <- c(0, 0.30)[s]
  r <- mclapply(sims[[s]], function(y) {
    g <- fit_gls(y, stg, dsn)
    o <- ols_t_mat(y, matrix(stg, 1), dsn)
    c(g_beta = g[["beta"]], g_se = g[["se"]], o_beta = o$beta, o_se = o$se, o_se3 = o$beta / o$t3)
  }, mc.cores = NCORES)
  r <- do.call(rbind, r); tq <- qt(0.975, 52)
  rbindlist(list(
    data.table(method = "OLS", beta_hat = r[, "o_beta"], se = r[, "o_se"]),
    data.table(method = "OLS_HC3", beta_hat = r[, "o_beta"], se = r[, "o_se3"]),
    data.table(method = "GLS_H59", beta_hat = r[, "g_beta"], se = r[, "g_se"])))[,
    .(true_beta = tb, n_draws = .N, n_failed = sum(!is.finite(beta_hat)), bias = mean(beta_hat - tb, na.rm = TRUE),
      coverage95 = mean(abs(beta_hat - tb) <= tq * se, na.rm = TRUE),
      reject05 = mean(2 * pt(abs(beta_hat / se), 52, lower.tail = FALSE) < 0.05, na.rm = TRUE)), by = method]
}))
wtsv(pc3, "PC3_known_truth_simulation.tsv")
g3cov <- pc3[method == "GLS_H59", min(coverage95)]
add_gate("G3", g3cov >= 0.92, sprintf("GLS 95%% CI coverage min over true beta 0 and 0.30 = %.4f (threshold 0.92)", g3cov))

# ------------------------------------------------------------------ PC4
crit <- apply(abs(Tn), 2, function(v) quantile(v, 1 - 0.05 / 117, na.rm = TRUE, type = 7))
wtsv(data.table(program = pid, null_abs_t_crit = crit), "null_critical_abs_t.tsv")
npc4 <- min(1000L, nrow(K))
pc4 <- rbindlist(lapply(seq_len(nrow(NAMED)), function(i) {
  d <- prog_data(pfit, NAMED$cell_type[i], NAMED$module[i], FIT_DATASETS)
  P <- perm_stage(d, K[seq_len(npc4), , drop = FALSE], donors62)
  cr <- crit[paste(NAMED$cell_type[i], NAMED$module[i], sep = "_")]
  rbindlist(lapply(c(0.20, 0.30, NAMED$printed_beta[i]), function(delta) {
    tt <- unlist(mclapply(seq_len(npc4), function(k) fit_gls(d$score + delta * P[k, ], P[k, ], d$dataset)[["t"]], mc.cores = NCORES))
    data.table(label = NAMED$label[i], delta = delta, n_draws = npc4, n_failed = sum(!is.finite(tt)),
               crit_abs_t = cr, detection = mean(abs(tt) > cr, na.rm = TRUE), median_t = median(tt, na.rm = TRUE))
  }))
}))
wtsv(pc4, "PC4_detection.tsv")

# ------------------------------------------------- E4 null (ambient corrected)
ambs <- fread(inp("ambient_scores"))
ambs[, module := as.integer(module)]
chk <- merge(p60[, .(cell_type, module, donor, score)], ambs[, .(cell_type, module, donor, raw, corrected)],
             by = c("cell_type", "module", "donor"), all.x = TRUE)
e4_maxdiff <- max(abs(chk$raw - chk$score))
e4_run <- !anyNA(chk$raw) && e4_maxdiff <= 1e-8
e4id <- data.table(rows = nrow(chk), raw_missing = sum(is.na(chk$raw)), max_abs_raw_minus_primary = e4_maxdiff,
                   corrected_missing = sum(is.na(chk$corrected)),
                   corrected_missing_tcells = chk[is.na(corrected) & cell_type == "tcells", .N], e4_runs = e4_run)
wtsv(e4id, "E4_identity_check.tsv")
msg("E4 identity", e4_maxdiff, "run", e4_run)
if (e4_run) {
  pcorr <- merge(chk[, .(cell_type, module, donor, score = corrected)], unique(p60[, .(donor, dataset, stage)]), by = "donor")
  pcorr <- pcorr[dataset %in% FIT_DATASETS & is.finite(score)]
  e4_one <- function(i) {
    d <- prog_data(pcorr, programs$cell_type[i], programs$module[i], FIT_DATASETS)
    if (nrow(d) < 8 || uniqueN(d$dataset) < 2) return(rep(NA_real_, nrow(K)))
    P <- perm_stage(d, K, donors62)
    vapply(seq_len(nrow(P)), function(k) fit_gls(d$score, P[k, ], d$dataset)[["t"]], 0)
  }
  msg("E4 null start")
  E4n <- do.call(cbind, mclapply(seq_len(nrow(programs)), e4_one, mc.cores = NCORES, mc.preschedule = FALSE))
  colnames(E4n) <- pid
  fwrite(as.data.table(E4n), file.path(OUT, "null_e4_corrected_gls_t.tsv.gz"), sep = "\t")
  msg("E4 null done")
}

wtsv(gate, "gates.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_stageA.txt"))
file.create(file.path(OUT, "STAGE_A_DONE"))
msg("stage A complete")
