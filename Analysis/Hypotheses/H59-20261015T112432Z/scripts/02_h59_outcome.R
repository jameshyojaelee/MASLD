# H59 stage B: observed statistics, permutation p, BH family, descriptives E1-E9,
# verdict by the prespec rule, RESULT_H59.md. Requires stage A outputs and
# CHECKED_BEFORE_OUTCOME.txt (written by the analyst after reading gates/controls).
RUN <- Sys.getenv("H59_RUN_DIR", "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H59-20261015T112432Z")
source(file.path(RUN, "scripts/h59_lib.R"))
NULLDATA <- Sys.getenv("H59_NULLDATA", "0") == "1"
ONAME <- Sys.getenv("H59_RESULT_SUBDIR", "outcome")
AOUT <- OUT
OUT <- file.path(AOUT, ONAME)
if (!file.exists(file.path(AOUT, "STAGE_A_DONE"))) stop("stage A not complete")
if (file.exists(file.path(AOUT, "STOP_RULE.txt"))) stop("stage A hit a stop rule")
if (!NULLDATA && !file.exists(file.path(RUN, "CHECKED_BEFORE_OUTCOME.txt"))) stop("CHECKED_BEFORE_OUTCOME.txt missing")
dir.create(OUT, recursive = FALSE)  # errors if it exists? no: warn; guard below
if (length(list.files(OUT))) stop("outcome dir not empty")
msg("stage B; OUT =", OUT, "NULLDATA =", NULLDATA)
hashes <- check_hashes()
dat <- load_data()
if (NULLDATA) {
  set.seed(1); e <- dat$elig
  e[, stage := sample(stage), by = dataset]
  m <- e[, .(donor, stage)]
  for (nm in c("p", "eq")) { x <- dat[[nm]]; x[, stage := NULL]; dat[[nm]] <- merge(x, m, by = "donor") }
}
programs <- dat$programs
pid <- paste(programs$cell_type, programs$module, sep = "_")
p62 <- dat$p; p60 <- p62[dataset != DROP_DATASET]; pfit <- p60[dataset %in% FIT_DATASETS]
eqfit <- dat$eq[dataset %in% FIT_DATASETS]
gate <- fread(file.path(AOUT, "gates.tsv"))
g1 <- fread(file.path(AOUT, "G1_min_n.tsv"))
diag <- fread(file.path(AOUT, "D2_D3_D4_D5_D6_diagnostics.tsv"))
nc2s <- fread(file.path(AOUT, "NC2_summary.tsv"))
pc4 <- fread(file.path(AOUT, "PC4_detection.tsv"))
pc3 <- fread(file.path(AOUT, "PC3_known_truth_simulation.tsv"))
nc1 <- fread(file.path(AOUT, "NC1_null_family.tsv"))
g2 <- fread(file.path(AOUT, "G2_mde_named.tsv"))
g4 <- fread(file.path(AOUT, "G4_G5_D6_identity_convergence.tsv"))
rnull <- function(f) { x <- as.matrix(fread(file.path(AOUT, f))); stopifnot(identical(colnames(x), pid)); x }
Tn <- rnull("null_gls_t.tsv.gz"); T3n <- rnull("null_wls_hc3_t3.tsv.gz")
O60 <- rnull("null_ols60_t.tsv.gz"); O62 <- rnull("null_ols62_t.tsv.gz")
Bdraw <- nrow(Tn)
gp <- function(g) gate[gate == g, pass]
G3 <- gp("G3")

perm_p <- function(obs, null) {
  vapply(seq_along(obs), function(j) {
    if (!is.finite(obs[j])) return(NA_real_)
    v <- abs(null[, j]); v <- v[is.finite(v)]
    (1 + sum(v >= abs(obs[j]) * (1 - 1e-12))) / (1 + length(v))
  }, 0)
}

# ------------------------------------------------------------ primary GLS
obs <- rbindlist(mclapply(seq_len(nrow(programs)), function(i) {
  d <- prog_data(pfit, programs$cell_type[i], programs$module[i], FIT_DATASETS)
  f <- fit_gls(d$score, d$stage, d$dataset)
  eqd <- prog_data(eqfit, programs$cell_type[i], programs$module[i], FIT_DATASETS)
  fe <- if (nrow(eqd)) fit_gls(eqd$score, eqd$stage, eqd$dataset) else c(beta = NA)
  s1 <- within_slope(d, "GSE202379"); s2 <- within_slope(d, "GSE244832")
  Q <- (s1[["slope"]] - s2[["slope"]])^2 / (s1[["se"]]^2 + s2[["se"]]^2)
  sxy <- sapply(FIT_DATASETS, function(ds) { x <- d[dataset == ds]; sum((x$stage - mean(x$stage)) * (x$score - mean(x$score))) })
  data.table(program = pid[i], cell_type = programs$cell_type[i], module = programs$module[i],
             module_name = programs$module_name[i], stability_median = programs$stability_median[i],
             n_fit = nrow(d), gls_beta = f[["beta"]], gls_se = f[["se"]], gls_t = f[["t"]], wls_hc3_t3 = f[["t3"]],
             sigma_GSE202379 = f[["s1"]], sigma_GSE244832 = f[["s2"]],
             gls_param_p = 2 * pt(abs(f[["t"]]), nrow(d) - 3, lower.tail = FALSE),
             gls_ci_low = f[["beta"]] - qt(0.975, nrow(d) - 3) * f[["se"]],
             gls_ci_high = f[["beta"]] + qt(0.975, nrow(d) - 3) * f[["se"]],
             equal_run_gls_beta = fe[["beta"]],
             slope_GSE202379 = s1[["slope"]], se_GSE202379 = s1[["se"]], p_GSE202379 = s1[["p"]], n_GSE202379 = s1[["n"]],
             slope_GSE244832 = s2[["slope"]], se_GSE244832 = s2[["se"]], p_GSE244832 = s2[["p"]], n_GSE244832 = s2[["n"]],
             cochran_Q = Q, cochran_Q_p = pchisq(Q, 1, lower.tail = FALSE), I2 = max(0, (Q - 1) / Q),
             ols_numerator_share_GSE202379 = sxy[1] / sum(sxy), ols_numerator_share_GSE244832 = sxy[2] / sum(sxy))
}, mc.cores = NCORES))
stopifnot(identical(obs$program, pid))
obs <- merge(obs, g1[, .(cell_type, module, g1_pass)], by = c("cell_type", "module"), sort = FALSE)
obs <- merge(obs, diag[, .(cell_type, module, info_gls_GSE202379, info_gls_GSE244832, info_ols_GSE202379,
                           info_ols_GSE244832, mde_lone_call, mde_alpha05, sd_ratio_244832_202379)],
             by = c("cell_type", "module"), sort = FALSE)
setorder(obs, cell_type, module); stopifnot(identical(obs$program, pid))
obs[, perm_p := perm_p(gls_t, Tn)]
obs[, perm_p_t3 := perm_p(wls_hc3_t3, T3n)]
obs[g1_pass == FALSE, `:=`(perm_p = 1, perm_p_t3 = 1)]
obs[, perm_q := bh117(perm_p)]
obs[, perm_q_t3 := bh117(perm_p_t3)]
mx <- apply(abs(Tn), 1, max, na.rm = TRUE)
obs[, maxT_fwer_p := vapply(abs(gls_t), function(a) if (is.finite(a)) (1 + sum(mx >= a * (1 - 1e-12))) / (1 + length(mx)) else NA_real_, 0)]
obs[, gls_param_q := bh117(gls_param_p)]
obs[, n_failed_null_draws := colSums(!is.finite(Tn))]
obs[, selected_H59 := is.finite(perm_q) & perm_q < 0.05 & is.finite(stability_median) & stability_median >= 0.5]
obs[, equal_run_sign_agree := is.finite(gls_beta) & is.finite(equal_run_gls_beta) & gls_beta != 0 & equal_run_gls_beta != 0 & sign(gls_beta) == sign(equal_run_gls_beta)]
obs[, robust_display_H59 := selected_H59 & is.finite(perm_q_t3) & perm_q_t3 < 0.05 & equal_run_sign_agree]

# ------------------------------------------------------------ E1 OLS 62 / 60
ols_obs <- function(scores, null, label) {
  r <- rbindlist(lapply(seq_len(nrow(programs)), function(i) {
    d <- prog_data(scores, programs$cell_type[i], programs$module[i])
    o <- ols_t_mat(d$score, matrix(d$stage, 1), d$dataset)
    data.table(program = pid[i], n = nrow(d), beta = o$beta, t = o$t, t3 = o$t3,
               p = 2 * pt(abs(o$t), o$rdf, lower.tail = FALSE), p_hc3 = 2 * pt(abs(o$t3), o$rdf, lower.tail = FALSE))
  }))
  r[, q := bh117(p)][, q_hc3 := bh117(p_hc3)]
  r[, perm_p := perm_p(t, null)][, perm_q := bh117(perm_p)]
  setnames(r, setdiff(names(r), "program"), paste0(label, "_", setdiff(names(r), "program")))
  r
}
e1 <- merge(ols_obs(p62, O62, "ols62"), ols_obs(p60, O60, "ols60"), by = "program", sort = FALSE)
obs <- merge(obs, e1, by = "program", sort = FALSE); setorder(obs, cell_type, module)

# ------------------------------------------------------------ E4 ambient
e4 <- fread(file.path(AOUT, "E4_identity_check.tsv"))
if (isTRUE(e4$e4_runs)) {
  E4n <- rnull("null_e4_corrected_gls_t.tsv.gz")
  ambs <- fread(inp("ambient_scores")); ambs[, module := as.integer(module)]
  pc <- merge(ambs[, .(cell_type, module, donor, score = corrected)], unique(pfit[, .(donor, dataset, stage)]), by = "donor")
  pc <- pc[is.finite(score)]
  e4o <- rbindlist(lapply(seq_len(nrow(programs)), function(i) {
    d <- prog_data(pc, programs$cell_type[i], programs$module[i], FIT_DATASETS)
    if (nrow(d) < 8 || uniqueN(d$dataset) < 2) return(data.table(program = pid[i], e4_n = nrow(d), e4_beta = NA_real_, e4_t = NA_real_))
    f <- fit_gls(d$score, d$stage, d$dataset)
    data.table(program = pid[i], e4_n = nrow(d), e4_beta = f[["beta"]], e4_t = f[["t"]])
  }))
  e4o[, e4_perm_p := perm_p(e4_t, E4n)][, e4_perm_q := bh117(e4_perm_p)]
  obs <- merge(obs, e4o, by = "program", sort = FALSE); setorder(obs, cell_type, module)
}

wtsv(obs, "results_primary_family.tsv")

# ------------------------------------------------------------ E7, E8
named_ids <- paste(NAMED$cell_type, NAMED$module, sep = "_")
e7_ids <- unique(c(named_ids, obs[selected_H59 == TRUE, program]))
e7 <- rbindlist(lapply(e7_ids, function(id) {
  i <- match(id, pid); d <- prog_data(pfit, programs$cell_type[i], programs$module[i], FIT_DATASETS)
  rbindlist(mclapply(seq_len(nrow(d)), function(r) {
    x <- d[-r]; f <- fit_gls(x$score, x$stage, x$dataset)
    data.table(program = id, deleted_donor = d$donor[r], dataset = d$dataset[r], stage = d$stage[r],
               beta = f[["beta"]], t = f[["t"]], param_p = 2 * pt(abs(f[["t"]]), nrow(x) - 3, lower.tail = FALSE))
  }, mc.cores = NCORES))
}))
wtsv(e7, "E7_single_donor_deletions.tsv")
e7s <- e7[, .(n_deletions = .N, beta_min = min(beta, na.rm = TRUE), beta_max = max(beta, na.rm = TRUE),
              n_param_p_above_lone_call = sum(param_p > 0.05 / 117, na.rm = TRUE)), by = program]
wtsv(e7s, "E7_summary.tsv")
e8 <- rbindlist(lapply(seq_len(nrow(NAMED)), function(i) {
  d <- prog_data(pfit, NAMED$cell_type[i], NAMED$module[i], FIT_DATASETS)
  df <- data.frame(y = d$score, stage = d$stage, dataset = factor(d$dataset), vg = interaction(d$dataset, d$stage, drop = TRUE))
  f <- tryCatch(gls(y ~ stage + dataset, data = df, weights = varIdent(form = ~ 1 | vg), method = "REML", control = GLS_CTRL),
                error = function(e) NULL)
  if (is.null(f)) return(data.table(label = NAMED$label[i], beta = NA_real_, param_p = NA_real_, converged = FALSE))
  tt <- summary(f)$tTable
  data.table(label = NAMED$label[i], beta = tt["stage", 1], param_p = tt["stage", 4], converged = TRUE)
}))
wtsv(e8, "E8_variance_by_dataset_and_stage.tsv")

# ------------------------------------------------------------ E9
e9rows <- rbindlist(lapply(setdiff(unique(programs$cell_type), "tcells"), function(ct) {
  x <- obs[cell_type == ct & is.finite(gls_beta)]
  rbind(x[gls_beta > 0][order(perm_q, perm_p, -abs(gls_t))][1:min(3, .N)][, side := "positive"],
        x[gls_beta < 0][order(perm_q, perm_p, -abs(gls_t))][1:min(3, .N)][, side := "negative"])
}))[, .(cell_type, module, module_name, side, gls_beta, perm_p, perm_q, selected_H59)]
wtsv(e9rows, "E9_fig4c_row_choice.tsv")
e9means <- p60[, .(mean_score = mean(score), n = .N), by = .(cell_type, module, stage)][order(cell_type, module, stage)]
wtsv(e9means, "E9_stage_means_60donors.tsv")

# ------------------------------------------------------------ verdict
gates_ok <- all(gate[gate %in% c("G0a", "G0b", "G1", "G4", "G5", "G6"), pass])
n_sel <- sum(obs$selected_H59)
verdict <- if (!gates_ok) "untestable" else if (n_sel >= 1) "pass" else "null"
named <- merge(NAMED, obs, by = c("cell_type", "module"))
named <- merge(named, pc4[abs(delta - 0.2) > 1e-9 & abs(delta - 0.3) > 1e-9, .(label, pc4_detection_printed = detection)], by = "label")
named <- merge(named, g2[, .(label, g2_pass)], by = "label")
named[, reading := fifelse(selected_H59, "selected",
                   fifelse(!g2_pass | pc4_detection_printed < 0.80, "not selected; cannot exclude the printed effect",
                           "not selected; printed effect would have been detected"))]
wtsv(named, "named_programs.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_stageB.txt"))

f6 <- function(x) formatC(x, digits = 6, format = "g")
sel <- obs[selected_H59 == TRUE]
L <- c(
  "# RESULT H59", "",
  if (NULLDATA) "DEBUG RUN ON STAGE-SHUFFLED DATA. NOT A RESULT." else NULL,
  sprintf("Run dir: %s. Prespec sha256 %s. Permutation draws: %d (seed %d).", RUN, PRESPEC_SHA, Bdraw, SEED), "",
  sprintf("**Verdict (computed by the prespec rule): %s.** %d of 117 programs selected (permutation BH q < 0.05 over 117 programs, n = 117 fixed, and stability_median >= 0.5); %d meet the robust-display rule.",
          toupper(verdict), n_sel, sum(obs$robust_display_H59)), "",
  "## Gates", "", "| Gate | Pass | Detail |", "|---|---|---|",
  sprintf("| %s | %s | %s |", gate$gate, gate$pass, gate$detail), "",
  "## Named programs (from the 117-program family)", "",
  "| Program | n | GLS beta (95% CI, model) | permutation p | BH q (117) | t3 permutation q | GSE202379 slope (SE) | GSE244832 slope (SE) | Cochran Q p | lone-call MDE | PC4 detection at printed beta | Reading |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %d | %.4f (%.4f to %.4f) | %s | %s | %s | %.4f (%.4f) | %.4f (%.4f) | %.3g | %.4f | %.3f | %s |",
          named$label, named$n_fit, named$gls_beta, named$gls_ci_low, named$gls_ci_high, f6(named$perm_p), f6(named$perm_q),
          f6(named$perm_q_t3), named$slope_GSE202379, named$se_GSE202379, named$slope_GSE244832, named$se_GSE244832,
          named$cochran_Q_p, named$mde_lone_call, named$pc4_detection_printed, named$reading), "",
  if (!G3) "G3 failed: model-based CIs must not be printed in any figure." else "G3 passed: model-based CIs may be shown.",
  "",
  "## Selected programs", "",
  if (nrow(sel)) c("| Program | name | n | GLS beta | perm q | t3 perm q | robust display | GSE202379 slope (SE) | GSE244832 slope (SE) | Q p |", "|---|---|---|---|---|---|---|---|---|---|",
                   sprintf("| %s | %s | %d | %.4f | %s | %s | %s | %.4f (%.4f) | %.4f (%.4f) | %.3g |", sel$program, sel$module_name, sel$n_fit,
                           sel$gls_beta, f6(sel$perm_q), f6(sel$perm_q_t3), sel$robust_display_H59, sel$slope_GSE202379, sel$se_GSE202379,
                           sel$slope_GSE244832, sel$se_GSE244832, sel$cochran_Q_p)) else "None.",
  "",
  "## Controls", "",
  sprintf("- NC1: %.4f of %d null draws gave at least one BH call (prediction <= 0.06).", nc1$frac_draws_any_call, nc1$n_draws),
  sprintf("- NC2: GLS parametric p < 0.05 under the permutation null, median rate %.4f; %d programs above 0.075 (anticonservative label: %s). Ductular injury OLS60 %.4f, OLS62 %.4f, GLS %.4f.",
          nc2s$gls_median, nc2s$n_gls_above_0.075, nc2s$gls_parametric_label_anticonservative, nc2s$ductular_ols60, nc2s$ductular_ols62, nc2s$ductular_gls),
  "- PC3 (known truth):", sprintf("  - %s, true beta %.2f: bias %.4f, coverage %.4f, rejection %.4f", pc3$method, pc3$true_beta, pc3$bias, pc3$coverage95, pc3$reject05),
  "- PC4 (detection):", sprintf("  - %s, delta %.6f: detection %.3f (critical |t| %.3f)", pc4$label, pc4$delta, pc4$detection, pc4$crit_abs_t), "",
  "## Descriptive (outside the family)", "",
  sprintf("- OLS 62 donors: %d calls at q < 0.05 (Ductular injury q %s, permutation q %s). OLS 60 donors: %d calls (Ductular injury q %s, permutation q %s).",
          sum(obs$ols62_q < .05), f6(obs[program == "hepatocytes_20", ols62_q]), f6(obs[program == "hepatocytes_20", ols62_perm_q]),
          sum(obs$ols60_q < .05), f6(obs[program == "hepatocytes_20", ols60_q]), f6(obs[program == "hepatocytes_20", ols60_perm_q])),
  sprintf("- GLS parametric BH q < 0.05: %d programs (parametric p label anticonservative: %s).", sum(obs$gls_param_q < .05), nc2s$gls_parametric_label_anticonservative),
  sprintf("- maxT family-wise p < 0.05: %d programs; smallest %s.", sum(obs$maxT_fwer_p < .05, na.rm = TRUE), f6(min(obs$maxT_fwer_p, na.rm = TRUE))),
  if ("e4_perm_q" %in% names(obs)) sprintf("- E4 ambient-corrected GLS: %d of 106 computable programs at permutation q < 0.05; Ductular injury beta %.4f q %s; Stromal ECM beta %.4f q %s.",
          sum(obs$e4_perm_q < .05, na.rm = TRUE), obs[program == "hepatocytes_20", e4_beta], f6(obs[program == "hepatocytes_20", e4_perm_q]),
          obs[program == "hepatocytes_8", e4_beta], f6(obs[program == "hepatocytes_8", e4_perm_q])) else "- E4 skipped (identity check failed).",
  sprintf("- E7 deletions: %s", paste(sprintf("%s beta %.4f to %.4f, %d of %d deletions with parametric p > 0.05/117", e7s$program, e7s$beta_min, e7s$beta_max, e7s$n_param_p_above_lone_call, e7s$n_deletions), collapse = "; ")),
  sprintf("- E8 variance by dataset x stage: %s", paste(sprintf("%s beta %.4f p %s", e8$label, e8$beta, f6(e8$param_p)), collapse = "; ")),
  sprintf("- Smallest permutation q over 117: %s (%s).", f6(min(obs$perm_q)), obs$program[which.min(obs$perm_q)]),
  "", "Units: donors (60-donor universe; fitted donors in GSE202379 and GSE244832). Stage is cross-sectional ordinal 0/1/2.",
  "Predicted (prespec): Ductular injury GLS beta 0.12 to 0.17, Stromal ECM 0.07 to 0.13, neither selected.")
writeLines(L, file.path(if (NULLDATA) OUT else RUN, "RESULT_H59.md"))
msg("verdict", verdict, "n_sel", n_sel)
