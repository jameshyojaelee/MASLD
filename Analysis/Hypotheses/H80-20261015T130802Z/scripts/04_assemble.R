# 04_assemble.R -- H80 job 4: diagnostics in prespecified order, stop rules S1-S7, primary statistics,
# window bootstrap (stream 3000), controls, readout, figures, RESULT_H80.md with the verdict computed here.
args0 <- commandArgs(trailingOnly = FALSE)
SCRIPTS <- dirname(normalizePath(sub("^--file=", "", args0[grep("^--file=", args0)])))
Sys.setenv(H80_SCRIPTS = SCRIPTS)
source(file.path(SCRIPTS, "h80_lib.R"))
suppressPackageStartupMessages(library(ggplot2))
A <- parse_args(list(setup = NA, arm_a = NA, arm_b = NA, readout = NA, out = NA, seed = "20260923",
                     result = "", n_boot = "2000", skip_preamble = "FALSE"))
H80 <- dirname(SCRIPTS); SETUP <- A$setup; OUT <- A$out
RESULT <- if (nzchar(A$result)) A$result else file.path(H80, "RESULT_H80.md")
NBOOT <- as.integer(A$n_boot)
stopifnot(as.integer(A$seed) == H80_SEED)
dir.create(file.path(OUT, "figures"), showWarnings = FALSE)
log_msg("H80 assembly start")

diag_lines <- character(0)
D <- function(...) { s <- paste0(...); diag_lines <<- c(diag_lines, s); log_msg(s) }
f3 <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "f"))
stop_rules <- list()
add_stop <- function(rule, text) stop_rules[[length(stop_rules) + 1L]] <<- data.table(rule = rule, text = text)

# ---- S1 again at assembly (inputs unchanged since setup) ----
pre <- if (identical(A$skip_preamble, "TRUE")) "PASS" else tryCatch(job_preamble(H80, SETUP, OUT, cores = 4L),
                                                                     error = function(e) paste0("ERROR:", conditionMessage(e)))
setup_status <- if (file.exists(file.path(SETUP, "STATUS"))) readLines(file.path(SETUP, "STATUS"))[1] else "MISSING"
job_state <- function(d) if (is.na(d) || !dir.exists(d)) "missing" else if (file.exists(file.path(d, "DONE"))) "done" else
  if (file.exists(file.path(d, "SKIPPED"))) "skipped" else "incomplete"
jobs <- c(arm_a = job_state(A$arm_a), arm_b = job_state(A$arm_b), readout = job_state(A$readout))
D("Job states: setup ", setup_status, "; Arm A ", jobs[["arm_a"]], "; Arm B + PC-1 ", jobs[["arm_b"]], "; readout ", jobs[["readout"]])

# ---- Diagnostic 1: hashes ----
hash_fail <- character(0)
for (d in c(SETUP, A$arm_a, A$arm_b, A$readout, OUT)) {
  f <- file.path(d, "input_hash_check.tsv")
  if (!is.na(d) && file.exists(f)) {
    h <- fread(f)
    D("D1 ", basename(d), ": ", sum(h$pass), " of ", nrow(h), " input hashes match the prespec")
    if (!all(h$pass)) hash_fail <- c(hash_fail, paste0(basename(d), ":", h[pass == FALSE, role]))
  }
}
for (d in c(SETUP, A$arm_a, A$arm_b, A$readout)) {
  f <- file.path(d, "opened_files_sha256.tsv")
  if (!is.na(d) && file.exists(f)) { h <- fread(f); D("D1 ", basename(d), ": sha256 of ", nrow(h), " opened files written (", round(sum(h$bytes) / 1e9, 1), " GB)") }
}
if (length(hash_fail) > 0 || startsWith(setup_status, "STOP:S1") || identical(pre, "STOP:S1")) {
  add_stop("S1", paste("input hash differs:", paste(c(hash_fail, setup_status), collapse = "; ")))
}
# ---- Diagnostic 2: G1 ----
sst <- if (file.exists(file.path(SETUP, "setup_status.tsv"))) fread(file.path(SETUP, "setup_status.tsv"), colClasses = "character") else data.table()
if (nrow(sst)) for (i in seq_len(nrow(sst))) D("D2 ", sst$check[i], " = ", sst$value[i], " (expected ", sst$expected[i], ") ", sst$pass[i])
if (startsWith(setup_status, "STOP:S3")) add_stop("S3", setup_status)
if (startsWith(setup_status, "STOP:S2") || startsWith(setup_status, "RESTART")) add_stop("S2", setup_status)
if (!setup_status %in% c("PASS") && length(stop_rules) == 0) add_stop("setup", paste("setup did not pass:", setup_status))

W <- if (file.exists(file.path(SETUP, "windows.tsv"))) fread(file.path(SETUP, "windows.tsv")) else data.table()
P1 <- if (file.exists(file.path(SETUP, "pc1_windows.tsv"))) fread(file.path(SETUP, "pc1_windows.tsv")) else data.table()
if (nrow(W)) D("D2 windows.tsv: ", nrow(W), " windows; ", paste(W[, .N, by = stratum][, paste0(stratum, "=", N)], collapse = " "),
               "; |E| median ", median(W$nE), "; |E and G| median ", median(W$nEG_plan))

blocked <- length(stop_rules) > 0
readA <- function(f) { p <- file.path(A$arm_a, f); if (file.exists(p)) fread(p) else NULL }
readB <- function(f) { p <- file.path(A$arm_b, f); if (file.exists(p)) fread(p) else NULL }
if (!blocked && jobs[["arm_a"]] != "done") { add_stop("G2", paste("Arm A job state:", jobs[["arm_a"]])); blocked <- TRUE }

STRATA <- c("EUR", "EUR_Seval", "EUR_Snores", "EAS", "AFR", "AMR", "SAS")
in_stratum <- function(st, s) if (s == "EUR") st %in% c("EUR_Seval", "EUR_Snores") else st == s
ERR_STATES <- c("fit_error", "coloc_error")
res <- list()

if (!blocked) {
  cellsA <- readA("cells.tsv.gz"); fitsA <- readA("fits.tsv.gz"); repsA <- readA("reps.tsv.gz")
  dgA <- readA("window_diagnostics.tsv"); wsA <- readA("window_status.tsv")
  cellsB <- readB("armB_cells.tsv.gz"); fitsB <- readB("armB_fits.tsv.gz"); repsB <- readB("armB_reps.tsv.gz"); dgB <- readB("armB_window_diagnostics.tsv")
  cellsP <- readB("pc1_cells.tsv.gz"); fitsP <- readB("pc1_fits.tsv.gz"); repsP <- readB("pc1_reps.tsv.gz"); dgP <- readB("pc1_window_diagnostics.tsv")
  wsB <- readB("window_status.tsv")

  # ---- Diagnostic 3: RNG ----
  rng_check <- function(reps, fits, cells, tag) {
    if (is.null(reps)) return(list(ok = NA, text = paste(tag, "not run")))
    reps <- copy(reps)
    reps[, dkey := if ("c_id" %in% names(reps)) paste(c_id, epsE_hash) else epsE_hash]
    dup <- reps[, .(n_dup = sum(duplicated(dkey))), by = row]
    n_dup_win <- sum(dup$n_dup > 0)
    n_stream_moved <- sum(!reps$stream_unchanged_by_fits)
    n_kind <- sum(reps$kind_after != "L'Ecuyer-CMRG")
    n_unrestored <- sum(!fits$restored, na.rm = TRUE) + sum(cells$coloc_restored %in% FALSE)
    drew <- fits[, .(n_fits = .N, n_drew = sum(drew_rng, na.rm = TRUE)), by = .(side, lambda, causal_type)]
    list(ok = n_dup_win == 0 && n_stream_moved == 0 && n_kind == 0 && n_unrestored == 0, drew = drew,
         text = paste0(tag, ": windows with duplicate replicates ", n_dup_win, "; replicates whose stream moved during fits ",
                       n_stream_moved, "; RNG kind not L'Ecuyer-CMRG after a fit ", n_kind, "; wrapped calls not restored ",
                       n_unrestored, "; fits whose purity step drew random numbers ", sum(drew$n_drew), " of ", sum(drew$n_fits)))
  }
  r3 <- list(A = rng_check(repsA, fitsA, cellsA, "Arm A"), B = rng_check(repsB, fitsB, cellsB, "Arm B"),
             P = rng_check(repsP, fitsP, cellsP, "PC-1"))
  for (x in r3) D("D3 ", x$text)
  fwrite(r3$A$drew, file.path(OUT, "D3_purity_rng_draws_armA.tsv"), sep = "\t")
  if (any(vapply(r3, function(x) isFALSE(x$ok), logical(1)))) add_stop("S6", "diagnostic 3: duplicate replicates or a moved/reset stream")

  # ---- Diagnostic 4 ----
  d4 <- file.path(SETUP, "D4_notconverged_summary_genes.tsv")
  if (file.exists(d4)) { x <- fread(d4); D("D4 genes with a not-converged last summary row (positive/tested): ", nrow(x), "; rds converged ", sum(x$rds_converged)) }
  # ---- Diagnostic 5 ----
  D("D5 Arm A windows: ", nrow(dgA), "; LD NA entries total ", sum(dgA$na_E) + sum(dgA$na_G), "; min eigenvalue of R + 1e-3 I (E side) ",
    signif(min(dgA$min_eig_E), 3), "; eigenvalues clipped to 0: E ", sum(dgA$n_clip_E), ", G ", sum(dgA$n_clip_G),
    "; E SNP ids identical to rds column names in ", sum(dgA$ids_E_identical_rds), " of ", nrow(dgA))
  # ---- Diagnostic 6: NC-3 ----
  dgA[, nc3_pass := nc3_T_pass & nc3_causal_pass]
  D("D6 NC-3 (Arm A): |T_w| <= 4.5 in ", sum(dgA$nc3_T_pass), " of ", nrow(dgA), " windows (max |T_w| ", f3(max(abs(dgA$T_w))),
    "); mean causal noise within +/- 4/sqrt(10) in ", sum(dgA$nc3_causal_pass), " of ", nrow(dgA),
    "; effective df min ", round(min(dgA$eff_df_E)), ", windows below 100 df ", sum(dgA$eff_df_E < 100))
  if (!is.null(dgB)) D("D6 NC-3 (Arm B, descriptive): max |T_w| matched ", f3(max(abs(dgB$T_w_matched))), ", mismatched ", f3(max(abs(dgB$T_w_mismatched))))
  if (!all(dgA$nc3_pass)) add_stop("S6", paste0("NC-3 failed in ", sum(!dgA$nc3_pass), " Arm A windows (rows ",
                                                 paste(dgA[nc3_pass == FALSE, row], collapse = ","), ")"))
  # ---- Diagnostic 7 ----
  poolA <- repsA[, .(min_pool = min(pool_size), n_na = sum(is.na(k_id))), by = row]
  D("D7 distinct-candidate pool: windows with a replicate below 20 SNPs ", sum(poolA$n_na > 0), "; median pool ",
    median(repsA$pool_size))
  # ---- Diagnostic 8 ----
  errs <- rbind(fitsA[!is.na(error), .(arm = "A", error)], if (!is.null(fitsB)) fitsB[!is.na(error), .(arm = "B", error)],
                if (!is.null(fitsP)) fitsP[!is.na(error), .(arm = "PC1", error)])
  D("D8 fits: Arm A ", nrow(fitsA), " (median ", f3(median(fitsA$time_s), 1), " s, max ", f3(max(fitsA$time_s), 1), " s); errors ", nrow(errs),
    if (nrow(errs)) paste0(" [", paste(errs[, .N, by = .(arm, error)][, paste0(arm, ":", error, "=", N)], collapse = "; "), "]") else "",
    "; peak RSS per window max ", f3(max(dgA$peak_rss_gb), 1), " GB; window failures Arm A ", sum(!wsA$ok),
    if (!is.null(wsB)) paste0(", Arm B/PC-1 ", sum(!wsB$ok)) else "")
  if (sum(!wsA$ok)) D("D8 Arm A window failures: ", paste(wsA[ok == FALSE, paste0(row, ":", error)], collapse = "; "))
  if (!is.null(wsB) && sum(!wsB$ok)) D("D8 Arm B/PC-1 failures: ", paste(wsB[ok == FALSE, paste0(task, ":", error)], collapse = "; "))
  blocked <- any(vapply(stop_rules, function(x) x$rule %in% c("S1", "S2", "S3", "S6"), logical(1)))
}

# ---- Diagnostic 9: PC-2 and PC-3 (from setup) ----
pc2 <- if (file.exists(file.path(SETUP, "pc2_results.tsv"))) fread(file.path(SETUP, "pc2_results.tsv")) else NULL
pc3 <- if (file.exists(file.path(SETUP, "pc3_results.tsv"))) fread(file.path(SETUP, "pc3_results.tsv")) else NULL
if (!is.null(pc2)) {
  pc2[, pass := ids_identical & max_pip_diff < 1e-4 & max_alpha_diff < 1e-4 & cs_det_stored == cs_det_rerun]
  D("D9 PC-2: ", sum(pc2$pass), " of ", nrow(pc2), " pass; max |PIP diff| ", signif(max(pc2$max_pip_diff, na.rm = TRUE), 3),
    "; max |alpha diff| ", signif(max(pc2$max_alpha_diff, na.rm = TRUE), 3), "; deterministic CS count equal ",
    sum(pc2$cs_det_stored == pc2$cs_det_rerun), "/", nrow(pc2), "; stored sets$cs count equal ", sum(pc2$cs_raw_stored == pc2$cs_raw_rerun), "/", nrow(pc2))
}
if (!is.null(pc3)) {
  for (pn in unique(pc3$panel)) D("D9 PC-3 (", pn, "): ", sum(pc3[panel == pn, pass]), " of ", nrow(pc3[panel == pn]),
                                   " pairs reproduce PP.H4 and PP.H3 within 1e-6; max |diff| PP.H4 ", signif(max(pc3[panel == pn, d4], na.rm = TRUE), 3),
                                   ", PP.H3 ", signif(max(pc3[panel == pn, d3], na.rm = TRUE), 3))
}

# ======================== primary statistics (only when not stopped) ========================
boot_idx <- list()
lam_cross <- function(lams, rec, thr) {
  j <- which(rec >= thr)[1]
  if (is.na(j)) return(list(v = max(lams), cens = "high"))
  if (j == 1L) return(list(v = min(lams), cens = "low"))
  v <- lams[j - 1] + (thr - rec[j - 1]) / (rec[j] - rec[j - 1]) * (lams[j] - lams[j - 1])
  list(v = v, cens = "none")
}
fmt_lam <- function(x, lo_edge, hi_edge) if (x$cens == "low") paste0("< ", lo_edge) else if (x$cens == "high") paste0("> ", hi_edge) else f3(x$v, 2)
ci_lam <- function(vals, cens, lo_edge, hi_edge) {
  q <- quantile(vals, c(0.025, 0.975), names = FALSE, type = 7)
  lo <- if (mean(cens == "low") > 0.025) paste0("< ", lo_edge) else f3(q[1], 2)
  hi <- if (mean(cens == "high") > 0.025) paste0("> ", hi_edge) else f3(q[2], 2)
  list(lo = lo, hi = hi, n_low = sum(cens == "low"), n_high = sum(cens == "high"), q = q)
}
pci <- function(M) apply(M, 2, function(v) quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 7))

if (!blocked) {
  # ---- G2 completion ----
  unit_fit_ok <- fitsA[, .(fit_ok = all(is.na(error))), by = .(row, rep)]
  unit_cell_ok <- cellsA[, .(cell_ok = !any(state %in% ERR_STATES)), by = .(row, rep)]
  units <- merge(unit_fit_ok, unit_cell_ok, by = c("row", "rep"))
  units[, ok := fit_ok & cell_ok]
  units <- merge(units, W[, .(row, stratum)], by = "row")
  g2 <- rbindlist(lapply(STRATA, function(s) {
    planned <- N_REP * sum(in_stratum(W$stratum, s))
    done <- sum(units[in_stratum(stratum, s), ok])
    data.table(stratum = s, planned_units = planned, complete_units = done, frac = done / planned, testable = done / planned >= 0.8)
  }))
  fwrite(g2, file.path(OUT, "G2_completion.tsv"), sep = "\t")
  for (i in seq_len(nrow(g2))) D("G2 ", g2$stratum[i], ": ", g2$complete_units[i], " of ", g2$planned_units[i], " units complete (",
                                  f3(g2$frac[i]), ")", if (!g2$testable[i]) " UNTESTABLE (S4)" else "")
  for (s in g2[testable == FALSE, stratum]) add_stop("S4", paste0("stratum ", s, " completed ", f3(g2[stratum == s, frac]), " of planned units"))

  # ---- per-window, per-cell rates (errored units excluded; distinct cells excluded where not applicable) ----
  cA <- merge(cellsA, units[, .(row, rep, ok)], by = c("row", "rep"))
  cA <- cA[ok == TRUE & state != "distinct_not_applicable"]
  cA[, cell := paste0(cell_set, "_", lam_e, "_", lam_g)]
  cA <- merge(cA, W[, .(row, stratum)], by = "row")
  wr <- cA[, .(recall = mean(state == "positive"), h3 = mean(!is.na(pp3) & pp3 > 0.8), n_rep = .N), by = .(row, cell, cell_set, lam_e, lam_g, causal)]
  wr <- merge(wr, unique(cA[, .(row, stratum)]), by = "row")
  cell_ids <- unique(wr[, .(cell, cell_set, lam_e, lam_g, causal)])
  setorder(cell_ids, cell_set, lam_e, lam_g)

  # ---- bootstrap indices, stream 3000, fixed order ----
  streams <- readRDS(file.path(SETUP, "streams.rds"))
  set_stream(streams, 3000L)
  win_sets <- list()
  for (s in STRATA) win_sets[[s]] <- sort(unique(wr[in_stratum(stratum, s), row]))
  b_rows <- if (!is.null(cellsB)) sort(unique(cellsB$row)) else integer(0)
  p_eur <- if (!is.null(cellsP)) sort(unique(merge(cellsP, P1[, .(row, ancestry)], by = "row")[ancestry == "EUR", row])) else integer(0)
  p_all <- if (!is.null(cellsP)) sort(unique(cellsP$row)) else integer(0)
  for (s in c(STRATA, "ArmB", "PC1_EUR", "PC1_pooled")) {
    ws <- switch(s, ArmB = b_rows, PC1_EUR = p_eur, PC1_pooled = p_all, win_sets[[s]])
    n <- length(ws)
    boot_idx[[s]] <- list(rows = ws, idx = if (n > 0) t(vapply(seq_len(NBOOT), function(b) sample.int(n, n, replace = TRUE), integer(n))) else NULL)
  }
  saveRDS(boot_idx, file.path(OUT, "bootstrap_indices.rds"))

  mat_for <- function(s, value) {
    ws <- boot_idx[[s]]$rows
    M <- matrix(NA_real_, length(ws), nrow(cell_ids), dimnames = list(ws, cell_ids$cell))
    x <- wr[row %in% ws]
    M[cbind(match(x$row, ws), match(x$cell, cell_ids$cell))] <- x[[value]]
    M
  }
  boot_means <- function(M, idx) t(vapply(seq_len(nrow(idx)), function(b) colMeans(M[idx[b, ], , drop = FALSE], na.rm = TRUE), numeric(ncol(M))))
  tab <- list(); lam_tab <- list(); boot_store <- list()
  E_GRID <- c(3, 4, 5, 6, 7, 8, 10, 15); G_GRID <- c(3, 4, 5, 6, 8, 10, 15)
  e_cells <- paste0("E_", E_GRID, "_10"); g_cells <- ifelse(G_GRID == 10, "E_10_10", paste0("G_10_", G_GRID))
  for (s in STRATA) {
    if (!g2[stratum == s, testable] || length(boot_idx[[s]]$rows) == 0) next
    Mr <- mat_for(s, "recall"); Mh <- mat_for(s, "h3")
    Br <- boot_means(Mr, boot_idx[[s]]$idx); Bh <- boot_means(Mh, boot_idx[[s]]$idx)
    boot_store[[s]] <- Br
    est_r <- colMeans(Mr, na.rm = TRUE); est_h <- colMeans(Mh, na.rm = TRUE)
    cr <- pci(Br); ch <- pci(Bh)
    st_comp <- cA[in_stratum(stratum, s), .N, by = .(cell, state)]
    st_comp <- dcast(st_comp, cell ~ state, value.var = "N", fill = 0)
    tb <- data.table(stratum = s, cell = cell_ids$cell, cell_set = cell_ids$cell_set, lam_e = cell_ids$lam_e, lam_g = cell_ids$lam_g,
                     causal = cell_ids$causal, n_windows = colSums(!is.na(Mr)),
                     rate_pp4_gt05 = est_r, ci_lo = cr[1, ], ci_hi = cr[2, ], rate_pp3_gt08 = est_h, h3_ci_lo = ch[1, ], h3_ci_hi = ch[2, ])
    tab[[s]] <- merge(tb, st_comp, by = "cell", all.x = TRUE)
    for (sw in c("E", "G")) {
      grid <- if (sw == "E") E_GRID else G_GRID; cc <- if (sw == "E") e_cells else g_cells
      for (thr in c(0.5, 0.8)) {
        pt <- lam_cross(grid, est_r[cc], thr)
        bs <- lapply(seq_len(nrow(Br)), function(b) lam_cross(grid, Br[b, cc], thr))
        vals <- vapply(bs, `[[`, numeric(1), "v"); cens <- vapply(bs, `[[`, character(1), "cens")
        ci <- ci_lam(vals, cens, min(grid), max(grid))
        lam_tab[[length(lam_tab) + 1L]] <- data.table(stratum = s, sweep = sw, threshold = thr, estimate = pt$v, censored = pt$cens,
                                                      printed = fmt_lam(pt, min(grid), max(grid)), ci_lo = ci$lo, ci_hi = ci$hi,
                                                      boot_censored_low = ci$n_low, boot_censored_high = ci$n_high)
      }
    }
  }
  tabA <- rbindlist(tab, fill = TRUE); lamA <- rbindlist(lam_tab)
  fwrite(tabA, file.path(OUT, "armA_cells_by_stratum.tsv"), sep = "\t")
  fwrite(lamA, file.path(OUT, "armA_lambda50_80.tsv"), sep = "\t")
  getr <- function(s, cl, col = "rate_pp4_gt05") { x <- tabA[stratum == s & cell == cl][[col]]; if (length(x)) x else NA_real_ }

  # ---- Arm B ----
  armB <- NULL
  if (!is.null(cellsB) && length(b_rows) > 0) {
    unitB <- merge(fitsB[, .(fit_ok = all(is.na(error))), by = .(row, rep)],
                   cellsB[, .(cell_ok = !any(state %in% ERR_STATES)), by = .(row, rep)], by = c("row", "rep"))
    cB <- merge(cellsB, unitB[fit_ok & cell_ok, .(row, rep)], by = c("row", "rep"))
    cB <- cB[state != "distinct_not_applicable"]
    cB[, cell := paste0(cell_set, "_", lam_e, "_", lam_g)]
    wb <- cB[, .(rate = mean(state == "positive")), by = .(row, cell, generator)]
    wb <- dcast(wb, row + cell ~ generator, value.var = "rate")
    wb[, diff := mismatched - matched]
    cellsBid <- sort(unique(wb$cell))
    rowsB <- boot_idx[["ArmB"]]$rows
    mk <- function(v) { M <- matrix(NA_real_, length(rowsB), length(cellsBid), dimnames = list(rowsB, cellsBid))
      M[cbind(match(wb$row, rowsB), match(wb$cell, cellsBid))] <- wb[[v]]; M }
    out <- list()
    for (v in c("matched", "mismatched", "diff")) {
      M <- mk(v); Bm <- boot_means(M, boot_idx[["ArmB"]]$idx); ci <- pci(Bm)
      out[[v]] <- data.table(cell = cellsBid, quantity = v, estimate = colMeans(M, na.rm = TRUE), ci_lo = ci[1, ], ci_hi = ci[2, ],
                             n_windows = colSums(!is.na(M)))
    }
    armB <- rbindlist(out)
    D("Arm B: ", length(rowsB), " windows; units complete ", nrow(unitB[fit_ok & cell_ok]), " of ", N_REP * sum(W$armB))
    fwrite(armB, file.path(OUT, "armB_cells.tsv"), sep = "\t")
  }

  # ---- PC-1 ----
  pc1 <- NULL
  if (!is.null(cellsP) && length(p_all) > 0) {
    unitP <- merge(fitsP[, .(fit_ok = all(is.na(error))), by = .(row, rep)],
                   cellsP[, .(cell_ok = !any(state %in% ERR_STATES)), by = .(row, rep)], by = c("row", "rep"))
    cP <- merge(cellsP, unitP[fit_ok & cell_ok, .(row, rep)], by = c("row", "rep"))
    wp <- cP[, .(recall = mean(state == "positive")), by = row]
    pc1 <- list()
    for (s in c("PC1_EUR", "PC1_pooled")) {
      rr <- boot_idx[[s]]$rows
      v <- wp$recall[match(rr, wp$row)]
      bm <- vapply(seq_len(NBOOT), function(b) mean(v[boot_idx[[s]]$idx[b, ]], na.rm = TRUE), numeric(1))
      q <- quantile(bm, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      pc1[[s]] <- data.table(set = s, n_windows = sum(!is.na(v)), recall = mean(v, na.rm = TRUE), ci_lo = q[1], ci_hi = q[2])
    }
    pc1 <- rbindlist(pc1)
    pc1_anc <- merge(wp, P1[, .(row, ancestry)], by = "row")[, .(n_windows = .N, recall = mean(recall)), by = ancestry]
    fwrite(pc1, file.path(OUT, "pc1_recall.tsv"), sep = "\t"); fwrite(pc1_anc, file.path(OUT, "pc1_recall_by_ancestry.tsv"), sep = "\t")
    pstate <- cP[, .N, by = state]
    D("PC-1: EUR window-mean recall ", f3(pc1[set == "PC1_EUR", recall]), " (95% CI ", f3(pc1[set == "PC1_EUR", ci_lo]), " to ",
      f3(pc1[set == "PC1_EUR", ci_hi]), ", ", pc1[set == "PC1_EUR", n_windows], " windows); pooled ", f3(pc1[set == "PC1_pooled", recall]),
      " (", pc1[set == "PC1_pooled", n_windows], " windows); states ", paste(pstate[, paste0(state, "=", N)], collapse = " "))
    if (pc1[set == "PC1_EUR", recall] < 0.80) add_stop("S5", paste0("PC-1 EUR recall ", f3(pc1[set == "PC1_EUR", recall]), " < 0.80"))
  } else {
    add_stop("S5", "PC-1 did not run; the curves cannot be reported as calibration of production")
  }
  # ---- PC-4, NC-1, NC-2 ----
  pc4 <- c(E15 = getr("EUR", "E_15_10"), G15 = getr("EUR", "G_10_15"))
  D("PC-4: EUR recall at (15,10) ", f3(pc4[["E15"]]), ", at (10,15) ", f3(pc4[["G15"]]), " (predicted >= 0.95)")
  f0 <- merge(fitsA[side == "eqtl" & lambda == 0 & is.na(error)], W[, .(row, stratum)], by = "row")
  f0[, cs := converged & n_cs >= 1]
  nc1_cs <- f0[stratum %in% c("EUR_Seval", "EUR_Snores"), .(f = mean(cs)), by = row][, mean(f)]
  nc1_pp4 <- getr("EUR", "E_0_10")
  D("NC-1: EUR eQTL credible set formed at lambda_e = 0 in ", f3(nc1_cs), " (window mean; predicted <= 0.05); PP.H4 > 0.5 in ", f3(nc1_pp4), " (predicted <= 0.02)")
  nc2 <- c(D6 = getr("EUR", "D_6_10"), D10 = getr("EUR", "D_10_10"), D20 = getr("EUR", "D_20_20"))
  D("NC-2: EUR false-shared rate (6,10) ", f3(nc2[["D6"]]), ", (10,10) ", f3(nc2[["D10"]]), ", (20,20) ", f3(nc2[["D20"]]), " (predicted <= 0.05 each)")

  # ---- Readout ----
  readout <- NULL
  readout_block <- any(vapply(stop_rules, function(x) x$rule %in% c("S1", "S2", "S3", "S5", "S6"), logical(1))) ||
    !g2[stratum == "EUR_Snores", testable] || jobs[["readout"]] != "done"
  sn <- merge(fitsA[side == "eqtl" & is.na(error)], W[, .(row, stratum, nE)], by = "row")[stratum == "EUR_Snores"]
  sn[, cs := converged & n_cs >= 1]
  build_map <- function(x) {
    L0 <- quantile(x[lambda == 0, lead_absz], 0.95, names = FALSE, type = 7)
    pcs <- x[, .(P_CS = mean(cs), med_lead = median(lead_absz), n = .N), by = lambda][order(lambda)]
    j <- which(pcs$P_CS >= 0.5)[1]
    if (is.na(j) || j == 1L) return(list(L0 = L0, pcs = pcs, lam50 = NA_real_, Lz50 = NA_real_))
    t <- (0.5 - pcs$P_CS[j - 1]) / (pcs$P_CS[j] - pcs$P_CS[j - 1])
    lam50 <- pcs$lambda[j - 1] + t * (pcs$lambda[j] - pcs$lambda[j - 1])
    Lz50 <- pcs$med_lead[j - 1] + t * (pcs$med_lead[j] - pcs$med_lead[j - 1])
    list(L0 = L0, pcs = pcs, lam50 = lam50, Lz50 = Lz50)
  }
  map <- build_map(sn)
  pc5 <- NA_real_
  fwrite(map$pcs, file.path(OUT, "readout_P_CS_by_lambda.tsv"), sep = "\t")
  D("Readout map (EUR S-nores, ", uniqueN(sn$row), " windows): L0 = ", f3(map$L0, 2), "; lambda_CS50 = ", f3(map$lam50, 2), "; Lz50 = ", f3(map$Lz50, 2),
    "; P_CS by lambda ", paste(map$pcs[, paste0(lambda, ":", f3(P_CS))], collapse = " "))
  if (!readout_block && is.finite(map$Lz50)) {
    LZ <- fread(file.path(A$readout, "lead_z.tsv"))
    cs0 <- LZ[susie_state_t12 == "no_susie_result" & fit == "cs0"]
    csg <- LZ[fit == "cs_ge1"]
    pc5 <- mean(csg$lead_absz > map$Lz50, na.rm = TRUE)
    D("PC-5: ", nrow(csg), " genes with an eQTL credible set; fraction with lead |z| > Lz50 = ", f3(pc5), " (predicted >= 0.85)")
    if (pc5 < 0.85) add_stop("S7", paste0("PC-5 ", f3(pc5), " < 0.85: readout categories not reported"))
  }
  if (!readout_block && is.finite(map$Lz50) && isTRUE(pc5 >= 0.85)) {
    cat_of <- function(z, L0, Lz50) fifelse(z <= L0, "A", fifelse(z <= Lz50, "B", "C"))
    cs0[, cat := cat_of(lead_absz, map$L0, map$Lz50)]
    cs0[, cat_nodup := cat_of(lead_absz_nodup, map$L0, map$Lz50)]
    cuts <- quantile(unique(sn[, .(row, nE)])$nE, c(1 / 3, 2 / 3), names = FALSE, type = 7)
    sn[, terc := findInterval(nE, cuts, left.open = TRUE) + 1L]
    tmaps <- lapply(1:3, function(k) build_map(sn[terc == k]))
    cs0[, terc := findInterval(nE, cuts, left.open = TRUE) + 1L]
    cs0[, cat_tercile := cat_of(lead_absz, vapply(tmaps, `[[`, numeric(1), "L0")[terc], vapply(tmaps, `[[`, numeric(1), "Lz50")[terc])]
    readout <- rbindlist(lapply(c("cat", "cat_nodup", "cat_tercile"), function(k)
      cs0[, .(n = .N), by = c(k)][, .(version = k, category = get(k), n, frac = n / nrow(cs0))]))
    readout <- readout[order(version, category)]
    readout_meta <- data.table(L0 = map$L0, lambda_CS50 = map$lam50, Lz50 = map$Lz50, n_cs0 = nrow(cs0), n_cs_ge1 = nrow(csg), pc5 = pc5,
                               n_cs0_with_dup_pos = sum(cs0$n_dup_pos > 0, na.rm = TRUE),
                               tercile_cuts = paste(round(cuts), collapse = ","),
                               tercile_L0 = paste(f3(vapply(tmaps, `[[`, numeric(1), "L0"), 2), collapse = ","),
                               tercile_Lz50 = paste(f3(vapply(tmaps, `[[`, numeric(1), "Lz50"), 2), collapse = ","))
    fwrite(readout, file.path(OUT, "readout_categories.tsv"), sep = "\t"); fwrite(readout_meta, file.path(OUT, "readout_map.tsv"), sep = "\t")
    fwrite(cs0[, .(ensembl, chr, nE, lead_absz, lead_absz_nodup, n_dup_pos, cat, cat_nodup, cat_tercile)], file.path(OUT, "readout_genes_cs0.tsv"), sep = "\t")
    for (i in seq_len(nrow(readout))) D("Readout ", readout$version[i], " ", readout$category[i], ": ", readout$n[i], " (", f3(readout$frac[i]), ")")
  } else {
    D("Readout categories not reported: ", if (readout_block) "blocked by a stop rule, an untestable EUR S-nores stratum, or a missing readout job"
      else if (!is.finite(map$Lz50)) "lambda_CS50 not inside the grid" else "PC-5 below 0.85 (S7)")
  }

  # ---- figures (PDF, 6 pt, no bold, black text, control gray for lambda = 0) ----
  pdf.options(useDingbats = FALSE)
  th <- theme_classic(base_size = 6) + theme(text = element_text(size = 6, colour = "black"), axis.text = element_text(size = 6, colour = "black"),
                                             strip.text = element_text(size = 6), strip.background = element_blank(),
                                             legend.text = element_text(size = 6), legend.title = element_text(size = 6))
  curve_dt <- function(cells_vec, grid, xname) {
    rbindlist(lapply(names(boot_store), function(s) {
      tb <- tabA[stratum == s]
      data.table(stratum = s, x = grid, recall = tb$rate_pp4_gt05[match(cells_vec, tb$cell)],
                 lo = tb$ci_lo[match(cells_vec, tb$cell)], hi = tb$ci_hi[match(cells_vec, tb$cell)])
    }))
  }
  ce <- curve_dt(e_cells, E_GRID); ce0 <- tabA[cell == "E_0_10", .(stratum, x = 0, recall = rate_pp4_gt05, lo = ci_lo, hi = ci_hi)]
  ce[, stratum := factor(stratum, STRATA)]; ce0[, stratum := factor(stratum, STRATA)]
  p <- ggplot(ce, aes(x, recall)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#C9265E", alpha = 0.2) +
    geom_line(colour = "#C9265E", linewidth = 0.4) + geom_point(colour = "#C9265E", size = 0.8) +
    geom_pointrange(data = ce0, aes(ymin = lo, ymax = hi), colour = "#9E9E9E", size = 0.15) +
    geom_hline(yintercept = c(0.5, 0.8), linetype = "dashed", linewidth = 0.2, colour = "#9E9E9E") +
    facet_wrap(~stratum, nrow = 2) + scale_x_continuous(breaks = c(0, E_GRID)) + coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Simulated causal eQTL z (lambda_e); GWAS z fixed at 10", y = "Fraction called colocalized (PP.H4 > 0.5)") + th
  ggsave(file.path(OUT, "figures", "h80_esweep_recall.pdf"), p, width = 5, height = 2.6, device = grDevices::pdf, useDingbats = FALSE)
  cg <- curve_dt(g_cells, G_GRID); cg[, stratum := factor(stratum, STRATA)]
  p <- ggplot(cg, aes(x, recall)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#1565C0", alpha = 0.2) +
    geom_line(colour = "#1565C0", linewidth = 0.4) + geom_point(colour = "#1565C0", size = 0.8) +
    geom_hline(yintercept = c(0.5, 0.8), linetype = "dashed", linewidth = 0.2, colour = "#9E9E9E") +
    facet_wrap(~stratum, nrow = 2) + scale_x_continuous(breaks = G_GRID) + coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Simulated causal GWAS z (lambda_g); eQTL z fixed at 10", y = "Fraction called colocalized (PP.H4 > 0.5)") + th
  ggsave(file.path(OUT, "figures", "h80_gsweep_recall.pdf"), p, width = 5, height = 2.6, device = grDevices::pdf, useDingbats = FALSE)
  comp <- cA[stratum %in% c("EUR_Seval", "EUR_Snores"), .N, by = .(cell, state)]
  comp[, frac := N / sum(N), by = cell]
  comp[, cell := factor(cell, cell_ids$cell)]
  p <- ggplot(comp, aes(state, cell, fill = frac)) + geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = formatC(frac, digits = 2, format = "f")), size = 6 / .pt, colour = "black") +
    scale_fill_gradient(low = "#F2F2F2", high = "#9E9E9E", limits = c(0, 1), name = "Fraction") +
    labs(x = "Gene-pair state (production rules)", y = "Cell (set_lambda_e_lambda_g), EUR windows") + th +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(OUT, "figures", "h80_state_composition_EUR.pdf"), p, width = 4.5, height = 3.4, device = grDevices::pdf, useDingbats = FALSE)
  if (!is.null(readout)) {
    hz <- rbind(cs0[, .(set = "No eQTL credible set (13,704 frame genes)", z = lead_absz)], csg[, .(set = "eQTL credible set (4,979 genes)", z = lead_absz)])
    p <- ggplot(hz, aes(z)) + geom_histogram(bins = 80, fill = "#9E9E9E", colour = NA) +
      geom_vline(xintercept = map$L0, linetype = "dashed", linewidth = 0.3) + geom_vline(xintercept = map$Lz50, linetype = "solid", linewidth = 0.3) +
      facet_wrap(~set, ncol = 1, scales = "free_y") + labs(x = "Observed lead eQTL |z| over the fit's SNPs", y = "Genes") + th
    ggsave(file.path(OUT, "figures", "h80_readout_lead_z.pdf"), p, width = 3.5, height = 2.6, device = grDevices::pdf, useDingbats = FALSE)
  }
  if (!is.null(armB)) {
    bd <- armB[quantity == "diff"]
    p <- ggplot(bd, aes(estimate, cell)) + geom_vline(xintercept = 0, colour = "#9E9E9E", linewidth = 0.3) +
      geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), size = 0.15) +
      labs(x = "Mismatched minus matched LD, rate of PP.H4 > 0.5 (120 EUR windows)", y = "Cell") + th
    ggsave(file.path(OUT, "figures", "h80_armB_difference.pdf"), p, width = 3.5, height = 1.8, device = grDevices::pdf, useDingbats = FALSE)
  }
  for (f in list.files(file.path(OUT, "figures"), pattern = "pdf$", full.names = TRUE)) {
    txt <- system2("pdftotext", c("-layout", shQuote(f), "-"), stdout = TRUE)
    writeLines(txt, sub("pdf$", "pdftotext.txt", f))
  }
}

# ======================== verdict, computed from the rules ========================
SR <- if (length(stop_rules)) rbindlist(stop_rules) else data.table(rule = character(0), text = character(0))
fwrite(SR, file.path(OUT, "stop_rules_triggered.tsv"), sep = "\t")
writeLines(diag_lines, file.path(OUT, "diagnostics.txt"))
fwrite(data.table(order = seq_along(diag_lines), line = diag_lines), file.path(OUT, "diagnostics.tsv"), sep = "\t")
hard <- SR[rule %in% c("S1", "S2", "S3", "S6", "G2", "setup")]
lam_line <- function(s, sw, thr) { x <- lamA[stratum == s & sweep == sw & threshold == thr]; if (!nrow(x)) "NA" else paste0(x$printed, " (95% CI ", x$ci_lo, " to ", x$ci_hi, ")") }
pred <- NULL
if (nrow(hard)) {
  verdict <- paste0("STOPPED by ", paste(unique(hard$rule), collapse = ", "), ": ", paste(hard$text, collapse = "; "), ". No simulated number is used.")
} else if ("S5" %in% SR$rule) {
  verdict <- paste0("FAILED CONTROL (S5): ", SR[rule == "S5", text][1], ". The simulated curves are not reported as calibration of production.")
} else {
  r3 <- getr("EUR", "E_3_10"); r15 <- getr("EUR", "E_15_10")
  informative <- isTRUE(r3 <= 0.20) && isTRUE(r15 >= 0.80)
  verdict <- paste0("CALIBRATION PASSES (PC-1 to PC-3 and NC-3 pass; S1-S6 not triggered). EUR lambda50_e = ", lam_line("EUR", "E", 0.5),
                    "; lambda80_e = ", lam_line("EUR", "E", 0.8), "; recall at lambda_e 3 = ", f3(r3), ", at 15 = ", f3(r15), ". ",
                    if (informative) "The E-sweep moves from <= 0.20 to >= 0.80 inside lambda_e 3 to 15 (informative calibration)."
                    else if (isTRUE(r15 < 0.80)) "Recall stays below 0.80 at lambda_e = 15: the pipeline misses strong shared signals (pipeline limit)."
                    else "Recall is above 0.20 at lambda_e = 3: the pipeline calls very weak shared signals.",
                    if ("S7" %in% SR$rule) " Readout categories not reported (S7)." else "",
                    if ("S4" %in% SR$rule) paste0(" Untestable strata (S4): ", paste(SR[rule == "S4", text], collapse = "; "), ".") else "")
  fs_max <- max(nc2, na.rm = TRUE)
  armB_drop <- if (!is.null(armB)) armB[quantity == "diff" & cell %in% c("E_8_10", "E_15_10"), min(estimate)] else NA_real_
  if (isTRUE(fs_max > 0.05)) verdict <- paste0(verdict, " False-shared rate above 0.05 in an EUR distinct cell (",
                                               paste(paste0(c("(6,10)", "(10,10)", "(20,20)"), " ", f3(nc2)), collapse = ", "), "); quoted beside the positive count, memberships unchanged.")
  if (isTRUE(armB_drop < -0.10)) verdict <- paste0(verdict, " Arm B recall drop larger than 0.10 (", f3(armB_drop),
                                                   "): pipeline limit under LD-reference mismatch; positive SuSiE calls are conservative by this amount.")
  pred <- data.table(prediction = c(
    "EUR recall <= 0.20 at lambda_e = 3", "EUR recall >= 0.80 at lambda_e = 15", "EUR lambda50_e in [4, 6]", "EUR lambda80_e in [5, 8]",
    "EUR lambda50_g <= lambda50_e", "EUR false-distinct (PP.H3 > 0.8) <= 0.10 at lambda_e 8, 10, 15",
    "Non-EUR recall at lambda_e = 8 lower than EUR by 0 to 0.15 (each stratum)",
    "Arm B recall change in [-0.10, 0] at lambda_e 8 and 15", "Arm B false-shared (mismatched) <= 0.05",
    "Readout: category A >= 60% and C <= 10%", "PC-1 recall >= 0.90 (EUR and pooled)", "PC-4 recall >= 0.95 at (15,10) and (10,15)",
    "NC-1 eQTL CS <= 0.05 and PP.H4 > 0.5 <= 0.02", "NC-2 false-shared <= 0.05 in every EUR distinct cell"))
  le <- lamA[stratum == "EUR" & sweep == "E"]; lg <- lamA[stratum == "EUR" & sweep == "G"]
  r8 <- getr("EUR", "E_8_10")
  nonEUR <- vapply(c("EAS", "AFR", "AMR", "SAS"), function(s) r8 - getr(s, "E_8_10"), numeric(1))
  obs <- c(f3(r3), f3(r15), le[threshold == 0.5, printed], le[threshold == 0.8, printed],
           paste0(lg[threshold == 0.5, printed], " vs ", le[threshold == 0.5, printed]),
           paste(f3(c(getr("EUR", "E_8_10", "rate_pp3_gt08"), getr("EUR", "E_10_10", "rate_pp3_gt08"), getr("EUR", "E_15_10", "rate_pp3_gt08"))), collapse = ", "),
           paste(paste0(names(nonEUR), " ", f3(nonEUR)), collapse = "; "),
           if (!is.null(armB)) paste(f3(armB[quantity == "diff" & cell %in% c("E_8_10", "E_15_10"), estimate]), collapse = ", ") else "not run",
           if (!is.null(armB)) f3(armB[quantity == "mismatched" & cell == "D_10_10", estimate]) else "not run",
           if (!is.null(readout)) paste0("A ", f3(sum(readout[version == "cat" & category == "A", frac])), ", C ", f3(sum(readout[version == "cat" & category == "C", frac]))) else "not reported",
           if (!is.null(pc1)) paste0("EUR ", f3(pc1[set == "PC1_EUR", recall]), ", pooled ", f3(pc1[set == "PC1_pooled", recall])) else "not run",
           paste0(f3(pc4[["E15"]]), ", ", f3(pc4[["G15"]])), paste0(f3(nc1_cs), ", ", f3(nc1_pp4)), paste(f3(nc2), collapse = ", "))
  ev <- function(x) if (length(x) == 0 || is.na(x)) "not evaluable" else if (x) "met" else "not met"
  # a <= b for interpolated lambdas that may be censored at the grid edges ("low" = below 3, "high" = above 15)
  cmp_le <- function(a, b) {
    if (a$censored == "none" && b$censored == "none") return(a$estimate <= b$estimate)
    if (a$censored == "low" && b$censored != "low") return(TRUE)
    if (b$censored == "high" && a$censored != "high") return(TRUE)
    if (a$censored == "high" && b$censored != "high") return(FALSE)
    if (b$censored == "low" && a$censored != "low") return(FALSE)
    NA
  }
  met <- c(ev(r3 <= 0.20), ev(r15 >= 0.80),
           ev(le[threshold == 0.5, censored == "none" & estimate >= 4 & estimate <= 6]),
           ev(le[threshold == 0.8, censored == "none" & estimate >= 5 & estimate <= 8]),
           ev(cmp_le(lg[threshold == 0.5], le[threshold == 0.5])),
           ev(all(c(getr("EUR", "E_8_10", "rate_pp3_gt08"), getr("EUR", "E_10_10", "rate_pp3_gt08"), getr("EUR", "E_15_10", "rate_pp3_gt08")) <= 0.10)),
           ev(all(nonEUR >= 0 & nonEUR <= 0.15, na.rm = TRUE)),
           if (!is.null(armB)) ev(all(armB[quantity == "diff" & cell %in% c("E_8_10", "E_15_10"), estimate >= -0.10 & estimate <= 0])) else "not run",
           if (!is.null(armB)) ev(armB[quantity == "mismatched" & cell == "D_10_10", estimate] <= 0.05) else "not run",
           if (!is.null(readout)) ev(sum(readout[version == "cat" & category == "A", frac]) >= 0.60 && sum(readout[version == "cat" & category == "C", frac]) <= 0.10) else "not reported",
           if (!is.null(pc1)) ev(all(pc1$recall >= 0.90)) else "not run",
           ev(all(pc4 >= 0.95)), ev(nc1_cs <= 0.05 && nc1_pp4 <= 0.02), ev(all(nc2 <= 0.05)))
  pred[, observed := obs][, result := met]
  fwrite(pred, file.path(OUT, "predictions_check.tsv"), sep = "\t")
}
log_msg("VERDICT: ", verdict)

# ---- proposed text (returned for approval; never written to docs/) ----
prop <- character(0)
if (!nrow(hard) && !("S5" %in% SR$rule)) {
  fig2 <- readLines(file.path(ROOT, "docs/manuscript/draft/fig2.md"))
  old <- paste(fig2[89:91], collapse = "\n")
  x <- lamA[stratum == "EUR" & sweep == "E" & threshold == 0.5]
  pve <- if (x$censored == "none") paste0(", equivalent to about ", formatC(100 * x$estimate^2 / EQTL_N, digits = 1, format = "f"), "% of expression variance at n = 1,183") else ""
  new <- paste0("Second, the absence of genetic support does not render a disease-state program\nunimportant, because colocalization interrogates inherited regulation rather than\n",
                "tissue activity and its sensitivity is limited: in simulations at 200 real tier-1/2\nwindows run with the pipeline's own settings, a simulated single shared causal variant\n",
                "was called colocalized in half of simulations (averaged over windows) only when its\neQTL z reached ", x$printed, " (95% CI ", x$ci_lo, " to ", x$ci_hi, ")", pve,
                " (Supplementary Note H80). We consequently")
  prop <- c(prop, "### docs/manuscript/draft/fig2.md lines 89-91", "", "old_text:", "", "```", old, "```", "", "new_text:", "", "```", new, "```", "")
  meth <- readLines(file.path(ROOT, "docs/manuscript/METHODS.md"))
  old_m <- meth[566]
  add_m <- paste0(" Sensitivity of this pipeline was calibrated by simulation (Supplementary Note H80): at 200 sampled tier-1/2 gene by GWAS-ancestry windows, ",
                  "z-scores for one causal variant shared by the eQTL and GWAS sides (or two distinct variants with r2 < 0.1) were drawn from the SuSiE-RSS model ",
                  "z ~ N(lambda R[, c], R + 10^-3 I) using each window's own LD matrices, 10 replicates per strength, and passed through the production eQTL and GWAS SuSiE fits ",
                  "and coloc.susie with the production settings; recall is the window-averaged fraction of replicates with SuSiE PP.H4 > 0.5, with 95% window-bootstrap intervals.")
  prop <- c(prop, "### docs/manuscript/METHODS.md line 566 (one sentence appended at the end of the paragraph)", "", "old_text:", "", "```", old_m, "```", "",
            "new_text:", "", "```", paste0(old_m, add_m), "```", "")
  if (!is.null(readout) && !("S7" %in% SR$rule)) {
    n <- function(k) format(sum(readout[version == "cat" & category == k, n]), big.mark = ",")
    prop <- c(prop, "### docs/manuscript/05_figure_legends.md, Figure 2B legend (optional line; user decision)", "", "new_text (added sentence):", "", "```",
              paste0("Of 13,704 genes without an eQTL credible set, ", n("A"), " were indistinguishable from no eQTL signal at their window size, ", n("B"),
                     " carried a signal below the pipeline's detection level, and ", n("C"), " exceeded it."), "```", "")
  }
} else {
  prop <- c(prop, "No number from this run enters the text; `docs/manuscript/draft/fig2.md:89-91` stays qualitative.")
}

# ---- RESULT_H80.md ----
L <- c("# RESULT H80: production SuSiE-COLOC calibration at real tier-1/2 LD windows", "",
       paste0("Written by `scripts/04_assemble.R` at ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), ". Prespec `PRESPEC_H80.md` sha256 ", H80_PRESPEC_SHA, "."), "",
       "## Verdict (computed by the prespecified rules)", "", verdict, "",
       "## Stop rules", "", if (nrow(SR)) paste0("- ", SR$rule, ": ", SR$text) else "- none triggered", "",
       "## Diagnostics (prespecified order)", "", paste0("- ", diag_lines), "")
if (!blocked && exists("lamA")) {
  L <- c(L, "## Primary statistics", "",
         "Unit: LD window (gene x GWAS ancestry). Replicates inside a window do not add n. Intervals are 95% window-bootstrap percentile intervals (2,000 resamples, stream 3000), one per cell and stratum, not adjusted, and none is read as a test. No p values. No testing family.", "",
         "### lambda50 and lambda80", "", "| Stratum | Sweep | Threshold | Estimate | 95% CI | Resamples censored low / high |", "|---|---|---|---|---|---|",
         lamA[, sprintf("| %s | %s | %s | %s | %s to %s | %d / %d |", stratum, sweep, threshold, printed, ci_lo, ci_hi, boot_censored_low, boot_censored_high)], "",
         "### Arm A cells, EUR pooled", "", "| Cell | Windows | PP.H4 > 0.5 | 95% CI | PP.H3 > 0.8 | 95% CI |", "|---|---|---|---|---|---|",
         tabA[stratum == "EUR", sprintf("| %s | %d | %s | %s to %s | %s | %s to %s |", cell, n_windows, f3(rate_pp4_gt05), f3(ci_lo), f3(ci_hi), f3(rate_pp3_gt08), f3(h3_ci_lo), f3(h3_ci_hi))], "",
         "All strata and the state composition per cell: `run_assemble/armA_cells_by_stratum.tsv`.", "")
  if (!is.null(armB)) L <- c(L, "### Arm B (1000 Genomes EUR generator, PolyFun fit) against the matched twin", "", "| Cell | Quantity | Estimate | 95% CI | Windows |", "|---|---|---|---|---|",
                             armB[, sprintf("| %s | %s | %s | %s to %s | %d |", cell, quantity, f3(estimate), f3(ci_lo), f3(ci_hi), n_windows)], "")
  if (!is.null(pred)) L <- c(L, "### Predictions written before the run", "", "| Prediction | Observed | Result |", "|---|---|---|",
                             pred[, sprintf("| %s | %s | %s |", prediction, observed, result)], "")
}
L <- c(L, "## Proposed text changes (for user approval; nothing was written to docs/)", "", prop,
       "## Files", "", paste0("- setup: `", SETUP, "`"), paste0("- Arm A: `", A$arm_a, "`"), paste0("- Arm B and PC-1: `", A$arm_b, "`"),
       paste0("- readout: `", A$readout, "`"), paste0("- assembly, tables, figures: `", OUT, "`"), "- deviations: `DEVIATIONS.md`", "")
writeLines(L, RESULT)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_R.txt"))
writeLines("DONE", file.path(OUT, "DONE"))
log_msg("RESULT written: ", RESULT)
