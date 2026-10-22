#!/usr/bin/env Rscript
# H81 summary: diagnostics first (prespec "Prespecified diagnostics" 1-12), stop rules S5-S8, then outputs O1-O5
# with replicate-cluster bootstrap (2,000 resamples of the replicates, seed 20260923 + 1000) and Clopper-Pearson
# intervals, thresholds T-O1 to T-O5, controls PC-A to PC-D and NC-1 to NC-4, the verdict by the prespec's
# Reading rules, three supplementary panels, and RESULT_H81.md.
# Unit: a spiked gene inside a replicate is a calibration draw, not a donor; replicates are permutation-plus-
# thinning draws on one fixed set of 774 units. No p value is reported; there is no test family.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
options(warn = 1, datatable.showProgress = FALSE)

ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
H81 <- file.path(ROOT, "Analysis/Hypotheses/H81-20261015T125233Z")
PRESPEC_SHA <- "ae52c1949dda4f549c6f35c8ce7170815a5a1534e7beaed016aaedc030f46875"
DEG <- file.path(ROOT, "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(name, default = NULL) { i <- match(paste0("--", name), args); if (is.na(i)) default else args[i + 1L] }
out <- getarg("out"); bseed <- as.integer(getarg("seed", "20261923")); B <- as.integer(getarg("B", "2000"))
debug <- "--debug" %in% args
gdir <- getarg("guard-dir", file.path(out, "guard"))
stopifnot(!is.null(out), dir.exists(out))
sdir <- file.path(out, getarg("summary-name", "summary"))
if (dir.exists(sdir)) stop("refusing to overwrite ", sdir)
dir.create(sdir); dir.create(file.path(sdir, "figures"))
logf <- file.path(sdir, "log_summary.txt")
say <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = logf, append = TRUE) }
sha256_file <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE)[1], "[[:space:]]+")[[1]][1]
h <- sha256_file(file.path(H81, "PRESPEC_H81.md"))
if (!identical(h, PRESPEC_SHA)) stop("PRESPEC_H81.md sha256 mismatch: ", h)
say("H81 summary. prespec sha256 OK. out=%s debug=%s B=%d bootstrap seed=%d", out, debug, B, bseed)
f3 <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 3))
f4 <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 4))
pct <- function(x) ifelse(is.na(x), "NA", formatC(100 * x, format = "f", digits = 1))

DIAG <- list()
dg <- function(no, item, value) DIAG[[length(DIAG) + 1L]] <<- data.table(diagnostic = no, item = item, value = as.character(value))
COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")
TRUE_LABELS <- "GSE126848:23/7;GSE130970:6/65;GSE135251:10/197;GSE162694:31/109;GSE213621:68/258"

# ---------------------------------------------------------------- guard (diagnostics 1, 2)
gst <- fread(file.path(gdir, "guard_status.tsv")); gck <- fread(file.path(gdir, "guard_checks.tsv"))
g1h <- fread(file.path(gdir, "g1_input_hashes.tsv"))
trig <- function(r) isTRUE(gst[rule == r, triggered])
for (i in seq_len(nrow(g1h))) dg(1, g1h$id[i], sprintf("%s match=%s", g1h$observed_sha256[i], g1h$match[i]))
for (i in seq_len(nrow(gck))) dg(2, paste(gck$gate[i], gck$check[i]), sprintf("value=%s target=%s pass=%s", gck$value[i], gck$target[i], gck$pass[i]))
if (!"gate_pass" %in% names(gst)) gst[, gate_pass := c(FALSE, NA, NA, NA)]
gate_pass <- setNames(gst$gate_pass, c("G1", "G-A", "G-B", "G-C"))
S <- list(S1 = trig("S1"), S2 = trig("S2"), S3 = trig("S3"), S4 = trig("S4"), S5 = NA, S6 = NA, S7 = NA)
say("guard: S1 %s S2 %s S3 %s S4 %s; gates G1 %s G-A %s G-B %s G-C %s", S$S1, S$S2, S$S3, S$S4,
    gate_pass["G1"], gate_pass["G-A"], gate_pass["G-B"], gate_pass["G-C"])

stopped <- NULL
if (S$S1) stopped <- "S1" else if (S$S2) stopped <- "S2"

# ---------------------------------------------------------------- replicates
rdir <- file.path(out, "replicates")
done <- sort(as.integer(sub("rep_(\\d+)\\.done", "\\1", list.files(rdir, pattern = "^rep_\\d+\\.done$"))))
errs <- list.files(rdir, pattern = "\\.error$")
REPS <- done; R <- length(REPS)
dg(11, "completed replicates", R)
for (e in errs) dg(11, paste("failure", e), paste(readLines(file.path(rdir, e)), collapse = " "))
if (is.null(stopped)) {
  S$S7 <- R < 45
  if (S$S7 && !debug) stopped <- "S7"
}
D <- if (R > 0) rbindlist(lapply(REPS, function(r) fread(file.path(rdir, sprintf("rep_%03d_genes.tsv.gz", r))))) else data.table()
RS <- if (R > 0) rbindlist(lapply(REPS, function(r) fread(file.path(rdir, sprintf("rep_%03d_summary.tsv", r))))) else data.table()
say("replicates read: %d (errors %d); gene rows %d", R, length(errs), nrow(D))

if (R > 0) {
  # diagnostic 3: permuted label table and label-change fraction
  for (i in seq_len(nrow(RS))) dg(3, sprintf("rep %d labels / frac changed", RS$rep[i]), paste(RS$labels[i], "|", RS$frac_label_changed[i]))
  lab_ok <- if (debug) uniqueN(RS$labels) == 1L else all(RS$labels == TRUE_LABELS)
  dg(3, "all permuted label tables equal the true margins", lab_ok)
  fc <- rbindlist(lapply(strsplit(RS$frac_label_changed, ";"), function(x) {
    y <- tstrsplit(x, ":"); data.table(cohort = y[[1]], frac = as.numeric(y[[2]])) }))
  for (co in unique(fc$cohort)) dg(3, paste("median fraction of labels changed", co), f3(median(fc[cohort == co, frac])))
  # diagnostic 4: thinning check
  D[, ra := fifelse(is.finite(realized_aligned), realized_aligned, NA_real_)]
  th <- D[, .(n = .N, median_aligned = median(ra, na.rm = TRUE), frac_positive = mean(ra > 0, na.rm = TRUE),
              n_nonfinite = sum(is.na(ra))), by = .(set, decile)][order(set, decile)]
  fwrite(th, file.path(sdir, "diag4_thinning_by_set_decile.tsv"), sep = "\t")
  for (i in seq_len(nrow(th))) dg(4, sprintf("%s decile %d median aligned change (n, frac>0)", th$set[i], th$decile[i]),
                                  sprintf("%s (%d, %s)", f4(th$median_aligned[i]), th$n[i], f3(th$frac_positive[i])))
  s6_med <- D[set == "U1.0" & decile >= 5, median(ra, na.rm = TRUE)]
  s6_dec <- th[set == "U1.0" & decile >= 5, median_aligned]
  dg(4, "S6: median aligned change U1.0 deciles 5-10 pooled (per-decile range)", sprintf("%s (%s to %s)", f4(s6_med), f4(min(s6_dec)), f4(max(s6_dec))))
  if (is.null(stopped)) { S$S6 <- !(s6_med >= 0.9 && s6_med <= 1.1); if (S$S6 && !debug) stopped <- "S6" }
  nc4_pos <- D[set == "U1.5" & decile >= 5, mean(ra > 0, na.rm = TRUE)]
  nc4_u0 <- D[set == "U0", median(ra, na.rm = TRUE)]
  dg(4, "NC-4: fraction of U1.5 deciles 5-10 with positive aligned change", f4(nc4_pos))
  dg(4, "NC-4: median aligned change U0", f4(nc4_u0))
  # diagnostic 5: universe and filtered-out spiked genes
  dg(5, "pooled universe per replicate (min/median/max)", sprintf("%d/%d/%d", min(RS$universe), as.integer(median(RS$universe)), max(RS$universe)))
  un <- D[, .(n = .N, frac_untested = mean(!tested)), by = .(set, decile)][order(set, decile)]
  fwrite(un, file.path(sdir, "diag5_untested_by_set_decile.tsv"), sep = "\t")
  for (i in seq_len(nrow(un))) dg(5, sprintf("%s decile %d fraction removed by filterByExpr", un$set[i], un$decile[i]), f4(un$frac_untested[i]))
  # diagnostic 6: attenuation
  at <- D[tested == TRUE & delta > 0, .(median_ratio = median(logFC * dir / delta)), by = .(set, decile)][order(set, decile)]
  fwrite(at, file.path(sdir, "diag6_attenuation_by_set_decile.tsv"), sep = "\t")
  for (i in seq_len(nrow(at))) dg(6, sprintf("%s decile %d median pooled logFC / intended delta", at$set[i], at$decile[i]), f3(at$median_ratio[i]))
  # diagnostic 7: quality weights
  qa <- fread(file.path(gdir, "quality_weight_range_adopted.tsv"))
  dg(7, "adopted quality weights min/median/max", sprintf("%s/%s/%s", f3(qa$adopted[1]), f3(qa$adopted[2]), f3(qa$adopted[3])))
  dg(7, "replicate quality weight min (range over replicates)", sprintf("%s to %s", f3(min(RS$qw_min)), f3(max(RS$qw_min))))
  dg(7, "replicate quality weight max (range over replicates)", sprintf("%s to %s", f3(min(RS$qw_max)), f3(max(RS$qw_max))))
  # diagnostic 8: unspiked background
  med_bg <- median(RS$supported_unspiked)
  for (i in seq_len(nrow(RS))) dg(8, sprintf("rep %d supported / tested_negative among unspiked (of %d)", RS$rep[i], RS$unspiked_in_universe[i]),
                                  sprintf("%d / %d", RS$supported_unspiked[i], RS$tested_negative_unspiked[i]))
  dg(8, "median supported among unspiked (S5 stop > 105; NC-2 < 20)", med_bg)
  if (is.null(stopped)) { S$S5 <- med_bg > 105; if (S$S5 && !debug) stopped <- "S5" }
  # diagnostic 9: realized tau
  tau <- D[tested == TRUE & set %in% c("H1.0", "H0.5", "U1.0", "U0"), .(median_tau = median(sqrt(pmax(tau2, 0)))), by = set]
  for (i in seq_len(nrow(tau))) dg(9, paste("median realized DerSimonian-Laird tau,", tau$set[i]), f3(tau$median_tau[i]))
  # diagnostic 12: BH-effective alpha
  dg(12, "median BH-effective alpha of the equivalence rule (range)", sprintf("%s (%s to %s)", f4(median(RS$eff_alpha)), f4(min(RS$eff_alpha)), f4(max(RS$eff_alpha))))
  dg(12, "median ratio to 0.0177 (adopted per-cohort fits)", f3(median(RS$eff_alpha) / (0.05 * 8548 / 24115)))
  dg(12, "replicates with non-finite per-cohort estimates (genes, max)", max(RS$eq_nonfinite))
  # extra: wrong-sign supported calls among spiked genes with delta > 0
  ws <- D[supported == TRUE & delta > 0 & type == "U", mean(sign(logFC) != dir)]
  dg(6, "fraction of supported U-set spiked genes with the wrong sign", f4(ws))
}

# ---------------------------------------------------------------- bootstrap machinery
CNT <- NULL
if (R > 0) {
  set.seed(bseed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  CNT <- matrix(0, B, R)
  for (b in seq_len(B)) CNT[b, ] <- tabulate(sample.int(R, R, replace = TRUE), R)
}
cp_int <- function(x, n) c(if (x == 0) 0 else qbeta(0.025, x, n - x + 1), if (x == n) 1 else qbeta(0.975, x + 1, n - x))
# rate over the gene-replicates in `rows` (logical over D) of event `ev`, by cell
EST <- list()
est <- function(output, quantity, rows, ev, set_label, by = character(), grid = NULL, conditional = FALSE) {
  X <- D[rows, c("rep", by), with = FALSE]; X[, ev_ := ev[rows]]
  if (is.null(grid)) grid <- if (length(by)) unique(X[, ..by]) else data.table(dummy_ = 1)
  if (!length(by)) X[, dummy_ := 1]
  byv <- if (length(by)) by else "dummy_"
  A <- X[, .(num = sum(ev_), den = .N), by = c(byv, "rep")]
  o <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
    key <- grid[i]
    a <- A[key, on = byv, nomatch = NULL]
    num <- numeric(R); den <- numeric(R)
    num[match(a$rep, REPS)] <- a$num; den[match(a$rep, REPS)] <- a$den
    n <- sum(den); x <- sum(num)
    measured <- if (conditional) n >= 50 else n > 0
    if (n > 0) {
      rb <- as.vector(CNT %*% num) / as.vector(CNT %*% den)
      q <- quantile(rb[is.finite(rb)], c(0.025, 0.975), type = 7, names = FALSE)
      ci <- cp_int(x, n)
    } else { q <- c(NA, NA); ci <- c(NA, NA) }
    cbind(key, data.table(n_den = n, n_num = x, rate = if (n > 0) x / n else NA_real_,
                          boot_lo = q[1], boot_hi = q[2], cp_lo = ci[1], cp_hi = ci[2], measured = measured,
                          n_reps_with_den = sum(den > 0)))
  }), fill = TRUE)
  o[, `:=`(output = output, quantity = quantity, set = set_label, conditional = conditional)]
  if ("dummy_" %in% names(o)) o[, dummy_ := NULL]
  if (!"decile" %in% names(o)) o[, decile := NA_integer_]
  if (!"spiked_cohort" %in% names(o)) o[, spiked_cohort := NA_character_]
  EST[[length(EST) + 1L]] <<- o
  invisible(o)
}
DEC <- data.table(decile = 1:10)
COH <- data.table(spiked_cohort = COHORTS)
CXD <- CJ(spiked_cohort = COHORTS, decile = 1:10)

# ---------------------------------------------------------------- diagnostic 10 (conditioning counts), before outputs
if (R > 0) {
  for (s in c("S1.0", "S2.0")) {
    cc <- D[set == s & supported == TRUE, .N, by = spiked_cohort]
    for (co in COHORTS) dg(10, sprintf("O4 %s supported gene-replicates spiked in %s", s, co), cc[spiked_cohort == co, if (.N) N else 0L])
    cd <- D[set == s & supported == TRUE, .N, by = decile]
    for (d in 1:10) dg(10, sprintf("O4 %s supported gene-replicates decile %d", s, d), cd[decile == d, if (.N) N else 0L])
    cxd <- D[set == s & supported == TRUE, .N, by = .(spiked_cohort, decile)]
    dg(10, sprintf("O4 %s cohort-by-decile cells with >= 50 supported", s), sum(cxd$N >= 50))
  }
  for (s in c("U0.5", "U0.75", "U1.0", "U1.5", "H1.0", "H0.5")) {
    cd <- D[set == s & supported == TRUE, .N, by = decile]
    dg(10, sprintf("O5 %s supported gene-replicates overall", s), D[set == s & supported == TRUE, .N])
    for (d in 1:10) dg(10, sprintf("O5 %s supported gene-replicates decile %d", s, d), cd[decile == d, if (.N) N else 0L])
  }
}
DG <- rbindlist(DIAG)
fwrite(DG, file.path(sdir, "diagnostics.tsv"), sep = "\t")
say("---- DIAGNOSTICS (printed before any output) ----")
for (i in seq_len(nrow(DG))) if (!grepl("^rep \\d+", DG$item[i]) && DG$diagnostic[i] != 1)
  say("D%s | %s | %s", DG$diagnostic[i], DG$item[i], DG$value[i])
say("diagnostic 1: %d of %d input sha256 match", sum(g1h$match), nrow(g1h))
say("stop rules: %s", paste(names(S), unlist(lapply(S, function(x) if (is.na(x)) "NA" else as.character(x))), sep = "=", collapse = " "))

# ---------------------------------------------------------------- outputs O1-O5
OUTS <- NULL
if (is.null(stopped) && R > 0) {
  say("---- OUTPUTS ----")
  # O1 recall and O2 false supported
  for (s in c("U0.5", "U0.75", "U1.0", "U1.5", "H1.0", "H0.5")) {
    est("O1", "P(supported)", D$set == s, D$supported, s)
    est("O1", "P(supported)", D$set == s, D$supported, s, by = "decile", grid = DEC)
  }
  est("O2", "P(supported)", D$set == "U0", D$supported, "U0")
  est("O2", "P(supported)", D$set == "U0", D$supported, "U0", by = "decile", grid = DEC)
  # O3 false tested_negative (L09-01 rule), size, power; D2 margin 0.25
  if (!S$S4) for (s in c("U0.75", "U1.0", "U0.5", "H0.5", "H1.0", "U0", "U1.5")) {
    est("O3", "P(tested_negative)", D$set == s, D$tested_negative, s)
    est("O3", "P(tested_negative)", D$set == s, D$tested_negative, s, by = "decile", grid = DEC)
    est("D2", "P(tested_negative at margin 0.25)", D$set == s, D$tested_negative_m025, s)
  }
  # O4 one-cohort effects
  for (s in c("S1.0", "S2.0")) {
    rs <- D$set == s
    rsup <- rs & D$supported
    for (g in list(list(by = character(), grid = NULL), list(by = "spiked_cohort", grid = COH),
                   list(by = "decile", grid = DEC), list(by = c("spiked_cohort", "decile"), grid = CXD))) {
      est("O4", "P(supported)", rs, D$supported, s, by = g$by, grid = g$grid)
      if (!S$S4) est("O4", "P(tested_negative)", rs, D$tested_negative, s, by = g$by, grid = g$grid)
      if (!S$S3) {
        est("O4", "P(not SD-any | supported)", rsup, !D$sd_any %in% TRUE, s, by = g$by, grid = g$grid, conditional = TRUE)
        est("O4", "P(supported and not SD-any)", rs, D$supported & !D$sd_any %in% TRUE, s, by = g$by, grid = g$grid)
        est("O4", "P(SD-own | supported)", rsup, D$sd_own %in% TRUE, s, by = g$by, grid = g$grid, conditional = TRUE)
        est("D3", "P(SD-sig | supported)", rsup, D$sd_sig %in% TRUE, s, by = g$by, grid = g$grid, conditional = TRUE)
        est("D3", "P(SD-any | supported)", rsup, D$sd_any %in% TRUE, s, by = g$by, grid = g$grid, conditional = TRUE)
      }
    }
    # D4 detection in the spiked cohort's own fit
    D[, own_det := FALSE]
    iS <- which(D$type == "S")
    fdr_m <- as.matrix(D[iS, paste0("coh_FDR_", COHORTS), with = FALSE]); lfc_m <- as.matrix(D[iS, paste0("coh_logFC_", COHORTS), with = FALSE])
    ci <- cbind(seq_along(iS), match(D$spiked_cohort[iS], COHORTS))
    D[iS, own_det := !is.na(fdr_m[ci]) & fdr_m[ci] < 0.05 & abs(lfc_m[ci]) > 0.5]
    est("D4", "P(FDR < 0.05 and |logFC| > 0.5 in the spiked cohort's own fit)", rs, D$own_det, s, by = "spiked_cohort", grid = COH)
  }
  # O5 false source-dependence on uniform and heterogeneous effects
  if (!S$S3) for (s in c("U0.5", "U0.75", "U1.0", "U1.5", "H1.0", "H0.5")) {
    rsup <- D$set == s & D$supported
    est("O5", "P(SD-any | supported)", rsup, D$sd_any %in% TRUE, s, conditional = TRUE)
    est("O5", "P(SD-any | supported)", rsup, D$sd_any %in% TRUE, s, by = "decile", grid = DEC, conditional = TRUE)
    est("D3", "P(SD-sig | supported)", rsup, D$sd_sig %in% TRUE, s, conditional = TRUE)
  }
  OUTS <- rbindlist(EST, fill = TRUE)
  setcolorder(OUTS, c("output", "quantity", "set", "spiked_cohort", "decile", "conditional", "n_den", "n_num", "rate",
                      "boot_lo", "boot_hi", "cp_lo", "cp_hi", "measured", "n_reps_with_den"))
  OUTS[, rate_reported := fifelse(conditional & !measured, NA_real_, rate)]
  fwrite(OUTS, file.path(sdir, "outputs_O1_O5.tsv"), sep = "\t", na = "NA")
}
g <- function(o, q, s, dec = NA, coh = NA, col = "rate") {
  x <- OUTS[output == o & quantity == q & set == s]
  x <- if (length(dec) == 1L && is.na(dec)) x[is.na(decile)] else x[decile %in% dec]
  x <- if (length(coh) == 1L && is.na(coh)) x[is.na(spiked_cohort)] else x[spiked_cohort %in% coh]
  if (!nrow(x)) return(NA_real_)
  x <- x[order(decile, spiked_cohort)]
  x[[col]]
}

# ---------------------------------------------------------------- thresholds, controls, predictions, verdict
CHK <- list()
ck <- function(name, kind, prediction, observed, met) CHK[[length(CHK) + 1L]] <<- data.table(
  name = name, kind = kind, prediction = prediction, observed = observed, met = met)
verdict <- NULL; reading <- character(); key <- list(); dl <- NULL; o3u <- NULL
if (!is.null(stopped)) {
  verdict <- sprintf("STOPPED by %s", stopped)
  reading <- c(reading, sprintf("Stopped (%s): no number from the run is used.", stopped))
  if (stopped == "S7") { verdict <- "UNTESTABLE (S7)"; reading <- sprintf("Untestable (S7): %d of 50 replicates completed; all outputs untestable.", R) }
} else if (!is.null(OUTS)) {
  # T-O2
  o2_hi <- c(g("O2", "P(supported)", "U0", col = "boot_hi"), g("O2", "P(supported)", "U0", dec = 1:10, col = "boot_hi"))
  t_o2 <- all(o2_hi <= 0.05)
  ck("T-O2", "threshold", "false supported rate <= 0.05 overall and in every decile (upper bootstrap bound)",
     sprintf("overall %s [%s, %s]; max decile upper bound %s", f4(g("O2", "P(supported)", "U0")), f4(g("O2", "P(supported)", "U0", col = "boot_lo")),
             f4(o2_hi[1]), f4(max(o2_hi[-1]))), t_o2)
  # recall by decile, detection limit
  lv <- c(0.5, 0.75, 1.0, 1.5); lsets <- c("U0.5", "U0.75", "U1.0", "U1.5")
  rec <- sapply(lsets, function(s) g("O1", "P(supported)", s, dec = 1:10))
  dl <- rbindlist(lapply(1:10, function(d) {
    r <- rec[d, ]; k <- which(r >= 0.8)[1]
    lim <- if (is.na(k)) "above 1.5" else if (k == 1) "at or below 0.5" else
      formatC(lv[k - 1] + (0.8 - r[k - 1]) / (r[k] - r[k - 1]) * (lv[k] - lv[k - 1]), format = "f", digits = 3)
    data.table(decile = d, recall_U0.5 = r[1], recall_U0.75 = r[2], recall_U1.0 = r[3], recall_U1.5 = r[4], detection_limit_log2FC = lim,
               U1.0_lower = g("O1", "P(supported)", "U1.0", dec = d, col = "boot_lo"),
               T_O1_adequate_0.8 = g("O1", "P(supported)", "U1.0", dec = d, col = "boot_lo") >= 0.8,
               T_O1_unlikely_0.95 = g("O1", "P(supported)", "U1.0", dec = d, col = "boot_lo") >= 0.95)
  }))
  fwrite(dl, file.path(sdir, "O1_detection_limit_by_decile.tsv"), sep = "\t")
  ck("T-O1", "threshold", "recall at delta 1.0 >= 0.8 (adequate) / >= 0.95 (missing call makes 2-fold unlikely) per decile, lower bound",
     sprintf("deciles adequate: %s; deciles >= 0.95: %s", paste(dl[T_O1_adequate_0.8 == TRUE, decile], collapse = ","),
             paste(dl[T_O1_unlikely_0.95 == TRUE, decile], collapse = ",")), NA)
  # PC-A
  pc_a <- all(gate_pass)
  ck("PC-A", "positive control", "G-A, G-B and G-C reproduce", sprintf("G-A %s, G-B %s, G-C %s", gate_pass["G-A"], gate_pass["G-B"], gate_pass["G-C"]), pc_a)
  # PC-B (binding part: U1.0 deciles 3-10 >= 0.95, point estimate)
  pcb_main <- all(rec[3:10, "U1.0"] >= 0.95)
  ck("PC-B", "positive control", "recall at U1.0 >= 0.95 in deciles 3-10", paste(f3(rec[3:10, "U1.0"]), collapse = ","), pcb_main)
  ck("PC-B(ii)", "positive control prediction", "recall at U1.5 >= 0.99 in every decile", paste(f3(rec[, "U1.5"]), collapse = ","), all(rec[, "U1.5"] >= 0.99))
  ck("PC-B(iii)", "positive control prediction", "recall at U0.5 between 0.35 and 0.65 in deciles 5-10", paste(f3(rec[5:10, "U0.5"]), collapse = ","),
     all(rec[5:10, "U0.5"] >= 0.35 & rec[5:10, "U0.5"] <= 0.65))
  # NC-1, NC-2, NC-4
  o2 <- g("O2", "P(supported)", "U0")
  ck("NC-1", "negative control", "U0 false supported rate <= 0.005 overall", f4(o2), o2 <= 0.005)
  ck("NC-2", "negative control", "median unspiked supported per replicate < 20 (S5 at 105)", med_bg, med_bg < 20)
  ck("NC-4", "negative control", ">= 0.99 of U1.5 deciles 5-10 positive aligned change; U0 median within 0.02 of 0",
     sprintf("%s; %s", f4(nc4_pos), f4(nc4_u0)), nc4_pos >= 0.99 && abs(nc4_u0) <= 0.02)
  if (!S$S4) {
    tn <- function(s, d = NA, col = "rate") g("O3", "P(tested_negative)", s, dec = d, col = col)
    ck("NC-3", "negative control", "U1.5 false tested_negative < 0.005", f4(tn("U1.5")), tn("U1.5") < 0.005)
    p_u0 <- tn("U0"); p_u0d <- tn("U0", 1:10)
    rho <- suppressWarnings(cor(1:10, p_u0d, method = "spearman"))
    ck("PC-C", "positive control", "P(tested_negative | U0) between 0.4 and 0.95 overall and rising with decile",
       sprintf("overall %s; Spearman rho with decile %s; decile 1 %s, decile 10 %s", f3(p_u0), f3(rho), f3(p_u0d[1]), f3(p_u0d[10])),
       p_u0 >= 0.4 && p_u0 <= 0.95 && isTRUE(rho > 0))
    # T-O3 per decile and overall
    o3u <- rbindlist(lapply(c(NA, 1:10), function(d) {
      hi <- sapply(c("U0.75", "U1.0", "H1.0", "U0.5", "H0.5"), function(s) tn(s, d, "boot_hi"))
      data.table(decile = d, U0.75_hi = hi[1], U1.0_hi = hi[2], H1.0_hi = hi[3], U0.5_hi = hi[4], H0.5_hi = hi[5],
                 usable = all(hi <= 0.05))
    }))
    fwrite(o3u, file.path(sdir, "T_O3_usable_by_decile.tsv"), sep = "\t")
    ck("T-O3", "threshold", "L09-01 rule usable in a decile: upper bound <= 0.05 at U0.75, U1.0, H1.0, U0.5 and H0.5",
       sprintf("overall usable %s; usable deciles: %s", o3u[is.na(decile), usable], paste(o3u[!is.na(decile) & usable == TRUE, decile], collapse = ",")),
       o3u[is.na(decile), usable])
    # O3 predictions
    u75d <- tn("U0.75", 1:10)
    ck("P-O3a", "prediction", "U0.75 false tested_negative 0.005-0.05 overall, highest in deciles 1 and 2",
       sprintf("%s; deciles 1-2 %s, max of deciles 3-10 %s", f4(tn("U0.75")), paste(f4(u75d[1:2]), collapse = ","), f4(max(u75d[3:10]))),
       tn("U0.75") >= 0.005 && tn("U0.75") <= 0.05 && min(u75d[1:2]) > max(u75d[3:10]))
    ck("P-O3b", "prediction", "U1.0 false tested_negative < 0.01", f4(tn("U1.0")), tn("U1.0") < 0.01)
    ck("P-O3c", "prediction", "U0.5 size 0.01-0.04 overall, every decile <= 0.06", sprintf("%s; max decile %s", f4(tn("U0.5")), f4(max(tn("U0.5", 1:10)))),
       tn("U0.5") >= 0.01 && tn("U0.5") <= 0.04 && max(tn("U0.5", 1:10)) <= 0.06)
    ck("P-O3d", "prediction", "H1.0 above U1.0", sprintf("%s vs %s", f4(tn("H1.0")), f4(tn("U1.0"))), tn("H1.0") > tn("U1.0"))
    ck("P-O3e", "prediction", "H0.5 size 0.005-0.05 and at least U0.5", sprintf("%s vs %s", f4(tn("H0.5")), f4(tn("U0.5"))),
       tn("H0.5") >= 0.005 && tn("H0.5") <= 0.05 && tn("H0.5") >= tn("U0.5"))
  }
  # PC-D, T-O4, O4 predictions
  if (!S$S3) {
    pd_sup <- g("O4", "P(supported)", "S2.0", coh = "GSE213621")
    pd_own <- g("O4", "P(SD-own | supported)", "S2.0", coh = "GSE213621")
    ck("PC-D", "positive control", "S2.0 in GSE213621: P(supported) >= 0.5 and P(SD-own | supported) >= 0.8",
       sprintf("%s; %s", f3(pd_sup), f3(pd_own)), isTRUE(pd_sup >= 0.5 && pd_own >= 0.8))
    to4 <- g("O4", "P(not SD-any | supported)", "S2.0")
    to4_cells <- OUTS[output == "O4" & quantity == "P(not SD-any | supported)" & set == "S2.0" & !is.na(spiked_cohort) & is.na(decile)]
    crossed <- to4_cells[measured == TRUE & rate > 0.20, spiked_cohort]
    ck("T-O4", "threshold (triggers on point estimate)", "S2.0 pooled P(not SD-any | supported) > 0.20 means LOCO retention cannot label source dependence",
       sprintf("%s [%s, %s]; cohorts crossing (cells with >= 50 supported): %s; contributing cohorts: %s", f3(to4),
               f3(g("O4", "P(not SD-any | supported)", "S2.0", col = "boot_lo")), f3(g("O4", "P(not SD-any | supported)", "S2.0", col = "boot_hi")),
               if (length(crossed)) paste(crossed, collapse = ",") else "none",
               paste(sprintf("%s:%d", to4_cells$spiked_cohort, to4_cells$n_den), collapse = ",")),
       isTRUE(to4 > 0.20))
    s1sup <- sapply(setdiff(COHORTS, "GSE213621"), function(co) g("O4", "P(supported)", "S1.0", coh = co))
    ck("P-O4a", "prediction", "S1.0 P(supported) < 0.1 in each of the four smaller cohorts", paste(names(s1sup), f3(s1sup), sep = ":", collapse = ","), all(s1sup < 0.1))
    ck("P-O4b", "prediction", "S2.0 P(supported and not SD-any) < 0.10 pooled", f3(g("O4", "P(supported and not SD-any)", "S2.0")),
       isTRUE(g("O4", "P(supported and not SD-any)", "S2.0") < 0.10))
    ck("P-O4c", "prediction", "S2.0 P(not SD-any | supported) <= 0.20 pooled", f3(to4), isTRUE(to4 <= 0.20))
    if (!S$S4) {
      s1tn <- sapply(c("GSE126848", "GSE130970"), function(co) g("O4", "P(tested_negative)", "S1.0", coh = co))
      ck("P-O4d", "prediction", "S1.0 P(tested_negative) > 0.2 when spiked in GSE126848 or GSE130970", paste(names(s1tn), f3(s1tn), sep = ":", collapse = ","), all(s1tn > 0.2))
    }
    # T-O5 and O5 predictions
    o5 <- function(s, col = "rate") g("O5", "P(SD-any | supported)", s, col = col)
    to5 <- o5("U1.0")
    ck("T-O5", "threshold (triggers on point estimate)", "U1.0 P(SD-any | supported) > 0.20 means LOCO loss is not evidence of source dependence",
       sprintf("U1.0 %s [%s, %s]; U1.5 %s [%s, %s]", f3(to5), f3(o5("U1.0", "boot_lo")), f3(o5("U1.0", "boot_hi")), f3(o5("U1.5")),
               f3(o5("U1.5", "boot_lo")), f3(o5("U1.5", "boot_hi"))), isTRUE(to5 > 0.20))
    ck("P-O5a", "prediction", "U0.5 SD-any rate near 0.5 (band 0.4-0.6, set by the implementer before the run)", f3(o5("U0.5")), isTRUE(abs(o5("U0.5") - 0.5) <= 0.1))
    ck("P-O5b", "prediction", "U0.75 SD-any rate 0.05-0.30", f3(o5("U0.75")), isTRUE(o5("U0.75") >= 0.05 && o5("U0.75") <= 0.30))
    ck("P-O5c", "prediction", "U1.0 and U1.5 SD-any rate < 0.05", sprintf("%s; %s", f3(o5("U1.0")), f3(o5("U1.5"))), isTRUE(o5("U1.0") < 0.05 && o5("U1.5") < 0.05))
    ck("P-O5d", "prediction", "H1.0 SD-any rate above U1.0", sprintf("%s vs %s", f3(o5("H1.0")), f3(o5("U1.0"))), isTRUE(o5("H1.0") > o5("U1.0")))
  }
  CK <- rbindlist(CHK)
  fwrite(CK, file.path(sdir, "controls_thresholds_predictions.tsv"), sep = "\t")

  # ---- verdict by the Reading rules
  if (!t_o2) {
    verdict <- "T-O2 NOT MET: the recipe calls unchanged genes on permuted labels more often than its BH level allows; no recall or detection-limit sentence"
  } else if (pc_a && pcb_main) {
    verdict <- "PASS (T-O2 met, PC-A and PC-B pass, no stop rule)"
  } else {
    verdict <- sprintf("NOT PASS: T-O2 met but %s", paste(c(if (!pc_a) sprintf("PC-A failed (%s)", paste(names(gate_pass)[!gate_pass], collapse = ",")),
                                                            if (!pcb_main) "PC-B failed in deciles 3-10 at U1.0 (simulation or recipe copy wrong; no recall-based sentence)"), collapse = "; "))
  }
  key$verdict <- verdict
  # reading sentences
  u1 <- rec[, "U1.0"]
  if (t_o2 && pcb_main) {
    lim_hi <- dl[decile == 10, detection_limit_log2FC]; lim_lo <- dl[decile == 1, detection_limit_log2FC]
    reading <- c(reading, sprintf(paste0("Resource sentence: \"In %d replicates that added known log2 fold changes to the 774-unit counts after permuting ",
      "disease labels within cohort, the pooled fit called a uniform 1.0 log2FC change in %s to %s%% of genes per expression decile and called %s%% of ",
      "unchanged genes; the smallest change detected with 80%% recall ranged from %s log2FC in the most highly expressed decile to %s in the least.\""),
      R, pct(min(u1)), pct(max(u1)), pct(o2), lim_hi, lim_lo))
    lowd <- dl[T_O1_adequate_0.8 == FALSE, decile]
    if (length(lowd)) reading <- c(reading, sprintf(paste0("O1 below 0.8 at delta 1.0 (lower bound) in deciles %s: a missing bulk call there does not exclude a ",
      "2-fold disease-versus-control change; those genes are below the detection limit rather than unchanged."), paste(lowd, collapse = ",")))
  }
  if (!t_o2) reading <- c(reading, sprintf("False supported rate by decile (upper bound): %s. The canonical rule is not calibrated on these 774 units.",
                                           paste(f4(o2_hi[-1]), collapse = ",")))
  if (t_o2 && med_bg >= 20 && med_bg <= 105) reading <- c(reading, sprintf("NC-2 missed (median unspiked supported %s): every sentence using these outputs states this background count.", med_bg))
  if (!S$S4 && !is.null(o3u)) {
    ok_d <- o3u[!is.na(decile) & usable == TRUE, decile]; bad_d <- o3u[!is.na(decile) & usable == FALSE, decile]
    ea <- sprintf("median BH-effective alpha %s against 0.0177 on the adopted per-cohort fits (ratio %s)", f4(median(RS$eff_alpha)),
                  f3(median(RS$eff_alpha) / (0.05 * 8548 / 24115)))
    if (length(ok_d)) reading <- c(reading, sprintf("O3 within T-O3 in deciles %s: the L09-01 rule called a 0.75 log2FC change tested_negative in %s of genes (%s). The Catalog still shows 0 strict tested_negative calls; adoption is a user decision.",
      paste(ok_d, collapse = ","), paste(sprintf("decile %d: %s%%", ok_d, pct(tn("U0.75", ok_d))), collapse = ", "), ea))
    if (length(bad_d)) reading <- c(reading, sprintf("O3 above T-O3 in deciles %s: the equivalence rule cannot certify absence of a 0.75 log2FC change there; any future adoption must leave them indeterminate (%s).",
      paste(bad_d, collapse = ","), ea))
    if (p_u0 < 0.2) reading <- c(reading, sprintf("PC-C: power on true nulls %s < 0.2; the rule has too little power on these cohorts to be worth adopting.", f3(p_u0)))
  }
  if (S$S4) reading <- c(reading, "O3 untestable (S4: G-C failed).")
  if (S$S3) reading <- c(reading, "O4 and O5 untestable (S3: G-B failed).")
  if (!S$S3) {
    if (isTRUE(to4 > 0.20)) reading <- c(reading, sprintf("O4 above T-O4: LOCO retention may not be described as a test of cohort dependence; a one-cohort effect can pass as a cohort-robust call in %s%% of cases (S2.0, pooled over cohorts).", pct(to4)))
    else reading <- c(reading, sprintf("O4 at or below T-O4: in these simulations a supported effect confined to one cohort lost its call in at least one LOCO refit in %s%% of cases (S2.0, denominator from %s). This does not show that LOCO retention proves cohort independence for real genes.",
      pct(1 - to4), paste(sprintf("%s:%d", to4_cells$spiked_cohort, to4_cells$n_den), collapse = ", ")))
    if (length(crossed)) reading <- c(reading, sprintf("T-O4 crossed in cohort cells: %s.", paste(crossed, collapse = ",")))
    if (isTRUE(to5 > 0.20)) reading <- c(reading, sprintf("O5 above T-O5 at U1.0 (%s): the 828 calls not kept in all five LOCO refits may not be called source-dependent; loss in a LOCO refit is common for uniform effects of that size.", f3(to5)))
    else if (isTRUE(o5("U1.5") <= 0.20)) reading <- c(reading, sprintf("O5 at or below T-O5 at U1.0 and U1.5: a uniform effect of that size kept its call in all five LOCO refits in %s%% (U1.0) and %s%% (U1.5) of cases; a LOCO loss for a call of that size is uncommon under a uniform effect. No call is relabelled.",
      pct(1 - to5), pct(1 - o5("U1.5"))))
    if (isTRUE(to5 <= 0.20) && isTRUE(o5("U1.5") <= 0.20) && (isTRUE(o5("U0.5") > 0.20) || isTRUE(o5("U0.75") > 0.20)))
      reading <- c(reading, sprintf("O5 high only at U0.5 (%s) and U0.75 (%s): LOCO losses concentrate in calls near the 0.5 floor.", f3(o5("U0.5")), f3(o5("U0.75"))))
  }
}

# ---------------------------------------------------------------- D5 lookup, figures
if (!is.null(dl)) {
  st <- fread(DEG)
  dec <- st[, .(gene, symbol, AveExpr, logFC, padj)][order(gene)]
  dec[, decile := as.integer(cut(frank(AveExpr, ties.method = "first"), breaks = quantile(seq_len(.N), seq(0, 1, 0.1), type = 7, names = FALSE),
                                 include.lowest = TRUE))]
  d5 <- dec[!(padj < 0.05 & abs(logFC) > 0.5)][dl[, .(decile, detection_limit_log2FC)], on = "decile"][order(gene)]
  fwrite(d5[, .(gene, symbol, AveExpr, decile, detection_limit_log2FC)], file.path(sdir, "D5_detection_limit_lookup_not_supported_genes.tsv.gz"), sep = "\t")
  edges <- dec[, .(lo = min(AveExpr), hi = max(AveExpr)), by = decile][order(decile)]
  xlab <- sprintf("%d\n%.1f\n%.1f", edges$decile, edges$lo, edges$hi)
  dodge <- position_dodge(width = 0.45)
  th <- theme_classic(base_size = 6) + theme(text = element_text(size = 6, colour = "black", face = "plain"),
    axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 6, colour = "black"),
    legend.text = element_text(size = 6, colour = "black"), legend.title = element_text(size = 6, colour = "black"),
    axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25), legend.key.size = unit(6, "pt"))
  cols_u <- c("U0.5" = "#9ccfee", "U0.75" = "#5d9fd3", "U1.0" = "#2f6aa8", "U1.5" = "#1b2f5b")
  pa <- OUTS[output == "O1" & !is.na(decile) & set %in% names(cols_u)]
  pa[, set := factor(set, levels = names(cols_u))]
  p1 <- ggplot(pa, aes(decile, rate, colour = set, shape = set, group = set)) +
    geom_line(linewidth = 0.4, position = dodge) + geom_errorbar(aes(ymin = boot_lo, ymax = boot_hi), width = 0.3, linewidth = 0.25, position = dodge) +
    geom_point(size = 0.9, position = dodge) +
    scale_x_continuous(breaks = 1:10, labels = xlab) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_colour_manual(values = cols_u, name = "Uniform log2FC") + scale_shape_manual(values = c(16, 17, 15, 18), name = "Uniform log2FC") +
    labs(x = "Expression decile (adopted-fit AveExpr minimum, maximum)", y = "Recall: P(supported)") + th
  ggsave(file.path(sdir, "figures", "H81_panel_a_recall_by_decile.pdf"), p1, width = 4.0, height = 2.2, units = "in", device = grDevices::pdf, useDingbats = FALSE)
  message("Panel a: recall of the canonical bulk rule (BH q < 0.05 and |log2FC| > 0.5) for uniform spiked log2FC by expression decile; bars are 2.5-97.5 percentiles of 2,000 replicate-cluster bootstrap resamples over ", R, " replicates.")
  if (!S$S4) {
    cols_b <- c("U0.5" = "#9ccfee", "U0.75" = "#5d9fd3", "U1.0" = "#2f6aa8", "H0.5" = "#e37faf", "H1.0" = "#7c256f")
    pb <- OUTS[output == "O3" & !is.na(decile) & set %in% names(cols_b)]
    pb[, set := factor(set, levels = names(cols_b))]
    p2 <- ggplot(pb, aes(decile, rate, colour = set, shape = set, linetype = set, group = set)) +
      geom_hline(yintercept = 0.05, linewidth = 0.25, linetype = "dashed", colour = "black") +
      geom_line(linewidth = 0.4, position = dodge) +
      geom_errorbar(aes(ymin = boot_lo, ymax = boot_hi), width = 0.3, linewidth = 0.25, linetype = "solid", position = dodge) +
      geom_point(size = 0.9, position = dodge) +
      scale_x_continuous(breaks = 1:10, labels = xlab) +
      scale_colour_manual(values = cols_b, name = "Spiked set") + scale_shape_manual(values = c(16, 17, 15, 1, 2), name = "Spiked set") +
      scale_linetype_manual(values = c("solid", "solid", "solid", "dotted", "dotted"), name = "Spiked set") +
      labs(x = "Expression decile (adopted-fit AveExpr minimum, maximum)", y = "P(tested_negative), L09-01 rule") + th
    ggsave(file.path(sdir, "figures", "H81_panel_b_false_tested_negative_by_decile.pdf"), p2, width = 4.0, height = 2.2, units = "in", device = grDevices::pdf, useDingbats = FALSE)
    message("Panel b: rate at which the unadopted L09-01 equivalence rule (random effects over five cohorts, TOST at margin 0.5, BH over the replicate universe) called a spiked gene tested_negative; U = uniform, H = between-cohort tau 0.306; dashed line 0.05.")
  }
  if (!S$S3) {
    pc <- OUTS[output == "O5" & is.na(decile)]
    pc[, x := as.numeric(sub("^[UH]", "", set))]; pc[, kind := ifelse(startsWith(set, "U"), "Uniform", "Heterogeneous (tau 0.306)")]
    pc[kind != "Uniform", x := x + 0.04]
    ps <- OUTS[output == "D3" & quantity == "P(SD-any | supported)" & set %in% c("S1.0", "S2.0") & !is.na(spiked_cohort) & is.na(decile) & measured == TRUE]
    ps[, x := as.numeric(sub("^S", "", set))]
    cols_c <- c("GSE126848" = "#40b499", "GSE130970" = "#4aa2c2", "GSE135251" = "#6470d0", "GSE162694" = "#d980dc", "GSE213621" = "#ea868d")
    p3 <- ggplot() +
      geom_line(data = pc[kind == "Uniform"], aes(x, rate), linewidth = 0.4, colour = "black") +
      geom_errorbar(data = pc, aes(x, ymin = boot_lo, ymax = boot_hi), width = 0.04, linewidth = 0.25, colour = "black") +
      geom_point(data = pc, aes(x, rate, shape = kind), size = 1.1, colour = "black") +
      geom_errorbar(data = ps, aes(x + 0.03 * (match(spiked_cohort, names(cols_c)) - 3), ymin = boot_lo, ymax = boot_hi, colour = spiked_cohort), width = 0.03, linewidth = 0.25) +
      geom_point(data = ps, aes(x + 0.03 * (match(spiked_cohort, names(cols_c)) - 3), rate, colour = spiked_cohort), shape = 18, size = 1.3) +
      scale_shape_manual(values = c("Uniform" = 16, "Heterogeneous (tau 0.306)" = 1), name = "All-cohort effect") +
      scale_colour_manual(values = cols_c, name = "One-cohort effect (S1.0, S2.0)") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + scale_x_continuous(breaks = c(0.5, 0.75, 1, 1.5, 2)) +
      labs(x = "Spiked log2FC", y = "P(SD-any | supported)") + th
    ggsave(file.path(sdir, "figures", "H81_panel_c_sd_any_by_delta.pdf"), p3, width = 3.6, height = 2.2, units = "in", device = grDevices::pdf, useDingbats = FALSE)
    message("Panel c: fraction of supported spiked genes that lose the canonical call in at least one leave-one-cohort-out refit (SD-any); one-cohort cells with fewer than 50 supported gene-replicates are not plotted (S8).")
  }
}

# ---------------------------------------------------------------- RESULT_H81.md
status_word <- if (!is.null(stopped)) (if (stopped == "S7") "untestable" else "stopped_by_rule") else "completed"
L <- c(sprintf("# RESULT H81%s", if (debug) " (DEBUG RUN: reduced data, non-primary seeds; not a result)" else ""), "",
       sprintf("Run: `%s`. Prespec sha256 `%s` (verified). Replicates completed: %d of 50. Bootstrap: %d replicate-cluster resamples, seed %d.", out, PRESPEC_SHA, R, B, bseed), "",
       sprintf("**Verdict (computed by the prespec Reading rules):** %s", verdict), "",
       "## Gates and stop rules", "",
       sprintf("- G1 inputs: %d of %d sha256 match.", sum(g1h$match), nrow(g1h)),
       sprintf("- G-A %s, G-B %s, G-C %s (details in `guard/guard_checks.tsv`).", gate_pass["G-A"], gate_pass["G-B"], gate_pass["G-C"]),
       sprintf("- Stop rules: %s.", paste(names(S), unlist(lapply(S, function(x) if (is.na(x)) "not evaluated" else if (x) "TRIGGERED" else "not triggered")), sep = " ", collapse = "; ")),
       if (R > 0) sprintf("- S5 background: median supported among unspiked genes %s (stop above 105). S6 thinning: median aligned change U1.0 deciles 5-10 %s (stop outside 0.9-1.1).", med_bg, f4(s6_med)) else NULL,
       "")
if (!is.null(OUTS)) {
  esc <- function(z) gsub("|", "\\|", z, fixed = TRUE)
  tab <- function(x) c("| output | quantity | set | cohort | decile | n | rate | bootstrap 95% | Clopper-Pearson 95% |", "|---|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %s | %s | %s | %s | %d | %s | %s-%s | %s-%s |", x$output, esc(x$quantity), x$set, ifelse(is.na(x$spiked_cohort), "pooled", x$spiked_cohort),
            ifelse(is.na(x$decile), "all", x$decile), as.integer(x$n_den), ifelse(x$conditional & !x$measured, "not measured (S8)", f4(x$rate)),
            f4(x$boot_lo), f4(x$boot_hi), f4(x$cp_lo), f4(x$cp_hi)))
  L <- c(L, "## Headline outputs (overall)", "", tab(OUTS[is.na(decile) & is.na(spiked_cohort) & output %in% c("O1", "O2", "O3", "O5")]), "",
         "## O4 one-cohort effects (pooled and by spiked cohort)", "", tab(OUTS[output %in% c("O4", "D4") & is.na(decile)]), "",
         "## Detection limit by decile (O1, T-O1)", "", "| decile | recall U0.5 | U0.75 | U1.0 | U1.5 | U1.0 lower bound | detection limit (log2FC, recall 0.8) |", "|---|---|---|---|---|---|---|",
         sprintf("| %d | %s | %s | %s | %s | %s | %s |", dl$decile, f3(dl$recall_U0.5), f3(dl$recall_U0.75), f3(dl$recall_U1.0), f3(dl$recall_U1.5), f3(dl$U1.0_lower), dl$detection_limit_log2FC), "")
  if (!is.null(o3u)) L <- c(L, "## T-O3 by decile (upper bootstrap bounds of P(tested_negative))", "", "| decile | U0.75 | U1.0 | H1.0 | U0.5 | H0.5 | usable |", "|---|---|---|---|---|---|---|",
         sprintf("| %s | %s | %s | %s | %s | %s | %s |", ifelse(is.na(o3u$decile), "all", o3u$decile), f4(o3u$U0.75_hi), f4(o3u$U1.0_hi), f4(o3u$H1.0_hi), f4(o3u$U0.5_hi), f4(o3u$H0.5_hi), o3u$usable), "")
  L <- c(L, "## Thresholds, controls and predictions", "", "| name | kind | prediction | observed | met |", "|---|---|---|---|---|",
         sprintf("| %s | %s | %s | %s | %s |", CK$name, CK$kind, esc(CK$prediction), esc(CK$observed), CK$met), "")
}
L <- c(L, "## Reading", "", paste("-", reading), "",
       "## Files", "", "- `summary/diagnostics.tsv` (diagnostics 1-12, written before the outputs)", "- `summary/outputs_O1_O5.tsv` (every cell, overall, by decile, by cohort)",
       "- `summary/controls_thresholds_predictions.tsv`, `summary/O1_detection_limit_by_decile.tsv`, `summary/T_O3_usable_by_decile.tsv`",
       "- `summary/D5_detection_limit_lookup_not_supported_genes.tsv.gz` (lookup only; no gene relabelled)", "- `summary/figures/H81_panel_[abc]_*.pdf`", "",
       "Unit: spiked genes inside replicates are calibration draws on one fixed set of 774 units, not donors. Every rate describes the recipe on these units and generalizes to nothing else. No p value is reported.")
writeLines(L, file.path(sdir, "RESULT_H81.md"))
rf <- file.path(out, "RESULT_H81.md")
if (!file.exists(rf)) writeLines(L, rf)
writeLines(c(sprintf("status\t%s", status_word), sprintf("verdict\t%s", verdict)), file.path(sdir, "status.tsv"))
writeLines(capture.output(sessionInfo()), file.path(sdir, "sessionInfo_summary.txt"))
say("VERDICT: %s", verdict)
say("DONE summary (%s)", status_word)
