#!/usr/bin/env Rscript
# Code review U7 (H81): independent summary of the stored spike-in replicate tables.
# Written from METHODS.md "Spike-in calibration of bulk evidence states", not from 02_h81_summary.R.
# Rates are pooled over gene-replicates (untested spiked genes count as not supported / not tested_negative).
# Bootstrap: resample the R replicates with replacement B times; ratio of summed numerators to summed
# denominators; 2.5 and 97.5 percentiles. Run once per bootstrap seed given on the command line.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
rdir <- args[1]; outf <- args[2]; seeds <- as.integer(strsplit(args[3], ",")[[1]]); B <- 2000L
fs <- sort(list.files(rdir, pattern = "^rep_\\d+_genes\\.tsv\\.gz$", full.names = TRUE))
D <- rbindlist(lapply(fs, fread))
reps <- sort(unique(D$rep)); R <- length(reps)
D[, ri := match(rep, reps)]
D[is.na(supported), supported := FALSE]; D[is.na(tested_negative), tested_negative := FALSE]
cat("replicates", R, "gene rows", nrow(D), "\n")

rate_cell <- function(rows, ev, CNT) {
  num <- tabulate(D$ri[rows & ev], R); den <- tabulate(D$ri[rows], R)
  rb <- as.vector(CNT %*% num) / as.vector(CNT %*% den)
  q <- quantile(rb[is.finite(rb)], c(0.025, 0.975), names = FALSE)
  c(n = sum(den), rate = sum(num) / sum(den), lo = q[1], hi = q[2])
}
res <- list()
for (s in seeds) {
  set.seed(s, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  CNT <- t(vapply(seq_len(B), function(b) tabulate(sample.int(R, R, replace = TRUE), R), numeric(R)))
  add <- function(what, set, dec, v) res[[length(res) + 1L]] <<- data.table(seed = s, what = what, set = set, decile = dec,
                                                                              n = v[["n"]], rate = v[["rate"]], lo = v[["lo"]], hi = v[["hi"]])
  for (st in c("U0", "U0.5", "U0.75", "U1.0", "U1.5", "H1.0", "H0.5")) {
    add("P(supported)", st, 0L, rate_cell(D$set == st, D$supported, CNT))
    add("P(tested_negative)", st, 0L, rate_cell(D$set == st, D$tested_negative, CNT))
    for (d in 1:10) {
      add("P(supported)", st, d, rate_cell(D$set == st & D$decile == d, D$supported, CNT))
      add("P(tested_negative)", st, d, rate_cell(D$set == st & D$decile == d, D$tested_negative, CNT))
    }
  }
  sup <- D$supported %in% TRUE
  add("P(not SD-any | supported)", "S2.0", 0L, rate_cell(D$set == "S2.0" & sup, !(D$sd_any %in% TRUE), CNT))
  add("P(SD-any | supported)", "U1.0", 0L, rate_cell(D$set == "U1.0" & sup, D$sd_any %in% TRUE, CNT))
  add("P(SD-any | supported)", "U1.5", 0L, rate_cell(D$set == "U1.5" & sup, D$sd_any %in% TRUE, CNT))
}
res <- rbindlist(res)
fwrite(res, outf, sep = "\t")
# derived: deciles with U1.0 lower bound >= 0.98 / 0.8; equivalence-rule usable deciles; detection limit
for (s in seeds) {
  x <- res[seed == s]
  u1 <- x[what == "P(supported)" & set == "U1.0" & decile > 0][order(decile)]
  tn <- x[what == "P(tested_negative)" & set %in% c("U0.75", "U1.0", "H1.0", "U0.5", "H0.5") & decile > 0,
          .(usable = all(hi <= 0.05)), by = decile][order(decile)]
  dl <- vapply(1:10, function(d) {
    r <- x[what == "P(supported)" & set %in% c("U0.5", "U0.75", "U1.0", "U1.5") & decile == d][order(match(set, c("U0.5", "U0.75", "U1.0", "U1.5")))]
    dd <- c(0.5, 0.75, 1.0, 1.5); k <- which(r$rate >= 0.8)[1]
    if (is.na(k)) Inf else if (k == 1) 0.5 else dd[k - 1] + (0.8 - r$rate[k - 1]) / (r$rate[k] - r$rate[k - 1]) * (dd[k] - dd[k - 1])
  }, 0)
  cat(sprintf("seed %d | U1.0 recall d1 %.3f d2 %.3f | U1.0 lower bounds d1-10: %s | min lower d3-10 %.4f | deciles lower>=0.8: %s\n",
              s, u1$rate[1], u1$rate[2], paste(sprintf("%.3f", u1$lo), collapse = ","), min(u1$lo[3:10]),
              paste(u1$decile[u1$lo >= 0.8], collapse = ",")))
  cat(sprintf("seed %d | equivalence usable deciles: %s | detection limit d1-10: %s\n", s,
              paste(tn$decile[tn$usable], collapse = ","), paste(sprintf("%.3f", dl), collapse = ",")))
  o4 <- x[what == "P(not SD-any | supported)"]; o5 <- x[what == "P(SD-any | supported)" & set == "U1.0"]
  u0 <- x[what == "P(supported)" & set == "U0" & decile == 0]
  cat(sprintf("seed %d | S2.0 lost in >=1 LOCO %.4f (n %d) [%.4f, %.4f] | U1.0 kept in all five %.4f (n %d) | U0 supported %.4f [%.4f, %.4f]\n",
              s, 1 - o4$rate, o4$n, 1 - o4$hi, 1 - o4$lo, 1 - o5$rate, o5$n, u0$rate, u0$lo, u0$hi))
}
