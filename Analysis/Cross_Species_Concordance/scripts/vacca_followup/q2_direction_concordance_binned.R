#!/usr/bin/env Rscript
# ============================================================================
# q2_direction_concordance_binned.R   (Vacca follow-up, Q2)
#
# QUESTION: Is the human-side direction agreement between OUR disease signature
# and Vacca's progression (Severe-vs-Mild) signature a real biological signal,
# or a thresholding artifact of the arbitrary LFC/padj cutoffs that define the
# Conserved_Core?
#
# APPROACH: take ALL shared genes (no significance cutoff at all). For each gene
# compute sign(our_human_logFC) vs sign(vacca_severe_vs_mild_logFC). Bin genes by
# |effect size| and compute direction-concordance per bin. If concordance rises
# MONOTONICALLY with effect size and is already > 50% even with no cutoff, the
# 172/1323 overlap is the strong-effect tail of a continuous, real concordance
# gradient — not an artifact of where we drew the LFC line.
#
# Three effect-size axes are reported (robustness):
#   axis A: |our human logFC|                         (our magnitude drives the bin)
#   axis B: |vacca severe-vs-mild logFC|              (their magnitude drives the bin)
#   axis C: min(|our|, |vacca|)  "both-strong"         (conservative joint magnitude)
#
# Baseline: 50% (coin-flip). Per-bin one-sided binomial test vs 0.5.
# Trend: Cochran-Armitage-style logistic slope of match~bin_index + Spearman of
#        per-bin concordance vs bin rank.
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
INT_I <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUT   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# ── Vacca S4: human Severe-vs-Mild progression logFC (mean of UCAM/VCU + EPoS) ──
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"),
                               sheet = "Table S4"))
setnames(s4, c("GeneSymbol", "Gene_used_in_DSEA(0:No/1:Yes)"),
         c("gene", "in_dsea"), skip_absent = TRUE)
hcols <- c("L2FC_UCAM/VCU: Severe vs Mild", "L2FC_EPoS: Severe vs Mild")
s4[, (hcols) := lapply(.SD, as.numeric), .SDcols = hcols]
s4[, vacca_lfc := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = hcols]
s4 <- s4[!is.na(gene) & gene != ""][!duplicated(gene)]
s4 <- s4[is.finite(vacca_lfc)]
cat(sprintf("Vacca S4: %d genes with finite Severe-vs-Mild logFC.\n", nrow(s4)))

# ── OUR human disease-vs-control signature (canonical), by symbol ─────────────
hum <- fread(file.path(INT_I, "canonical_deg_results.csv"))
hum <- hum[!is.na(symbol) & symbol != "" & is.finite(logFC)]
hum <- hum[order(padj)][!duplicated(symbol)]   # keep most-significant per symbol
hum <- hum[, .(gene = symbol, our_lfc = logFC, our_padj = padj)]
cat(sprintf("Our human DEG: %d unique-symbol genes with finite logFC.\n", nrow(hum)))

# ── Merge on shared symbols, NO significance cutoff ──────────────────────────
m <- merge(hum, s4[, .(gene, vacca_lfc, in_dsea)], by = "gene")
m <- m[abs(our_lfc) > 0 & abs(vacca_lfc) > 0]   # drop exact-zero (sign undefined)
m[, match := sign(our_lfc) == sign(vacca_lfc)]
cat(sprintf("Shared genes (no cutoff): %d ; overall direction-concordance = %.1f%%\n",
            nrow(m), 100 * mean(m$match)))

binom_p <- function(x, n) if (n == 0) NA_real_ else
  binom.test(sum(x), n, 0.5, alternative = "greater")$p.value

# ── Binning helper: decile bins of an effect-size axis ───────────────────────
bin_curve <- function(dt, axis_col, axis_label, nbin = 10) {
  d <- copy(dt)
  d[, esize := get(axis_col)]
  # decile bins by the effect-size axis
  qs <- unique(quantile(d$esize, probs = seq(0, 1, length.out = nbin + 1), na.rm = TRUE))
  d[, bin := cut(esize, breaks = qs, include.lowest = TRUE, labels = FALSE)]
  d <- d[!is.na(bin)]
  curve <- d[, .(n = .N,
                 esize_lo = min(esize), esize_hi = max(esize),
                 esize_median = median(esize),
                 n_match = sum(match),
                 concordance_pct = 100 * mean(match),
                 binom_p_gt50 = binom_p(match, .N)),
             by = bin][order(bin)]
  curve[, `:=`(axis = axis_label, bin_rank = seq_len(.N))]
  # monotonic-trend stats over genes (not bins)
  rho_bin <- suppressWarnings(cor(curve$bin_rank, curve$concordance_pct, method = "spearman"))
  fit <- glm(match ~ esize, data = d, family = binomial())
  slope <- coef(summary(fit))["esize", ]
  list(curve = curve, rho_bin = rho_bin,
       logit_slope = unname(slope["Estimate"]),
       logit_p = unname(slope["Pr(>|z|)"]))
}

axes <- list(
  list(col = "abs_our",   lab = "abs_our_logFC"),
  list(col = "abs_vacca", lab = "abs_vacca_logFC"),
  list(col = "abs_min",   lab = "min_abs_both"))
m[, `:=`(abs_our = abs(our_lfc), abs_vacca = abs(vacca_lfc),
         abs_min = pmin(abs(our_lfc), abs(vacca_lfc)))]

all_curves <- list(); trend <- list()
for (a in axes) {
  r <- bin_curve(m, a$col, a$lab, nbin = 10)
  all_curves[[a$lab]] <- r$curve
  trend[[a$lab]] <- data.table(axis = a$lab, rho_bin_vs_rank = r$rho_bin,
                               logit_slope = r$logit_slope, logit_p = r$logit_p)
  cat(sprintf("\n=== Axis: %s ===  (Spearman conc-vs-binrank = %.3f ; logit slope = %.3f, p = %.2e)\n",
              a$lab, r$rho_bin, r$logit_slope, r$logit_p))
  print(r$curve[, .(bin = bin_rank, n, esize_lo = round(esize_lo,2),
                    esize_hi = round(esize_hi,2), n_match,
                    concordance_pct = round(concordance_pct,1),
                    binom_p = signif(binom_p_gt50,2))])
}

curve_dt <- rbindlist(all_curves)
trend_dt <- rbindlist(trend)

# ── Also: simple fixed |logFC| thresholds (interpretable curve) ──────────────
thr <- c(0, 0.25, 0.5, 1.0, 1.5, 2.0)
fixed <- rbindlist(lapply(thr, function(t) {
  sub <- m[abs_min >= t]
  data.table(min_abs_both_threshold = t, n = nrow(sub),
             concordance_pct = 100 * mean(sub$match),
             binom_p_gt50 = binom_p(sub$match, nrow(sub)))
}))
cat("\n=== Concordance vs increasing |min(both)| threshold (cumulative) ===\n")
print(fixed[, .(min_abs_both_threshold, n, concordance_pct = round(concordance_pct,1),
                binom_p = signif(binom_p_gt50,2))])

# in-DSEA (Vacca 951) tail check — are strong-effect genes the DSEA genes?
cat(sprintf("\nGenes with min(|both|) >= 1.0: %d ; fraction in Vacca-951 DSEA sig = %.1f%%\n",
            nrow(m[abs_min >= 1]), 100 * mean(m[abs_min >= 1]$in_dsea == 1, na.rm = TRUE)))
cat(sprintf("Genes with min(|both|) <  0.25: %d ; fraction in Vacca-951 DSEA sig = %.1f%%\n",
            nrow(m[abs_min < 0.25]), 100 * mean(m[abs_min < 0.25]$in_dsea == 1, na.rm = TRUE)))

# ── Write outputs ────────────────────────────────────────────────────────────
fwrite(curve_dt, file.path(OUT, "q2_direction_concordance_binned.csv"))
fwrite(trend_dt, file.path(OUT, "q2_direction_concordance_trend.csv"))
fwrite(fixed,    file.path(OUT, "q2_direction_concordance_threshold.csv"))

cat("\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("HEADLINE: %d shared genes, NO cutoff. Overall dir-concordance = %.1f%%.\n",
            nrow(m), 100 * mean(m$match)))
prim <- all_curves[["min_abs_both"]]
cat(sprintf("  min_abs_both axis: bin1 (weakest) = %.1f%% -> bin10 (strongest) = %.1f%%.\n",
            prim[bin_rank == 1]$concordance_pct, prim[bin_rank == max(bin_rank)]$concordance_pct))
cat(sprintf("  Monotone trend: Spearman(conc, bin-rank) = %.3f ; logistic slope p = %.2e.\n",
            trend_dt[axis == "min_abs_both"]$rho_bin_vs_rank,
            trend_dt[axis == "min_abs_both"]$logit_p))
cat(strrep("=", 78), "\n")
