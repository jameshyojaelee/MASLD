#!/usr/bin/env Rscript
# ============================================================================
# q2_rrho_concordance.R  (Vacca benchmark, Q2)
#
# WHY: the prior 172/1323 conserved-core overlap with Vacca's 951 signature is
# bottlenecked by ARBITRARY human+mouse LFC/padj cutoffs on BOTH sides. A
# cutoff-FREE rank-rank concordance answers: do our continuous cross-species
# concordance ranking and Vacca's continuous human-progression ranking agree,
# and WHERE does the agreement concentrate (top-up / top-down)?
#
# Two continuous rankings (no thresholds applied to either):
#  (A) OURS  = signed cross-species concordance score
#              = translatability_score * sign(mean_h_lfc)
#       translatability_score is our continuous [0,1] gene-level concordance stat
#       (encodes #diets concordant + module preservation + magnitude); sign by
#       human disease direction so "up in disease & cross-species concordant"
#       ranks high-positive and "down & concordant" ranks high-negative.
#  (B) VACCA = human Severe-vs-Mild progression L2FC, mean(UCAM/VCU + EPoS)
#       (the published human late-progression effect; continuous, all genes).
#
# RRHO2 / RRHO are NOT installed -> manual rank-rank hypergeometric (phyper)
# over a sliding grid of top-k thresholds, in all 4 quadrants:
#   uu = top of A (up) vs top of B (up)        [concordant up]
#   dd = bottom of A (down) vs bottom of B (down) [concordant down]
#   ud = top A up vs bottom B down  (discordant)
#   du = bottom A down vs top B up   (discordant)
# Peak -log10(p) over the grid, sided so concordant overlap is positive.
# Hypergeometric tail: phyper(overlap-1, K, N-K, n, lower.tail=FALSE).
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })
BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CONC  <- file.path(BASE, "Analysis/Cross_Species_Concordance/results")
VACCA <- file.path(BASE, "data/external/vacca_2024")
OUT   <- file.path(CONC, "vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
strip <- function(x) gsub("\\..*", "", x)

# ── (A) OUR continuous signed concordance ranking ───────────────────────────
u <- fread(file.path(CONC, "concordance_atlas_unified.csv"))
u <- u[!is.na(translatability_score) & !is.na(mean_h_lfc) & !is.na(human_symbol) &
       human_symbol != ""]
# signed concordance: magnitude = translatability_score, direction = human LFC.
# genes with mean_h_lfc exactly 0 get sign 0 -> tiny jitter by sign of score_gene
u[, sgn := sign(mean_h_lfc)]
u[sgn == 0, sgn := 1]
u[, our_stat := translatability_score * sgn]
ourdt <- u[, .(symbol = toupper(human_symbol), our_stat)]
# collapse dup symbols by max |stat|
ourdt <- ourdt[order(-abs(our_stat))][!duplicated(symbol)]

# ── (B) VACCA human Severe-vs-Mild progression L2FC ─────────────────────────
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"),
                               sheet = "Table S4"))
prog_cols <- c("L2FC_UCAM/VCU: Severe vs Mild", "L2FC_EPoS: Severe vs Mild")
for (cc in prog_cols) suppressWarnings(s4[, (cc) := as.numeric(get(cc))])
s4[, vacca_stat := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = prog_cols]
s4[, symbol := toupper(GeneSymbol)]
vac <- s4[is.finite(vacca_stat) & !is.na(symbol) & symbol != "",
          .(symbol, vacca_stat)]
vac <- vac[order(-abs(vacca_stat))][!duplicated(symbol)]

# ── common gene universe (cutoff-free: ALL shared genes) ────────────────────
m <- merge(ourdt, vac, by = "symbol")
N <- nrow(m)
cat(sprintf("Common gene universe (both rankings, no cutoffs): N = %d\n", N))
cat(sprintf("OUR signed stat range [%.3f, %.3f]; VACCA L2FC range [%.3f, %.3f]\n",
            min(m$our_stat), max(m$our_stat), min(m$vacca_stat), max(m$vacca_stat)))

# rank: descending for "up" side, ascending for "down" side
m[, rA_up   := frank(-our_stat,   ties.method = "first")]
m[, rA_down := frank( our_stat,   ties.method = "first")]
m[, rB_up   := frank(-vacca_stat, ties.method = "first")]
m[, rB_down := frank( vacca_stat, ties.method = "first")]

# ── sliding-grid rank-rank hypergeometric ───────────────────────────────────
# grid step ~ N/100 (standard RRHO stratification level)
step <- max(1L, floor(N / 100))
ks   <- seq(step, N - step, by = step)
hgrid <- CJ(kA = ks, kB = ks)

# precompute membership matrices for speed
topset <- function(rankcol, k) m$symbol[m[[rankcol]] <= k]
# overlap counter via rank thresholds
overlap_n <- function(rA, rB, kA, kB) sum(m[[rA]] <= kA & m[[rB]] <= kB)

hyper_logp <- function(ov, kA, kB) {
  # P(X >= ov) for drawing kB from N with kA "successes"
  -phyper(ov - 1L, kA, N - kA, kB, lower.tail = FALSE, log.p = TRUE) / log(10)
}

run_quadrant <- function(rA, rB, label) {
  res <- hgrid[, {
    ov <- overlap_n(rA, rB, kA, kB)
    exp_ov <- kA * kB / N
    lp <- hyper_logp(ov, kA, kB)
    .(overlap = ov, expected = exp_ov, neglog10p = lp,
      enrich = ifelse(exp_ov > 0, ov / exp_ov, NA_real_))
  }, by = .(kA, kB)]
  res[, quadrant := label]
  res
}

res_uu <- run_quadrant("rA_up",   "rB_up",   "up_up")     # concordant up
res_dd <- run_quadrant("rA_down", "rB_down", "down_down") # concordant down
res_ud <- run_quadrant("rA_up",   "rB_down", "up_down")   # discordant
res_du <- run_quadrant("rA_down", "rB_up",   "down_up")   # discordant
allres <- rbindlist(list(res_uu, res_dd, res_ud, res_du))

# peak per quadrant
peaks <- allres[, .SD[which.max(neglog10p)], by = quadrant][
  order(-neglog10p)]
cat("\n=== PEAK rank-rank hypergeometric per quadrant ===\n")
print(peaks[, .(quadrant, kA, kB, overlap, expected = round(expected,1),
                enrich = round(enrich,2), neglog10p = round(neglog10p,1))])

# global signed peak: concordant (uu/dd) positive, discordant (ud/du) negative
conc_peak <- allres[quadrant %in% c("up_up","down_down")][which.max(neglog10p)]
disc_peak <- allres[quadrant %in% c("up_down","down_up")][which.max(neglog10p)]
cat(sprintf("\nCONCORDANT peak  -log10p = %.1f  (%s, kA=%d kB=%d, overlap=%d vs expected %.1f, %.2fx)\n",
            conc_peak$neglog10p, conc_peak$quadrant, conc_peak$kA, conc_peak$kB,
            conc_peak$overlap, conc_peak$expected, conc_peak$enrich))
cat(sprintf("DISCORDANT peak  -log10p = %.1f  (%s, kA=%d kB=%d, overlap=%d vs expected %.1f, %.2fx)\n",
            disc_peak$neglog10p, disc_peak$quadrant, disc_peak$kA, disc_peak$kB,
            disc_peak$overlap, disc_peak$expected, disc_peak$enrich))

# global continuous concordance for context (Spearman over full common universe)
sp <- cor(m$our_stat, m$vacca_stat, method = "spearman")
cat(sprintf("\nGlobal Spearman (signed our_stat vs Vacca Severe-vs-Mild L2FC, all %d genes) = %.3f\n",
            N, sp))

# ── write results ───────────────────────────────────────────────────────────
fwrite(allres, file.path(OUT, "q2_rrho_concordance_grid.csv"))
summ <- data.table(
  metric = c("N_common_genes",
             "concordant_peak_neglog10p", "concordant_peak_quadrant",
             "concordant_peak_kA", "concordant_peak_kB",
             "concordant_peak_overlap", "concordant_peak_expected",
             "concordant_peak_enrichment",
             "discordant_peak_neglog10p", "discordant_peak_quadrant",
             "up_up_peak_neglog10p", "down_down_peak_neglog10p",
             "up_down_peak_neglog10p", "down_up_peak_neglog10p",
             "global_spearman"),
  value  = c(N,
             round(conc_peak$neglog10p,2), conc_peak$quadrant,
             conc_peak$kA, conc_peak$kB,
             conc_peak$overlap, round(conc_peak$expected,2),
             round(conc_peak$enrich,3),
             round(disc_peak$neglog10p,2), disc_peak$quadrant,
             round(peaks[quadrant=="up_up", neglog10p],2),
             round(peaks[quadrant=="down_down", neglog10p],2),
             round(peaks[quadrant=="up_down", neglog10p],2),
             round(peaks[quadrant=="down_up", neglog10p],2),
             round(sp,3)))
fwrite(summ, file.path(OUT, "q2_rrho_concordance.csv"))
cat("\nWrote ", file.path(OUT, "q2_rrho_concordance.csv"), "\n", sep="")
cat("Wrote ", file.path(OUT, "q2_rrho_concordance_grid.csv"), "\n", sep="")
