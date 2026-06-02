#!/usr/bin/env Rscript
# 319_nichenet_bulk_reverse_validation.R
#
# Analysis B2 (v1) — Reverse validation of scRNA LIANA LR predictions in bulk.
#
# Spec: docs/superpowers/specs/2026-04-17-cell-type-resolved-masld-biology-design.md
#
# Strategy:
#   1. Load LIANA differential LR pairs (score_masld, score_control, score_diff).
#   2. Load bulk dream DE (logFC + padj by gene symbol).
#   3. For each LIANA-predicted MASLD-enriched LR pair:
#        - Get ligand bulk logFC + padj
#        - Get receptor bulk logFC + padj
#        - Define concordance: sign(score_diff) == sign(ligand_lfc) AND
#                              sign(score_diff) == sign(receptor_lfc)
#   4. Compute:
#        - Overall concordance rate (vs. expected 0.25 under null)
#        - Concordance stratified by source-target cell-type pair
#        - fgsea of LIANA-ranked LR genes in bulk t-stat ranking
#   5. Produce per-LR validation table + per-cell-type concordance summary.
#
# v2 (needs install): NicheNet v2 ligand-target fgsea + deconv-weighted LR
# bulk co-expression with permutation null.
#
# Env: rnaseq
# Outputs: Analysis/SingleCell/results_gpu_v2/ccc/

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
LIANA   <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data")
INT_RES <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
OUTDIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PADJ_BULK <- 0.10   # bulk sig threshold for call-level flag
SCORE_DIFF_MIN <- 0.10  # LIANA score_diff threshold for "strongly directional"

message("[1] Loading LIANA differential LR pairs...")
liana <- fread(file.path(LIANA, "liana_differential_interactions.csv"))
message(sprintf("  LR pairs: %d", nrow(liana)))
message(sprintf("  Cell type pairs: %d",
                length(unique(paste(liana$source, liana$target, sep = "->")))))

message("[2] Loading bulk dream DE (symbol-indexed)...")
bulk <- fread(file.path(INT_RES, "dream_results_ashr.csv"),
              select = c("gene","symbol","logFC","t","padj"))
bulk <- bulk[!is.na(symbol) & symbol != ""]
bulk <- bulk[order(-abs(t))][!duplicated(symbol)]  # keep highest |t| per symbol

message("[3] Annotating LIANA pairs with bulk ligand/receptor DE...")
# Join ligand
setnames(bulk, c("logFC","t","padj"), c("lig_lfc","lig_t","lig_padj"))
liana <- merge(liana,
               bulk[, .(symbol, lig_lfc, lig_t, lig_padj)],
               by.x = "ligand_complex", by.y = "symbol", all.x = TRUE)

# Re-name + join receptor
setnames(bulk, c("lig_lfc","lig_t","lig_padj"), c("rec_lfc","rec_t","rec_padj"))
liana <- merge(liana,
               bulk[, .(symbol, rec_lfc, rec_t, rec_padj)],
               by.x = "receptor_complex", by.y = "symbol", all.x = TRUE)

n_both_found <- sum(!is.na(liana$lig_lfc) & !is.na(liana$rec_lfc))
message(sprintf("  LR pairs with both ligand + receptor in bulk: %d (%.1f%%)",
                n_both_found, 100 * n_both_found / nrow(liana)))

message("[4] Computing concordance per LR pair...")
liana[, score_sign := sign(score_diff)]
liana[, lig_concordant := !is.na(lig_lfc) & sign(lig_lfc) == score_sign]
liana[, rec_concordant := !is.na(rec_lfc) & sign(rec_lfc) == score_sign]
liana[, both_concordant := lig_concordant & rec_concordant]

# Strong pairs: score_diff beyond threshold
strong <- liana[abs(score_diff) >= SCORE_DIFF_MIN & !is.na(lig_lfc) & !is.na(rec_lfc)]
overall_concordance <- mean(strong$both_concordant, na.rm = TRUE)
n_strong <- nrow(strong)

# ---------------------------------------------------------------------------
# T1.13 2026-04-22: replace theoretical 0.25 null with empirical shuffle null
# ---------------------------------------------------------------------------
# Team 3 §4 C1 flagged that `P(both concordant) = 0.25` assumes balanced bulk
# LFC signs. Bulk dream DE has a directional bias (more up- than down-regulated
# at |LFC| threshold), so the theoretical null underestimates the random rate.
# Replace with a permutation null that preserves (a) the marginal bulk LFC
# distribution, (b) the score_diff sign distribution, by shuffling bulk gene
# labels onto LR pairs. Shuffle ligand and receptor independently to reflect
# the expected rate under random LR pairing.
# ---------------------------------------------------------------------------
compute_rate <- function(score_sign_vec, lig_lfc_vec, rec_lfc_vec) {
  lig_ok <- !is.na(lig_lfc_vec) & sign(lig_lfc_vec) == score_sign_vec
  rec_ok <- !is.na(rec_lfc_vec) & sign(rec_lfc_vec) == score_sign_vec
  mean(lig_ok & rec_ok, na.rm = TRUE)
}

set.seed(42)
N_SHUFFLES <- 1000L
bulk_lfc_pool <- bulk$bulk_lfc  # renamed after step [6] below — but not yet
# The `bulk` object was renamed to rec_* after the join; reconstruct a clean
# pool of bulk logFCs usable for sampling. `strong$lig_lfc` and `strong$rec_lfc`
# are already the matched bulk LFC values for observed LR pairs; sample from
# those in aggregate to preserve the empirical distribution.
.t113_lfc_pool <- c(strong$lig_lfc, strong$rec_lfc)
.t113_lfc_pool <- .t113_lfc_pool[!is.na(.t113_lfc_pool)]

if (n_strong > 10 && length(.t113_lfc_pool) > 0) {
  null_rates <- replicate(N_SHUFFLES, {
    # Shuffle bulk LFC independently for ligand and receptor slots
    lig_sh <- sample(.t113_lfc_pool, n_strong, replace = TRUE)
    rec_sh <- sample(.t113_lfc_pool, n_strong, replace = TRUE)
    compute_rate(strong$score_sign, lig_sh, rec_sh)
  })
  empirical_null_mean <- mean(null_rates)
  empirical_null_sd   <- sd(null_rates)
  empirical_p         <- (sum(null_rates >= overall_concordance) + 1) / (N_SHUFFLES + 1)
  # Retain theoretical 0.25 test for comparison (legacy)
  binom_p <- binom.test(sum(strong$both_concordant), n_strong,
                        p = 0.25, alternative = "greater")$p.value
} else {
  null_rates          <- NA_real_
  empirical_null_mean <- NA_real_
  empirical_null_sd   <- NA_real_
  empirical_p         <- NA_real_
  binom_p             <- NA_real_
}

message(sprintf("  Strong LIANA pairs (|score_diff| >= %.2f): %d", SCORE_DIFF_MIN, n_strong))
message(sprintf("  Both-concordant rate: %.1f%%", 100 * overall_concordance))
message(sprintf("  Empirical null (%d shuffles): mean=%.3f, sd=%.3f, emp_p=%.3g",
                N_SHUFFLES, empirical_null_mean, empirical_null_sd, empirical_p))
message(sprintf("  Legacy theoretical null (p=0.25): binom_p=%.3g", binom_p))

# Save null distribution for provenance
if (!all(is.na(null_rates))) {
  fwrite(data.table(shuffle_idx = seq_len(N_SHUFFLES),
                    null_both_concordant_rate = null_rates),
         file.path(OUTDIR, "liana_bulk_reverse_null.csv"))
}
# --- end T1.13 empirical null block ---

message("[5] Per cell-type-pair concordance...")
ct_pair <- strong[, .(
  n_pairs = .N,
  frac_lig_concordant = mean(lig_concordant),
  frac_rec_concordant = mean(rec_concordant),
  frac_both_concordant = mean(both_concordant)
), by = .(source, target)][order(-n_pairs)]

message("[6] fgsea: LIANA MASLD-ranked ligands + receptors enriched in bulk t-stat?")
# Rank bulk genes by signed t-stat
setnames(bulk, c("rec_lfc","rec_t","rec_padj"), c("bulk_lfc","bulk_t","bulk_padj"))
ranks <- setNames(bulk$bulk_t, bulk$symbol)
ranks <- ranks[!is.na(ranks)]

# Build gene sets from LIANA: top MASLD-enriched, top Control-enriched (by score_diff)
liana_masld_top_lr <- unique(c(
  liana[score_diff > 0.2 & !is.na(lig_lfc), ligand_complex],
  liana[score_diff > 0.2 & !is.na(rec_lfc), receptor_complex]))
liana_ctrl_top_lr <- unique(c(
  liana[score_diff < -0.2 & !is.na(lig_lfc), ligand_complex],
  liana[score_diff < -0.2 & !is.na(rec_lfc), receptor_complex]))

gene_sets <- list(
  "LIANA_MASLD_up_LR"   = intersect(liana_masld_top_lr, names(ranks)),
  "LIANA_Control_up_LR" = intersect(liana_ctrl_top_lr, names(ranks))
)
message(sprintf("  Gene sets sizes: MASLD-up=%d, Control-up=%d",
                length(gene_sets[["LIANA_MASLD_up_LR"]]),
                length(gene_sets[["LIANA_Control_up_LR"]])))

set.seed(42)
fres <- if (all(sapply(gene_sets, length) > 0)) {
  fgsea(pathways = gene_sets, stats = ranks, eps = 0, nPermSimple = 2000)
} else data.table()

message("[7] Writing outputs...")
fwrite(liana, file.path(OUTDIR, "liana_bulk_concordance_perLR.csv"))
fwrite(ct_pair, file.path(OUTDIR, "liana_bulk_concordance_by_ct_pair.csv"))
if (nrow(fres) > 0) fwrite(fres[, !"leadingEdge"], file.path(OUTDIR, "liana_fgsea_in_bulk.csv"))

summary_lines <- c(
  sprintf("Total LIANA differential LR pairs: %d", nrow(liana)),
  sprintf("Pairs with both ligand + receptor in bulk: %d", n_both_found),
  sprintf("Strong directional LIANA pairs (|score_diff|>=%.2f): %d",
          SCORE_DIFF_MIN, n_strong),
  sprintf("Both-concordant rate: %.1f%%", 100 * overall_concordance),
  sprintf("Empirical null (1000 LR-shuffles): mean=%.3f (sd=%.3f); emp_p=%.3g",
          empirical_null_mean, empirical_null_sd, empirical_p),
  sprintf("Legacy theoretical 0.25 null: binom_p=%.3g", binom_p),
  "",
  "Per cell-type-pair concordance (top 20 by n):",
  capture.output(print(ct_pair[1:20], nrows = 20)),
  "",
  if (nrow(fres) > 0) c("fgsea of LIANA-ranked LR sets in bulk t-stat:",
                         capture.output(print(fres[, !"leadingEdge"], nrows = 10))) else "fgsea skipped (empty gene sets)"
)
writeLines(summary_lines, file.path(OUTDIR, "liana_bulk_concordance_summary.txt"))
writeLines(summary_lines)

message("Done. Outputs in: ", OUTDIR)
