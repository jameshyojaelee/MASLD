#!/usr/bin/env Rscript
# =============================================================================
# B9 Driver score weight perturbation (V5 verification)
# The driver score is a weighted sum of 7 components. Perturb each weight ±20%
# (0.8× and 1.2×), recompute driver_score, re-rank, and track THRB's rank.
# Also include uniform weights (1/7 each).
#
# Pass: THRB stays in top 1% at F3→F4 across all perturbations.
# Output: verification/controls/b9_weight_perturbation.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
DRIVER_CSV <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/driver_scores.csv")
OUT_CSV <- file.path(BASE, "docs/manuscript/verification/controls/b9_weight_perturbation.csv")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

cat("=== B9 weight perturbation ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ---- Load ----
cat("Loading driver_scores.csv ...\n")
ds <- fread(DRIVER_CSV)
cat("Rows:", nrow(ds), "  cols:", ncol(ds), "\n")

# Component columns (from Script 126 lines 548-557)
comp_cols <- c(
  "s1_transition_specificity",
  "s1_is_transition_deg",
  "s2_genetic_causality",
  "s3_transport_cost",
  "s4_bifurcation_divergence",
  "s5_celltype_resolution",
  "s6_cross_species",
  "s7_communication_rewiring"
)
# Note: 8 component columns, but manuscript claim says "7 features".
# Original weights in Script 126:
orig_weights <- c(
  s1_transition_specificity = 0.18,
  s1_is_transition_deg      = 0.13,
  s2_genetic_causality      = 0.18,
  s3_transport_cost         = 0.09,
  s4_bifurcation_divergence = 0.13,
  s5_celltype_resolution    = 0.09,
  s6_cross_species          = 0.10,
  s7_communication_rewiring = 0.10
)

# Coerce component columns to numeric with NA→0
for (c in comp_cols) {
  ds[[c]] <- as.numeric(ds[[c]])
  ds[[c]][is.na(ds[[c]])] <- 0
}

# Transitions to analyse (focus on F3_to_F4 per claim, report all)
transitions <- unique(ds$transition)
cat("Transitions:", paste(transitions, collapse=", "), "\n")

# Function: given weight vector, compute score, rank, return THRB rank
compute_ranks <- function(ds_sub, w, gene_target = "THRB") {
  # w: named vector matching comp_cols
  score <- Reduce("+", Map(function(c) ds_sub[[c]] * w[[c]], comp_cols))
  # Rank: 1 = best (highest score)
  rk <- rank(-score, ties.method = "average")
  n_total <- length(score)
  # THRB
  thrb_idx <- which(ds_sub$gene_symbol == gene_target)
  if (length(thrb_idx) == 0) {
    return(data.table(thrb_rank = NA_real_, thrb_pct = NA_real_, n_total = n_total))
  }
  data.table(
    thrb_rank = rk[thrb_idx[1]],
    thrb_pct = 100 * rk[thrb_idx[1]] / n_total,
    thrb_score = score[thrb_idx[1]],
    n_total = n_total
  )
}

# ---- Build perturbation list ----
perturbations <- list()

# 1. Original
perturbations[["original"]] <- orig_weights

# 2. Each weight ±20% (with others held at original, renormalized)
for (c in comp_cols) {
  for (factor in c(0.8, 1.2)) {
    w <- orig_weights
    w[[c]] <- w[[c]] * factor
    # normalize to sum = 1 (matches original sum)
    w <- w * (sum(orig_weights) / sum(w))
    label <- sprintf("%s_x%.1f", c, factor)
    perturbations[[label]] <- w
  }
}

# 3. All +20% and all -20% (should give same ranks as original after renorm,
# include as sanity check)
perturbations[["uniform"]] <- setNames(rep(1/length(comp_cols), length(comp_cols)),
                                         comp_cols)

# 4. Genetic-upweight (+50%) and DEG-status-downweight (-50%) — adversarial
w_gen <- orig_weights
w_gen[["s2_genetic_causality"]] <- w_gen[["s2_genetic_causality"]] * 1.5
w_gen <- w_gen * (sum(orig_weights) / sum(w_gen))
perturbations[["genetic_upweight_1.5x"]] <- w_gen

w_deg <- orig_weights
w_deg[["s1_is_transition_deg"]] <- w_deg[["s1_is_transition_deg"]] * 0.5
w_deg <- w_deg * (sum(orig_weights) / sum(w_deg))
perturbations[["deg_status_downweight_0.5x"]] <- w_deg

cat("Total perturbations:", length(perturbations), "\n")

# ---- Loop over transitions and perturbations ----
all_rows <- list()
for (trans in transitions) {
  ds_sub <- ds[transition == trans]
  for (pname in names(perturbations)) {
    w <- perturbations[[pname]]
    r <- compute_ranks(ds_sub, w)
    r[, perturbation := pname]
    r[, transition := trans]
    r[, weight_vec := paste(sprintf("%s=%.3f", names(w), w), collapse = "; ")]
    all_rows[[length(all_rows) + 1]] <- r
  }
}
result <- rbindlist(all_rows)
setcolorder(result, c("transition", "perturbation", "thrb_rank", "thrb_pct",
                      "thrb_score", "n_total", "weight_vec"))
# Top-1% flag
result[, thrb_in_top_1pct := thrb_pct <= 1.0]

cat("\n=== F3_to_F4 results ===\n")
print(result[transition == "F3_to_F4"][order(perturbation)],
      topn = length(perturbations))

cat("\nAcross all transitions, fraction where THRB is top 1%:\n")
print(result[, .(frac_top1 = mean(thrb_in_top_1pct, na.rm=TRUE),
                 n = .N), by = perturbation])

fwrite(result, OUT_CSV)
cat("\nSaved:", OUT_CSV, "\n")
cat("End:", format(Sys.time()), "\n")
