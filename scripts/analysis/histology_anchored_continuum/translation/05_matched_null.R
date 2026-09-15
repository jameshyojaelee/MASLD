#!/usr/bin/env Rscript

# Size- and expression-decile-matched random gene-set null for the 117 fixed
# programs, on the continuum axes.
#
# THE QUESTION. The extension reports 93 of 117 programs supported under the
# two-axis stage/sex rule. That count has no denominator. Bulk liver expression
# is broadly correlated with the continuum axis, so a random gene set of the
# same size is NOT expected to score zero. Without the matched-random rate the
# 93 is uninterpretable.
#
# This is not a new method. The matched-draw and fast-scoring machinery is
# reused verbatim from scripts/manuscript/program_observability_map/
# t3_vocabulary_lib.R, which already ran this design on the fibrosis and NAS
# axes. Only the outcome, the substrate and the model change: here the axis is
# the continuum score, the cohorts are the extension's two evaluation cohorts,
# and the model/BH/support rule are tr_fit()/tr_meta()/tr_bh() from
# lib_translation.R, unchanged, so the null regenerates the statistic end to end
# through exactly the pipeline that produced the 93.
#
# SCORING IDENTITY. The canonical program score is a WEIGHTED mean of
# within-cohort gene z-scores. A random set has no weights. Comparing a weighted
# observed score against unweighted random scores would confound "better genes"
# with "has weights", so the primary comparison scores BOTH observed and random
# sets unweighted. The observed programs are additionally scored with their
# native weights so the canonical count can be read off the same run.
#
# LANGUAGE. Fibrosis stage and sex are recorded cross-sectional covariates.
# Nothing here orders donors or stages, and nothing here refits or reselects a
# program: the 117 memberships and weights are read, never learned.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(Matrix)
  library(matrixStats)
  library(parallel)
})

source("scripts/analysis/histology_anchored_continuum/lib_continuum.R")
source("scripts/analysis/histology_anchored_continuum/translation/lib_translation.R")
source("scripts/manuscript/program_observability_map/t3_vocabulary_lib.R")

args <- commandArgs(TRUE)
stopifnot(length(args) >= 1L)
root <- args[[1]]
n_draws <- if (length(args) >= 2L) as.integer(args[[2]]) else 200L
mc_cores <- if (length(args) >= 3L) as.integer(args[[3]]) else 4L
out <- file.path(root, "matched_null")
stopifnot(!dir.exists(out))
dir.create(out, recursive = TRUE)
NULL_SEED <- 20260906L
set.seed(NULL_SEED)
setDTthreads(mc_cores)

EVAL_COHORTS <- c("GSE162694", "GSE213621")
AXES <- c("fixed_projection", "signature_pc1")

# ---------------------------------------------------------------------------
# Substrate, built exactly as 04_score_programs_and_validate.R builds it.
# ---------------------------------------------------------------------------
message("[1/7] Expression substrate")
dge <- readRDS(file.path("RNA-seq/results/manuscript_release/candidates",
  "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION",
  "arms/F_five/results/integration/merged_dge.rds"))
meta <- fread(file.path("figures/candidates/pi-figure-redesign-2026-08-13-v3/analysis",
  "stage_extensions/five_cohort_sample_manifest.tsv"), na.strings = c("", "NA"))
annotation <- fread(file.path("RNA-seq/results/manuscript_release/candidates",
  "resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE",
  "frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz"),
  select = c("gene_id", "gene_name"))
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)[, meta$sample_id, drop = FALSE]
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm, dge); gc()

meta <- meta[dataset %in% EVAL_COHORTS]
symbol_matrix <- symbol_matrix[, meta$sample_id, drop = FALSE]

# ---------------------------------------------------------------------------
# Fixed programs and the signature exclusion, reproduced from the same sources.
# ---------------------------------------------------------------------------
message("[2/7] Fixed program memberships")
pr <- "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
registry <- fread(file.path(pr, "program_registry_v2.tsv"))
membership_raw <- fread(file.path(pr, "program_membership_v2.tsv"))
assert_true(nrow(registry) == 117L, "Program family is not 117")

signature <- fread(file.path("RNA-seq/results/histology_anchored_continuum/candidates",
  "hac-continuum-20260818T024923Z/reproduction/signature_genes.tsv"))
annotation[, gene_id_base := base_gene_id(gene_id)]
sig_sym <- unique(toupper(annotation[match(base_gene_id(signature$gene_id_base),
                                           gene_id_base), gene_name]))
sig_sym <- sig_sym[!is.na(sig_sym)]
assert_true(length(sig_sym) > 100L, "Signature symbol mapping collapsed")

membership <- membership_raw[!is.na(mapped_symbol) & mapped_symbol != "",
  .(weight = sum(as.numeric(original_l1_weight))),
  by = .(program_uid, gene_symbol = toupper(mapped_symbol))]
membership_use <- membership[!gene_symbol %in% sig_sym]

# ---------------------------------------------------------------------------
# Universe. A gene may enter a draw only where it could have entered a real
# program score: present, and non-constant in BOTH evaluation cohorts.
# Signature genes are excluded exactly as they are from program scoring.
# ---------------------------------------------------------------------------
message("[3/7] Draw universe and expression deciles")
cohort_cols <- lapply(EVAL_COHORTS, function(co) which(meta$dataset == co))
names(cohort_cols) <- EVAL_COHORTS
nonconstant <- Reduce(intersect, lapply(cohort_cols, function(j) {
  sdv <- matrixStats::rowSds(symbol_matrix[, j, drop = FALSE])
  rownames(symbol_matrix)[is.finite(sdv) & sdv > 0]
}))
universe <- sort(setdiff(nonconstant, sig_sym))
assert_true(length(universe) > 5000L, "Draw universe collapsed")
bins <- t3_expression_deciles(symbol_matrix[universe, , drop = FALSE])

observed_sets <- split(membership_use$gene_symbol, membership_use$program_uid)
observed_sets <- lapply(observed_sets, function(g) intersect(g, universe))
observed_sets <- observed_sets[registry$program_uid]
names(observed_sets) <- registry$program_uid
set_sizes <- lengths(observed_sets)
assert_true(all(set_sizes > 0L), "A program has no member gene inside the universe")
message("  universe ", length(universe), " genes; program sizes ",
        min(set_sizes), "-", max(set_sizes), " (median ", median(set_sizes), ")")

# ---------------------------------------------------------------------------
# Matched draws. Seeded per program, so one program reruns identically.
# ---------------------------------------------------------------------------
message("[4/7] ", n_draws, " matched random draws per program")
draws <- t3_matched_random_sets(observed_sets, universe, bins, n_draws, NULL_SEED)

# ---------------------------------------------------------------------------
# Scoring. Every set, observed and random, through the same unweighted rule.
# ---------------------------------------------------------------------------
message("[5/7] Scoring ", 117L * (n_draws + 1L), " gene sets")
zlist <- t3_cohort_z(symbol_matrix, meta, EVAL_COHORTS, universe)

obs_idx <- lapply(observed_sets, function(g) sort(match(g, universe)))
flat <- vector("list", 117L * (n_draws + 1L))
k <- 0L
key_prog <- character(length(flat)); key_draw <- integer(length(flat))
for (p in seq_along(obs_idx)) {
  k <- k + 1L; flat[[k]] <- obs_idx[[p]]
  key_prog[k] <- names(obs_idx)[p]; key_draw[k] <- 0L
}
for (p in seq_along(draws)) {
  for (r in seq_len(n_draws)) {
    k <- k + 1L; flat[[k]] <- draws[[p]][[r]]
    key_prog[k] <- names(draws)[p]; key_draw[k] <- r
  }
}
names(flat) <- paste0(key_prog, "|", key_draw)
scores <- t3_score_sets_fast(zlist, flat)

# Native weighted scores for the observed programs only, so the canonical count
# is reproduced in the same run rather than quoted from elsewhere.
weighted_scores <- matrix(NA_real_, nrow = 117L, ncol = ncol(scores),
                          dimnames = list(registry$program_uid, colnames(scores)))
for (co in EVAL_COHORTS) {
  j <- match(colnames(zlist[[co]]), colnames(scores))
  Z <- zlist[[co]]
  for (uid in registry$program_uid) {
    mm <- membership_use[program_uid == uid & gene_symbol %in% universe]
    w <- mm$weight / sum(mm$weight)
    weighted_scores[uid, j] <- standardize_vector(
      as.numeric(crossprod(w, Z[match(mm$gene_symbol, universe), , drop = FALSE])))
  }
}

# The prespecified signature-retention gate. One program
# (hotspot_fibroblasts_8584ba834d31e608, retained L1 0.720) falls below it and
# carries no score in the extension; blanking it here is what makes the observed
# count reproduce exactly. A random set has no signature-retention property, so
# no gate applies to the null. The asymmetry is one program in the null's
# favour, which is conservative for the observed claim.
testability <- fread(file.path("RNA-seq/results/histology_anchored_continuum/candidates",
  "hac-continuum-20260818T024923Z/programs/program_testability.tsv"))
untestable <- testability[dataset %in% EVAL_COHORTS & testable == FALSE,
                          unique(program_uid)]
message("  gating ", length(untestable), " untestable program(s) from the observed arm")
for (uid in untestable) {
  weighted_scores[uid, ] <- NA_real_
  scores[which(key_prog == uid & key_draw == 0L), ] <- NA_real_
}

# ---------------------------------------------------------------------------
# Fitting. tr_fit / tr_meta / tr_bh unchanged, on the extension's own design.
# ---------------------------------------------------------------------------
message("[6/7] Fitting")
axes <- fread(file.path("RNA-seq/results/histology_anchored_continuum/candidates",
  "hac-continuum-20260818T024923Z/projection/score_registry.tsv"))
base <- merge(meta[, .(sample_id, dataset, fibrosis_stage, inferred_sex)],
  axes[, .(sample_id, fixed_projection_raw, signature_pc1_raw)], by = "sample_id")

# One synthetic vocabulary = the 117 matched sets at draw r. BH is re-applied
# inside each vocabulary, so the family size is 117 in the null exactly as it is
# in the observed analysis.
support_for <- function(score_matrix, tag) {
  rows <- list()
  for (co in EVAL_COHORTS) {
    d0 <- base[dataset == co]
    for (axis in AXES) {
      d0[, axis_raw := get(paste0(axis, "_raw"))]
      keep <- is.finite(d0$axis_raw) & !is.na(d0$fibrosis_stage) & !is.na(d0$inferred_sex)
      dd <- d0[keep]
      dd[, axis_z := as.numeric(scale(axis_raw))]
      cols <- match(dd$sample_id, colnames(score_matrix))
      for (i in seq_len(nrow(score_matrix))) {
        z <- copy(dd); z[, outcome_z := score_matrix[i, cols]]
        z <- z[is.finite(outcome_z)]
        r <- tr_fit(z, FALSE)
        r[, `:=`(feature = rownames(score_matrix)[i], dataset = co, axis_id = axis)]
        rows[[length(rows) + 1L]] <- r
      }
    }
  }
  f <- rbindlist(rows)
  f[, q_value := tr_bh(p_value, n = .N), by = .(dataset, axis_id)]
  m <- f[, tr_meta(.SD), by = .(feature, axis_id)]
  m[, q_value := tr_bh(p_value, n = .N), by = axis_id]
  m[, supported := is.finite(q_value) & q_value < .05 & direction_consistent]
  dirs <- f[, .(all_four_directions = all(is.finite(beta)) & uniqueN(sign(beta)) == 1L,
                four_cohort_axis_BH = all(is.finite(q_value) & q_value < .05)),
            by = feature]
  s <- m[, .(both_axes_supported = all(supported) & .N == 2L), by = feature]
  s <- merge(s, dirs, by = "feature")
  s[, both_axes_supported := both_axes_supported & all_four_directions]
  s[, arm := tag]
  s[]
}

observed_unweighted <- support_for(scores[key_draw == 0L, , drop = FALSE], "observed_unweighted")
observed_weighted <- support_for(weighted_scores, "observed_weighted")

null_rows <- parallel::mclapply(seq_len(n_draws), function(r) {
  sel <- which(key_draw == r)
  sm <- scores[sel, , drop = FALSE]
  rownames(sm) <- key_prog[sel]
  s <- support_for(sm, paste0("null_", r))
  data.table(draw = r, n_supported = sum(s$both_axes_supported),
             n_each_cohort_BH = sum(s$both_axes_supported & s$four_cohort_axis_BH))
}, mc.cores = mc_cores)
bad <- vapply(null_rows, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) fail("A null draw failed: ", sum(bad), " of ", n_draws)
null_counts <- rbindlist(null_rows)

# ---------------------------------------------------------------------------
message("[7/7] Writing")
obs_n <- sum(observed_unweighted$both_axes_supported)
obs_w <- sum(observed_weighted$both_axes_supported)
obs_strict <- sum(observed_unweighted$both_axes_supported &
                    observed_unweighted$four_cohort_axis_BH)
summary <- data.table(
  metric = c("observed_supported_weighted_canonical", "observed_supported_unweighted",
             "observed_each_cohort_BH_unweighted", "null_draws",
             "null_mean_supported", "null_sd_supported", "null_median_supported",
             "null_q025_supported", "null_q975_supported", "null_max_supported",
             "n_null_draws_ge_observed", "empirical_p_supported",
             "null_mean_each_cohort_BH", "empirical_p_each_cohort_BH"),
  value = c(obs_w, obs_n, obs_strict, n_draws,
            mean(null_counts$n_supported), sd(null_counts$n_supported),
            median(null_counts$n_supported),
            quantile(null_counts$n_supported, .025, names = FALSE),
            quantile(null_counts$n_supported, .975, names = FALSE),
            max(null_counts$n_supported),
            sum(null_counts$n_supported >= obs_n),
            t3_empirical_p_ge(obs_n, null_counts$n_supported),
            mean(null_counts$n_each_cohort_BH),
            t3_empirical_p_ge(obs_strict, null_counts$n_each_cohort_BH)))
fwrite(summary, file.path(out, "matched_null_summary.tsv"), sep = "\t")
fwrite(null_counts, file.path(out, "matched_null_draw_counts.tsv"), sep = "\t")
fwrite(rbind(observed_unweighted, observed_weighted),
       file.path(out, "matched_null_observed_support.tsv"), sep = "\t")
fwrite(data.table(program_uid = names(set_sizes), n_genes_in_universe = as.integer(set_sizes)),
       file.path(out, "matched_null_set_sizes.tsv"), sep = "\t")
writeLines(c(
  paste0("seed=", NULL_SEED), paste0("n_draws=", n_draws),
  paste0("universe_genes=", length(universe)),
  paste0("untestable_programs_gated_from_observed=", length(untestable)),
  "matching=size_and_mean_expression_decile",
  "scoring_primary=unweighted_mean_of_within_cohort_gene_z",
  "scoring_secondary=native_l1_weights_observed_programs_only",
  "model=tr_fit_stage_sex_no_composition",
  "BH_family=117_reapplied_inside_each_synthetic_vocabulary",
  "support_rule=both_axes_meta_BH_q<0.05_and_all_four_cohort_axis_directions",
  "cohorts=GSE162694;GSE213621",
  "signature_genes_excluded_from_universe_and_programs=true"),
  file.path(out, "matched_null_contract.txt"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(summary)
