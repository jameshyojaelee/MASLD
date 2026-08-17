#!/usr/bin/env Rscript

# SYS-D: molecular-histologic discordance.
#
# THE QUESTION
# Which donors are molecularly advanced relative to what their recorded histology
# says, and which are the reverse? Histology is the reference standard for MASLD
# and it is also a semi-quantitative read with known inter-observer variability.
# Quantifying where it disagrees with molecular state says something the stage
# association cannot.
#
# WHY THIS CONSTRUCTION AND NOT A PREDICTOR
# This repo has already retired supervised staging twice over. A random 500-gene
# set predicted F>=3 as well as a curated one (AUROC 0.860 vs 0.881, p = 0.733),
# and the cross-cohort F-stage QWK fell from 0.638 to 0.225 once same-donor
# train/test contamination was removed. So nothing here is fit to predict stage.
# The discordance is a residual from the same mutually adjusted model the axis
# map already fits, and no ordering is constructed at any point.
#
# THE CIRCULARITY THAT HAD TO BE DESIGNED OUT
# A residual is only interpretable once projected onto a direction. Projecting a
# donor's residual onto a fibrosis effect estimated from that donor's own cohort
# would be circular. Directions are therefore leave-one-cohort-out: the direction
# applied to a GSE130970 donor is meta-analysed from the other three cohorts.
#
# CALIBRATION
# Program residuals are correlated because programs share genes, so a donor-level
# summary has no closed-form null. The null permutes fibrosis and NAS jointly
# within cohort and re-runs the entire procedure, directions included, which
# preserves both the program correlation structure and the fibrosis-NAS
# correlation. No number here is reportable without it.

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "analysis_lib.R"))

set.seed(seed)
n_null <- as.integer(Sys.getenv("MASLD_DISCORDANCE_PERMUTATIONS", "1000"))
output_root <- file.path(workstream_root, "discordance")
if (dir.exists(output_root)) fail("Refusing to overwrite: ", output_root)
if (!file.exists(nonholdout_dge)) fail("Run partition_inputs.R first")

message("[1/6] Loading sealed non-holdout input")
dge <- readRDS(nonholdout_dge)
meta <- as.data.table(readRDS(nonholdout_meta))
assert_true(!holdout_cohort %in% meta$dataset, "HOLDOUT LEAK")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
annotation <- fread(gene_annotation, select = c("gene_id", "gene_name"))
symbol_matrix <- collapse_symbols(logcpm, rownames(logcpm), annotation)
rm(logcpm, dge); gc()

registry <- fread(program_registry)
membership <- fread(program_membership)
program_scores <- score_programs(symbol_matrix, meta, membership, registry,
                                 program_weight_coverage)
testable <- program_scores$coverage[testable == TRUE, feature_id]
S <- program_scores$primary[testable, , drop = FALSE]
rm(symbol_matrix); gc()

discovery_meta <- meta[dataset %in% discovery_cohorts &
                         !is.na(fibrosis_stage) & !is.na(nas_score)]
assert_true(nrow(discovery_meta) == 469L, "Discovery census drift")

# Biological-unit control, stated honestly. The sealed metadata carries no
# participant column, so participant_id falls back to sample_id exactly as the
# v8 driver does. That makes a "one sample per participant" assertion trivially
# true and it is NOT evidence against pseudoreplication. The substantive
# evidence is external: no donor-pairing table exists for any of the four
# discovery cohorts, whereas one exists for every dataset that does have repeat
# sampling (GSE244832, GSE202379, GSE185477, GSE136103 on the single-cell side,
# and GSE193066 via its participant crosswalk on the bulk side). GSE193066 is
# the only bulk cohort with repeat biopsies and it is the sealed holdout.
if (!"participant_id" %in% names(meta)) {
  message("      note: no participant column in sealed metadata; ",
          "falling back to sample_id (see comment above)")
  meta[, participant_id := sample_id]
  discovery_meta[, participant_id := sample_id]
}
pairing_files <- file.path(project_root, "data", discovery_cohorts,
                           "metadata", "donor_pairing.csv")
assert_true(!any(file.exists(pairing_files)),
            paste0("A discovery cohort has a donor-pairing table, so repeat ",
                   "sampling is possible and the residual model needs a donor term: ",
                   paste(pairing_files[file.exists(pairing_files)], collapse = ", ")))
holdout_participants <- fread(gse193066_crosswalk)[, unique(participant_token)]
assert_true(!any(discovery_meta$sample_id %in% holdout_participants),
            "HOLDOUT LEAK: a discovery sample matches a GSE193066 participant token")

# ------------------------------------------------------------------ procedure

# One pass over the four cohorts, matrix form. Returns per-cohort studentized
# residuals plus the per-cohort fibrosis coefficient and its HC3 standard error,
# which is everything the leave-one-cohort-out direction needs.
# The residual model must be SATURATED in stage, and the first version of this
# script was not. With a linear fibrosis term the non-linear part of the stage
# response stays in the residual, and because the permutation null destroys that
# structure the leftover trend scores as significant. The symptoms were
# unmistakable: every "advanced" donor sat at F0-F2 and every "early" donor at
# F2-F4, the discordance rate by stage was U-shaped (5.9, 1.6, 2.7, 5.3, 16.1
# percent), mean z trended monotonically with stage, and sd(z) was 1.29 rather
# than 1. The axis map had already shown 23 fibrosis programs whose shape a
# linear term misses, so this was predictable in hindsight.
#
# Fibrosis therefore enters as a full 5-level factor (every level is populated in
# every discovery cohort; the smallest cell is n = 2). NAS is coarsened to the
# standard not-NASH / borderline / NASH cut because its raw 0-8 levels are too
# thin to saturate at n = 76. Discordance now means deviating from donors with
# the SAME recorded histology, which is the question worth asking.
#
# Note the categorical block is used only to remove structure, never tested, so
# the anticonservative HC3 categorical Wald behaviour measured elsewhere in this
# workstream does not apply here.
nas_bin <- function(x) cut(x, breaks = c(-Inf, 2, 4, Inf),
                           labels = c("le2", "3to4", "ge5"), right = TRUE)

cohort_fits <- function(md) {
  out <- list()
  for (cohort in discovery_cohorts) {
    idx <- which(md$dataset == cohort)
    d <- copy(md[idx])
    d[, `:=`(fibrosis_factor = factor(fibrosis_stage),
             nas_factor = nas_bin(nas_score),
             fibrosis_centered = fibrosis_stage - mean(fibrosis_stage),
             sex_final = factor(sex_final))]
    d[, (c("fibrosis_factor", "nas_factor", "sex_final")) :=
        lapply(.SD, droplevels), .SDcols = c("fibrosis_factor", "nas_factor", "sex_final")]
    X_resid <- stats::model.matrix(~ sex_final + fibrosis_factor + nas_factor, data = d)
    Y <- S[, d$sample_id, drop = FALSE]
    fit_resid <- hc3_fit_matrix(Y, X_resid)
    if (!isTRUE(fit_resid$estimable)) {
      fail("Saturated residual design not estimable in ", cohort, ": ",
           fit_resid$failure_reason)
    }
    # The projection direction stays the linear fibrosis effect: it is the axis
    # being projected onto, and it is estimated leave-one-cohort-out.
    X_dir <- stats::model.matrix(
      ~ sex_final + fibrosis_centered + nas_factor, data = d)
    fit_dir <- hc3_fit_matrix(Y, X_dir)
    if (!isTRUE(fit_dir$estimable)) fail("Direction design not estimable in ", cohort)
    out[[cohort]] <- list(
      residuals = fit_resid$studentized_residuals,
      beta = fit_dir$coefficients["fibrosis_centered", ],
      se = fit_dir$se["fibrosis_centered", ],
      samples = d$sample_id,
      stage = d$fibrosis_stage
    )
  }
  out
}

# Fixed-effect inverse-variance direction from the three training cohorts. A
# direction only needs a sign and a magnitude, so this deliberately avoids the
# REML fit: with k = 3 the between-cohort variance is barely identified, and the
# permutation null is re-run through this same estimator anyway.
loco_direction <- function(fits, held) {
  train <- setdiff(names(fits), held)
  b <- do.call(cbind, lapply(train, function(k) fits[[k]]$beta))
  s <- do.call(cbind, lapply(train, function(k) fits[[k]]$se))
  w <- 1 / s^2
  w[!is.finite(w)] <- 0
  ok <- rowSums(w) > 0
  out <- rep(NA_real_, nrow(b))
  names(out) <- rownames(b)
  out[ok] <- rowSums((b * w)[ok, , drop = FALSE]) / rowSums(w[ok, , drop = FALSE])
  out
}

discordance_once <- function(md) {
  fits <- cohort_fits(md)
  rows <- lapply(discovery_cohorts, function(held) {
    dirs <- loco_direction(fits, held)
    dirs <- dirs[is.finite(dirs)]
    if (!length(dirs)) return(NULL)
    d <- donor_discordance(fits[[held]]$residuals, dirs, weights = abs(dirs))
    d[, cohort := held][]
  })
  rbindlist(Filter(Negate(is.null), rows))
}

message("[2/6] Observed discordance")
observed <- discordance_once(discovery_meta)
observed <- merge(observed,
                  discovery_meta[, .(sample_id, participant_id, dataset,
                                     fibrosis_stage, nas_score, sex_final)],
                  by = "sample_id", sort = FALSE)

message("[3/6] ", n_null, " joint within-cohort permutations")
null_mat <- matrix(NA_real_, nrow = n_null, ncol = nrow(observed),
                   dimnames = list(NULL, observed$sample_id))
for (b in seq_len(n_null)) {
  md <- copy(discovery_meta)
  md[, .perm := sample.int(.N), by = dataset]
  md[, `:=`(fibrosis_stage = fibrosis_stage[.perm], nas_score = nas_score[.perm]),
     by = dataset]
  md[, .perm := NULL]
  nb <- discordance_once(md)
  null_mat[b, nb$sample_id] <- nb$discordance
  if (b %% 100L == 0L) message("      ", b, " / ", n_null)
}

message("[4/6] Calibrating")
observed[, `:=`(
  null_sd = apply(null_mat[, sample_id, drop = FALSE], 2L, stats::sd, na.rm = TRUE),
  null_mean = colMeans(null_mat[, sample_id, drop = FALSE], na.rm = TRUE)
)]
observed[, z_discordance := (discordance - null_mean) / null_sd]
observed[, empirical_p := vapply(seq_len(.N), function(i) {
  nulls <- null_mat[, sample_id[i]]
  nulls <- nulls[is.finite(nulls)]
  if (!length(nulls)) return(NA_real_)
  (1 + sum(abs(nulls - null_mean[i]) >= abs(discordance[i] - null_mean[i]))) /
    (length(nulls) + 1)
}, numeric(1))]
observed[, q_value := p.adjust(empirical_p, method = "BH")]
observed[, discordance_class := fcase(
  q_value < 0.05 & z_discordance > 0, "molecularly_advanced",
  q_value < 0.05 & z_discordance < 0, "molecularly_early",
  default = "concordant"
)]

# Misspecification guard. Under a residual model saturated in stage there must be
# no systematic trend of discordance with stage: any monotone drift means the
# model is leaving stage structure in the residual, which is precisely the
# artifact that invalidated the first version of this analysis. Spearman rho of
# z against stage is reported and gated rather than left for a reader to notice.
stage_trend <- observed[, .(mean_z = mean(z_discordance), n = .N), by = fibrosis_stage][
  order(fibrosis_stage)]
trend_rho <- suppressWarnings(
  stats::cor(observed$z_discordance, observed$fibrosis_stage, method = "spearman"))
trend_p <- suppressWarnings(
  stats::cor.test(observed$z_discordance, observed$fibrosis_stage,
                  method = "spearman")$p.value)
sign_stage_segregation <- observed[discordance_class != "concordant", {
  adv <- fibrosis_stage[discordance_class == "molecularly_advanced"]
  early <- fibrosis_stage[discordance_class == "molecularly_early"]
  .(separable = length(adv) > 0 && length(early) > 0 && max(adv) <= min(early))
}]$separable
message("      stage trend in discordance: rho = ", signif(trend_rho, 3),
        ", p = ", signif(trend_p, 3),
        " | sd(z) = ", signif(stats::sd(observed$z_discordance), 3))

# Procedure-level calibration: how many donors would be called discordant under
# the null? A calibrated procedure returns close to zero after BH.
null_calls <- vapply(seq_len(min(n_null, 200L)), function(b) {
  z <- (null_mat[b, ] - observed$null_mean) / observed$null_sd
  p <- 2 * stats::pnorm(-abs(z))
  sum(p.adjust(p, method = "BH") < 0.05, na.rm = TRUE)
}, numeric(1))

# PRIMARY STATISTIC: global excess dispersion.
#
# Per-donor calls are resolution-limited. The smallest achievable empirical p is
# 1/(B+1), while BH across 469 donors needs 0.05/469 = 1.07e-4 for the most
# extreme donor, so B must exceed about 470 / 0.05 before any donor can be called
# at all. At B = 1,000 the answer was structurally zero and said nothing.
#
# The dispersion of z across donors has no such ceiling: if donors deviate from
# their same-histology peers more than the design's own null permits, sd(z)
# exceeds 1. That is the question SYS-D exists to answer, it is well powered at
# this n, and it does not require naming an individual donor.
observed_sd <- stats::sd(observed$z_discordance, na.rm = TRUE)
null_sd_dist <- apply(null_mat, 1L, function(r) {
  stats::sd((r - observed$null_mean) / observed$null_sd, na.rm = TRUE)
})
null_sd_dist <- null_sd_dist[is.finite(null_sd_dist)]
dispersion_p <- (1 + sum(null_sd_dist >= observed_sd)) / (length(null_sd_dist) + 1)
per_donor_resolvable <- (1 / (n_null + 1)) < (0.05 / nrow(observed))
message("      global dispersion: observed sd(z) ", signif(observed_sd, 4),
        " vs null mean ", signif(mean(null_sd_dist), 4),
        ", p = ", signif(dispersion_p, 3),
        " | per-donor calls resolvable: ", per_donor_resolvable)

message("[5/6] Writing")
tmp <- atomic_dir(output_root)
on.exit(if (dir.exists(tmp)) unlink(tmp, recursive = TRUE), add = TRUE)
write_tsv(observed, file.path(tmp, "donor_discordance.tsv"))
write_tsv(stage_trend, file.path(tmp, "stage_trend_diagnostic.tsv"))
write_tsv(data.table(
  metric = c("n_donors", "n_programs_used", "n_permutations",
             "n_molecularly_advanced", "n_molecularly_early", "n_concordant",
             "null_calls_mean", "null_calls_max",
             "stage_trend_rho", "stage_trend_p", "sd_z",
             "sign_stage_segregation", "misspecification_suspected",
             "dispersion_observed_sd", "dispersion_null_mean_sd",
             "dispersion_null_p95_sd", "dispersion_empirical_p",
             "per_donor_resolvable", "empirical_p_floor", "bh_threshold_needed"),
  value = c(nrow(observed), observed$n_programs_used[[1L]], n_null,
            sum(observed$discordance_class == "molecularly_advanced"),
            sum(observed$discordance_class == "molecularly_early"),
            sum(observed$discordance_class == "concordant"),
            mean(null_calls), max(null_calls),
            trend_rho, trend_p, stats::sd(observed$z_discordance),
            as.numeric(isTRUE(sign_stage_segregation)),
            as.numeric(isTRUE(trend_p < 0.05) || isTRUE(sign_stage_segregation)),
            observed_sd, mean(null_sd_dist),
            as.numeric(stats::quantile(null_sd_dist, 0.95)), dispersion_p,
            as.numeric(per_donor_resolvable), 1 / (n_null + 1),
            0.05 / nrow(observed))
), file.path(tmp, "discordance_summary.tsv"))
saveRDS(null_mat, file.path(tmp, "discordance_null.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(tmp, "sessionInfo.txt"))
inputs <- c(nonholdout_dge, nonholdout_meta, program_registry, program_membership)
write_tsv(data.table(path = inputs, sha256 = vapply(inputs, sha256_file, character(1))),
          file.path(tmp, "input_manifest.tsv"))
jsonlite::write_json(list(
  state = "DISCORDANCE_COMPLETE_HOLDOUT_UNOPENED",
  workstream_id = workstream_id, n_donors = nrow(observed),
  n_permutations = n_null,
  null_false_call_rate_mean = mean(null_calls),
  holdout_accessed = FALSE
), file.path(tmp, "DISCORDANCE_READY.json"), pretty = TRUE, auto_unbox = TRUE)

message("[6/6] Publishing")
publish_dir(tmp, output_root)
Sys.chmod(list.files(output_root, full.names = TRUE), mode = "0440")

cat("\n=== DISCORDANCE ===\n")
print(observed[, .N, by = discordance_class])
cat("\nnull false-call rate: mean ", round(mean(null_calls), 2),
    ", max ", max(null_calls), " of ", nrow(observed), " donors\n", sep = "")
cat("SYS_D_COMPLETE\t", output_root, "\n", sep = "")
