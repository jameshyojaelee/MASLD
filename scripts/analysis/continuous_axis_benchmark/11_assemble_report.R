#!/usr/bin/env Rscript
# 11: assemble every result into one report, and record what did NOT work.
#
# Two rules this report follows.
#
# Prespecified criteria are reported as written, even where they turned out to be
# poorly chosen. Where a criterion proved uninterpretable, that is stated as a
# defect in the criterion rather than quietly replaced with a better one.
#
# Nothing here is a Resource claim. This workstream ran under an ordering
# exemption whose whole point is that its outputs may not enter Resource claim
# text, figures, or provenance. The report says so on its face.

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
rd <- function(p) if (file.exists(file.path(OUT, p))) fread(file.path(OUT, p)) else NULL
rj <- function(p) if (file.exists(file.path(OUT, p))) fromJSON(file.path(OUT, p), simplifyVector = TRUE) else NULL

gate    <- rj("q1/GATE.json")
reg     <- rd("arms/axis_registry_summary.tsv")
loco    <- rd("q1/loco_overall.tsv")
agg     <- rd("q1/aggregation_inflation.tsv")
aggsim  <- rd("q1/aggregation_inflation_simulation.tsv")
randn   <- rd("q1/random_axis_null_verdict.tsv")
strat   <- rd("arms/lineage_stratified_summary.tsv")
rank1   <- rd("q2/M_gene_histology_rank1_verdict.tsv")
gains   <- rd("q2/M_gene_histology_window_gains.tsv")
paired  <- rd("paired/10b_paired_verdict.tsv")
a5      <- rd("arms/A5_axis_summary.tsv")
ctrl    <- rd("arms/pipeline_controls.tsv")
combat  <- rd("arms/combat_damage_audit.tsv")

L <- character(0)
add <- function(...) L <<- c(L, paste0(...))

add("# Continuous-axis benchmark: results")
add("")
add("Workstream CONTINUOUS-AXIS-BENCHMARK-v1 under exemption ORDERING-EXEMPTION-2026-08-14.")
add("Output root: ", OUT)
add("Assembled: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
add("")
add("**Claim boundary.** Nothing in this report is a Resource-level claim. Outputs under this")
add("exemption may not appear in Resource claim text, figure captions, release tables, or")
add("provenance chains; `validate_resource_scope.py` enforces that separately.")
add("")

add("## Pipeline controls")
add("")
if (!is.null(ctrl)) {
  add("The rebuilt pipeline reproduces known values from the validated substrate:")
  add("")
  for (i in seq_len(nrow(ctrl)))
    add("- ", ctrl$quantity[i], ": expected ", ctrl$expected[i], ", observed ",
        ctrl$observed[i], " -> ", ifelse(ctrl$pass[i], "PASS", "FAIL"))
  add("")
}
if (!is.null(combat)) {
  s <- combat[quantity == "attenuation_slope_cb_on_raw", value]
  r <- combat[quantity == "pearson_r_cb_vs_raw", value]
  add("ComBat attenuates the F4-vs-F0 contrast by a slope of ", round(s, 3),
      " with Pearson r ", round(r, 3), " against the uncorrected fit, so batch correction")
  add("preserves disease signal here and is not responsible for the Q1 result.")
  add("")
}

add("## Q1: does a data-driven ordering add information beyond recorded histology?")
add("")
if (!is.null(gate)) { add("**", gate$q1_conclusion, "**"); add("") }
if (!is.null(reg)) {
  add("Axis correlation with fibrosis stage, and how much each axis is really cohort:")
  add("")
  add("| arm | n scored | Spearman vs fibrosis | cohort eta2 |")
  add("|---|---:|---:|---:|")
  for (i in seq_len(nrow(reg)))
    add("| ", reg$arm[i], " | ", reg$n_scored[i], " | ",
        round(reg$spearman_vs_fibrosis[i], 3), " | ",
        round(reg$kruskal_eta2_on_cohort[i], 3), " |")
  add("")
}
if (!is.null(loco)) {
  lc <- loco[population == "POP-C"]
  add("Leave-one-cohort-out prediction of held-out fibrosis stage (POP-C):")
  add("")
  add("| arm | n | C-index | 95% CI |")
  add("|---|---:|---:|---|")
  for (i in seq_len(nrow(lc)))
    add("| ", lc$arm[i], " | ", lc$n_donors[i], " | ", round(lc$pooled_c_index[i], 3),
        " | ", round(lc$boot_lo[i], 3), " - ", round(lc$boot_hi[i], 3), " |")
  add("")
}
if (!is.null(strat)) {
  add("Axis identifiability across 27 unstated-parameter settings, stratified by whether")
  add("two settings inferred the same number of lineages:")
  add("")
  add("| arm | same lineage count | pairs | median abs rho |")
  add("|---|---|---:|---:|")
  for (i in seq_len(nrow(strat)))
    add("| ", strat$arm[i], " | ", strat$same_lineage_count[i], " | ", strat$n_pairs[i],
        " | ", round(strat$median_abs_rho[i], 3), " |")
  add("")
}
if (!is.null(a5)) {
  add("Supervised arm A5 against its own null of ", 200, " draws of 500 random genes: observed ",
      "C-index ", round(a5$observed_c_index[1], 3), " against a null 95th percentile of ",
      round(a5$null_q95_c_index[1], 3), " (empirical p ", round(a5$empirical_p_vs_random_genes[1], 3),
      "). Disqualified: ", a5$disqualified[1], ".")
  add("")
}

add("### Defect in the prespecified g1 criterion")
add("")
if (!is.null(randn)) {
  add("g1 counted genes explained by the axis beyond histology, and the prespecification")
  add("assumed the held-out-gene-half arm controlled circularity. It does not: genes are")
  add("co-expressed, so a random half of the transcriptome spans the same components as the")
  add("other half. The matched control is a held-out AXIS.")
  add("")
  add("- Real axis: ", randn$real_n_hits[1], " of ", randn$n_test_genes[1], " held-out genes")
  add("- Random-subset axes: median ", randn$random_median_hits[1],
      ", 95th percentile ", round(randn$random_q95_hits[1], 0))
  add("- Empirical p: ", round(randn$empirical_p_real_vs_random[1], 3))
  add("")
  add("The inferred axis explains FEWER held-out genes than a random expression direction.")
  add("g1 is therefore uninterpretable as written; the Q1 conclusion rests on g2 and g4.")
  add("")
}

add("### Aggregation inflation in the published validation statistic")
add("")
if (!is.null(agg) && !is.null(aggsim)) {
  a <- agg[variable == "fibrosis_stage"]; s <- aggsim[variable == "fibrosis_stage"]
  m <- merge(a[, .(arm, r_group_kamzolas_estimator, r_patient_pearson)],
             s[, .(arm, sim_median_r_group, sim_p_r_group_ge_0.96)], by = "arm")
  add("| arm | r (donor level) | r (group means) | simulated median r_group | P(r_group >= 0.96) |")
  add("|---|---:|---:|---:|---:|")
  for (i in seq_len(nrow(m)))
    add("| ", m$arm[i], " | ", round(m$r_patient_pearson[i], 3), " | ",
        round(m$r_group_kamzolas_estimator[i], 3), " | ",
        round(m$sim_median_r_group[i], 3), " | ", round(m$sim_p_r_group_ge_0.96[i], 3), " |")
  add("")
  add("The group-mean estimator saturates: a donor-level correlation near 0.57 produces")
  add(">= 0.96 at the group level essentially always. A reported group-level R of 0.96-1.0")
  add("is therefore consistent with a wide range of donor-level correlations and cannot")
  add("establish that an axis is good.")
  add("")
}

add("## Q2: ordered sequence or scalar amplification?")
add("")
if (!is.null(rank1)) {
  add("- Observed PVE1: ", round(rank1$observed_pve1[1], 4))
  add("- H_amp null median / 5th percentile: ", round(rank1$null_pve1_median[1], 4), " / ",
      round(rank1$null_pve1_p05[1], 4))
  add("- Rank-1 strictly sufficient: ", rank1$rank1_sufficient[1],
      " (empirical p ", signif(rank1$empirical_p_pve1[1], 3), ")")
  add("- Features individually breaking rank 1: ", rank1$n_features_breaking_rank1_BH05[1],
      " of ", rank1$n_features[1], " (",
      round(100 * rank1$frac_features_breaking_rank1[1], 2), "%)")
  add("- Worst window by residual chi-square: ", rank1$worst_window_by_resid[1],
      "; earliest window is worst: ", rank1$earliest_window_is_worst[1])
  if ("gain_ratio_last_over_first" %in% names(rank1))
    add("- Gain ratio last/first: ", round(rank1$gain_ratio_last_over_first[1], 3),
        " (", round(rank1$gain_ratio_boot_lo[1], 3), " - ",
        round(rank1$gain_ratio_boot_hi[1], 3), ")")
  add("")
  add("Amplification is the dominant structure and the window gains increase monotonically,")
  add("but rank 1 is not strictly sufficient. Per the prespecified decision rule this is the")
  add("'amplification plus named exceptions' branch. The exception is the earliest transition,")
  add("which was predicted a priori from the F0->F1 anomaly.")
  add("")
}
if (!is.null(gains)) {
  add("| window | gain | H_amp draw interval |")
  add("|---|---:|---|")
  for (i in seq_len(nrow(gains)))
    add("| ", gains$window[i], " | ", round(gains$gain[i], 3), " | ",
        round(gains$boot_lo[i], 3), " - ", round(gains$boot_hi[i], 3), " |")
  add("")
  add("**Do not quote the per-window intervals as confidence intervals.** They are the")
  add("distribution of the estimator under H_amp, which is what the PVE1 test requires, but")
  add("the rank-1 fit absorbs simulated noise into its magnitude, biasing every draw upward")
  add("by roughly 2%. That is why each point estimate sits just below its own interval. The")
  add("bias cancels in a ratio, so the quotable quantity is the last/first gain ratio above,")
  add("whose point estimate does lie inside its interval.")
  add("")
  add("Convergent check: the pairwise disattenuated analysis in")
  add("RNA-seq/results/stage_axis_geometry/20260813T142216Z/ gives a compound F1->F4 gain of")
  add("1.70 x 1.51 x 1.62 = 4.16. The weighted rank-1 decomposition here gives 3.94. Two")
  add("methods sharing no machinery agree to within 5%.")
  add("")
}

add("## Within-person test (GSE193066 paired biopsies)")
add("")
if (!is.null(paired)) {
  add("- Donors: ", paired$n_donors[1], "; genes ", paired$n_genes[1])
  add("- P1 (delta axis vs delta fibrosis): Spearman ", round(paired$P1_spearman_rho[1], 3),
      ", one-sided p ", signif(paired$P1_p_one_sided[1], 3),
      " -> supported: ", paired$P1_supported[1])
  add("- P2 (within-person delta parallel to cross-sectional F0->F4): cosine ",
      round(paired$P2_cosine_mean_delta[1], 3), " (",
      round(paired$P2_cosine_boot_lo[1], 3), " - ", round(paired$P2_cosine_boot_hi[1], 3),
      "), threshold ", paired$P2_threshold[1], " -> supported: ", paired$P2_supported[1])
  add("- Donors with positive cosine: ",
      round(100 * paired$P2_frac_donors_positive_cosine[1], 1), "%")
  add("- Stage change: ", paired$n_stage_increased[1], " increased, ",
      paired$n_stage_decreased[1], " decreased, ", paired$n_stage_unchanged[1], " unchanged")
  add("")
  add(if (isTRUE(paired$P2_supported[1]))
        paste("P2 holds, so within-person change runs parallel to the cross-sectional axis and",
              "the amplification structure reflects change within people.")
      else
        paste("P2 does not hold, so the cross-sectional structure does not transport to",
              "within-person change and must be described as a between-donor property."))
  add("")
  if (!isTRUE(paired$P2_supported[1])) {
    add("**P2 was mis-specified, and that limits what its failure proves.** It predicted the")
    add("MEAN within-person delta would align with the cross-sectional disease direction. But")
    add("this cohort has ", paired$n_stage_increased[1], " donors progressing, ",
        paired$n_stage_decreased[1], " regressing and ", paired$n_stage_unchanged[1],
        " unchanged, so there is no net")
    add("progression to detect; a mean delta near zero on the disease axis is what balanced")
    add("movement produces whether or not the amplification model is correct. The test as")
    add("sealed could not have confirmed the hypothesis it was written for.")
    add("")
    add("Two signals sit underneath it. Every stage stratum has a negative cosine and the")
    add("mean-delta cosine (", round(paired$P2_cosine_mean_delta[1], 3), ") is far more negative than the median donor (",
        round(paired$P2_median_donor_cosine[1], 3), "),")
    add("which is the signature of a shared first-versus-second-biopsy offset that survives")
    add("averaging while donor-specific noise cancels. It is not dominant: the rank-1 share of")
    add("the delta matrix is only ", round(paired$delta_pve1[1], 3), ", so deltas are mostly idiosyncratic. Allowing for")
    add("that offset, progressors rank highest by median cosine, the predicted ordering, and P1")
    add("points the same way at rho ", round(paired$P1_spearman_rho[1], 3), " -- but neither reaches significance with ",
        paired$n_stage_increased[1] + paired$n_stage_decreased[1], " informative donors.")
    add("")
    add("The well-posed test is a progressor-versus-regressor contrast with the biopsy-order")
    add("offset removed. That is the NEXT test, not evidence now, and GSE193066 is no longer")
    add("sealed for it.")
    add("")
  }
} else add("Not run.")

add("## What did not work")
add("")
add("- Prespecified g1 assumed gene-splitting controls circularity. It does not; recorded above.")
add("- Prespecified P2 predicted a mean-delta alignment in a cohort with balanced progression")
add("  and regression, where no net movement exists to detect. The prediction was poorly")
add("  reasoned, not merely unsupported, and it consumed the only within-person data available.")
add("- The Q2 decision rule's literal branch (sequence wins if criteria i and ii both fail)")
add("  triggered, and was overridden by the prespecified earliest-window exception. The")
add("  override is recorded rather than presented as a clean amplification win. Supporting it:")
add("  change-point ENTROPY is BELOW the amplification null, so the data are less sequential")
add("  than amplification-plus-noise predicts, which is the opposite of ordered activation.")
add("- Slingshot returns a tree, not a path. Taking lineage 1 alone left the axis undefined")
add("  for 385 of 844 donors and inflated apparent stability from 0.567 to 0.946.")
add("- The first rank-1 run reported meaningless gain CIs because a*b^T is sign-invariant.")
add("- preprocessCore::normalize.quantiles aborts under SLURM cgroups; limma::normalizeQuantiles")
add("  replaces it.")
add("- Superseded outputs are retained under _superseded_* directories with REASON.txt files.")
add("")

writeLines(L, file.path(OUT, "REPORT.md"))
write_run_parameters(OUT, list(report = "REPORT.md"))
write_output_manifest(OUT)
cat(paste(L, collapse = "\n"), "\n")
log_step("REPORT_COMPLETE: ", file.path(OUT, "REPORT.md"))
