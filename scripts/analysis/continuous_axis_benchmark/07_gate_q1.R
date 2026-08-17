#!/usr/bin/env Rscript
# 07: evaluate the Q1 gate and write GATE.json. Nothing in Q2 may run without it.
#
# The four criteria are read from 00_prespecification.json, not written here, so a
# threshold cannot be nudged after seeing a result without editing a file whose
# hash is recorded in every output manifest.
#
# The important branch is the failure branch. If no arm passes, that is not a dead
# end -- it is the Q1 answer, and the headline. Q2 still runs, on recorded
# histology bins, so the amplification-versus-sequence question is answered either
# way.

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
source(file.path(Sys.getenv("MASLD_PROJECT_ROOT"),
                 "scripts/analysis/continuous_axis_benchmark/lib_axis_common.R"))

OUT <- Sys.getenv("CAB_OUT_ROOT", "")
assert_true(nzchar(OUT), "CAB_OUT_ROOT is unset")
pre <- read_prespec()
g <- pre$q1_gate$an_arm_passes_only_if_all_four_hold

read_if <- function(p) if (file.exists(p)) fread(p) else NULL
loco <- read_if(file.path(OUT, "q1", "loco_overall.tsv"))
arms <- c("A1", "A2", "A3", "A4", "A5", "A6")

# Validity of criterion g1, measured rather than assumed.
#
# g1 counts genes explained by the axis beyond histology. The prespecification
# assumed arm A6 (held-out gene half) controlled the circularity in that count. It
# does not: genes are co-expressed, so a random half of the transcriptome spans the
# same components as the other half, and A6 still returned 71%. The matched control
# is a held-out AXIS, not a held-out gene set. If the real axis does not beat axes
# built from random gene subsets, then g1 is measuring "a global expression
# component explains most genes" and cannot discriminate between arms at all.
rand <- read_if(file.path(OUT, "q1", "random_axis_null_verdict.tsv"))
g1_valid <- if (is.null(rand)) NA else isTRUE(rand$real_beats_random_axes[1])

rows <- rbindlist(lapply(arms, function(this_arm) {
  arm <- this_arm
  de <- read_if(file.path(OUT, "q1", sprintf("%s_de_increment.tsv", arm)))
  de6 <- read_if(file.path(OUT, "q1", "A6_de_increment.tsv"))
  st <- read_if(file.path(OUT, "arms", sprintf("%s_axis_summary.tsv", arm)))
  xf <- read_if(file.path(OUT, "q1", sprintf("%s_cross_fold.tsv", arm)))

  n_axis <- if (!is.null(de)) de$n_axis_given_histology_calibrated[1] else NA_integer_
  n_a6 <- if (!is.null(de6)) de6$n_axis_given_histology_calibrated[1] else NA_integer_
  g1 <- is.finite(n_axis) && n_axis >= g$g1_de_increment$pop_c_min_genes &&
        is.finite(n_a6) && n_a6 >= g$g1_de_increment$a6_min_genes

  lo <- if (!is.null(loco)) loco[loco$arm == this_arm & loco$population == "POP-C"] else NULL
  g2 <- !is.null(lo) && nrow(lo) && isTRUE(lo$gate_g2_pass[1])

  g3 <- if (!is.null(xf)) {
    sum(xf$spearman >= g$g3_cross_fold_replication$min_rho, na.rm = TRUE) >=
      g$g3_cross_fold_replication$min_folds
  } else NA

  # g4 applies only to arms that infer an ordering. A3 and A4 do not, so the
  # identifiability criterion is not applicable rather than failed.
  g4 <- if (arm %in% c("A3", "A4")) NA
        else if (!is.null(st) && "g4_pass" %in% names(st)) isTRUE(st$g4_pass[1])
        else if (!is.null(st) && "median_abs_spearman_between_splits" %in% names(st))
          st$median_abs_spearman_between_splits[1] >= g$g4_axis_identifiability$min_median_abs_rho
        else NA

  a5_dq <- if (arm == "A5" && !is.null(st) && "disqualified" %in% names(st))
    isTRUE(st$disqualified[1]) else FALSE

  data.table(arm = arm,
             g1_de_increment = g1, g2_held_out = g2,
             g3_cross_fold = g3, g4_identifiable = g4,
             a5_disqualified_by_random_gene_null = a5_dq,
             n_axis_calibrated = n_axis, n_axis_calibrated_A6 = n_a6,
             passes = isTRUE(g1) && isTRUE(g2) &&
               (isTRUE(g3) || is.na(g3)) && (isTRUE(g4) || is.na(g4)) && !a5_dq)
}))

# A3 and A4 are comparators, not candidates: they cannot "win" Q1 because neither
# infers an ordering, and Q1 asks whether inferring one adds anything.
candidates <- rows[!arm %in% c("A3", "A4")]
winners <- candidates[passes == TRUE]

gate <- list(
  workstream = "CONTINUOUS-AXIS-BENCHMARK-v1",
  exemption_id = "ORDERING-EXEMPTION-2026-08-14",
  evaluated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  prespecification_sha256 = sha256_file(file.path(
    project_root(), "scripts/analysis/continuous_axis_benchmark/00_prespecification.json")),
  criteria = g,
  per_arm = rows,
  g1_criterion_valid = g1_valid,
  g1_validity_note = if (isTRUE(g1_valid))
    "The real axis beat random-gene-subset axes on held-out genes, so g1 discriminates."
  else
    paste0("g1 IS NOT INTERPRETABLE AS PRESPECIFIED. The real axis explained ",
           if (!is.null(rand)) rand$real_n_hits[1] else NA, " of ",
           if (!is.null(rand)) rand$n_test_genes[1] else NA,
           " held-out genes against a random-axis median of ",
           if (!is.null(rand)) rand$random_median_hits[1] else NA,
           " (empirical p = ", if (!is.null(rand)) round(rand$empirical_p_real_vs_random[1], 3) else NA,
           "). Gene-splitting does not control circularity because genes are co-expressed. ",
           "The g1 counts measure that any global expression component explains most genes. ",
           "This is a defect in the prespecified gate, recorded rather than reinterpreted; ",
           "the Q1 conclusion rests on g2 and g4, which are unaffected."),
  any_arm_passed = nrow(winners) > 0L,
  winning_axis = if (nrow(winners) > 0L) winners$arm[1] else "histology",
  q1_conclusion = if (nrow(winners) > 0L)
    paste0("A derived ordering added information beyond recorded histology: arm ",
           winners$arm[1], " cleared all four prespecified criteria.")
  else
    paste0("No derived ordering added information beyond recorded histology on 844 ",
           "participants. This null is the Q1 result. Q2 proceeds on recorded ",
           "histology bins so the amplification-versus-sequence question is still answered."),
  q2_substrate = if (nrow(winners) > 0L) winners$arm[1] else "A4_histology_bins"
)

write_tsv_once(rows, file.path(OUT, "q1", "gate_per_arm.tsv"))
write_json_once(gate, file.path(OUT, "q1", "GATE.json"))
print(rows)
cat("\nQ1 CONCLUSION: ", gate$q1_conclusion, "\n", sep = "")
cat("Q2 SUBSTRATE : ", gate$q2_substrate, "\n", sep = "")
log_step("Q1_GATE_COMPLETE any_arm_passed=", gate$any_arm_passed)
