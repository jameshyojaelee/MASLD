#!/usr/bin/env Rscript

# Key message: this candidate tests whether externally oriented continuous scores
# add within-stage molecular resolution; it never changes Figure 3 automatically.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- dirname(normalizePath(script_arg))
source(file.path(script_dir, "lib_continuum.R"))

pre <- read_prespec()
out <- out_root()
gate_dir <- file.path(out, "decision")

paths <- c(
  prepared_input_manifest = file.path(out, "inputs", "input_manifest.tsv"),
  reproduction_summary = file.path(out, "reproduction", "reproduction_summary.tsv"),
  signature_coverage = file.path(out, "unsupervised", "coverage.tsv"),
  histology_anchors = file.path(out, "programs", "axis_histology_anchors.tsv"),
  co_primary_agreement = file.path(out, "programs", "co_primary_agreement.tsv"),
  focal_program_gate = file.path(gate_dir, "focal_program_gate.tsv"),
  paired_tests = file.path(out, "paired", "paired_tests.tsv"),
  stable_visit_offset = file.path(out, "paired", "stable_visit_offset_sensitivity.tsv"),
  gate_decisions = file.path(gate_dir, "gate_decisions.tsv"),
  promotion_decision = file.path(gate_dir, "promotion_decision.tsv"),
  composition_state = file.path(out, "programs", "composition_sensitivity_state.tsv"),
  figure_manifest = file.path(out, "figures", "figure_manifest.tsv"),
  gate_input_checksums = file.path(gate_dir, "input_checksums.tsv"),
  gate_output_checksums = file.path(gate_dir, "output_checksums.tsv")
)
absent <- names(paths)[!file.exists(paths)]
assert_true(!length(absent), paste0(
  "Report inputs are absent: ", paste(absent, collapse = ", ")
))
paths <- vapply(paths, normalizePath, character(1L), mustWork = TRUE)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "fg", flag = "#"), "NA")
}

fmt_p <- function(x) {
  ifelse(!is.finite(x), "NA",
         ifelse(x < 1e-4, formatC(x, digits = 3L, format = "e"),
                formatC(x, digits = 4L, format = "f")))
}

escape_md <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- gsub("\\|", "\\\\|", x)
  gsub("[\r\n]+", " ", x)
}

markdown_table <- function(tab) {
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  assert_true(ncol(tab) > 0L, "Cannot render a zero-column Markdown table")
  header <- paste0("| ", paste(escape_md(names(tab)), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(tab)), collapse = " | "), " |")
  if (!nrow(tab)) return(c(header, separator))
  body <- apply(tab, 1L, function(row) {
    paste0("| ", paste(escape_md(row), collapse = " | "), " |")
  })
  c(header, separator, body)
}

prepared_inputs <- fread(paths[["prepared_input_manifest"]], na.strings = c("", "NA"))
assert_true(all(c("input_id", "path", "bytes", "sha256") %in% names(prepared_inputs)),
            "Prepared input manifest schema drift")
reproduction <- fread(paths[["reproduction_summary"]], na.strings = c("", "NA"))
coverage <- fread(paths[["signature_coverage"]], na.strings = c("", "NA"))
anchors <- fread(paths[["histology_anchors"]], na.strings = c("", "NA"))
agreement <- fread(paths[["co_primary_agreement"]], na.strings = c("", "NA"))
focal <- fread(paths[["focal_program_gate"]], na.strings = c("", "NA"))
paired <- fread(paths[["paired_tests"]], na.strings = c("", "NA"))
stable_offset <- fread(paths[["stable_visit_offset"]], na.strings = c("", "NA"))
gates <- fread(paths[["gate_decisions"]], na.strings = c("", "NA"))
decision <- fread(paths[["promotion_decision"]], na.strings = c("", "NA"))
composition <- fread(paths[["composition_state"]], na.strings = c("", "NA"))
figures <- fread(paths[["figure_manifest"]], na.strings = c("", "NA"))
gate_inputs <- fread(paths[["gate_input_checksums"]], na.strings = c("", "NA"))
gate_outputs <- fread(paths[["gate_output_checksums"]], na.strings = c("", "NA"))

assert_true(nrow(decision) == 1L, "Promotion decision must contain one row")
assert_true(identical(as.character(decision$promotion_status),
                      "candidate_only_not_promoted"),
            "Report refuses a decision that claims automatic promotion")
assert_true(!isTRUE(decision$automatically_promoted),
            "Report refuses automatic promotion")
assert_true(!isTRUE(decision$figure_3_modified),
            "Report refuses a Figure 3 mutation")
assert_true(nrow(composition) == 1L && grepl("^deferred_", composition$state),
            "Composition sensitivity is not recorded as deferred")
assert_true(nrow(figures) == 4L &&
              all(c("figure_id", "role", "path", "sha256") %in% names(figures)),
            "Candidate figure manifest schema drift")

co_primary <- as.character(pre$axes$co_primary)
independent <- as.character(pre$cohorts$independent_primary)
all_resource_cohorts <- c(as.character(pre$cohorts$source_overlap), independent)

released <- reproduction[metric == "released_pc1_spearman"]
assert_true(nrow(released) == 1L, "Released-PC1 result is missing")
released_rho <- as.numeric(released$value)
operational_sign <- reproduction[metric == "operational_pc1_multiplier"]
assert_true(nrow(operational_sign) == 1L &&
              as.numeric(operational_sign$value) ==
                pre$orientation$operational_continuum_multiplier,
            "Released operational-PC1 sign evidence is missing")

coverage_show <- coverage[match(all_resource_cohorts, dataset), .(
  Cohort = dataset,
  Participants = n_participants,
  `Signature genes` = paste0(n_signature_used, "/", n_signature_total),
  Coverage = paste0(fmt_num(100 * coverage_fraction, 4L), "%"),
  Pass = coverage_fraction >= pre$signature$minimum_coverage_fraction
)]

agreement_show <- agreement[
  dataset %in% independent,
  .(
    Cohort = dataset,
    Donors = n_donors,
    `Spearman rho` = fmt_num(spearman_rho),
    `Threshold` = paste0(">=", fmt_num(
      pre$promotion_gates$co_primary_agreement_spearman_minimum, 2L
    )),
    Pass = spearman_rho >= pre$promotion_gates$co_primary_agreement_spearman_minimum
  )
]
agreement_show[, order__ := match(Cohort, independent)]
setorder(agreement_show, order__)
agreement_show[, order__ := NULL]

anchor_show <- anchors[
  dataset %in% independent & axis_id %in% co_primary,
  .(
    Score = axis_id,
    Cohort = dataset,
    Donors = n_donors,
    `Fibrosis Spearman rho` = fmt_num(spearman_rho),
    `Bootstrap 95% CI` = paste0(
      "[", fmt_num(bootstrap_ci_low), ", ", fmt_num(bootstrap_ci_high), "]"
    ),
    Pass = spearman_rho > 0 & bootstrap_ci_low > 0
  )
]
anchor_show[, `:=`(
  score_order__ = match(Score, co_primary),
  cohort_order__ = match(Cohort, independent)
)]
setorder(anchor_show, score_order__, cohort_order__)
anchor_show[, c("score_order__", "cohort_order__") := NULL]

focal_show <- focal[, .(
  Program = focal_program,
  `Signature-retained L1` = paste0(fmt_num(100 * minimum_retained_l1_fraction, 4L), "%"),
  `Observed mapped L1` = paste0(fmt_num(100 * minimum_observed_l1_fraction, 4L), "%"),
  `Minimum genes` = minimum_observed_genes,
  `Signature PC1 beta` = fmt_num(signature_pc1_beta),
  `Signature PC1 BH q` = fmt_p(signature_pc1_q_value),
  `Fixed projection beta` = fmt_num(fixed_projection_beta),
  `Fixed projection BH q` = fmt_p(fixed_projection_q_value),
  `Cross-score sign concordant` = cross_axis_meta_beta_sign_concordant,
  Pass = program_gate_pass
)]

paired_show <- paired[
  endpoint == "delta_fibrosis" & axis_id %in% co_primary,
  .(
    Score = axis_id,
    Donors = n_donors,
    `Delta continuum vs delta fibrosis rho` = fmt_num(spearman_rho),
    `Permutation p` = fmt_p(permutation_p_value),
    `Holm p` = fmt_p(holm_p_value),
    Pass = spearman_rho > 0 &
      holm_p_value < pre$promotion_gates$paired_holm_maximum
  )
]
paired_show[, order__ := match(Score, co_primary)]
setorder(paired_show, order__)
paired_show[, order__ := NULL]

offset_show <- stable_offset[, .(
  Score = axis_id,
  `Stable-pair median offset` = fmt_num(stable_visit_offset),
  `Adjusted rho` = fmt_num(adjusted_spearman_rho),
  `Adjusted permutation p` = fmt_p(adjusted_permutation_p_value),
  `Rank invariant` = algebraically_rank_invariant
)]

family_show <- gates[gate_level == "family", .(
  Gate = gate_family,
  `Passing components` = paste0(observed_primary, "/", observed_secondary),
  Pass = pass
)]

figure_show <- figures[, .(
  Panel = figure_id,
  Role = role,
  SHA256 = sha256,
  Path = path
)]

assert_true(all(c("overall_pass", "disposition", "figure_route",
                  "explicit_approval_required") %in% names(decision)),
            "Promotion decision renderer contract drift")
assert_true(decision$figure_route %in% c("figure3_proposal", "figure_s3_only"),
            "Unknown figure route")
assert_true(isTRUE(decision$explicit_approval_required),
            "Explicit approval must remain required")
overall_pass <- isTRUE(decision$overall_pass)
decision_sentence <- if (overall_pass) {
  paste(
    "All prespecified gates passed. This result is not a promotion.",
    "The candidate remains unpromoted until explicit approval and a successful",
    "scripts/manuscript/validate_resource_scope.py run."
  )
} else {
  paste(
    "At least one prespecified gate failed. The analysis remains a Figure S3",
    "candidate, and Figure 3 remains unchanged."
  )
}

input_show <- gate_inputs[, .(
  Direction = "input",
  Artifact = artifact_role,
  Digest = sha256,
  Path = path
)]
source_show <- prepared_inputs[, .(
  Direction = "frozen_source",
  Artifact = input_id,
  Digest = sha256,
  Path = path
)]
output_show <- gate_outputs[, .(
  Direction = "output",
  Artifact = artifact_role,
  Digest = sha256,
  Path = path
)]
provenance_show <- rbindlist(
  list(source_show, input_show, output_show), use.names = TRUE
)

lines <- c(
  "# Histology-Anchored Molecular Continuum candidate report",
  "",
  paste0("Status: **", decision$promotion_status, "**. ", decision_sentence),
  "",
  "## Prespecified claim",
  "",
  paste0("> ", pre$scientific_claim),
  "",
  paste(
    "This is a donor-level, cross-sectional analysis. Recorded fibrosis stage is",
    "the clinical anchor. The continuum is a complementary molecular layer, not a",
    "replacement for histology and not evidence of longitudinal patient progression."
  ),
  "",
  "## Scope decision",
  "",
  paste(
    "The prior custom multi-lineage Slingshot benchmark is preserved as an alternative",
    "sensitivity. It did not implement the competitor release's operational downstream",
    "PC1, so it cannot support a conclusion that histology wins over the released",
    "data-driven continuum method."
  ),
  "",
  paste(
    "No existing Figure 3 file, NMF definition, Hotspot weight, Gene Catalog evidence",
    "class, hero gene, or biomarker claim was changed. Composition adjustment remains",
    "deferred until the synchronized five-cohort composition release is promoted."
  ),
  "",
  "## Gate decision",
  "",
  markdown_table(family_show),
  "",
  paste0(
    "Current placement: `", decision$current_placement, "`. Recommended placement: `",
    decision$recommended_placement, "`. Decision disposition: `", decision$disposition,
    "`. Figure route: `", decision$figure_route,
    "`. Explicit approval required: `TRUE`. Automatic promotion: `FALSE`."
  ),
  "",
  "## Exact released-PC1 reproduction",
  "",
  paste0(
    "The donor-level reproduced PC1 versus released PC1 Spearman rho was ",
    fmt_num(released_rho, 8L), "; absolute rho was ", fmt_num(abs(released_rho), 8L),
    " against the prespecified >=",
    fmt_num(pre$promotion_gates$reproduction_abs_spearman_minimum, 2L), " gate."
  ),
  "",
  paste(
    "The frozen reproduction uses log2(count+1), quantile normalization,",
    "ComBat batch=Dataset x Sex with NAS protected, and PCA. The released downstream",
    "ordering is the supplied descending raw-PC1 order, equivalent to increasing",
    "-PC1, not extracted Slingshot pseudotime."
  ),
  "",
  paste(
    "The raw released CSV coordinate is retained for exact reproduction. The frozen",
    "v1.0 downstream source explicitly negates that coordinate, and its RF code calls",
    "inverse PC1 the real trajectory. Progression-facing",
    "loadings use -raw PC1, consistent with the paper's low=early and high=advanced",
    "definition. This direction was fixed from discovery source code; target histology",
    "was not used to orient any validation axis."
  ),
  "",
  "## Signature coverage",
  "",
  markdown_table(coverage_show),
  "",
  "All coverage fractions use the frozen 145-gene denominator and GENCODE v49 Ensembl base IDs.",
  "",
  "## Independent-cohort score agreement",
  "",
  markdown_table(agreement_show),
  "",
  "GSE162694 and GSE213621 are the competitor-independent primary validation cohorts.",
  "",
  "## Independent-cohort fibrosis calibration",
  "",
  markdown_table(anchor_show),
  "",
  paste0(
    "Confidence intervals are unstratified donor-bootstrap intervals with ",
    pre$resampling$bootstrap_replicates, " replicates. Positive rho and a lower bound",
    " above zero were required for both co-primary scores in both cohorts."
  ),
  "",
  "## Prespecified focal Hotspot programs",
  "",
  markdown_table(focal_show),
  "",
  paste0(
    "Every 145-signature gene was excluded before scoring. The 80% rule uses original",
    " L1 mass remaining after that exclusion, matching the frozen 95.7% and 80.5%",
    " focal audits. Observed mapped L1 is reported separately; original relative",
    " weights are renormalized over observed genes for scoring. Program meta-analysis",
    " used fibrosis-stage and inferred-sex",
    " adjustment within each independent cohort, fixed-effect cohort synthesis, and",
    " BH adjustment with n=", pre$programs$family_size,
    ". The focal gate required the same program to pass for both co-primary scores,",
    " concordant cohort directions for each score, and concordant meta-effect signs",
    " across scores."
  ),
  "",
  "## Paired within-person replication",
  "",
  markdown_table(paired_show),
  "",
  markdown_table(offset_show),
  "",
  paste0(
    "All 54 QC-valid GSE193066 participant pairs, including stable-fibrosis pairs,",
    " were retained. Tests were one-sided permutation Spearman tests with ",
    pre$resampling$permutation_replicates,
    " permutations; the two co-primary delta-fibrosis tests form the Holm family. NAS",
    " change is a sensitivity. The stable-pair median visit offset is reported as",
    " an algebraically rank-invariant calibration audit, not independent evidence."
  ),
  "",
  "## Cohort roles",
  "",
  paste0(
    "Source-overlap replication: ",
    paste(as.character(pre$cohorts$source_overlap), collapse = ", "), "."
  ),
  paste0(
    "Competitor-independent primary validation: ", paste(independent, collapse = ", "), "."
  ),
  paste0("Paired within-person replication: ", pre$cohorts$paired_replication, "."),
  "",
  "## Candidate-only figures",
  "",
  markdown_table(figure_show),
  "",
  paste(
    "These PDFs exist only inside the timestamped candidate. Passing gates route",
    "two panels to a Figure 3 proposal; failing gates route them to Figure S3.",
    "Neither route edits the current manuscript figure without explicit approval."
  ),
  "",
  "## Input and gate-output provenance",
  "",
  markdown_table(provenance_show),
  "",
  "The candidate root is:",
  "",
  paste0("`", out, "`"),
  ""
)

report_dir <- file.path(out, "report")
ensure_new_dir(report_dir)
report_path <- file.path(report_dir, "candidate_report.md")
assert_true(!file.exists(report_path), paste0("Refusing to overwrite file: ", report_path))
writeLines(lines, report_path, useBytes = TRUE)
write_session_info(file.path(report_dir, "sessionInfo.txt"))

report_artifacts <- c(
  candidate_report = report_path,
  report_session_info = file.path(report_dir, "sessionInfo.txt")
)
report_provenance <- rbindlist(list(
  prepared_inputs[, .(
    artifact_class = "frozen_source",
    artifact_role = input_id,
    path,
    size_bytes = bytes,
    sha256
  )],
  gate_inputs[, .(artifact_class = "gate_input", artifact_role, path, size_bytes, sha256)],
  gate_outputs[, .(artifact_class = "gate_output", artifact_role, path, size_bytes, sha256)],
  figures[, .(
    artifact_class = "candidate_figure",
    artifact_role = figure_id,
    path,
    size_bytes = bytes,
    sha256
  )],
  data.table(
    artifact_class = "report_output",
    artifact_role = names(report_artifacts),
    path = unname(report_artifacts),
    size_bytes = as.numeric(file.size(unname(report_artifacts))),
    sha256 = vapply(unname(report_artifacts), sha256_file, character(1L))
  )
), use.names = TRUE)
write_tsv_once(report_provenance, file.path(report_dir, "provenance.tsv"))

message(
  "CANDIDATE_REPORT_COMPLETE: ", report_path,
  "; all_prespecified_gates_pass=", overall_pass,
  "; promotion_status=candidate_only_not_promoted"
)
