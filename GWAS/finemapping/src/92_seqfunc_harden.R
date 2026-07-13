#!/usr/bin/env Rscript

# Additive hardening of current seqfunc artifacts. This script never overwrites
# the historical result tables; it emits calibrated/reconciled replacements and
# an explicit apply-only firewall audit.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

ROOT <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
SF <- file.path(ROOT, "GWAS/finemapping/results/seqfunc")
OUT <- file.path(SF, "hardening")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

num <- function(x) suppressWarnings(as.numeric(x))
bool <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "YES", "1")

# ---------------------------------------------------------------- Decima repair
decima_in <- file.path(SF, "decima_celltype.tsv")
if (!file.exists(decima_in)) stop("Missing ", decima_in)
d <- fread(decima_in)
abs_cols <- grep("^abs_", names(d), value = TRUE)
if (length(abs_cols) < 2) stop("Decima table lacks compartment absolute-effect columns")
for (z in abs_cols) set(d, j = z, value = num(d[[z]]))
d[, argmax_abs_effect := num(argmax_abs_effect)]

second_max <- function(x) {
  x <- sort(x[is.finite(x)], decreasing = TRUE)
  if (length(x) < 2) return(0)
  x[2]
}
d[, second_abs_effect := apply(.SD, 1, second_max), .SDcols = abs_cols]
d[, argmax_margin := fifelse(
  second_abs_effect > 0,
  argmax_abs_effect / second_abs_effect,
  Inf
)]
d[, confident_call := argmax_abs_effect >= 0.005 & argmax_margin >= 1.2]
d[, is_coding_variant := grepl("coding_effector", var_set, fixed = TRUE)]
d[, regulatory_celltype_eligible := confident_call & !is_coding_variant]
d[, celltype_call_status := fifelse(
  is_coding_variant,
  "ineligible_coding_variant_expression_VEP_not_mechanism_localizer",
  fifelse(
    !confident_call,
    "below_effect_or_margin_floor",
    "coarse_regulatory_hypothesis"
  )
)]
d[, interpretation_scope := "coarse_hypothesis_only_pending_primary_liver_model"]

fwrite(
  d,
  file.path(OUT, "decima_celltype_calibrated.tsv"),
  sep = "\t", na = ""
)

reg <- d[is_coding_variant == FALSE]
eligible <- reg[regulatory_celltype_eligible == TRUE]
decima_summary <- list(
  status = "reconciled_from_existing_decima_effect_table",
  historical_input = "GWAS/finemapping/results/seqfunc/decima_celltype.tsv",
  historical_input_predates_confidence_columns = TRUE,
  confidence_contract = list(
    absolute_effect_floor = 0.005,
    top1_top2_margin_floor = 1.2
  ),
  n_all = nrow(d),
  n_coding_ineligible = sum(d$is_coding_variant),
  n_regulatory = nrow(reg),
  n_regulatory_confident = nrow(eligible),
  regulatory_confident_distribution = as.list(table(eligible$argmax_celltype)),
  inference = paste(
    "Only confident noncoding variants receive a coarse cell-type hypothesis.",
    "Coding variants are explicitly ineligible because an expression VEP cannot",
    "localize a protein-altering mechanism. Adult primary-liver ChromBPNet models",
    "are required before cell-type localization is treated as more than a hypothesis."
  )
)
write_json(
  decima_summary,
  file.path(OUT, "decima_celltype_calibrated_summary.json"),
  pretty = TRUE, auto_unbox = TRUE
)

# -------------------------------------------------------------- firewall audit
consumer_candidates <- c(
  file.path(ROOT, "RNA-seq/27a_assemble_evidence_atlas.R"),
  file.path(ROOT, "RNA-seq/27b_benchmark_presets.R"),
  file.path(ROOT, "RNA-seq/46d_convergence_evidence.R"),
  file.path(ROOT, "RNA-seq/75_integrate_causal_overhaul.R"),
  file.path(ROOT, "RNA-seq/217_integrate_stratified_causal.R")
)
consumer_candidates <- consumer_candidates[file.exists(consumer_candidates)]
forbidden <- c(
  "seqfunc", "borzoi", "alphagenome", "chrombpnet", "decima_celltype",
  "mpra_benchmark", "coding_hardening_v2", "direction_hardened"
)
audit <- rbindlist(lapply(consumer_candidates, function(path) {
  txt <- readLines(path, warn = FALSE)
  hits <- which(vapply(txt, function(line) {
    any(vapply(forbidden, function(p) grepl(p, line, ignore.case = TRUE), logical(1)))
  }, logical(1)))
  if (!length(hits)) {
    return(data.table(
      file = sub(paste0("^", ROOT, "/"), "", path),
      status = "PASS", line = NA_integer_, text = ""
    ))
  }
  data.table(
    file = sub(paste0("^", ROOT, "/"), "", path),
    status = "FAIL", line = hits, text = trimws(txt[hits])
  )
}))
fwrite(audit, file.path(OUT, "apply_only_firewall_audit.tsv"), sep = "\t", na = "")

firewall_pass <- nrow(audit[status == "FAIL"]) == 0
contract <- list(
  status = if (firewall_pass) "PASS" else "FAIL",
  apply_only_firewall = firewall_pass,
  canonical_consumers_checked = nrow(audit),
  prohibited_behaviors = c(
    "No seqfunc field enters the atlas or convergence score.",
    "No AlphaGenome output is used to train or calibrate another model.",
    "No per-locus eQTL direction is emitted from a failed benchmark.",
    "QTL modalities are not counted as independent convergence channels."
  ),
  supervised_direction_script_80 = "disabled_retired_ToS_noncompliant",
  calibrated_decima_output = "hardening/decima_celltype_calibrated.tsv"
)
write_json(contract, file.path(OUT, "seqfunc_contract.json"),
           pretty = TRUE, auto_unbox = TRUE)

readme <- c(
  "# Seqfunc hardening outputs",
  "",
  "This directory contains additive corrected contracts derived from the historical",
  "seqfunc artifacts. Historical outputs are preserved for provenance.",
  "",
  "## Decima",
  "",
  paste0("- Existing variants: ", nrow(d), "."),
  paste0("- Coding variants explicitly ineligible for mechanism localization: ",
         sum(d$is_coding_variant), "."),
  paste0("- Confident noncoding coarse hypotheses: ", nrow(eligible), "."),
  "- Confidence floor: |effect| >= 0.005 and top1/top2 margin >= 1.2.",
  "- These remain coarse hypotheses pending adult primary-liver accessibility models.",
  "",
  "## Apply-only firewall",
  "",
  paste0("- Status: **", if (firewall_pass) "PASS" else "FAIL", "**."),
  "- Canonical atlas/convergence consumers were scanned for seqfunc dependencies.",
  "- The retired supervised direction script is disabled because it fit models on",
  "  AlphaGenome outputs; zero-shot Scripts 84 and 85 are the valid replacements."
)
writeLines(readme, file.path(OUT, "README.md"))

message("[92 harden] Decima calibrated noncoding calls: ", nrow(eligible),
        "/", nrow(reg))
message("[92 harden] Apply-only firewall: ", if (firewall_pass) "PASS" else "FAIL")
